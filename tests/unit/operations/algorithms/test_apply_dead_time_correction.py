"""Unit tests for the ApplyDeadTimeCorrection Mantid algorithm."""

import os

import numpy as np
import pytest
from mantid.api import IEventWorkspace
from mantid.simpleapi import (
    CloneWorkspace,
    ConvertUnits,
    CreateSampleWorkspace,
    CreateWorkspace,
    DeleteWorkspace,
    LoadEventNexus,
    mtd,
)

from lr_reduction.legacy import event_reduction
from lr_reduction.legacy.reduction_template_reader import ReductionParameters
from lr_reduction.legacy.utils import amend_config
from lr_reduction.operations.algorithms.apply_dead_time_correction import ApplyDeadTimeCorrection

# Events of the sample workspace span this TOF range.
TOF_MIN, TOF_MAX = 8000.0, 56000.0

OUTPUT_WORKSPACE = "dead_time_corrected"


@pytest.fixture
def workspaces():
    """Collects workspace names to delete on teardown, starting with the output workspace."""
    names = [OUTPUT_WORKSPACE]
    yield names
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


@pytest.fixture
def events(workspaces):
    name = mtd.unique_hidden_name()
    CreateSampleWorkspace(
        WorkspaceType="Event",
        Function="Flat background",
        NumBanks=1,
        BankPixelWidth=2,
        XMin=TOF_MIN,
        XMax=TOF_MAX,
        BinWidth=TOF_MAX - TOF_MIN,
        OutputWorkspace=name,
    )
    workspaces.append(name)
    return name


@pytest.fixture
def make_curve(workspaces):
    """Factory for a correction curve with four bins spanning ``[x_min, x_max]``."""

    def make(factors=(1.0, 1.1, 1.2, 1.3), x_min=0.0, x_max=60000.0, unit="TOF", spectra=1):
        name = mtd.unique_hidden_name()
        edges = np.linspace(x_min, x_max, len(factors) + 1)
        CreateWorkspace(
            DataX=np.tile(edges, spectra),
            DataY=np.tile(factors, spectra),
            NSpec=spectra,
            UnitX=unit,
            OutputWorkspace=name,
        )
        workspaces.append(name)
        return name

    return make


def make_algorithm(input_workspace, curve):
    algorithm = ApplyDeadTimeCorrection()
    algorithm.initialize()
    algorithm.setProperty("InputWorkspace", input_workspace)
    algorithm.setProperty("DeadTimeCorrectionWorkspace", curve)
    algorithm.setProperty("OutputWorkspace", OUTPUT_WORKSPACE)
    return algorithm


def validate(input_workspace, curve):
    return make_algorithm(input_workspace, curve).validateInputs()


def run_algorithm(input_workspace, curve):
    algorithm = make_algorithm(input_workspace, curve)
    algorithm.execute()
    return mtd[OUTPUT_WORKSPACE]


def test_valid_inputs_pass(events, make_curve):
    assert validate(events, make_curve()) == {}


def test_zero_factors_are_allowed(events, make_curve):
    """Bins above the threshold are set to zero by the curve-computing algorithm."""
    assert validate(events, make_curve(factors=(1.0, 1.1, 0.0, 0.0))) == {}


def test_input_not_in_tof_is_rejected(events, make_curve):
    ConvertUnits(InputWorkspace=events, Target="Wavelength", OutputWorkspace=events)

    issues = validate(events, make_curve())

    assert "TOF" in issues["InputWorkspace"]


def test_curve_not_in_tof_is_rejected(events, make_curve):
    issues = validate(events, make_curve(unit="Wavelength"))

    assert "TOF" in issues["DeadTimeCorrectionWorkspace"]


def test_curve_with_several_spectra_is_rejected(events, make_curve):
    issues = validate(events, make_curve(spectra=2))

    assert "one spectrum" in issues["DeadTimeCorrectionWorkspace"]


@pytest.mark.parametrize(("x_min", "x_max"), [(TOF_MIN + 1000.0, 60000.0), (0.0, TOF_MAX - 1000.0)])
def test_curve_not_covering_the_events_is_rejected(events, make_curve, x_min, x_max):
    issues = validate(events, make_curve(x_min=x_min, x_max=x_max))

    assert "does not cover" in issues["DeadTimeCorrectionWorkspace"]


def test_negative_factor_is_rejected(events, make_curve):
    issues = validate(events, make_curve(factors=(1.0, -0.1, 1.2, 1.3)))

    assert "negative" in issues["DeadTimeCorrectionWorkspace"]


@pytest.mark.parametrize("bad", [np.nan, np.inf])
def test_undefined_factor_is_rejected(events, make_curve, bad):
    issues = validate(events, make_curve(factors=(1.0, bad, 1.2, 1.3)))

    assert "undefined" in issues["DeadTimeCorrectionWorkspace"]


def test_each_event_is_weighted_by_the_factor_of_its_curve_bin(events, make_curve):
    factors = np.array([1.0, 1.1, 0.0, 1.3])
    curve = make_curve(factors=factors)

    corrected = run_algorithm(events, curve)

    assert isinstance(corrected, IEventWorkspace)
    assert corrected.getNumberEvents() == mtd[events].getNumberEvents()
    np.testing.assert_array_equal(corrected.readX(0), mtd[events].readX(0))
    edges = mtd[curve].readX(0)
    for index in range(corrected.getNumberHistograms()):
        event_list = corrected.getSpectrum(index)
        expected = factors[np.searchsorted(edges, event_list.getTofs(), side="right") - 1]
        np.testing.assert_allclose(event_list.getWeights(), expected, rtol=1e-6)


def test_events_in_a_zero_bin_get_zero_weight_and_are_kept(events, make_curve):
    corrected = run_algorithm(events, make_curve(factors=(0.0, 0.0, 0.0, 0.0)))

    assert corrected.getNumberEvents() == mtd[events].getNumberEvents()
    for index in range(corrected.getNumberHistograms()):
        np.testing.assert_array_equal(corrected.getSpectrum(index).getWeights(), 0.0)


def test_input_is_unchanged(events, make_curve):
    """Rebinning sorts the input's events by TOF, so they are compared regardless of order."""

    def snapshot(workspace):
        result = []
        for index in range(workspace.getNumberHistograms()):
            spectrum = workspace.getSpectrum(index)
            order = np.argsort(spectrum.getTofs(), kind="stable")
            result.append((spectrum.getTofs()[order], spectrum.getWeights()[order], workspace.readX(index).copy()))
        return result

    before = snapshot(mtd[events])

    run_algorithm(events, make_curve())

    for (tofs, weights, x), (tofs_after, weights_after, x_after) in zip(before, snapshot(mtd[events]), strict=True):
        np.testing.assert_array_equal(tofs_after, tofs)
        np.testing.assert_array_equal(weights_after, weights)
        np.testing.assert_array_equal(x_after, x)


def test_coarse_input_against_a_fine_curve(events, make_curve):
    """Four spectra with one bin against a 190-bin curve: plain ``Multiply`` refuses this."""
    factors = np.linspace(1.0, 1.5, 190)

    corrected = run_algorithm(events, make_curve(factors=factors))

    assert corrected.getNumberHistograms() == 4
    assert corrected.blocksize() == 1
    np.testing.assert_array_equal(corrected.readX(0), mtd[events].readX(0))
    assert corrected.getNumberEvents() == mtd[events].getNumberEvents()


def test_history_records_the_algorithm_and_its_inputs(events, make_curve):
    curve = make_curve()

    corrected = run_algorithm(events, curve)

    history = corrected.getHistory()
    last = history.getAlgorithmHistory(history.size() - 1)
    assert last.name() == "ApplyDeadTimeCorrection"
    assert last.getPropertyValue("InputWorkspace") == events
    assert last.getPropertyValue("DeadTimeCorrectionWorkspace") == curve


def test_invalid_curve_stops_the_algorithm(events, make_curve):
    with pytest.raises(RuntimeError, match="negative"):
        run_algorithm(events, make_curve(factors=(1.0, -0.1, 1.2, 1.3)))


def total_weight(workspace):
    return sum(workspace.getSpectrum(index).getWeights().sum() for index in range(workspace.getNumberHistograms()))


@pytest.mark.datarepo
def test_total_weight_matches_the_legacy_correction(nexus_dir, workspaces):
    run = "dead_time_input"
    LoadEventNexus(Filename=os.path.join(nexus_dir, "REF_L_198409.nxs.h5"), OutputWorkspace=run)
    workspaces.append(run)
    legacy = CloneWorkspace(InputWorkspace=run, OutputWorkspace="dead_time_legacy")
    workspaces.append(legacy.name())
    # The legacy code loads the run's rejected events by run number, and leaves the curve it multiplied
    # by in "corr_ws" alongside its intermediate workspaces.
    with amend_config(new_config={"default.instrument": "REF_L"}, data_dir=nexus_dir):
        legacy = event_reduction.apply_dead_time_correction(legacy, ReductionParameters())
    workspaces.extend(["corr_ws", "corr", "error_ws", "_errors", "_ws_sc"])

    corrected = run_algorithm(run, "corr_ws")

    assert total_weight(corrected) != pytest.approx(total_weight(mtd[run]))
    assert total_weight(corrected) == pytest.approx(total_weight(legacy), rel=1e-9)
