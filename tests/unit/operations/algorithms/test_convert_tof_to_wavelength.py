"""Unit tests for the ConvertTofToWavelength Mantid algorithm."""

import os

import numpy as np
import pytest
from mantid.api import AlgorithmManager, IEventWorkspace
from mantid.kernel import DateAndTime, PhysicalConstants
from mantid.simpleapi import (
    CloneWorkspace,
    ConvertToEventWorkspace,
    ConvertUnits,
    CreateSampleWorkspace,
    CreateWorkspace,
    DeleteWorkspace,
    LoadEventNexus,
    LoadInstrument,
    ModeratorTzero,
    MoveInstrumentComponent,
    SetInstrumentParameter,
    mtd,
)

from lr_reduction.models.config import EmissionTimeConfig
from lr_reduction.operations.algorithms.convert_tof_to_wavelength import (
    FORMULA_PARAMETER,
    T0_LOG,
    T1_LOG,
    ConvertTofToWavelength,
    emission_time_formula,
)
from lr_reduction.operations.emission_time import (
    CHOPPER_MODERATOR_LOG,
    CHOPPER_MULTIPLIER_LOG,
    CHOPPER_OFFSET_LOG,
    resolve_emission_coefficients,
)
from lr_reduction.utils.sample_logs import SampleLogs

# Realistic coefficients: t1 is a few us/A against m_n * L / h ~ 3460 us/A on REF_L.
T0 = 10.0  # us
T1 = 2.0  # us / A

TOF_RANGE = (5000.0, 60000.0)  # us, about 1.4-17 Angstrom on REF_L
EVENTS_PER_PIXEL = 200

# The design's tolerance for the exact formula, met only with ModeratorTzero's Niter=2 and tolTOF=0.1.
EXACT_RTOL = 1e-8
# The documented gap between the two-pass iteration and the exact answer.
TWO_PASS_RTOL = 1e-5


def closed_form(tof_us, path_m, t0=T0, t1=T1):
    """lambda = h (TOF - t0) / (m_n L + h t1), in Angstrom, with Mantid's own constants."""
    h, m_n = PhysicalConstants.h, PhysicalConstants.NeutronMass
    return h * (tof_us - t0) * 1e-6 / (m_n * path_m + h * t1 * 1e-6 * 1e10) * 1e10


def two_pass(tof_us, path_m, t0=T0, t1=T1):
    """The two-pass iteration of requirement §2.1.3, as the production workflow computes it."""
    h, m_n = PhysicalConstants.h, PhysicalConstants.NeutronMass
    us_per_angstrom = m_n * path_m / h * 1e-10 * 1e6
    wavelength_0 = tof_us / us_per_angstrom
    return (tof_us - (t0 + t1 * wavelength_0)) / us_per_angstrom


def flight_path(workspace, index):
    spectrum_info = workspace.spectrumInfo()
    return spectrum_info.l1() + spectrum_info.l2(index)


def assert_matches(output, reference_input, oracle, rtol):
    """Every output event's wavelength matches `oracle(tof, L)` of the same event in `reference_input`."""
    assert reference_input.getNumberEvents() > 0, "comparison would pass vacuously on an empty workspace"
    for index in range(reference_input.getNumberHistograms()):
        expected = oracle(reference_input.getSpectrum(index).getTofs(), flight_path(reference_input, index))
        np.testing.assert_allclose(output.getSpectrum(index).getTofs(), expected, rtol=rtol, atol=0)


def run(input_workspace, output_name=None, **properties):
    """Run the algorithm directly and return the output workspace."""
    algorithm = ConvertTofToWavelength()
    algorithm.initialize()
    algorithm.setChild(True)
    algorithm.setProperty("InputWorkspace", input_workspace)
    for name, value in properties.items():
        algorithm.setProperty(name, value)
    algorithm.setPropertyValue("OutputWorkspace", output_name or mtd.unique_hidden_name())
    algorithm.execute()
    return algorithm.getProperty("OutputWorkspace").value


@pytest.fixture
def workspaces():
    """Tracks workspaces created by a test and deletes them on teardown."""
    names = []
    yield names
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


@pytest.fixture
def event_workspace(workspaces):
    """Factory for a small four-pixel event workspace in TOF carrying the real REF_L instrument."""

    def make():
        name = mtd.unique_hidden_name()
        workspaces.append(name)
        workspace = CreateSampleWorkspace(
            WorkspaceType="Event",
            NumBanks=1,
            BankPixelWidth=2,
            NumEvents=0,
            XMin=TOF_RANGE[0],
            XMax=TOF_RANGE[1],
            BinWidth=100.0,
            OutputWorkspace=name,
        )
        # Evenly spaced events over the whole band, so every pixel is populated and the
        # comparison covers short and long wavelengths alike.
        pulse_time = DateAndTime("2026-01-01T00:00:00")
        for index in range(workspace.getNumberHistograms()):
            events = workspace.getSpectrum(index)
            for tof in np.linspace(*TOF_RANGE, EVENTS_PER_PIXEL):
                events.addEventQuickly(float(tof), pulse_time)
        # Keep the sample workspace's spectrum map: the four pixels map to REF_L detector IDs,
        # so each spectrum gets a real pixel position and flight path.
        LoadInstrument(Workspace=name, InstrumentName="REF_L", RewriteSpectraMap=False)
        return mtd[name]

    return make


def test_matches_exact_formula_per_pixel(event_workspace):
    workspace = event_workspace()
    output = run(workspace, T0=T0, T1=T1)
    assert_matches(output, workspace, closed_form, EXACT_RTOL)


def test_single_iteration_misses_exact_formula(event_workspace):
    # Pins the iteration settings: one pass (the two-pass result) is not good enough for the
    # tolerance the algorithm promises, so ModeratorTzero must be told to iterate.
    workspace = event_workspace()
    copy = mtd.unique_hidden_name()
    CloneWorkspace(InputWorkspace=workspace, OutputWorkspace=copy)
    SetInstrumentParameter(
        Workspace=copy, ParameterName=FORMULA_PARAMETER, ParameterType="String", Value=emission_time_formula(T0, T1)
    )
    ModeratorTzero(InputWorkspace=copy, EMode="Elastic", Niter=1, OutputWorkspace=copy)
    once = ConvertUnits(InputWorkspace=copy, Target="Wavelength", EMode="Elastic", OutputWorkspace=copy)
    try:
        with pytest.raises(AssertionError):
            assert_matches(once, workspace, closed_form, EXACT_RTOL)
    finally:
        DeleteWorkspace(copy)


def test_zero_coefficients_is_plain_convert_units(event_workspace, workspaces):
    workspace = event_workspace()
    output = run(workspace, T0=0.0, T1=0.0)

    plain_name = mtd.unique_hidden_name()
    workspaces.append(plain_name)
    plain = ConvertUnits(InputWorkspace=workspace, Target="Wavelength", EMode="Elastic", OutputWorkspace=plain_name)
    assert plain.getNumberEvents() > 0
    for index in range(workspace.getNumberHistograms()):
        np.testing.assert_array_equal(output.getSpectrum(index).getTofs(), plain.getSpectrum(index).getTofs())
    assert not output.getInstrument().hasParameter(FORMULA_PARAMETER)


def test_input_is_unchanged(event_workspace):
    workspace = event_workspace()
    assert workspace.getNumberEvents() > 0
    tofs_before = [workspace.getSpectrum(i).getTofs().copy() for i in range(workspace.getNumberHistograms())]

    run(workspace, T0=T0, T1=T1)

    assert workspace.getAxis(0).getUnit().unitID() == "TOF"
    assert not workspace.getInstrument().hasParameter(FORMULA_PARAMETER)
    for index, before in enumerate(tofs_before):
        np.testing.assert_array_equal(workspace.getSpectrum(index).getTofs(), before)


def test_output_is_event_workspace_in_wavelength_carrying_formula(event_workspace):
    workspace = event_workspace()
    output = run(workspace, T0=T0, T1=T1)

    assert isinstance(output, IEventWorkspace)
    assert output.getAxis(0).getUnit().unitID() == "Wavelength"
    assert output.getNumberEvents() == workspace.getNumberEvents()
    assert list(output.getInstrument().getStringParameter(FORMULA_PARAMETER)) == [emission_time_formula(T0, T1)]


def test_geometry_reaches_conversion_through_instrument(event_workspace):
    # Moving the moderator changes L1, and with it every wavelength, exactly as the formula
    # predicts for the new flight path: the algorithm takes geometry from nowhere else.
    workspace = event_workspace()
    l1_before = workspace.spectrumInfo().l1()
    MoveInstrumentComponent(Workspace=workspace, ComponentName="moderator", Z=-15.0, RelativePosition=False)
    assert workspace.spectrumInfo().l1() != pytest.approx(l1_before)

    output = run(workspace, T0=T0, T1=T1)
    assert_matches(output, workspace, closed_form, EXACT_RTOL)


def test_sample_logs_record_applied_coefficients(event_workspace):
    output = run(event_workspace(), T0=T0, T1=T1)
    logs = SampleLogs(output)
    assert logs[T0_LOG] == T0
    assert logs[T1_LOG] == T1


def test_agrees_with_two_pass_calculation(event_workspace):
    workspace = event_workspace()
    output = run(workspace, T0=T0, T1=T1)
    assert_matches(output, workspace, two_pass, TWO_PASS_RTOL)


def test_history_records_algorithm_and_coefficients(event_workspace, workspaces):
    # Run through the factory as a non-child, the way the operations invoke registered
    # algorithms, so the call is recorded in the output's history.
    workspace = event_workspace()
    output_name = mtd.unique_hidden_name()
    workspaces.append(output_name)
    algorithm = AlgorithmManager.create("ConvertTofToWavelength")
    algorithm.setProperty("InputWorkspace", workspace.name())
    algorithm.setProperty("T0", T0)
    algorithm.setProperty("T1", T1)
    algorithm.setProperty("OutputWorkspace", output_name)
    algorithm.execute()

    last = mtd[output_name].getHistory().lastAlgorithm()
    assert last.name() == "ConvertTofToWavelength"
    assert float(last.getPropertyValue("T0")) == T0
    assert float(last.getPropertyValue("T1")) == T1


def test_rejects_input_not_in_tof(event_workspace, workspaces):
    workspace = event_workspace()
    name = mtd.unique_hidden_name()
    workspaces.append(name)
    in_wavelength = ConvertUnits(InputWorkspace=workspace, Target="Wavelength", OutputWorkspace=name)
    with pytest.raises(RuntimeError, match="InputWorkspace: X axis must be in TOF"):
        run(in_wavelength)


def test_rejects_workspace_without_source_or_sample(workspaces):
    histogram_name, event_name = mtd.unique_hidden_name(), mtd.unique_hidden_name()
    workspaces.extend([histogram_name, event_name])
    CreateWorkspace(DataX=[1000.0, 2000.0, 3000.0], DataY=[1.0, 1.0], UnitX="TOF", OutputWorkspace=histogram_name)
    no_instrument = ConvertToEventWorkspace(InputWorkspace=histogram_name, OutputWorkspace=event_name)
    with pytest.raises(RuntimeError, match="InputWorkspace: instrument has no source or sample position"):
        run(no_instrument)


@pytest.mark.parametrize("name", ["T0", "T1"])
@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_rejects_non_finite_coefficient(event_workspace, name, value):
    with pytest.raises(RuntimeError, match=f"{name}: must be a finite number"):
        run(event_workspace(), **{name: value})


def test_formula_round_trips_negative_and_small_coefficients(event_workspace):
    # The formula is text: a negative slope and an exponent-notation offset must survive it.
    workspace = event_workspace()
    output = run(workspace, T0=1e-5, T1=-1.5)
    assert_matches(output, workspace, lambda tof, path: closed_form(tof, path, t0=1e-5, t1=-1.5), EXACT_RTOL)


@pytest.fixture(scope="module")
def real_run(nexus_dir):
    """A real REF_L run with chopper emission-time logs injected: no run in the test-data
    repository records them, so the values here stand in for what the DAS writes."""
    name = mtd.unique_hidden_name()
    LoadEventNexus(Filename=os.path.join(nexus_dir, "REF_L_198409.nxs.h5"), OutputWorkspace=name)
    logs = SampleLogs(name)
    logs.insert(CHOPPER_MODERATOR_LOG, [1])
    logs.insert(CHOPPER_OFFSET_LOG, [T0], unit="microsecond")
    logs.insert(CHOPPER_MULTIPLIER_LOG, [T1], unit="microsecond/Angstrom")
    yield mtd[name]
    DeleteWorkspace(name)


@pytest.mark.datarepo
def test_real_run_agrees_with_two_pass_from_chopper_logs(real_run):
    t0, t1 = resolve_emission_coefficients(EmissionTimeConfig(), real_run)
    assert (t0, t1) == (T0, T1)

    output = run(real_run, T0=t0, T1=t1)

    assert output.getNumberEvents() == real_run.getNumberEvents()
    assert SampleLogs(output)["gd_prtn_chrg"] == SampleLogs(real_run)["gd_prtn_chrg"]
    np.testing.assert_array_equal(output.getRun()["proton_charge"].value, real_run.getRun()["proton_charge"].value)
    # Every pixel with events, not a sample of them.
    populated = [i for i in range(real_run.getNumberHistograms()) if real_run.getSpectrum(i).getNumberEvents()]
    assert populated
    for index in populated:
        expected = two_pass(real_run.getSpectrum(index).getTofs(), flight_path(real_run, index), t0=t0, t1=t1)
        np.testing.assert_allclose(output.getSpectrum(index).getTofs(), expected, rtol=TWO_PASS_RTOL, atol=0)
