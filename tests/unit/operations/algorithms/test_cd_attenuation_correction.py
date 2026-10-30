"""Unit tests for the CdAttenuationCorrection Mantid algorithm."""

import os

import numpy as np
import pytest
from mantid.api import IEventWorkspace
from mantid.kernel import Int32TimeSeriesProperty
from mantid.simpleapi import (
    ConvertUnits,
    CreateSampleWorkspace,
    CreateWorkspace,
    DeleteWorkspace,
    LoadEventNexus,
    mtd,
)

from lr_reduction.operations.algorithms.cd_attenuation_correction import CdAttenuationCorrection
from lr_reduction.properties.cd_attenuation import (
    CD_FOILS,
    MENU_LOG,
    STAGE_LOGS,
    _load_cd_attenuation_data,
    get_transmission_workspace,
)

# Event weights are stored as float32.
RTOL = 1e-6

OUTPUT_WORKSPACE = "cd_attenuation_corrected"

FULL_STACK_THICKNESS = sum(CD_FOILS) * 1e-4


def _delete_transmission_workspaces():
    for name in mtd.getObjectNames():
        if name.startswith("transmission_"):
            DeleteWorkspace(name)


@pytest.fixture(autouse=True)
def clean_workspaces():
    """Transmission workspaces are cached in the process-global analysis data service."""
    _delete_transmission_workspaces()
    yield
    _delete_transmission_workspaces()
    if mtd.doesExist(OUTPUT_WORKSPACE):
        DeleteWorkspace(OUTPUT_WORKSPACE)


def add_int_series(ws, name, values):
    """Record `values` as a time series log, as the DAS writes the actuator logs."""
    series = Int32TimeSeriesProperty(name)
    for second, value in enumerate(values):
        series.addValue(f"2020-01-01T00:00:0{second}", value)
    mtd[ws].getRun().addProperty(name, series, True)


@pytest.fixture
def make_direct_beam():
    """Factory for a small four-pixel event workspace in wavelength (~2.1-14.8 Angstrom).

    ``menu`` is recorded as the Menu log and ``stages`` as the older per-stage logs, each unless
    it is None; deletes each workspace on teardown.
    """
    names = []

    def make(menu=None, stages=None):
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(
            WorkspaceType="Event",
            Function="Flat background",
            NumBanks=1,
            BankPixelWidth=2,
            XMin=8000.0,
            XMax=56000.0,
            BinWidth=48000.0,
            OutputWorkspace=name,
        )
        ConvertUnits(InputWorkspace=name, Target="Wavelength", OutputWorkspace=name)
        if menu is not None:
            add_int_series(name, MENU_LOG, [menu])
        if stages is not None:
            for log, state in zip(STAGE_LOGS, stages, strict=True):
                add_int_series(name, log, [state])
        names.append(name)
        return name

    yield make
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


def run_algorithm(input_workspace, **properties):
    """Run the algorithm on ``input_workspace`` and return the output workspace."""
    algorithm = CdAttenuationCorrection()
    algorithm.PyInit()
    algorithm.setProperty("InputWorkspace", input_workspace)
    algorithm.setProperty("OutputWorkspace", OUTPUT_WORKSPACE)
    for name, value in properties.items():
        algorithm.setProperty(name, value)
    algorithm.PyExec()
    return algorithm.getProperty("OutputWorkspace").value


def expected_scale(wavelengths, cd_thickness):
    """1 / T for each wavelength, read off the transmission histogram's bins."""
    _, mu = _load_cd_attenuation_data()
    bin_edges = get_transmission_workspace(cd_thickness).readX(0)
    bins = np.searchsorted(bin_edges, wavelengths, side="right") - 1
    return np.exp(mu[bins] * cd_thickness)


def assert_corrected(corrected, cd_thickness):
    """Every event's weight and error equal 1 / T at its wavelength, for unit-weight input events."""
    for index in range(corrected.getNumberHistograms()):
        event_list = corrected.getSpectrum(index)
        scale = expected_scale(event_list.getTofs(), cd_thickness)
        np.testing.assert_allclose(event_list.getWeights(), scale, rtol=RTOL)
        np.testing.assert_allclose(event_list.getWeightErrors(), scale, rtol=RTOL)


def test_algorithm_metadata():
    algorithm = CdAttenuationCorrection()

    assert algorithm.name() == "CdAttenuationCorrection"
    assert algorithm.category() == "Reflectometry\\SNS"
    assert algorithm.version() == 1


@pytest.mark.parametrize("logs", [{"menu": 5}, {"stages": [1, 0, 1, 0]}], ids=["menu", "stage_logs"])
def test_correction_uses_the_stages_in_the_beam(make_direct_beam, logs):
    """Stages 1 and 3 in the beam put 57.5 + 249.5 microns of Cd in its path."""
    direct_beam = make_direct_beam(**logs)

    corrected = run_algorithm(direct_beam)

    assert isinstance(corrected, IEventWorkspace)
    assert corrected.getNumberEvents() == mtd[direct_beam].getNumberEvents()
    assert_corrected(corrected, (57.5 + 249.5) * 1e-4)


def test_flip_atten_reads_zero_as_in_the_beam(make_direct_beam):
    """Menu code 14 flags stages 2-4, so with inverted polarity only stage 1 is in the beam."""
    corrected = run_algorithm(make_direct_beam(menu=14), FlipAtten=True)

    assert_corrected(corrected, 57.5e-4)


def test_histogram_input_is_rejected():
    """The correction is applied per event, so the input must be an event workspace."""
    histogram = CreateWorkspace(DataX=[1.0, 2.0], DataY=[1.0], UnitX="Wavelength", OutputWorkspace="cd_histogram")
    algorithm = CdAttenuationCorrection()
    algorithm.PyInit()

    try:
        with pytest.raises(ValueError):
            algorithm.setProperty("InputWorkspace", histogram)
    finally:
        DeleteWorkspace(histogram)


def _load_in_wavelength(nexus_dir, run_number):
    name = f"cd_attenuation_input_{run_number}"
    LoadEventNexus(Filename=os.path.join(nexus_dir, f"REF_L_{run_number}.nxs.h5"), OutputWorkspace=name)
    ConvertUnits(InputWorkspace=name, Target="Wavelength", OutputWorkspace=name)
    return name


@pytest.fixture(scope="module")
def menu_run(nexus_dir):
    """A REF_L run in wavelength (2.5-10.6 Angstrom) whose Menu log is code 3: stages 1 and 2 in the beam."""
    name = _load_in_wavelength(nexus_dir, 234266)
    yield name
    DeleteWorkspace(name)


@pytest.fixture(scope="module")
def full_stack_run(nexus_dir):
    """A REF_L run in wavelength (2.1-10.7 Angstrom) that predates the Menu log, with all four stages in the beam."""
    name = _load_in_wavelength(nexus_dir, 220063)
    yield name
    DeleteWorkspace(name)


@pytest.mark.datarepo
def test_correction_on_a_measured_run(menu_run):
    corrected = run_algorithm(menu_run)

    assert isinstance(corrected, IEventWorkspace)
    assert corrected.getNumberEvents() == mtd[menu_run].getNumberEvents()
    assert_corrected(corrected, (57.5 + 126.5) * 1e-4)


@pytest.mark.datarepo
def test_full_stack_on_a_measured_run_drops_overflowing_events(full_stack_run):
    """932.5 microns of Cd attenuates the longest wavelengths past what float32 weights can undo,
    so those events are dropped and the rest corrected."""
    corrected = run_algorithm(full_stack_run)

    assert 0 < corrected.getNumberEvents() < mtd[full_stack_run].getNumberEvents()
    assert_corrected(corrected, FULL_STACK_THICKNESS)
