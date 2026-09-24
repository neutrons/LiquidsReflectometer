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

from lr_reduction.exceptions import LogNotFoundError, LogTypeError, SampleLogsError
from lr_reduction.operations.algorithms.cd_attenuation_correction import CdAttenuationCorrection
from lr_reduction.properties.cd_attenuation import _load_cd_attenuation_data, get_transmission_workspace
from lr_reduction.utils.sample_logs import SampleLogs

# Event weights are stored as float32.
RTOL = 1e-6

OUTPUT_WORKSPACE = "cd_attenuation_corrected"


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


@pytest.fixture
def make_direct_beam():
    """Factory for a small four-pixel event workspace in wavelength (~2.1-14.8 Angstrom).

    ``atten`` is recorded as the vector ``Atten`` log unless it is None; deletes each workspace
    on teardown.
    """
    names = []

    def make(atten=None):
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
        if atten is not None:
            SampleLogs(name).insert("Atten", atten)
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


def test_correction_uses_the_stages_flagged_in_the_atten_log(make_direct_beam):
    """Stages 1 and 3 in the beam put 57.5 + 249.5 microns of Cd in its path."""
    direct_beam = make_direct_beam(atten=[1, 0, 1, 0])

    corrected = run_algorithm(direct_beam)

    assert isinstance(corrected, IEventWorkspace)
    assert corrected.getNumberEvents() == mtd[direct_beam].getNumberEvents()
    assert_corrected(corrected, (57.5 + 249.5) * 1e-4)


def test_flip_atten_reads_zero_as_in_the_beam(make_direct_beam):
    corrected = run_algorithm(make_direct_beam(atten=[0, 1, 1, 1]), FlipAtten=True)

    assert_corrected(corrected, 57.5e-4)


def test_missing_atten_log_raises(make_direct_beam):
    with pytest.raises(LogNotFoundError, match="Atten"):
        run_algorithm(make_direct_beam())


@pytest.mark.parametrize(
    ("atten", "match"),
    [
        ([1, 0, 1], "fewer than the 4"),
        ([1, 0, 2, 0], "0 or 1"),
        (5, "one flag per attenuator stage"),
    ],
)
def test_malformed_atten_log_raises(make_direct_beam, atten, match):
    with pytest.raises(SampleLogsError, match=match):
        run_algorithm(make_direct_beam(atten=atten))


def test_time_series_atten_log_raises(make_direct_beam):
    """A series' value is its list of entries over time, not one flag per stage."""
    direct_beam = make_direct_beam()
    series = Int32TimeSeriesProperty("Atten")
    for second, state in enumerate([1, 0, 1, 0]):
        series.addValue(f"2020-01-01T00:00:0{second}", state)
    mtd[direct_beam].getRun().addProperty("Atten", series, True)

    with pytest.raises(LogTypeError, match="time series"):
        run_algorithm(direct_beam)


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


@pytest.fixture(scope="module")
def event_workspace(nexus_dir):
    """A REF_L run in wavelength (12.9-17.1 Angstrom, inside the table), loaded once for this module."""
    name = "cd_attenuation_input"
    LoadEventNexus(Filename=os.path.join(nexus_dir, "REF_L_198409.nxs.h5"), OutputWorkspace=name)
    ConvertUnits(InputWorkspace=name, Target="Wavelength", OutputWorkspace=name)
    SampleLogs(name).insert("Atten", [1, 0, 1, 0])
    yield name
    DeleteWorkspace(name)


@pytest.mark.datarepo
def test_correction_on_a_measured_run(event_workspace):
    corrected = run_algorithm(event_workspace)

    assert isinstance(corrected, IEventWorkspace)
    assert corrected.getNumberEvents() == mtd[event_workspace].getNumberEvents()
    assert_corrected(corrected, (57.5 + 249.5) * 1e-4)
