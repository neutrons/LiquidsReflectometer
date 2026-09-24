"""Unit tests for the Cd attenuation transmission correction functions."""

import numpy as np
import pytest
from mantid.api import IEventWorkspace
from mantid.simpleapi import ConvertUnits, CreateSampleWorkspace, DeleteWorkspace, Rebin, mtd

from lr_reduction.exceptions import LrValidationError, SampleLogsError
from lr_reduction.properties import cd_attenuation
from lr_reduction.properties.cd_attenuation import (
    CD_FOILS,
    _load_cd_attenuation_data,
    apply_correction,
    compute_cd_thickness,
    get_transmission_workspace,
)

# Event weights are stored as float32.
RTOL = 1e-6

# Stages 1 and 3 in the beam: 57.5 + 249.5 microns.
MODERATE_THICKNESS = 0.0307
FULL_STACK_THICKNESS = sum(CD_FOILS) * 1e-4


def _delete_transmission_workspaces():
    for name in mtd.getObjectNames():
        if name.startswith("transmission_"):
            DeleteWorkspace(name)


@pytest.fixture(autouse=True)
def clean_transmission_cache():
    """Transmission workspaces are cached in the process-global analysis data service."""
    _delete_transmission_workspaces()
    yield
    _delete_transmission_workspaces()


@pytest.fixture
def loader_cache():
    """A cold loader cache, restored to cold afterwards so no test sees another's load."""
    _load_cd_attenuation_data.cache_clear()
    yield
    _load_cd_attenuation_data.cache_clear()


@pytest.fixture
def make_events():
    """Factory for small four-pixel event workspaces; deletes each one on teardown.

    With the default TOF bounds the events span roughly 2.1-14.8 Angstrom, inside the table.
    A single TOF bin reproduces a freshly loaded run, whose blocksize is 1.
    """
    names = []

    def make(tof_min=8000.0, tof_max=56000.0, wavelength=True):
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(
            WorkspaceType="Event",
            Function="Flat background",
            NumBanks=1,
            BankPixelWidth=2,
            XMin=tof_min,
            XMax=tof_max,
            BinWidth=tof_max - tof_min,
            OutputWorkspace=name,
        )
        if wavelength:
            ConvertUnits(InputWorkspace=name, Target="Wavelength", OutputWorkspace=name)
        names.append(name)
        return name

    yield make
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


@pytest.fixture
def output_name():
    name = mtd.unique_hidden_name()
    yield name
    if mtd.doesExist(name):
        DeleteWorkspace(name)


def expected_bin_edges():
    """Bin edges at the midpoints between table wavelengths, the outer two a half-step out."""
    wavelength, _ = _load_cd_attenuation_data()
    midpoints = (wavelength[:-1] + wavelength[1:]) / 2
    first = wavelength[0] - (wavelength[1] - wavelength[0]) / 2
    last = wavelength[-1] + (wavelength[-1] - wavelength[-2]) / 2
    return np.concatenate([[first], midpoints, [last]])


def expected_scale(wavelengths, cd_thickness):
    """1 / T for each wavelength, with μ taken from the table bin holding it."""
    _, mu = _load_cd_attenuation_data()
    bins = np.searchsorted(expected_bin_edges(), wavelengths, side="right") - 1
    return np.exp(mu[bins] * cd_thickness)


##############
### Loader ###
##############


@pytest.mark.usefixtures("loader_cache")
def test_table_is_returned_as_points():
    """The table comes back as one μ per wavelength, with no appended bin edge."""
    wavelength, mu = _load_cd_attenuation_data()

    assert len(wavelength) == len(mu) == 1000
    assert wavelength[0] == 1.5
    assert wavelength[-1] == 20.0
    assert np.all(np.diff(wavelength) > 0)


@pytest.mark.usefixtures("loader_cache")
def test_table_arrays_are_read_only():
    """Every caller shares the cached arrays, so none may write to them."""
    wavelength, mu = _load_cd_attenuation_data()

    with pytest.raises(ValueError, match="read-only"):
        wavelength[0] = 0.0
    with pytest.raises(ValueError, match="read-only"):
        mu[0] = 0.0


@pytest.mark.usefixtures("loader_cache")
def test_table_is_read_once():
    first = _load_cd_attenuation_data()
    second = _load_cd_attenuation_data()

    assert first is second
    assert _load_cd_attenuation_data.cache_info().hits == 1


@pytest.mark.usefixtures("loader_cache")
def test_missing_table_raises_and_is_not_cached(monkeypatch, tmp_path):
    """A failed read propagates numpy's error, and the next call reads the file again."""
    monkeypatch.setattr(cd_attenuation, "cd_attenuation_file", tmp_path / "missing.csv")
    with pytest.raises(FileNotFoundError):
        _load_cd_attenuation_data()

    monkeypatch.undo()
    wavelength, _ = _load_cd_attenuation_data()
    assert len(wavelength) == 1000


####################
### Cd thickness ###
####################


@pytest.mark.parametrize(
    ("atten", "thickness"),
    [
        ([0, 0, 0, 0], 0.0),
        ([1, 0, 0, 0], 0.00575),
        ([0, 1, 1, 0], 0.0376),
        ([1, 0, 1, 0], MODERATE_THICKNESS),
        ([1, 1, 1, 1], 0.09325),
    ],
)
def test_thickness_sums_the_stages_in_the_beam(atten, thickness):
    assert compute_cd_thickness(atten) == pytest.approx(thickness)


@pytest.mark.parametrize("atten", [[0, 0, 0, 0], [1, 0, 1, 0], [0, 1, 1, 1], [1, 1, 1, 1]])
def test_flip_inverts_every_stage(atten):
    """With inverted polarity, a stage flagged 0 is the one in the beam."""
    flipped = [1 - state for state in atten]

    assert compute_cd_thickness(atten, flip_atten=True) == pytest.approx(compute_cd_thickness(flipped))


@pytest.mark.parametrize("atten", [[1.0, 0.0, 1.0, 0.0], np.array([1, 0, 1, 0]), (1, 0, 1, 0)])
def test_thickness_accepts_any_flat_numeric_sequence(atten):
    assert compute_cd_thickness(atten) == pytest.approx(MODERATE_THICKNESS)


def test_entries_beyond_the_stages_are_ignored():
    assert compute_cd_thickness([1, 0, 1, 0, 1]) == pytest.approx(MODERATE_THICKNESS)


def test_short_log_raises_sample_logs_error():
    with pytest.raises(SampleLogsError, match="fewer than the 4") as excinfo:
        compute_cd_thickness([1, 0, 1])

    assert type(excinfo.value) is SampleLogsError


@pytest.mark.parametrize("bad_value", [0.5, 2, -1, np.nan])
def test_non_binary_flag_raises(bad_value):
    with pytest.raises(SampleLogsError, match="0 or 1"):
        compute_cd_thickness([1, 0, bad_value, 0])


@pytest.mark.parametrize("atten", [5, "1010", ["a", "b", "c", "d"], [[1, 0], [1, 0]], None])
def test_non_vector_log_raises(atten):
    """A log that is not a flat sequence of numbers stays within the package's exceptions."""
    with pytest.raises(SampleLogsError):
        compute_cd_thickness(atten)


####################
### Transmission ###
####################


def test_transmission_is_a_wavelength_histogram_centred_on_the_table():
    transmission = get_transmission_workspace(MODERATE_THICKNESS)

    assert transmission.getAxis(0).getUnit().unitID() == "Wavelength"
    assert transmission.isHistogramData()
    assert transmission.getNumberHistograms() == 1
    np.testing.assert_allclose(transmission.readX(0), expected_bin_edges(), rtol=1e-12)


def test_transmission_follows_beer_lambert_without_uncertainty():
    _, mu = _load_cd_attenuation_data()
    transmission = get_transmission_workspace(MODERATE_THICKNESS)

    np.testing.assert_allclose(transmission.readY(0), np.exp(-mu * MODERATE_THICKNESS), rtol=1e-12)
    np.testing.assert_array_equal(transmission.readE(0), 0.0)


def test_zero_thickness_transmits_everything():
    np.testing.assert_array_equal(get_transmission_workspace(0.0).readY(0), 1.0)


def test_transmission_is_reused_from_the_analysis_data_service():
    """A second call for the same thickness returns the stored workspace instead of rebuilding it."""
    transmission = get_transmission_workspace(MODERATE_THICKNESS)
    transmission.setY(0, np.full(transmission.blocksize(), -1.0))

    np.testing.assert_array_equal(get_transmission_workspace(MODERATE_THICKNESS).readY(0), -1.0)


def test_nearby_thicknesses_do_not_share_a_transmission():
    first = get_transmission_workspace(0.00574)
    second = get_transmission_workspace(0.00575)

    assert first.name() != second.name()
    assert not np.allclose(first.readY(0), second.readY(0), rtol=1e-6)


@pytest.mark.parametrize("cd_thickness", [-0.001, np.nan, np.inf])
def test_nonsensical_thickness_raises(cd_thickness):
    with pytest.raises(LrValidationError, match="non-negative"):
        get_transmission_workspace(cd_thickness)


##################
### Correction ###
##################


def test_correction_divides_each_event_by_its_transmission(make_events, output_name):
    """Each event's weight is divided by T at its wavelength, and its error by T as well."""
    events = make_events()

    corrected = apply_correction(events, MODERATE_THICKNESS, output_name)

    assert isinstance(corrected, IEventWorkspace)
    assert mtd.doesExist(output_name)
    assert corrected.getNumberEvents() == mtd[events].getNumberEvents()
    for index in range(corrected.getNumberHistograms()):
        event_list = corrected.getSpectrum(index)
        scale = expected_scale(event_list.getTofs(), MODERATE_THICKNESS)
        np.testing.assert_allclose(event_list.getWeights(), scale, rtol=RTOL)
        np.testing.assert_allclose(event_list.getWeightErrors(), scale, rtol=RTOL)


def test_correction_rejects_a_workspace_not_in_wavelength(make_events, output_name):
    """A single-bin TOF run would pass Mantid's own unit check and come back unscaled."""
    tof_events = make_events(wavelength=False)
    assert mtd[tof_events].blocksize() == 1

    with pytest.raises(LrValidationError, match="Wavelength"):
        apply_correction(tof_events, MODERATE_THICKNESS, output_name)


@pytest.mark.parametrize(
    ("tof_min", "tof_max"),
    [
        (1000.0, 56000.0),  # down to ~0.26 Angstrom
        (8000.0, 80000.0),  # up to ~21 Angstrom
    ],
)
def test_correction_rejects_events_beyond_the_table(make_events, output_name, tof_min, tof_max):
    """Mantid would leave events outside the transmission's bins unscaled."""
    events = make_events(tof_min=tof_min, tof_max=tof_max)

    with pytest.raises(LrValidationError, match="beyond the Cd attenuation table"):
        apply_correction(events, MODERATE_THICKNESS, output_name)


def test_correction_rejects_a_transmission_that_would_overflow_float32(make_events, output_name):
    """The full stack at ~15 Angstrom scales squared errors past the float32 range."""
    events = make_events()

    with pytest.raises(LrValidationError, match="float32"):
        apply_correction(events, FULL_STACK_THICKNESS, output_name)


def test_correction_leaves_mismatched_histogram_binning_to_mantid(make_events, output_name):
    histogram = Rebin(make_events(), Params="2.5,0.5,14.5", PreserveEvents=False, OutputWorkspace=output_name)

    with pytest.raises(ValueError, match="X arrays must match"):
        apply_correction(histogram, MODERATE_THICKNESS, output_name)
