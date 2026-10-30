"""Unit tests for the ConvertTofToWavelength Mantid algorithm."""

import os

import numpy as np
import pytest
from mantid.api import AlgorithmManager, IEventWorkspace
from mantid.kernel import DateAndTime, FloatTimeSeriesProperty, PhysicalConstants
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

from lr_reduction.operations.algorithms.convert_tof_to_wavelength import (
    FORMULA_PARAMETER,
    MULTIPLIER_LOG,
    OFFSET_LOG,
    T0_LOG,
    T1_LOG,
    ConvertTofToWavelength,
    emission_time_formula,
)
from lr_reduction.utils.sample_logs import SampleLogs

# Production chopper settings, as run 234266 records them in its ChopperOffset and ChopperMultiplier logs.
T0 = 114.0  # us
T1 = 29.5  # us / A

TOF_RANGE = (5000.0, 60000.0)  # us, about 1.4-17 Angstrom on REF_L
EVENTS_PER_PIXEL = 200
PULSE_TIME = DateAndTime("2026-01-01T00:00:00")

# The algorithm's promise against the exact formula.
EXACT_RTOL = 1e-8
# The legacy two-pass formula is not exact. Its relative gap is
#   t1 * (t0 + t1 * TOF / A) / (A * (TOF - t0)),   A = m_n * L / h ~ 3460 us/A on REF_L,
# about 2.4e-4 at the shortest TOF for the production coefficients.
TWO_PASS_RTOL = 5e-4


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


def populated(workspace):
    """Indices of the spectra holding events; most pixels of a real run hold none."""
    indices = [i for i in range(workspace.getNumberHistograms()) if workspace.getSpectrum(i).getNumberEvents()]
    assert indices, "comparison would pass vacuously on an empty workspace"
    return indices


def assert_matches(output, reference_input, oracle, rtol):
    """Every output event's wavelength matches `oracle(tof, L)` of the same event in `reference_input`."""
    for index in populated(reference_input):
        expected = oracle(reference_input.getSpectrum(index).getTofs(), flight_path(reference_input, index))
        np.testing.assert_allclose(output.getSpectrum(index).getTofs(), expected, rtol=rtol, atol=0)


def assert_plain_conversion(output, reference_input, workspaces):
    """`output` is exactly a plain ConvertUnits of `reference_input`, with no emission correction."""
    plain_name = mtd.unique_hidden_name()
    workspaces.append(plain_name)
    plain = ConvertUnits(
        InputWorkspace=reference_input, Target="Wavelength", EMode="Elastic", OutputWorkspace=plain_name
    )
    for index in populated(reference_input):
        np.testing.assert_array_equal(output.getSpectrum(index).getTofs(), plain.getSpectrum(index).getTofs())
    assert not output.getInstrument().hasParameter(FORMULA_PARAMETER)


def exact_for(t0, t1):
    """The closed form for coefficients other than the module defaults, as an `assert_matches` oracle."""
    return lambda tof, path: closed_form(tof, path, t0=t0, t1=t1)


def run(input_workspace, **properties):
    """Run the algorithm directly and return the output workspace."""
    algorithm = ConvertTofToWavelength()
    algorithm.initialize()
    algorithm.setChild(True)
    algorithm.setProperty("InputWorkspace", input_workspace)
    for name, value in properties.items():
        algorithm.setProperty(name, value)
    algorithm.setPropertyValue("OutputWorkspace", mtd.unique_hidden_name())
    algorithm.execute()
    return algorithm.getProperty("OutputWorkspace").value


def add_series_log(workspace, name, values):
    """Record `values` as a time series, the shape the DAS writes chopper PVs in."""
    series = FloatTimeSeriesProperty(name)
    for seconds, value in enumerate(values):
        series.addValue(PULSE_TIME + int(seconds * 1e9), value)
    workspace.mutableRun().addProperty(name, series, True)


@pytest.fixture
def workspaces():
    """Tracks workspaces created by a test and deletes them on teardown."""
    names = []
    yield names
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


@pytest.fixture(scope="module")
def real_run(nexus_dir):
    """A real REF_L run whose DAS recorded the chopper coefficient logs. Not to be modified."""
    name = mtd.unique_hidden_name()
    workspace = LoadEventNexus(Filename=os.path.join(nexus_dir, "REF_L_234266.nxs.h5"), OutputWorkspace=name)
    logs = SampleLogs(workspace)
    assert (logs[OFFSET_LOG], logs[MULTIPLIER_LOG]) == (T0, T1)
    yield mtd[name]
    DeleteWorkspace(name)


@pytest.fixture
def real_copy(real_run, workspaces):
    """A copy of `real_run` for tests that modify their input."""
    name = mtd.unique_hidden_name()
    workspaces.append(name)
    return CloneWorkspace(InputWorkspace=real_run, OutputWorkspace=name)


@pytest.fixture
def event_workspace(workspaces):
    """Factory for a small four-pixel event workspace in TOF carrying the real REF_L instrument,
    for conditions a real run cannot provide.

    `logs` maps chopper log names to values; by default there are none.
    """

    def make(logs=None):
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
        for index in range(workspace.getNumberHistograms()):
            events = workspace.getSpectrum(index)
            for tof in np.linspace(*TOF_RANGE, EVENTS_PER_PIXEL):
                events.addEventQuickly(float(tof), PULSE_TIME)
        # Keep the sample workspace's spectrum map: the four pixels map to REF_L detector IDs,
        # so each spectrum gets a real pixel position and flight path.
        LoadInstrument(Workspace=name, InstrumentName="REF_L", RewriteSpectraMap=False)
        for log_name, value in (logs or {}).items():
            add_series_log(workspace, log_name, [value])
        return mtd[name]

    return make


@pytest.mark.datarepo
def test_coefficients_default_to_chopper_logs(real_run):
    output = run(real_run)

    assert_matches(output, real_run, closed_form, EXACT_RTOL)
    logs = SampleLogs(output)
    assert (logs[T0_LOG], logs[T1_LOG]) == (T0, T1)
    assert output.getNumberEvents() == real_run.getNumberEvents()
    assert logs["gd_prtn_chrg"] == SampleLogs(real_run)["gd_prtn_chrg"]
    np.testing.assert_array_equal(output.getRun()["proton_charge"].value, real_run.getRun()["proton_charge"].value)


@pytest.mark.datarepo
def test_given_coefficients_override_logs(real_run):
    output = run(real_run, T0=50.0, T1=10.0)
    assert_matches(output, real_run, exact_for(50.0, 10.0), EXACT_RTOL)


@pytest.mark.datarepo
def test_each_coefficient_falls_back_to_its_own_log(real_run):
    output = run(real_run, T0=50.0)
    assert_matches(output, real_run, exact_for(50.0, T1), EXACT_RTOL)


def test_given_coefficients_need_no_logs(event_workspace):
    workspace = event_workspace()
    output = run(workspace, T0=T0, T1=T1)
    assert_matches(output, workspace, closed_form, EXACT_RTOL)


@pytest.mark.datarepo
def test_zero_coefficients_turn_correction_off(real_run, workspaces):
    # The logs are present, but explicit zeros override them and leave nothing to correct.
    output = run(real_run, T0=0.0, T1=0.0)
    assert_plain_conversion(output, real_run, workspaces)


@pytest.mark.datarepo
def test_flag_turns_correction_off_ignoring_logs(real_run, workspaces):
    output = run(real_run, ApplyEmissionTimeCorrection=False)

    assert_plain_conversion(output, real_run, workspaces)
    logs = SampleLogs(output)
    assert (logs[T0_LOG], logs[T1_LOG]) == (0.0, 0.0)


def test_flag_turns_correction_off_without_logs(event_workspace, workspaces):
    workspace = event_workspace()
    output = run(workspace, ApplyEmissionTimeCorrection=False)
    assert_plain_conversion(output, workspace, workspaces)


@pytest.mark.datarepo
def test_large_accepted_slope_matches_exact_formula(real_run):
    # Ten times the production slope: each iteration gains only ~12x, still within the cap.
    output = run(real_run, T1=300.0)
    assert_matches(output, real_run, exact_for(T0, 300.0), EXACT_RTOL)


@pytest.mark.datarepo
def test_rejects_slope_too_large_to_converge(real_run):
    # At this slope ten iterations would leave ~0.3% error.
    with pytest.raises(RuntimeError, match=r"T1: 2000.0 us/A is too large .*must not exceed"):
        run(real_run, T1=2000.0)


@pytest.mark.datarepo
@pytest.mark.parametrize("niter", [1, 2], ids=["mantid-default", "two-iterations"])
def test_looser_iteration_settings_miss_exact_formula(real_run, real_copy, niter):
    # Pins the iteration settings: at the production slope each iteration gains only ~120x, so
    # Mantid's defaults (and the two iterations once thought enough) fall short of 1e-8.
    copy = real_copy.name()
    SetInstrumentParameter(
        Workspace=copy, ParameterName=FORMULA_PARAMETER, ParameterType="String", Value=emission_time_formula(T0, T1)
    )
    ModeratorTzero(InputWorkspace=copy, EMode="Elastic", Niter=niter, OutputWorkspace=copy)
    loose = ConvertUnits(InputWorkspace=copy, Target="Wavelength", EMode="Elastic", OutputWorkspace=copy)
    with pytest.raises(AssertionError):
        assert_matches(loose, real_run, closed_form, EXACT_RTOL)


@pytest.mark.datarepo
def test_input_is_unchanged(real_run):
    indices = populated(real_run)
    tofs_before = [real_run.getSpectrum(i).getTofs().copy() for i in indices]

    run(real_run)

    assert real_run.getAxis(0).getUnit().unitID() == "TOF"
    assert not real_run.getInstrument().hasParameter(FORMULA_PARAMETER)
    for index, before in zip(indices, tofs_before, strict=True):
        np.testing.assert_array_equal(real_run.getSpectrum(index).getTofs(), before)


@pytest.mark.datarepo
def test_output_is_event_workspace_in_wavelength_carrying_formula(real_run):
    output = run(real_run)

    assert isinstance(output, IEventWorkspace)
    assert output.getAxis(0).getUnit().unitID() == "Wavelength"
    assert output.getNumberEvents() == real_run.getNumberEvents()
    assert list(output.getInstrument().getStringParameter(FORMULA_PARAMETER)) == [emission_time_formula(T0, T1)]


@pytest.mark.datarepo
def test_geometry_reaches_conversion_through_instrument(real_copy):
    # Moving the moderator changes L1, and with it every wavelength, exactly as the formula
    # predicts for the new flight path: the algorithm takes geometry from nowhere else.
    l1_before = real_copy.spectrumInfo().l1()
    MoveInstrumentComponent(Workspace=real_copy, ComponentName="moderator", Z=-15.0, RelativePosition=False)
    assert real_copy.spectrumInfo().l1() != pytest.approx(l1_before)

    output = run(real_copy)
    assert_matches(output, real_copy, closed_form, EXACT_RTOL)


@pytest.mark.datarepo
def test_agrees_with_two_pass_calculation(real_run):
    output = run(real_run)
    assert_matches(output, real_run, two_pass, TWO_PASS_RTOL)
    # The gap is real: the two-pass formula is not a stand-in for the exact one.
    with pytest.raises(AssertionError):
        assert_matches(output, real_run, two_pass, 1e-5)


@pytest.mark.datarepo
def test_formula_round_trips_negative_and_small_coefficients(real_run):
    # The formula is text: a negative slope and an exponent-notation offset must survive it.
    output = run(real_run, T0=1e-5, T1=-1.5)
    assert_matches(output, real_run, exact_for(1e-5, -1.5), EXACT_RTOL)


@pytest.mark.datarepo
def test_history_records_algorithm_and_coefficients(real_run, workspaces):
    # Run through the factory as a non-child, the way the operations invoke registered
    # algorithms, so the call is recorded in the output's history.
    output_name = mtd.unique_hidden_name()
    workspaces.append(output_name)
    algorithm = AlgorithmManager.create("ConvertTofToWavelength")
    algorithm.setProperty("InputWorkspace", real_run.name())
    algorithm.setProperty("T0", 50.0)
    algorithm.setProperty("T1", 10.0)
    algorithm.setProperty("OutputWorkspace", output_name)
    algorithm.execute()

    last = mtd[output_name].getHistory().lastAlgorithm()
    assert last.name() == "ConvertTofToWavelength"
    assert float(last.getPropertyValue("T0")) == 50.0
    assert float(last.getPropertyValue("T1")) == 10.0


@pytest.mark.datarepo
def test_rejects_input_not_in_tof(real_run, workspaces):
    name = mtd.unique_hidden_name()
    workspaces.append(name)
    in_wavelength = ConvertUnits(InputWorkspace=real_run, Target="Wavelength", OutputWorkspace=name)
    with pytest.raises(RuntimeError, match="InputWorkspace: X axis must be in TOF"):
        run(in_wavelength)


def test_rejects_workspace_without_source_or_sample(workspaces):
    histogram_name, event_name = mtd.unique_hidden_name(), mtd.unique_hidden_name()
    workspaces.extend([histogram_name, event_name])
    CreateWorkspace(DataX=[1000.0, 2000.0, 3000.0], DataY=[1.0, 1.0], UnitX="TOF", OutputWorkspace=histogram_name)
    no_instrument = ConvertToEventWorkspace(InputWorkspace=histogram_name, OutputWorkspace=event_name)
    with pytest.raises(RuntimeError, match="InputWorkspace: instrument has no source or sample position"):
        run(no_instrument, T0=T0, T1=T1)


@pytest.mark.parametrize(("name", "log_name"), [("T0", OFFSET_LOG), ("T1", MULTIPLIER_LOG)])
def test_rejects_missing_log_when_coefficient_not_given(event_workspace, name, log_name):
    with pytest.raises(RuntimeError, match=f"{name}: not set, and the run's log cannot supply it: .*{log_name}"):
        run(event_workspace())


def test_rejects_log_that_varies_across_run(event_workspace):
    workspace = event_workspace(logs={OFFSET_LOG: T0})
    add_series_log(workspace, MULTIPLIER_LOG, [T1, T1 + 1.0])
    with pytest.raises(RuntimeError, match="T1: not set, and the run's log cannot supply it: .*varies across the run"):
        run(workspace)


def test_rejects_non_finite_log(event_workspace):
    workspace = event_workspace(logs={OFFSET_LOG: np.inf, MULTIPLIER_LOG: T1})
    with pytest.raises(RuntimeError, match=f"T0: log '{OFFSET_LOG}' is not a finite number"):
        run(workspace)


@pytest.mark.datarepo
@pytest.mark.parametrize("name", ["T0", "T1"])
@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_rejects_non_finite_given_coefficient(real_run, name, value):
    with pytest.raises(RuntimeError, match=f"{name}: the given value is not a finite number"):
        run(real_run, **{name: value})
