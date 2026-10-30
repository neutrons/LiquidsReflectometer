"""Unit tests for the TrimDirectBeamTof Mantid algorithm and its bound functions."""

import numpy as np
import pytest
from mantid.api import EventType, IEventWorkspace
from mantid.kernel import DateAndTime, FloatTimeSeriesProperty
from mantid.simpleapi import (
    AddSampleLog,
    AddTimeSeriesLog,
    ConvertUnits,
    CreateSampleWorkspace,
    CreateWorkspace,
    DeleteWorkspace,
    mtd,
)

from lr_reduction.operations.algorithms.trim_direct_beam_tof import (
    CHOPPER2_PHASE_LOG,
    TrimDirectBeamTof,
    chopper_cut,
    dead_time_cut,
    intensity_cut,
)

OUTPUT_WORKSPACE = "trimmed_direct_beam"
RUN_NUMBER = 4242

# With these chopper constants the frame-overlap cut is the phase itself, in µs, so each test
# sets the chopper cut directly through the phase log.
IDENTITY_CHOPPER = {"ChopperCutCoefficients": [1.0, 0.0], "ChopperCutOffset": 0.0}

# A bright region of five events per 10 µs from 1000 to 9990 µs, then a sparse tail of single
# events every 200 µs up to 20000 µs. On the default 50 µs grid the bright bins hold 25 events
# and the tail bins at most 1, so an IntensityThreshold of 2 cuts the tail at 10000 µs.
BRIGHT_TOFS = np.repeat(np.arange(1000.0, 10000.0, 10.0), 5)
TAIL_TOFS = np.arange(10000.0, 20001.0, 200.0)
DEFAULT_TOFS = np.concatenate([BRIGHT_TOFS, TAIL_TOFS])
FLOOR = 10000.0

# Dead-time-corrected events carry weights: here the bright events from 8000 µs on weigh 0.05, so
# their bins sum to 1.25, below an IntensityThreshold of 2, though each still holds 25 events.
WEIGHTED_FLOOR = 8000.0
DEFAULT_WEIGHTS = np.where((DEFAULT_TOFS >= WEIGHTED_FLOOR) & (DEFAULT_TOFS < FLOOR), 0.05, 1.0)

# The correction curve's grid: 1000 µs bins over [0, 60000] µs.
CURVE_EDGES = np.arange(0.0, 60001.0, 1000.0)


@pytest.fixture
def workspaces():
    """Collects workspace names to delete on teardown, starting with the output workspace."""
    names = [OUTPUT_WORKSPACE]
    yield names
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


@pytest.fixture
def make_events(workspaces):
    """Factory for an event workspace holding `tofs`, dealt round-robin over four spectra, weighted by `weights`."""

    def make(tofs=DEFAULT_TOFS, phase=2000.0, weights=None):
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(
            WorkspaceType="Event",
            NumBanks=1,
            BankPixelWidth=2,
            NumEvents=0,
            XMin=0.0,
            XMax=60000.0,
            BinWidth=60000.0,
            OutputWorkspace=name,
        )
        workspace = mtd[name]
        spectra = [workspace.getSpectrum(index) for index in range(workspace.getNumberHistograms())]
        if weights is None:
            for index, tof in enumerate(tofs):
                spectra[index % len(spectra)].addEventQuickly(tof, DateAndTime(0))
        else:
            for spectrum in spectra:
                spectrum.switchTo(EventType.WEIGHTED)
            for index, (tof, weight) in enumerate(zip(tofs, weights, strict=True)):
                spectra[index % len(spectra)].addWeightedEventQuickly(tof, weight, weight**2, DateAndTime(0))
        AddSampleLog(Workspace=name, LogName="run_number", LogText=str(RUN_NUMBER), LogType="Number")
        if phase is not None:
            AddTimeSeriesLog(Workspace=name, Name=CHOPPER2_PHASE_LOG, Time="2026-01-01T00:00:00", Value=phase)
        workspaces.append(name)
        return name

    return make


@pytest.fixture
def make_curve(workspaces):
    """Factory for a correction curve on `CURVE_EDGES`, 1 everywhere except the `over` bins."""

    def make(over=(), value=2.0, edges=CURVE_EDGES):
        name = mtd.unique_hidden_name()
        factors = np.ones(len(edges) - 1)
        factors[list(over)] = value
        CreateWorkspace(DataX=edges, DataY=factors, NSpec=1, UnitX="TOF", OutputWorkspace=name)
        workspaces.append(name)
        return name

    return make


def make_algorithm(input_workspace, curve, **properties):
    algorithm = TrimDirectBeamTof()
    algorithm.initialize()
    # As through mantid.simpleapi, so a refusal in PyExec raises rather than only being logged.
    algorithm.setRethrows(True)
    algorithm.setProperty("InputWorkspace", input_workspace)
    algorithm.setProperty("CorrectionWorkspace", curve)
    algorithm.setProperty("DeadTimeThreshold", 1.25)
    algorithm.setProperty("IntensityThreshold", 2.0)
    for name, value in {**IDENTITY_CHOPPER, **properties}.items():
        algorithm.setProperty(name, value)
    algorithm.setProperty("OutputWorkspace", OUTPUT_WORKSPACE)
    return algorithm


def run_algorithm(input_workspace, curve, **properties):
    """Execute the algorithm, returning it for its output properties along with the trimmed workspace."""
    algorithm = make_algorithm(input_workspace, curve, **properties)
    algorithm.execute()
    return algorithm, mtd[OUTPUT_WORKSPACE]


def all_tofs(workspace):
    return np.sort(np.concatenate([workspace.getSpectrum(i).getTofs() for i in range(workspace.getNumberHistograms())]))


def all_weights(workspace):
    """The event weights in TOF order; weights here depend only on TOF, so ties need no ordering."""
    spectra = [workspace.getSpectrum(i) for i in range(workspace.getNumberHistograms())]
    tofs = np.concatenate([spectrum.getTofs() for spectrum in spectra])
    weights = np.concatenate([spectrum.getWeights() for spectrum in spectra])
    return weights[np.argsort(tofs, kind="stable")]


# Bound functions


def test_chopper_cut_with_the_prototype_constants():
    # (2.077 * 10000 - 16818) / 1000 + 1.0 = 4.952 ms
    assert chopper_cut(10000.0, (2.077, 16818.0), 1.0) == pytest.approx(4952.0)


EDGES = np.array([0.0, 10.0, 20.0, 30.0, 40.0])


@pytest.mark.parametrize(
    ("factors", "expected"),
    [
        pytest.param([2.0, 1.0, 1.0, 1.0], 10.0, id="over-at-start"),
        pytest.param([1.0, 2.0, 1.1, 1.0], 20.0, id="over-in-middle"),
        pytest.param([2.0, 1.0, 2.0, 1.0], 30.0, id="last-over-bin-wins"),
        pytest.param([1.0, 1.0, 1.0, 1.0], None, id="none-over"),
        pytest.param([1.25, 1.0, 1.0, 1.0], None, id="at-threshold-is-kept"),
    ],
)
def test_dead_time_cut(factors, expected):
    assert dead_time_cut(EDGES, factors, 1.25) == expected


@pytest.mark.parametrize(
    ("counts", "expected"),
    [
        pytest.param([5.0, 5.0, 1.0, 0.0], 20.0, id="low-tail"),
        pytest.param([5.0, 0.0, 5.0, 0.0], 30.0, id="interior-low-bins-are-kept"),
        pytest.param([2.0, 0.0, 0.0, 0.0], 10.0, id="at-threshold-is-kept"),
        pytest.param([1.0, 1.0, 1.0, 1.0], None, id="none-reaching"),
    ],
)
def test_intensity_cut(counts, expected):
    assert intensity_cut(EDGES, counts, 2.0) == expected


def prototype_trim_and_chop(tofs, intensity, dtc, chopper_tof, intensity_threshold, dead_time_threshold):
    """``Direct_Beam._trim_and_chop`` from the ``exp`` branch, ported onto the same index-aligned arrays.

    Returns the TOFs of the surviving bins.
    """
    kept = np.where(intensity >= intensity_threshold)[0]
    tofs, dtc = tofs[kept], dtc[kept]

    kept = np.where(tofs >= chopper_tof)[0]
    tofs, dtc = tofs[kept], dtc[kept]

    above = np.where(dtc > dead_time_threshold)[0]
    if above.size == 0:
        kept = np.where(dtc < dead_time_threshold)[0]
    else:
        kept = np.where(dtc[above[-1] + 1 :] < dead_time_threshold)[0] + above[-1] + 1
    return tofs[kept]


def test_combined_bounds_match_the_prototype_sequential_cuts():
    """The equivalence argument of the design: the three cuts compose into one range."""
    edges = np.arange(1000.0, 30001.0, 50.0)
    bin_starts = edges[:-1]
    # Rate-driven correction falling with TOF, with a spike above threshold past the chopper cut.
    dtc = 1.0 + 2.0 * np.exp(-(bin_starts - 1000.0) / 1000.0)
    dtc[(bin_starts >= 6000.0) & (bin_starts < 6200.0)] = 1.4
    # Intensity rising then falling into a tail below the floor.
    intensity = 100.0 * np.exp(-(((bin_starts - 12000.0) / 6000.0) ** 2))
    phase, chopper_constants, offset = 10000.0, (2.077, 16818.0), 1.0

    prototype = prototype_trim_and_chop(
        bin_starts,
        intensity,
        dtc,
        chopper_tof=((phase * 2.077 - 16818.0) / 1000 + offset) * 1000,
        intensity_threshold=1.0,
        dead_time_threshold=1.25,
    )

    tof_min = max(chopper_cut(phase, chopper_constants, offset), dead_time_cut(edges, dtc, 1.25))
    tof_max = intensity_cut(edges, intensity, 1.0)
    combined = bin_starts[(bin_starts >= tof_min) & (bin_starts < tof_max)]
    np.testing.assert_array_equal(combined, prototype)
    assert tof_min == 6200.0  # the dead time spike, not the chopper cut, is binding


# The algorithm


@pytest.mark.parametrize(
    ("phase", "over", "intensity_threshold", "expected"),
    [
        pytest.param(3000.0, (), 2.0, (3000.0, FLOOR), id="chopper-binding"),
        pytest.param(2000.0, (3,), 2.0, (4000.0, FLOOR), id="dead-time-binding"),
        # A floor of 1 keeps the whole tail: the histogram ends just past its last event, at 20000 µs.
        pytest.param(2000.0, (), 1.0, (2000.0, float(np.nextafter(20000.0, np.inf))), id="floor-binding"),
    ],
)
def test_surviving_range_is_the_combined_bounds(make_events, make_curve, phase, over, intensity_threshold, expected):
    events, curve = make_events(phase=phase), make_curve(over=over)

    algorithm, trimmed = run_algorithm(events, curve, IntensityThreshold=intensity_threshold)

    tof_min, tof_max = expected
    assert (algorithm.getProperty("TOFMin").value, algorithm.getProperty("TOFMax").value) == expected
    inside = DEFAULT_TOFS[(DEFAULT_TOFS >= tof_min) & (DEFAULT_TOFS < tof_max)]
    np.testing.assert_array_equal(all_tofs(trimmed), np.sort(inside))


def test_individual_bounds_are_reported(make_events, make_curve):
    algorithm, _ = run_algorithm(make_events(phase=2000.0), make_curve(over=(3,)))

    assert algorithm.getProperty("ChopperCut").value == 2000.0
    assert algorithm.getProperty("DeadTimeCut").value == 4000.0
    assert algorithm.getProperty("IntensityCut").value == FLOOR


def test_dead_time_cut_is_zero_without_over_threshold_bins(make_events, make_curve):
    algorithm, _ = run_algorithm(make_events(), make_curve())

    assert algorithm.getProperty("DeadTimeCut").value == 0.0
    assert algorithm.getProperty("TOFMin").value == algorithm.getProperty("ChopperCut").value


def test_events_at_tof_min_are_kept_and_at_tof_max_removed(make_events, make_curve):
    # Both bounds fall exactly on event TOFs: 3000 µs in the bright region, 10000 µs in the tail.
    _, trimmed = run_algorithm(make_events(phase=3000.0), make_curve())

    tofs = all_tofs(trimmed)
    assert tofs[0] == 3000.0
    assert tofs[-1] == 9990.0
    assert FLOOR not in tofs


def test_output_is_an_event_workspace_on_the_single_bin_range(make_events, make_curve):
    algorithm, trimmed = run_algorithm(make_events(), make_curve())

    assert isinstance(trimmed, IEventWorkspace)
    expected = [algorithm.getProperty("TOFMin").value, algorithm.getProperty("TOFMax").value]
    for index in range(trimmed.getNumberHistograms()):
        np.testing.assert_array_equal(trimmed.readX(index), expected)


def test_intensity_floor_sums_event_weights(make_events, make_curve):
    algorithm, trimmed = run_algorithm(make_events(weights=DEFAULT_WEIGHTS), make_curve())

    assert algorithm.getProperty("IntensityCut").value == WEIGHTED_FLOOR
    inside = DEFAULT_TOFS[(DEFAULT_TOFS >= 2000.0) & (DEFAULT_TOFS < WEIGHTED_FLOOR)]
    np.testing.assert_array_equal(all_tofs(trimmed), np.sort(inside))


def test_inputs_are_unchanged(make_events, make_curve):
    events, curve = make_events(weights=DEFAULT_WEIGHTS), make_curve(over=(3,))
    tofs_before, weights_before, x_before = all_tofs(mtd[events]), all_weights(mtd[events]), mtd[events].readX(0).copy()
    curve_x_before, curve_y_before = mtd[curve].readX(0).copy(), mtd[curve].readY(0).copy()

    run_algorithm(events, curve)

    np.testing.assert_array_equal(all_tofs(mtd[events]), tofs_before)
    np.testing.assert_array_equal(all_weights(mtd[events]), weights_before)
    np.testing.assert_array_equal(mtd[events].readX(0), x_before)
    np.testing.assert_array_equal(mtd[curve].readX(0), curve_x_before)
    np.testing.assert_array_equal(mtd[curve].readY(0), curve_y_before)


def test_no_intermediate_workspace_is_left_in_the_ads(make_events, make_curve):
    events, curve = make_events(), make_curve()
    before = set(mtd.getObjectNames())

    run_algorithm(events, curve)

    assert set(mtd.getObjectNames()) - before == {OUTPUT_WORKSPACE}


def test_history_records_the_algorithm_and_its_inputs(make_events, make_curve):
    events, curve = make_events(), make_curve()

    _, trimmed = run_algorithm(events, curve)

    history = trimmed.getHistory()
    last = history.getAlgorithmHistory(history.size() - 1)
    assert last.name() == "TrimDirectBeamTof"
    assert last.getPropertyValue("InputWorkspace") == events
    assert last.getPropertyValue("CorrectionWorkspace") == curve


# Validation and errors


def test_input_not_in_tof_is_rejected(make_events, make_curve):
    events = make_events()
    ConvertUnits(InputWorkspace=events, Target="Wavelength", OutputWorkspace=events)

    with pytest.raises(RuntimeError, match=f"Run {RUN_NUMBER}.*TOF"):
        run_algorithm(events, make_curve())


def test_input_without_events_is_rejected(make_events, make_curve):
    with pytest.raises(RuntimeError, match=f"Run {RUN_NUMBER}.*no events"):
        run_algorithm(make_events(tofs=[]), make_curve())


def test_invalid_curve_is_rejected(make_events, make_curve):
    with pytest.raises(RuntimeError, match="dead time correction curve.*negative"):
        run_algorithm(make_events(), make_curve(over=(3,), value=-1.0))


def test_curve_not_covering_the_events_is_rejected(make_events, make_curve):
    with pytest.raises(RuntimeError, match="does not cover"):
        run_algorithm(make_events(), make_curve(edges=np.arange(2000.0, 60001.0, 1000.0)))


def test_missing_chopper_phase_log_is_rejected(make_events, make_curve):
    with pytest.raises(RuntimeError, match=f"frame overlap.*no {CHOPPER2_PHASE_LOG}"):
        run_algorithm(make_events(phase=None), make_curve())


def test_empty_chopper_phase_log_is_rejected(make_events, make_curve):
    events = make_events(phase=None)
    mtd[events].getRun().addProperty(CHOPPER2_PHASE_LOG, FloatTimeSeriesProperty(CHOPPER2_PHASE_LOG), True)

    with pytest.raises(RuntimeError, match="has no entries"):
        run_algorithm(events, make_curve())


@pytest.mark.parametrize(
    ("name", "value"),
    [("DeadTimeThreshold", 0.99), ("IntensityThreshold", -0.1), ("IntensityTOFStep", 0.0)],
)
def test_out_of_bounds_thresholds_are_refused(name, value):
    algorithm = TrimDirectBeamTof()
    algorithm.initialize()

    with pytest.raises(ValueError, match=name):
        algorithm.setProperty(name, value)


@pytest.mark.parametrize("name", ["DeadTimeThreshold", "IntensityThreshold"])
def test_thresholds_are_mandatory(make_events, make_curve, name):
    """The calling operation chooses both thresholds per run; the algorithm supplies no default."""
    algorithm = TrimDirectBeamTof()
    algorithm.initialize()
    algorithm.setProperty("InputWorkspace", make_events())
    algorithm.setProperty("CorrectionWorkspace", make_curve())
    for threshold in {"DeadTimeThreshold", "IntensityThreshold"} - {name}:
        algorithm.setProperty(threshold, 1.25)
    algorithm.setProperty("OutputWorkspace", OUTPUT_WORKSPACE)

    assert algorithm.getProperty(name).isValid != ""
    with pytest.raises(RuntimeError, match="invalid Properties"):
        algorithm.execute()


def test_floor_above_every_bin_is_refused(make_events, make_curve):
    with pytest.raises(RuntimeError, match=r"IntensityThreshold=26.0 \(histogram maximum 25.0\)"):
        run_algorithm(make_events(), make_curve(), IntensityThreshold=26.0)


def test_empty_range_names_all_three_bounds(make_events, make_curve):
    message = "no TOF range survives.*chopper cut 15000.*dead time cut None.*intensity cut 10000"
    with pytest.raises(RuntimeError, match=message):
        run_algorithm(make_events(phase=15000.0), make_curve())


def test_range_without_events_is_refused(make_events, make_curve):
    # The floor ends the range at 1550 µs, after the last bright bin [1500, 1550), whose last events
    # are at 1540 µs, so the chopper cut at 1545 µs leaves a gap without events. The lone event at
    # 9000 µs is below the floor.
    events = make_events(tofs=np.concatenate([BRIGHT_TOFS[BRIGHT_TOFS < 1550.0], [9000.0]]), phase=1545.0)

    with pytest.raises(RuntimeError, match=r"no events remain in the TOF range \[1545.0, 1550.0\)"):
        run_algorithm(events, make_curve())
