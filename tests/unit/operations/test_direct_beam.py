from itertools import groupby

import numpy as np
import pytest
from mantid import simpleapi as api
from mantid.api import IEventWorkspace
from mantid.kernel import DateAndTime, FloatTimeSeriesProperty, config
from mantid.simpleapi import (
    AddSampleLog,
    AddTimeSeriesLog,
    CreateSampleWorkspace,
    DeleteWorkspace,
    SetInstrumentParameter,
    mtd,
)

from lr_reduction.exceptions import ConfigError, OperationInputError, SampleLogsError
from lr_reduction.models.config import ReductionConfig, RunParameters
from lr_reduction.models.run_data import RunData
from lr_reduction.operations.algorithms.convert_tof_to_wavelength import MULTIPLIER_LOG, OFFSET_LOG, T0_LOG, T1_LOG
from lr_reduction.operations.algorithms.trim_direct_beam_tof import CHOPPER2_PHASE_LOG
from lr_reduction.operations.direct_beam import DirectBeamCompositionOperation, cd_order, direct_beam_roi

MENU_LOG = "BL4B:Actuator:Menu"

# Attenuator menu codes (LSB-first stage flags) and the Cd thickness they put in the beam, in cm.
NO_CD = (0, 0.0)  # no stage in the beam
THIN_CD = (1, 57.5e-4)  # stage 1
THICK_CD = (5, (57.5 + 126.5 + 123.0) * 1e-4)  # stages 1 and 3


@pytest.fixture
def make_run():
    """Factory for direct-beam runs carrying an attenuator menu log, deleted afterwards."""
    names = []

    def _make(run_number, menu_code=THIN_CD[0], workspace_type="Event"):
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(
            WorkspaceType=workspace_type, NumBanks=1, BankPixelWidth=1, NumEvents=1, OutputWorkspace=name
        )
        if menu_code is not None:
            AddSampleLog(Workspace=name, LogName=MENU_LOG, LogText=str(menu_code), LogType="Number")
        names.append(name)
        return RunData(workspace=name, run_numbers=(run_number,))

    yield _make
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


# --- cd_order ----------------------------------------------------------------------------


def test_cd_order_sorts_by_ascending_thickness(make_run):
    thick = make_run(1, THICK_CD[0])
    none = make_run(2, NO_CD[0])
    thin = make_run(3, THIN_CD[0])

    ordered = cd_order([thick, none, thin], flip_atten=False)

    assert [run for run, _ in ordered] == [none, thin, thick]
    assert [thickness for _, thickness in ordered] == pytest.approx([NO_CD[1], THIN_CD[1], THICK_CD[1]])


def test_cd_order_keeps_input_order_for_equal_thickness(make_run):
    first = make_run(1)
    second = make_run(2)
    assert [run for run, _ in cd_order([first, second], flip_atten=False)] == [first, second]
    assert [run for run, _ in cd_order([second, first], flip_atten=False)] == [second, first]


def test_cd_order_applies_flip_atten(make_run):
    # Inverted polarity: menu code 15 (every flag set) means no stage is in the beam.
    all_flags = make_run(1, 15)
    no_flags = make_run(2, 0)
    ordered = cd_order([no_flags, all_flags], flip_atten=True)
    assert [run for run, _ in ordered] == [all_flags, no_flags]
    assert ordered[0][1] == 0.0


def test_cd_order_without_attenuator_log_raises(make_run):
    with pytest.raises(SampleLogsError):
        cd_order([make_run(1, menu_code=None)], flip_atten=False)


# --- direct_beam_roi ---------------------------------------------------------------------


def test_direct_beam_roi_pads_peak_and_uses_low_res():
    params = RunParameters(peak={"min": 140, "max": 150}, low_res={"min": 90, "max": 160}, peak_fit={"pad": 3})
    assert direct_beam_roi(params) == ((90, 160), (137, 153))


def test_direct_beam_roi_default_pad_without_peak_fit():
    params = RunParameters(peak={"min": 140, "max": 150}, low_res={"min": 90, "max": 160})
    assert direct_beam_roi(params) == ((90, 160), (139, 151))


@pytest.mark.parametrize(
    ("params", "match"),
    [
        (RunParameters(low_res={"min": 90, "max": 160}), "no `peak` range"),
        (RunParameters(peak={"min": 140, "max": 150}), "no `low_res` range"),
        (
            RunParameters(peak={"min": 1, "max": 10}, low_res={"min": 90, "max": 160}, peak_fit={"pad": 2}),
            "start the box at pixel -1",
        ),
    ],
    ids=["no peak", "no low_res", "pad below pixel 0"],
)
def test_direct_beam_roi_rejects_incomplete_or_out_of_range_box(params, match):
    with pytest.raises(ConfigError, match=match):
        direct_beam_roi(params)


# --- DirectBeamCompositionOperation.validate_input -----------------------------------------


def _config(run_numbers=(1, 2), corrections=None, **direct_beam):
    return ReductionConfig(
        direct_beams={"db_8mm": {"direct_beam_run_numbers": list(run_numbers), **direct_beam}},
        corrections=corrections or {},
    )


def _validate(runs, config=None, name="db_8mm"):
    DirectBeamCompositionOperation(runs, config or _config(), name).validate_input()


def _deleted(run):
    DeleteWorkspace(run.workspace)
    return run


def test_validate_input_accepts_complete_runs(make_run):
    _validate([make_run(1), make_run(2)])


def test_validate_input_rejects_unknown_direct_beam_name(make_run):
    with pytest.raises(OperationInputError, match="No direct beam named 'db_20mm'.*db_8mm"):
        _validate([make_run(1), make_run(2)], name="db_20mm")


def test_validate_input_rejects_empty_data():
    with pytest.raises(OperationInputError, match="no direct-beam runs"):
        _validate([])


@pytest.mark.parametrize(
    "given",
    [(1,), (1, 2, 3), (1, 3), (1, 1, 2)],
    ids=["missing", "extra", "different", "duplicate"],
)
def test_validate_input_rejects_runs_differing_from_configuration(make_run, given):
    with pytest.raises(OperationInputError, match=r"runs \[.*\] were given, but the configuration lists \[1, 2\]"):
        _validate([make_run(run_number) for run_number in given])


@pytest.mark.parametrize(
    ("unusable", "match"),
    [
        (lambda make_run: make_run(2, workspace_type="Histogram"), "run 2: .* is not an event workspace"),
        (lambda make_run: _deleted(make_run(2)), "run 2: no workspace named"),
        (lambda make_run: make_run(2, menu_code=None), "run 2: cannot read the attenuator log: No attenuator log"),
    ],
    ids=["histogram", "deleted", "no attenuator log"],
)
def test_validate_input_reports_every_unusable_run_together(make_run, unusable, match):
    with pytest.raises(OperationInputError, match=match) as raised:
        _validate([make_run(1, menu_code=None), unusable(make_run)])
    assert "run 1: cannot read the attenuator log" in str(raised.value)


@pytest.mark.parametrize(("apply", "warned"), [(True, True), (False, False)])
def test_validate_input_warns_about_missing_error_events_when_correcting_dead_time(
    make_run, monkeypatch, apply, warned
):
    warnings = []
    monkeypatch.setattr("lr_reduction.operations.direct_beam.logger.warning", warnings.append)
    _validate([make_run(1), make_run(2)], _config(corrections={"dead_time": {"apply": apply}}))
    assert bool(warnings) is warned
    if warned:
        assert [message.split(" has ")[0] for message in warnings] == ["run 1", "run 2"]


# --- DirectBeamCompositionOperation.process --------------------------------------------------

# Synthetic direct-beam runs on a 4 x 4 detector, in configuration order. Run 12 has no Cd, run 11
# thin Cd and run 13 thick Cd, so the Cd order is 12, 11, 13. Each run has its own TOF range,
# chopper-2 phase, arm angle and emission-time logs, so every per-run value is traceable.
SYNTHETIC_RUNS = {
    11: {"menu": THIN_CD[0], "tofs": (8000.0, 25000.0), "phase": 9000.0, "tthd": 0.010, "t0": 114.0, "t1": 29.5},
    12: {"menu": NO_CD[0], "tofs": (10000.0, 30000.0), "phase": 11000.0, "tthd": 0.012, "t0": 110.0, "t1": 29.0},
    13: {"menu": THICK_CD[0], "tofs": (15000.0, 38000.0), "phase": 16000.0, "tthd": 0.014, "t0": 118.0, "t1": 30.0},
}
CD_ORDER = (12, 11, 13)
THICKNESS_UM = {11: THIN_CD[1] * 1e4, 12: NO_CD[1] * 1e4, 13: THICK_CD[1] * 1e4}
DETECTOR_SIDE = 4
EVENTS_PER_PIXEL = 100
PULSES = 10000

#: The per-run steps 2-9, in the order the story fixes.
SEGMENT_STEPS = [
    "ExtractDetectorROI",
    "NormaliseByCurrent",
    "ApplyDeadTimeCorrection",
    "TrimDirectBeamTof",
    "ConvertTofToWavelength",
    "CdAttenuationCorrection",
    "SlitScanCorrection",
]


@pytest.fixture
def ads_lists_hidden_workspaces():
    """Makes `mtd.getObjectNames()` include hidden (``__``-prefixed) names."""
    option = "MantidOptions.InvisibleWorkspaces"
    previous = config[option]
    config[option] = "1"
    yield
    config[option] = previous


@pytest.fixture
def synthetic_runs(ads_lists_hidden_workspaces):  # noqa: ARG001
    """Factory for the synthetic direct-beam runs, carrying everything each step reads.

    `overrides` maps a run number to entries replacing those of its `SYNTHETIC_RUNS` spec. Run 11
    also has a rejected-events workspace. Deletes the runs on teardown.
    """
    names = []

    def _make(run_numbers=(11, 12, 13), overrides=None):
        runs = []
        for run_number in run_numbers:
            spec = {**SYNTHETIC_RUNS[run_number], **(overrides or {}).get(run_number, {})}
            name = mtd.unique_hidden_name()
            names.append(name)
            workspace = CreateSampleWorkspace(
                WorkspaceType="Event",
                NumBanks=1,
                BankPixelWidth=DETECTOR_SIDE,
                NumEvents=0,
                XMin=0.0,
                XMax=50000.0,
                BinWidth=50000.0,
                OutputWorkspace=name,
            )
            for index in range(workspace.getNumberHistograms()):
                events = workspace.getSpectrum(index)
                for tof in np.linspace(*spec["tofs"], EVENTS_PER_PIXEL):
                    events.addEventQuickly(float(tof), DateAndTime(0))
            for parameter in ("number-of-x-pixels", "number-of-y-pixels"):
                SetInstrumentParameter(
                    Workspace=name, ParameterName=parameter, ParameterType="Number", Value=str(float(DETECTOR_SIDE))
                )
            charge = FloatTimeSeriesProperty("proton_charge")
            start = DateAndTime("2026-01-01T00:00:00")
            for pulse in range(PULSES):
                charge.addValue(start + pulse * 16_666_667, 1.0e7)
            workspace.getRun().addProperty("proton_charge", charge, True)
            for log, value in (
                ("run_number", run_number),
                ("gd_prtn_chrg", 2.0),
                (MENU_LOG, spec["menu"]),
                (OFFSET_LOG, spec["t0"]),
                (MULTIPLIER_LOG, spec["t1"]),
            ):
                AddSampleLog(Workspace=name, LogName=log, LogText=str(value), LogType="Number")
            for log, value in ((CHOPPER2_PHASE_LOG, spec["phase"]), ("tthd", spec["tthd"])):
                AddTimeSeriesLog(Workspace=name, Name=log, Time="2026-01-01T00:00:00", Value=value)

            error_events = None
            if run_number == 11:
                error_events = mtd.unique_hidden_name()
                names.append(error_events)
                CreateSampleWorkspace(
                    WorkspaceType="Event", NumBanks=1, BankPixelWidth=1, NumEvents=50, OutputWorkspace=error_events
                )
            runs.append(RunData(workspace=name, run_numbers=(run_number,), error_events_workspace=error_events))
        return runs

    yield _make
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


def _composition_config(run_numbers=(11, 12, 13), plots=True, corrections=None, **direct_beam):
    tof_trim = {"chopper_cut": {"coefficients": [1.0, 0.0], "offset_ms": 0.0}, "intensity_threshold": 0.0}
    return ReductionConfig(
        diagnostics={"plots": plots},
        corrections=corrections or {},
        direct_beams={
            "db": {
                "direct_beam_run_numbers": list(run_numbers),
                "peak": {"min": 1, "max": 2},
                "low_res": {"min": 1, "max": 2},
                "tof_trim": tof_trim,
                **direct_beam,
            }
        },
    )


def _events(workspace):
    """Each spectrum's events as (TOF, weight, error) rows, sorted."""
    rows = []
    for index in range(workspace.getNumberHistograms()):
        spectrum = workspace.getSpectrum(index)
        table = np.column_stack([spectrum.getTofs(), spectrum.getWeights(), spectrum.getWeightErrors()])
        rows.append(table[np.lexsort(table.T[::-1])])
    return rows


def _summary(workspace):
    """What the composition must leave unchanged on a loaded run."""
    return (
        workspace.getNumberEvents(),
        workspace.getNumberHistograms(),
        workspace.getAxis(0).getUnit().unitID(),
        sorted(workspace.getRun().keys()),
    )


class _RecordingComposition(DirectBeamCompositionOperation):
    """Stands in for the steps whose algorithms do not exist yet, and records each working copy.

    The instrument settings stay those of the loaded runs, and the beam-peak fit returns
    ``100 + run / 10`` with error ``run / 1000``. Before `cleanup` deletes the working copies, the
    event workspaces among the names it releases, their emission-time logs and events are kept in
    `recorded`, by run.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.recorded = {}

    def _apply_instrument_settings(self, work):
        pass

    def _fit_beam_peak(self, work, peak_fit):  # noqa: ARG002
        run_number = mtd[work].getRunNumber()
        return 100.0 + run_number / 10, run_number / 1000

    def cleanup(self):
        for name in self._created:
            if name not in self._kept and mtd.doesExist(name) and isinstance(mtd[name], IEventWorkspace):
                workspace = mtd[name]
                self.recorded[workspace.getRunNumber()] = {
                    "emission_time": (
                        (workspace.getRun()[T0_LOG].value, workspace.getRun()[T1_LOG].value)
                        if workspace.getRun().hasProperty(T0_LOG)
                        else None
                    ),
                    "events": _events(workspace),
                }
        super().cleanup()


@pytest.fixture
def compose(synthetic_runs):
    """Run the recording composition; returns (result, operation). Deletes what it kept on teardown."""
    results = []

    def _compose(run_numbers=(11, 12, 13), overrides=None, **config):
        runs = synthetic_runs(run_numbers, overrides)
        operation = _RecordingComposition(runs, _composition_config(run_numbers, **config), "db")
        result = operation.execute()
        results.append(result)
        return result, operation

    yield _compose
    for result in results:
        for name in (result.workspace, result.weights_workspace):
            if name is not None and mtd.doesExist(name):
                DeleteWorkspace(name)


def _steps_by_run(result):
    """The composite's algorithm history before the merge, as ``(run, name, properties)`` per step.

    A step belongs to the run whose working copy or dead-time curve it names.
    """
    steps = []
    for history in mtd[result.workspace].getHistory().getAlgorithmHistories()[:-1]:
        properties = {prop.name(): prop.value() for prop in history.getProperties()}
        runs = {run for run in result.run_numbers if f"REF_L_{run}__" in " ".join(properties.values())}
        if runs:
            (run,) = runs
            steps.append((run, history.name(), properties))
    return steps


def _history(result, run_number):
    """The ``(name, properties)`` of each step on `run_number`, in order."""
    return [(name, properties) for run, name, properties in _steps_by_run(result) if run == run_number]


def _created_names(before):
    """Workspaces registered since `before`, leaving out the Cd transmission cache."""
    return {name for name in mtd.getObjectNames() if name not in before and not name.startswith("transmission_")}


def test_composite_is_a_single_bin_event_workspace_in_wavelength(compose):
    result, _ = compose()
    composite = mtd[result.workspace]

    assert isinstance(composite, IEventWorkspace)
    assert composite.getAxis(0).getUnit().unitID() == "Wavelength"
    low = min(segment.wavelength_band[0] for segment in result.segments)
    high = max(segment.wavelength_band[1] for segment in result.segments)
    np.testing.assert_allclose(composite.readX(0), [low, high])
    assert result.workspace.startswith("REF_L_12__") and result.workspace.endswith("_db")
    assert result.run_numbers == CD_ORDER


def test_segment_steps_run_in_cd_order_with_the_fixed_parameters(compose):
    result, _ = compose(flip_atten=False)

    assert mtd[result.workspace].getHistory().getAlgorithmHistories()[-1].name() == "MergeDirectBeamSegments"
    assert [run for run, _ in groupby(run for run, _, _ in _steps_by_run(result))] == list(CD_ORDER)
    for run_number in CD_ORDER:
        history = _history(result, run_number)
        names = [name for name, _ in history]
        assert [name for name in names if name in SEGMENT_STEPS] == SEGMENT_STEPS, run_number
        assert dict(history)["SingleReadoutDeadTimeCorrection"]["UseDeadTimeThreshold"] == "0"
        assert names.index("SingleReadoutDeadTimeCorrection") < names.index("ApplyDeadTimeCorrection")
        assert dict(history)["CdAttenuationCorrection"]["FlipAtten"] == "0"
        assert dict(history)["ConvertTofToWavelength"]["ApplyEmissionTimeCorrection"] == "1"


def test_thickest_cd_run_is_trimmed_with_the_lowest_wavelength_threshold(compose):
    result, _ = compose()
    thresholds = {segment.run_number: segment.dead_time_threshold for segment in result.segments}
    assert thresholds == {12: 1.25, 11: 1.25, 13: 1.5}


def test_trim_without_dead_time_correction_still_uses_the_curve(compose):
    result, _ = compose(corrections={"dead_time": {"apply": False}})
    names = [name for name, _ in _history(result, 12)]
    assert "ApplyDeadTimeCorrection" not in names
    assert names.index("SingleReadoutDeadTimeCorrection") < names.index("TrimDirectBeamTof")


def test_without_trim_the_segments_carry_no_trim_bounds(compose):
    result, _ = compose(tof_trim={"apply": False})
    names = [name for name, _ in _history(result, 12)]
    assert [name for name in names if name in SEGMENT_STEPS] == [
        step for step in SEGMENT_STEPS if step != "TrimDirectBeamTof"
    ]
    for segment in result.segments:
        assert segment.dead_time_threshold is None
        assert segment.tof_range is None
        assert (segment.chopper_cut, segment.dead_time_cut, segment.intensity_cut) == (None, None, None)


def test_loaded_runs_are_unchanged_and_only_the_results_are_left_registered(synthetic_runs):
    runs = synthetic_runs()
    summaries = [_summary(mtd[run.workspace]) for run in runs]
    before = set(mtd.getObjectNames())
    result = _RecordingComposition(runs, _composition_config(), "db").execute()
    try:
        assert [_summary(mtd[run.workspace]) for run in runs] == summaries
        assert _created_names(before) == {result.workspace, result.weights_workspace}
    finally:
        DeleteWorkspace(result.workspace)
        DeleteWorkspace(result.weights_workspace)


def test_metadata_logs_summarize_the_segments(compose):
    result, _ = compose()
    run = mtd[result.workspace].getRun()

    pixels = [100.0 + run_number / 10 for run_number in CD_ORDER]
    tthds = [SYNTHETIC_RUNS[run_number]["tthd"] for run_number in CD_ORDER]
    assert run.getLogData("db_pixel").value == pytest.approx(np.mean(pixels))
    assert run.getLogData("db_pixel_err").value == pytest.approx(np.std(pixels, ddof=1))
    assert run.getLogData("tthd").value == pytest.approx(np.mean(tthds))
    assert run.getLogData("tthd_err").value == pytest.approx(np.std(tthds, ddof=1))
    np.testing.assert_allclose(
        run.getLogData("Cd_attenuator").value, [THICKNESS_UM[run_number] for run_number in CD_ORDER]
    )
    assert run.getLogData("Cd_attenuator").units == "um"
    assert list(run.getLogData("db_runs").value) == list(CD_ORDER)


def test_segments_record_each_runs_values(compose):
    result, _ = compose()

    for segment in result.segments:
        spec = SYNTHETIC_RUNS[segment.run_number]
        assert segment.cd_thickness_um == pytest.approx(THICKNESS_UM[segment.run_number])
        assert segment.slit_scan_factor == 1.0
        assert segment.chopper_cut == spec["phase"]
        assert segment.tof_range[0] == spec["phase"]
        assert segment.dead_time_cut == 0.0
        assert segment.intensity_cut == pytest.approx(spec["tofs"][1])
        assert (segment.db_pixel, segment.db_pixel_err) == (100.0 + segment.run_number / 10, segment.run_number / 1000)
        assert segment.tthd == pytest.approx(spec["tthd"])
        low, high = segment.wavelength_band
        assert 0 < low < high
    assert result.gaps == ()


def test_gaps_follow_the_merge(compose):
    # Run 13 trimmed so the bands cannot meet: run 11 ends at 25000 µs, run 13 starts at 30000 µs.
    result, _ = compose(run_numbers=(11, 13), overrides={13: {"phase": 30000.0}})
    (thin, thick) = result.segments
    assert result.gaps == ((thin.wavelength_band[1], thick.wavelength_band[0]),)


@pytest.mark.parametrize("plots", [True, False])
def test_weights_workspace_only_in_diagnostic_mode(compose, plots):
    result, _ = compose(plots=plots)
    if plots:
        assert result.weights_workspace.endswith("_db_weights") and mtd.doesExist(result.weights_workspace)
    else:
        assert result.weights_workspace is None


@pytest.mark.parametrize(
    "failing_step",
    [
        "SingleReadoutDeadTimeCorrection",
        *SEGMENT_STEPS,
        "_apply_instrument_settings",
        "_fit_beam_peak",
        "MergeDirectBeamSegments",
        "_write_metadata",
    ],
)
def test_a_failure_at_any_step_leaves_nothing_registered(synthetic_runs, monkeypatch, failing_step):
    runs = synthetic_runs()
    operation = _RecordingComposition(runs, _composition_config(), "db")
    if failing_step.startswith("_"):

        def fail(*_args, **_kwargs):
            raise RuntimeError("injected")

        monkeypatch.setattr(operation, failing_step, fail)
    else:
        algorithm = getattr(api, failing_step)
        calls = []

        def algorithm_failing_on_the_last_run(**properties):
            calls.append(failing_step)
            if failing_step == "MergeDirectBeamSegments" or len(calls) == len(runs):
                raise RuntimeError("injected")
            return algorithm(**properties)

        monkeypatch.setattr(api, failing_step, algorithm_failing_on_the_last_run)

    before = set(mtd.getObjectNames())
    with pytest.raises(RuntimeError, match="injected"):
        operation.execute()
    assert _created_names(before) == set()


@pytest.mark.parametrize("step", ["ApplyInstrumentSettings", "FitBeamPeak"])
def test_missing_algorithms_stop_the_composition(synthetic_runs, step):
    runs = synthetic_runs()
    operation = DirectBeamCompositionOperation(runs, _composition_config(), "db")
    if step == "FitBeamPeak":
        operation._apply_instrument_settings = lambda _work: None
    before = set(mtd.getObjectNames())
    with pytest.raises(NotImplementedError, match=step):
        operation.execute()
    assert _created_names(before) == set()


def test_one_run_composite_equals_its_corrected_segment(compose):
    result, operation = compose(run_numbers=(11,))

    assert result.run_numbers == (11,)
    run = mtd[result.workspace].getRun()
    assert (run.getLogData("db_pixel_err").value, run.getLogData("tthd_err").value) == (0.0, 0.0)
    composite = _events(mtd[result.workspace])
    segment = operation.recorded[11]["events"]
    assert len(composite) == len(segment)
    for merged, corrected in zip(composite, segment, strict=True):
        np.testing.assert_array_equal(merged, corrected)


@pytest.mark.parametrize(
    ("emission_time", "expected"),
    [
        ({}, {run: (SYNTHETIC_RUNS[run]["t0"], SYNTHETIC_RUNS[run]["t1"]) for run in CD_ORDER}),
        ({"coefficients": [100.0, 20.0]}, dict.fromkeys(CD_ORDER, (100.0, 20.0))),
        ({"apply": False}, dict.fromkeys(CD_ORDER, (0.0, 0.0))),
    ],
    ids=["own chopper logs", "configured", "disabled"],
)
def test_emission_time_coefficients_of_each_run(compose, emission_time, expected):
    _, operation = compose(corrections={"emission_time": emission_time})
    assert {run: record["emission_time"] for run, record in operation.recorded.items()} == expected


def test_cleanup_keeps_the_results_and_can_be_repeated(compose):
    result, operation = compose()
    operation.cleanup()
    operation.cleanup()
    assert mtd.doesExist(result.workspace)
    assert mtd.doesExist(result.weights_workspace)
