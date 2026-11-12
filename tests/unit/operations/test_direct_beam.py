import pytest
from mantid.simpleapi import AddSampleLog, CreateSampleWorkspace, DeleteWorkspace, mtd

from lr_reduction.exceptions import ConfigError, OperationInputError, SampleLogsError
from lr_reduction.models.config import ReductionConfig, RunParameters
from lr_reduction.models.run_data import RunData
from lr_reduction.operations.algorithms.convert_tof_to_wavelength import MULTIPLIER_LOG, OFFSET_LOG
from lr_reduction.operations.algorithms.trim_direct_beam_tof import CHOPPER2_PHASE_LOG
from lr_reduction.operations.direct_beam import DirectBeamCompositionOperation, cd_order, direct_beam_roi
from lr_reduction.properties.cd_attenuation import STAGE_LOGS

MENU_LOG = "BL4B:Actuator:Menu"

# Attenuator menu codes (LSB-first stage flags) and the Cd thickness they put in the beam, in cm.
NO_CD = (0, 0.0)  # no stage in the beam
THIN_CD = (1, 57.5e-4)  # stage 1
THICK_CD = (5, (57.5 + 126.5 + 123.0) * 1e-4)  # stages 1 and 3


@pytest.fixture
def make_run():
    """Factory for direct-beam runs carrying an attenuator menu log, deleted afterwards."""
    names = []

    def _make(run_number, menu_code=None, logs=(), workspace_type="Event", unit="TOF"):
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(
            WorkspaceType=workspace_type, XUnit=unit, NumBanks=1, BankPixelWidth=1, NumEvents=1, OutputWorkspace=name
        )
        if menu_code is not None:
            AddSampleLog(Workspace=name, LogName=MENU_LOG, LogText=str(menu_code), LogType="Number")
        for log in logs:
            AddSampleLog(Workspace=name, LogName=log, LogText="1.0", LogType="Number")
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
    first = make_run(1, THIN_CD[0])
    second = make_run(2, THIN_CD[0])
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
        cd_order([make_run(1)], flip_atten=False)


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


# --- DirectBeamCompositionOperation --------------------------------------------------------


def test_operation_keeps_full_config_and_direct_beam_name(make_run):
    config = ReductionConfig(direct_beams={"db_8mm": {"direct_beam_run_numbers": [1]}})
    runs = [make_run(1, THIN_CD[0])]

    operation = DirectBeamCompositionOperation(runs, config, "db_8mm")

    assert operation.data == runs
    assert operation.config is config
    assert operation.direct_beam_name == "db_8mm"


# --- DirectBeamCompositionOperation.validate_input -----------------------------------------

#: The logs the default configuration's steps read beyond the attenuator log: the trim's
#: chopper-2 phase and the emission-time chopper offset and multiplier.
STEP_LOGS = (CHOPPER2_PHASE_LOG, OFFSET_LOG, MULTIPLIER_LOG)


def _config(run_numbers=(1, 2), corrections=None, **direct_beam):
    return ReductionConfig(
        direct_beams={"db_8mm": {"direct_beam_run_numbers": list(run_numbers), **direct_beam}},
        corrections=corrections or {},
    )


@pytest.fixture
def complete_run(make_run):
    """Factory for a direct-beam run carrying every log the default configuration reads."""

    def _make(run_number, **kwargs):
        kwargs.setdefault("menu_code", THIN_CD[0])
        kwargs.setdefault("logs", STEP_LOGS)
        return make_run(run_number, **kwargs)

    return _make


def _validate(runs, config=None, name="db_8mm"):
    DirectBeamCompositionOperation(runs, config or _config(), name).validate_input()


def test_validate_input_accepts_complete_runs(complete_run):
    _validate([complete_run(1), complete_run(2)])


def test_validate_input_rejects_unknown_direct_beam_name(complete_run):
    with pytest.raises(OperationInputError, match="No direct beam named 'db_20mm'.*db_8mm"):
        _validate([complete_run(1), complete_run(2)], name="db_20mm")


def test_validate_input_rejects_empty_data():
    with pytest.raises(OperationInputError, match="no direct-beam runs"):
        _validate([])


@pytest.mark.parametrize(
    "given",
    [(1,), (1, 2, 3), (1, 3), (1, 1, 2)],
    ids=["missing", "extra", "different", "duplicate"],
)
def test_validate_input_rejects_runs_differing_from_configuration(complete_run, given):
    with pytest.raises(OperationInputError, match=r"runs \[.*\] were given, but the configuration lists \[1, 2\]"):
        _validate([complete_run(run_number) for run_number in given])


def test_validate_input_rejects_workspace_in_wavelength(complete_run):
    with pytest.raises(OperationInputError, match="run 2: .* is in Wavelength, not TOF"):
        _validate([complete_run(1), complete_run(2, unit="Wavelength")])


def test_validate_input_rejects_histogram_workspace(complete_run):
    with pytest.raises(OperationInputError, match="run 2: .* is not an event workspace"):
        _validate([complete_run(1), complete_run(2, workspace_type="Histogram")])


def test_validate_input_rejects_deleted_workspace(complete_run):
    deleted = complete_run(2)
    DeleteWorkspace(deleted.workspace)
    with pytest.raises(OperationInputError, match="run 2: no workspace named"):
        _validate([complete_run(1), deleted])


def test_validate_input_rejects_run_without_attenuator_log(complete_run):
    with pytest.raises(OperationInputError, match="run 2: no attenuator log"):
        _validate([complete_run(1), complete_run(2, menu_code=None)])


def test_validate_input_accepts_per_stage_attenuator_logs(complete_run):
    _validate([complete_run(1), complete_run(2, menu_code=None, logs=STEP_LOGS + tuple(STAGE_LOGS))])


def test_validate_input_requires_chopper2_phase_log_only_when_trimming(complete_run):
    without_phase = (OFFSET_LOG, MULTIPLIER_LOG)
    with pytest.raises(OperationInputError, match=f"run 2: no {CHOPPER2_PHASE_LOG} log, which the TOF trim"):
        _validate([complete_run(1), complete_run(2, logs=without_phase)])
    _validate(
        [complete_run(1), complete_run(2, logs=without_phase)],
        _config(tof_trim={"apply": False}),
    )


@pytest.mark.parametrize("missing", [OFFSET_LOG, MULTIPLIER_LOG])
def test_validate_input_requires_chopper_emission_logs_without_coefficients(complete_run, missing):
    logs = tuple(log for log in STEP_LOGS if log != missing)
    with pytest.raises(OperationInputError, match=f"run 2: no {missing} log, which the emission-time correction"):
        _validate([complete_run(1), complete_run(2, logs=logs)])


@pytest.mark.parametrize(
    "emission_time",
    [{"coefficients": [100.0, 20.0]}, {"apply": False}],
    ids=["coefficients configured", "correction disabled"],
)
def test_validate_input_does_not_need_chopper_emission_logs(complete_run, emission_time):
    runs = [complete_run(1), complete_run(2, logs=(CHOPPER2_PHASE_LOG,))]
    _validate(runs, _config(corrections={"emission_time": emission_time}))


def test_validate_input_reports_every_issue_together(complete_run):
    with pytest.raises(OperationInputError) as raised:
        _validate([complete_run(1, unit="Wavelength"), complete_run(2, menu_code=None)])
    message = str(raised.value)
    assert "run 1: " in message and "not TOF" in message
    assert "run 2: no attenuator log" in message


@pytest.mark.parametrize(("apply", "warned"), [(True, True), (False, False)])
def test_validate_input_warns_about_missing_error_events_when_correcting_dead_time(
    complete_run, monkeypatch, apply, warned
):
    warnings = []
    monkeypatch.setattr("lr_reduction.operations.direct_beam.logger.warning", warnings.append)
    _validate([complete_run(1), complete_run(2)], _config(corrections={"dead_time": {"apply": apply}}))
    assert bool(warnings) is warned
    if warned:
        assert [message.split(" has ")[0] for message in warnings] == ["run 1", "run 2"]
