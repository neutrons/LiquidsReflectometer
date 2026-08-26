from pathlib import Path

from lr_reduction.api.manual import ManualRunSequence, ManualSingleRun, main, reduce_and_combine_runs, reduce_run
from lr_reduction.io import RunData
from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult


def _config(run_number: int) -> ReductionConfig:
    return ReductionConfig(
        direct_beams={"db": DirectBeamConfig(db_run_numbers=[11111])},
        runs={1: ReflectedRunConfig(sequence_number=1, direct_beam="db", run_number=run_number)},
    )


def test_reduce_run_executes_end_to_end(tmp_path, monkeypatch):
    """Proves the walking skeleton runs today, even with call_operations still a placeholder."""
    monkeypatch.setattr("lr_reduction.api.manual.ConfigLoader.load", lambda _self, _path: _config(12345))

    result = reduce_run(12345, tmp_path / "run_12345.yaml", output_dir=str(tmp_path))

    assert isinstance(result, ReductionResult)


def test_reduce_run_uses_manual_single_run(tmp_path, monkeypatch):
    monkeypatch.setattr("lr_reduction.api.manual.ConfigLoader.load", lambda _self, _path: _config(12345))

    result = ManualSingleRun(12345, tmp_path / "run_12345.yaml", output_dir=str(tmp_path)).execute()

    assert isinstance(result, ReductionResult)


def test_reduce_and_combine_runs_executes_end_to_end(tmp_path, monkeypatch):
    monkeypatch.setattr("lr_reduction.api.manual.ConfigLoader.load", lambda _self, _path: _config(54321))

    result = reduce_and_combine_runs([54321], tmp_path / "seq_999.yaml", output_dir=str(tmp_path))

    assert isinstance(result, CombinedReductionResult)


def test_reduce_and_combine_runs_uses_manual_run_sequence(tmp_path, monkeypatch):
    monkeypatch.setattr("lr_reduction.api.manual.ConfigLoader.load", lambda _self, _path: _config(54321))

    result = ManualRunSequence([54321], tmp_path / "seq_999.yaml", output_dir=str(tmp_path)).execute()

    assert isinstance(result, CombinedReductionResult)


def test_manual_run_sequence_loads_every_configured_run(tmp_path, monkeypatch):
    config = _config(54321)
    loaded_run_numbers = []

    monkeypatch.setattr("lr_reduction.api.manual.ConfigLoader.load", lambda _self, _path: config)

    def _capture_run_number(_self, run_number):
        loaded_run_numbers.append(run_number)
        return RunData()

    monkeypatch.setattr("lr_reduction.api.manual.RunLoader.load", _capture_run_number)

    ManualRunSequence([54321], tmp_path / "seq_999.yaml", output_dir=str(tmp_path)).execute()

    assert loaded_run_numbers == [54321, 11111]  # The reflected run and the direct beam run are both loaded


def test_main_run_subcommand_parses_and_dispatches(monkeypatch):
    captured = {}
    monkeypatch.setattr(
        "lr_reduction.api.manual.reduce_run",
        lambda *args, **kwargs: captured.update(args=args, kwargs=kwargs),
    )

    main(["run", "12345", "--configuration", "config.yaml"])

    assert captured["args"] == (12345, Path("config.yaml"))


def test_main_sequence_subcommand_parses_and_dispatches(monkeypatch):
    captured = {}
    monkeypatch.setattr(
        "lr_reduction.api.manual.reduce_and_combine_runs",
        lambda *args, **kwargs: captured.update(args=args, kwargs=kwargs),
    )

    main(
        [
            "sequence",
            "--run-numbers",
            "111",
            "112",
            "--configuration",
            "config.yaml",
            "--sequence-numbers",
            "1",
            "2",
        ]
    )

    assert captured["args"] == ([111, 112], Path("config.yaml"))
    assert captured["kwargs"] == {"sequence_numbers": [1, 2]}


def test_main_sequence_subcommand_sequence_numbers_is_optional(monkeypatch):
    captured = {}
    monkeypatch.setattr(
        "lr_reduction.api.manual.reduce_and_combine_runs",
        lambda *args, **kwargs: captured.update(args=args, kwargs=kwargs),
    )

    main(["sequence", "--run-numbers", "111", "112", "--configuration", "config.yaml"])

    assert captured["args"] == ([111, 112], Path("config.yaml"))
    assert captured["kwargs"] == {"sequence_numbers": None}
