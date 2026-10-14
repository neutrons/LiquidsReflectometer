from pathlib import Path

from lr_reduction.api._single_run import reduce_with_direct_beam
from lr_reduction.api.manual import ManualRunSequence, ManualSingleRun, main, reduce_and_combine_runs, reduce_run
from lr_reduction.io import RunLoader
from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult
from lr_reduction.types import ID
from lr_reduction.utils.sample_logs import SampleLogs


def _config(run_number: ID) -> ReductionConfig:
    return ReductionConfig(
        direct_beams={"db": DirectBeamConfig(run_numbers=[11111])},
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

    original_load = RunLoader.load

    def _capture_run_number(_self, run_number):
        loaded_run_numbers.append(run_number)
        return original_load(_self, run_number)

    monkeypatch.setattr("lr_reduction.api.manual.RunLoader.load", _capture_run_number)

    ManualRunSequence([54321], tmp_path / "seq_999.yaml", output_dir=str(tmp_path)).execute()

    assert loaded_run_numbers == [54321, 11111]  # The reflected run and the direct beam run are both loaded


def _shared_direct_beam_config() -> ReductionConfig:
    """Three reflected runs: sequence numbers 1 and 3 share the composite direct beam `db_a`."""
    return ReductionConfig(
        direct_beams={
            "db_a": DirectBeamConfig(run_numbers=[11111, 11112]),
            "db_b": DirectBeamConfig(run_numbers=[22222]),
        },
        runs={
            1: ReflectedRunConfig(sequence_number=1, direct_beam="db_a", run_number=101),
            2: ReflectedRunConfig(sequence_number=2, direct_beam="db_b", run_number=102),
            3: ReflectedRunConfig(sequence_number=3, direct_beam="db_a", run_number=103),
        },
    )


def _load_with_configured_sequence_numbers(monkeypatch, config: ReductionConfig) -> list[ID]:
    """Patch `RunLoader.load` to record each reflected run's configured sequence_number
    (the stub loader records 1 on every run); returns the list of loaded run numbers."""
    sequence_numbers = {run.run_number: run.sequence_number for run in config.runs.values()}
    loaded_run_numbers = []
    original_load = RunLoader.load

    def _load(_self, run_number):
        loaded_run_numbers.append(run_number)
        run = original_load(_self, run_number)
        if run_number in sequence_numbers:
            SampleLogs(run.workspace).insert("sequence_number", sequence_numbers[run_number])
        return run

    monkeypatch.setattr("lr_reduction.api.manual.RunLoader.load", _load)
    return loaded_run_numbers


def test_manual_run_sequence_loads_a_shared_direct_beam_once(tmp_path, monkeypatch):
    """Reflected runs referencing the same composite direct beam share one load of its runs."""
    config = _shared_direct_beam_config()
    loaded_run_numbers = _load_with_configured_sequence_numbers(monkeypatch, config)

    data = ManualRunSequence([101, 102, 103], tmp_path / "seq.yaml").load_data(config)

    assert loaded_run_numbers == [101, 11111, 11112, 102, 22222, 103]
    assert data[0].direct_beams is data[2].direct_beams


def test_manual_run_sequence_composes_a_shared_direct_beam_once(tmp_path, monkeypatch):
    """Reflected runs referencing the same composite direct beam share one composition of it."""
    config = _shared_direct_beam_config()
    _load_with_configured_sequence_numbers(monkeypatch, config)
    composed = []
    monkeypatch.setattr(
        "lr_reduction.api._single_run.DirectBeamCompositionOperation.process",
        lambda self: composed.append(self.config) or f"composite_{len(composed)}",
    )
    reduced_against = []
    original_reduce = reduce_with_direct_beam

    def _capture_composite(run_data, config, comp_db, sequence_number):
        reduced_against.append(comp_db)
        return original_reduce(run_data, config, comp_db, sequence_number)

    monkeypatch.setattr("lr_reduction.api.manual.reduce_with_direct_beam", _capture_composite)
    sequence = ManualRunSequence([101, 102, 103], tmp_path / "seq.yaml")

    sequence.call_operations(sequence.load_data(config), config)

    assert composed == [config.direct_beams["db_a"], config.direct_beams["db_b"]]
    assert reduced_against == ["composite_1", "composite_2", "composite_1"]


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
