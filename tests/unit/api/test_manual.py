from lr_reduction.api.manual import ManualRunSequence, ManualSingleRun, reduce_run, reduce_run_sequence
from lr_reduction.io import RunData
from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult


def _config(run_number: int) -> ReductionConfig:
    return ReductionConfig(
        direct_beams={"db": DirectBeamConfig(db_runs=[11111])},
        runs=[ReflectedRunConfig(sequence_number=1, direct_beam="db", run_number=run_number)],
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


def test_reduce_run_sequence_executes_end_to_end(tmp_path, monkeypatch):
    monkeypatch.setattr("lr_reduction.api.manual.ConfigLoader.load", lambda _self, _path: _config(54321))

    result = reduce_run_sequence(999, tmp_path / "seq_999.yaml", output_dir=str(tmp_path))

    assert isinstance(result, CombinedReductionResult)


def test_reduce_run_sequence_uses_manual_run_sequence(tmp_path, monkeypatch):
    monkeypatch.setattr("lr_reduction.api.manual.ConfigLoader.load", lambda _self, _path: _config(54321))

    result = ManualRunSequence(999, tmp_path / "seq_999.yaml", output_dir=str(tmp_path)).execute()

    assert isinstance(result, CombinedReductionResult)


def test_manual_run_sequence_loads_every_configured_run(tmp_path, monkeypatch):
    config = _config(54321)
    loaded_run_numbers = []

    monkeypatch.setattr("lr_reduction.api.manual.ConfigLoader.load", lambda _self, _path: config)

    def _capture_run_number(_self, run_number):
        loaded_run_numbers.append(run_number)
        return RunData()

    monkeypatch.setattr("lr_reduction.api.manual.RunLoader.load", _capture_run_number)

    ManualRunSequence(999, tmp_path / "seq_999.yaml", output_dir=str(tmp_path)).execute()

    assert loaded_run_numbers == [54321]
