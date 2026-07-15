from lr_reduction.api.manual import ManualRunSequence, ManualSingleRun, reduce_run, reduce_run_sequence
from lr_reduction.config.model import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.io import RunData
from lr_reduction.types import SequenceResult, SingleRunResult


def test_reduce_run_executes_end_to_end(tmp_path):
    """Proves the walking skeleton runs today, even with call_operations still a placeholder."""
    result = reduce_run(12345, str(tmp_path / "run_12345.yaml"), output_dir=str(tmp_path))
    assert isinstance(result, SingleRunResult)
    assert result.sequence_id == "12345"


def test_reduce_run_uses_manual_single_run(tmp_path):
    result = ManualSingleRun(12345, tmp_path / "run_12345.yaml", output_dir=str(tmp_path)).execute()
    assert isinstance(result, SingleRunResult)


def test_reduce_run_sequence_executes_end_to_end(tmp_path):
    result = reduce_run_sequence(999, str(tmp_path / "seq_999.yaml"), output_dir=str(tmp_path))
    assert isinstance(result, SequenceResult)
    assert result.sequence_id == "999"


def test_reduce_run_sequence_uses_manual_run_sequence(tmp_path):
    result = ManualRunSequence(999, tmp_path / "seq_999.yaml", output_dir=str(tmp_path)).execute()
    assert isinstance(result, SequenceResult)


def test_manual_run_sequence_converts_run_id_to_int(tmp_path, monkeypatch):
    config = ReductionConfig(
        sequence_id="999",
        direct_beams=[DirectBeamConfig(name="db", db_runs=["12345"])],
        reflected_runs=[ReflectedRunConfig(run_id="54321", direct_beam="db")],
    )
    loaded_run_numbers = []

    monkeypatch.setattr("lr_reduction.api.manual.ConfigLoader.load", lambda _self, _path: config)

    def _capture_run_number(_self, run_number):
        loaded_run_numbers.append(run_number)
        return RunData()

    monkeypatch.setattr("lr_reduction.api.manual.RunLoader.load", _capture_run_number)

    ManualRunSequence(999, tmp_path / "seq_999.yaml").execute()

    assert loaded_run_numbers == [54321]
