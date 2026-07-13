from lr_reduction.api.manual import ManualRunSequence, ManualSingleRun, reduce_run, reduce_run_sequence
from lr_reduction.types import SequenceResult, SingleRunResult


def test_reduce_run_executes_end_to_end(tmp_path):
    """Proves the walking skeleton runs today, even with call_operations still a placeholder."""
    result = reduce_run(12345, tmp_path / "run_12345.yaml", output_dir=str(tmp_path))
    assert isinstance(result, SingleRunResult)
    assert result.sequence_id == "12345"


def test_reduce_run_uses_manual_single_run(tmp_path):
    result = ManualSingleRun(12345, tmp_path / "run_12345.yaml", output_dir=str(tmp_path)).execute()
    assert isinstance(result, SingleRunResult)


def test_reduce_run_sequence_executes_end_to_end(tmp_path):
    result = reduce_run_sequence(999, tmp_path / "seq_999.yaml", output_dir=str(tmp_path))
    assert isinstance(result, SequenceResult)
    assert result.sequence_id == "999"


def test_reduce_run_sequence_uses_manual_run_sequence(tmp_path):
    result = ManualRunSequence(999, tmp_path / "seq_999.yaml", output_dir=str(tmp_path)).execute()
    assert isinstance(result, SequenceResult)
