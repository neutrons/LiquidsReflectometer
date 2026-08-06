from pathlib import Path

from lr_reduction.api._single_run import SingleRunReduction
from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig


def _config() -> ReductionConfig:
    return ReductionConfig(
        direct_beams={"db": DirectBeamConfig(db_runs=[11111])},
        runs=[ReflectedRunConfig(sequence_number=1, direct_beam="db", run_number=12345)],
    )


class _RecordingSingleRun(SingleRunReduction):
    def load_configuration(self) -> ReductionConfig:
        return _config()

    def load_data(self, _config):
        return None


def test_output_directory_defaults_to_cwd():
    assert _RecordingSingleRun().output_directory == Path(".")


def test_output_directory_reflects_override(tmp_path):
    assert _RecordingSingleRun(output_dir=str(tmp_path)).output_directory == tmp_path


def test_save_output_writes_orso_then_publishes(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(
        "lr_reduction.api._single_run.write_orso", lambda _result, output_dir: calls.append(("write", output_dir))
    )

    entrypoint = _RecordingSingleRun(output_dir=str(tmp_path))
    monkeypatch.setattr(entrypoint, "publish", lambda result: calls.append(("publish", result)))

    result = entrypoint.execute()

    assert calls == [("write", str(tmp_path)), ("publish", result)]


def test_default_publish_is_a_no_op():
    entrypoint = _RecordingSingleRun()
    assert entrypoint.publish(object()) is None
