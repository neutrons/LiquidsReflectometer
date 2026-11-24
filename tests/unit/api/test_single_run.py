from pathlib import Path

from lr_reduction.api._single_run import SingleRunReduction, reduce_one
from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.types import SingleReductionInput


def _config() -> ReductionConfig:
    return ReductionConfig(
        direct_beams={"db": DirectBeamConfig(direct_beam_run_numbers=[11111])},
        runs={1: ReflectedRunConfig(sequence_number=1, direct_beam="db", run_number=12345)},
    )


class _RecordingSingleRun(SingleRunReduction):
    def load_configuration(self) -> ReductionConfig:
        return _config()

    def load_data(self, _config):
        return None

    def call_operations(self, _data, _config):
        return object()


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


def test_reduce_one_builds_the_composite_from_the_full_config_and_its_name(monkeypatch):
    calls = {}
    composite = object()

    class _RecordingComposition:
        def __init__(self, data, config, direct_beam_name):
            calls["composition"] = (data, config, direct_beam_name)

        def execute(self):
            return composite

    class _RecordingSingleRunOperation:
        def __init__(self, data, config, composite_direct_beam, sequence_number):
            calls["single_run"] = (data, config, composite_direct_beam, sequence_number)

        def execute(self):
            return "result"

    monkeypatch.setattr("lr_reduction.api._single_run.DirectBeamCompositionOperation", _RecordingComposition)
    monkeypatch.setattr("lr_reduction.api._single_run.SingleRunReductionOperation", _RecordingSingleRunOperation)

    config = _config()
    reflected, direct_beams = object(), [object()]
    data = SingleReductionInput(run_data=reflected, direct_beams=direct_beams, direct_beam_name="db")

    assert reduce_one(data, config, 1) == "result"
    assert calls["composition"] == (direct_beams, config, "db")
    assert calls["single_run"] == (reflected, config, composite, 1)
