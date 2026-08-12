import pytest

from lr_reduction.api.interfaces import Entrypoint
from lr_reduction.models.config import ReductionConfig


def test_entrypoint_cannot_be_instantiated_directly():
    with pytest.raises(TypeError):
        Entrypoint()


def test_execute_calls_steps_in_order():
    calls = []

    class RecordingEntrypoint(Entrypoint):
        def load_configuration(self) -> ReductionConfig:
            calls.append("load_configuration")
            return "config"

        def load_data(self, config):
            calls.append(("load_data", config))
            return "data"

        def call_operations(self, config, data):
            calls.append(("call_operations", config, data))
            return "result"

        def save_output(self, result) -> None:
            calls.append(("save_output", result))

    result = RecordingEntrypoint().execute()

    assert result == "result"
    assert calls == [
        "load_configuration",
        ("load_data", "config"),
        ("call_operations", "config", "data"),
        ("save_output", "result"),
    ]
