import pytest

from lr_reduction.operations.interfaces import OperationInterface


class _Recording(OperationInterface[None, None, str]):
    """Records the order of its calls; fails where told to."""

    def __init__(self, fail_in=None):
        super().__init__(None, None)
        self.fail_in = fail_in
        self.calls = []

    def _step(self, name):
        self.calls.append(name)
        if name == self.fail_in:
            raise RuntimeError(f"{name} failed")

    def validate_input(self):
        self._step("validate_input")

    def process(self):
        self._step("process")
        return "result"

    def cleanup(self):
        self._step("cleanup")


def test_execute_cleans_up_after_processing():
    operation = _Recording()
    assert operation.execute() == "result"
    assert operation.calls == ["validate_input", "process", "cleanup"]


def test_execute_cleans_up_when_processing_fails():
    operation = _Recording(fail_in="process")
    with pytest.raises(RuntimeError, match="process failed"):
        operation.execute()
    assert operation.calls == ["validate_input", "process", "cleanup"]


def test_execute_neither_processes_nor_cleans_up_invalid_input():
    operation = _Recording(fail_in="validate_input")
    with pytest.raises(RuntimeError, match="validate_input failed"):
        operation.execute()
    assert operation.calls == ["validate_input"]
