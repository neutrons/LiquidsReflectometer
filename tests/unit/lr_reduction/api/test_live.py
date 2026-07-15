import pytest
from mantid.simpleapi import CreateSampleWorkspace

from lr_reduction.api.live import LiveEntrypoint, reduce_live
from lr_reduction.types import SingleRunResult


@pytest.fixture
def reflected_run():
    ws = CreateSampleWorkspace(WorkspaceType="Event")
    ws.getRun().addProperty("run_number", "54321", True)
    return ws


def test_reduce_live_executes_end_to_end(reflected_run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result = reduce_live(reflected_run)
    assert isinstance(result, SingleRunResult)
    assert result.sequence_id == "54321"


def test_reduce_live_uses_live_entrypoint(reflected_run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result = LiveEntrypoint(reflected_run).execute()
    assert isinstance(result, SingleRunResult)
