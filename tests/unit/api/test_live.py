import pytest
from mantid.simpleapi import CreateSampleWorkspace

from lr_reduction.api.live import LiveEntrypoint, reduce_live
from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult


@pytest.fixture
def reflected_run():
    ws = CreateSampleWorkspace(WorkspaceType="Event")
    ws.getRun().addProperty("run_number", "54321", True)
    return ws


@pytest.fixture(autouse=True)
def _config_loader(monkeypatch):
    config = ReductionConfig(
        direct_beams={"db": DirectBeamConfig(direct_beam_run_numbers=[11111])},
        runs={1: ReflectedRunConfig(sequence_number=1, direct_beam="db", run_number=54321)},
    )
    monkeypatch.setattr("lr_reduction.io.config_loader.ConfigLoader.load", lambda _self, _path: config)


def test_reduce_live_executes_end_to_end(reflected_run, tmp_path, monkeypatch):
    """Mirrors autoreduction (§6.4.4): reduces the run, then assembles from on-disk partials."""
    monkeypatch.chdir(tmp_path)

    result = reduce_live(reflected_run)

    assert isinstance(result, CombinedReductionResult)


def test_reduce_live_uses_live_entrypoint(reflected_run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    result = LiveEntrypoint(reflected_run).execute()

    assert isinstance(result, ReductionResult)
