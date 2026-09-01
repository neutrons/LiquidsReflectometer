import pytest
from mantid.kernel import Int32TimeSeriesProperty
from mantid.simpleapi import CreateSampleWorkspace

from lr_reduction.api.live import LiveEntrypoint, reduce_live
from lr_reduction.exceptions import LogNotFoundError
from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult


def _live_workspace(sequence_number: int | None = 1):
    """A live reflected-run workspace, with `sequence_number` recorded as the DAS records
    it: a time series that is constant across the run (`legacy/workflow.py` reads it as
    `.value[0]`). Pass None to omit the log entirely."""
    ws = CreateSampleWorkspace(WorkspaceType="Event")
    ws.getRun().addProperty("run_number", "54321", True)
    if sequence_number is not None:
        recorded = Int32TimeSeriesProperty("sequence_number")
        recorded.addValue("2020-01-01T00:00:00", sequence_number)
        recorded.addValue("2020-01-01T00:00:01", sequence_number)
        ws.getRun().addProperty("sequence_number", recorded, True)
    return ws


@pytest.fixture
def reflected_run():
    return _live_workspace()


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


def test_load_data_reads_sequence_number_from_the_sample_logs(reflected_run, tmp_path, monkeypatch):
    """The logs are the single source of truth for the metadata-derived value (§3.1.2),
    and it must arrive as a plain int: numpy scalars leak into ReductionResult and ORSO."""
    monkeypatch.chdir(tmp_path)
    entrypoint = LiveEntrypoint(reflected_run)

    entrypoint.execute()

    assert entrypoint.sequence_number == 1
    assert isinstance(entrypoint.sequence_number, int)


def test_load_data_requires_a_recorded_sequence_number(tmp_path, monkeypatch):
    """No silent default: a run whose sequence_number was never recorded would otherwise
    be reduced against whichever configuration entry the default happened to name."""
    monkeypatch.chdir(tmp_path)

    with pytest.raises(LogNotFoundError, match="sequence_number"):
        LiveEntrypoint(_live_workspace(sequence_number=None)).execute()
