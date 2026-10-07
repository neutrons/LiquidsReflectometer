import re
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest
from mantid.kernel import amend_config
from mantid.simpleapi import CreateSampleWorkspace, DeleteWorkspace, mtd

from lr_reduction.exceptions import RunNotFoundError
from lr_reduction.io.run_loader import RunLoader
from lr_reduction.utils.sample_logs import SampleLogs


@pytest.fixture
def loaded_workspaces(monkeypatch):
    """Patch `_load_single_workspace` to register a one-event run 198409 and its rejected
    events, recording each path it is asked to load."""
    loaded_paths = []
    created = []

    def _load_single_workspace(_self, path):
        loaded_paths.append(path)
        name, error_events_name = mtd.unique_hidden_name(), mtd.unique_hidden_name()
        for workspace in (name, error_events_name):
            CreateSampleWorkspace(
                WorkspaceType="Event", NumBanks=1, BankPixelWidth=1, NumEvents=1, OutputWorkspace=workspace
            )
            SampleLogs(workspace).insert("run_number", "198409")
            created.append(workspace)
        return name, error_events_name

    monkeypatch.setattr("lr_reduction.io.run_loader.RunLoader._load_single_workspace", _load_single_workspace)
    yield loaded_paths
    for workspace in created:
        DeleteWorkspace(workspace)


def test_load_from_path_builds_run_data_from_the_loaded_workspaces(loaded_workspaces, tmp_path):
    path = tmp_path / "REF_L_198409.nxs.h5"

    run = RunLoader().load_from_path(path)

    assert loaded_workspaces == [path]
    assert run.run_numbers == (198409,)
    assert run.error_events_workspace is not None
    assert run.source_paths == (path,)


def test_load_from_path_accepts_a_str(loaded_workspaces, tmp_path):
    path = tmp_path / "REF_L_198409.nxs.h5"

    run = RunLoader().load_from_path(str(path))

    assert loaded_workspaces == [path]
    assert run.source_paths == (path,)


def test_load_from_path_without_rejected_events(monkeypatch, tmp_path):
    name = mtd.unique_hidden_name()
    CreateSampleWorkspace(WorkspaceType="Event", NumBanks=1, BankPixelWidth=1, NumEvents=1, OutputWorkspace=name)
    SampleLogs(name).insert("run_number", "198409")
    monkeypatch.setattr(
        "lr_reduction.io.run_loader.RunLoader._load_single_workspace",
        lambda _self, _path: (name, None),
    )
    try:
        run = RunLoader().load_from_path(tmp_path / "REF_L_198409.nxs.h5")

        assert run.error_events_workspace is None
    finally:
        DeleteWorkspace(name)


def test_load_loads_the_resolved_path(loaded_workspaces, monkeypatch, tmp_path):
    path = tmp_path / "REF_L_198409.nxs.h5"
    monkeypatch.setattr("lr_reduction.io.run_loader.RunLoader.get_filepath_for_run", lambda _self, _run_number: path)

    run = RunLoader().load(198409)

    assert loaded_workspaces == [path]
    assert run.source_paths == (path,)


def test_load_raises_run_not_found_without_loading(loaded_workspaces, monkeypatch):
    def _get_filepath_for_run(_self, run_number):
        raise RunNotFoundError(f"No NeXus file found for run {run_number}")

    monkeypatch.setattr("lr_reduction.io.run_loader.RunLoader.get_filepath_for_run", _get_filepath_for_run)

    with pytest.raises(RunNotFoundError, match="12345"):
        RunLoader().load(12345)

    assert loaded_workspaces == []


def test_get_filepath_for_run_returns_the_file_mantid_finds(monkeypatch, tmp_path):
    found = tmp_path / "REF_L_12345.nxs.h5"
    searched = []

    def _find_runs(hint):
        searched.append(hint)
        return [str(found)]

    monkeypatch.setattr("lr_reduction.io.run_loader.FileFinder", SimpleNamespace(findRuns=_find_runs))

    assert RunLoader().get_filepath_for_run(12345) == found
    assert searched == ["REF_L_12345"]


def test_get_filepath_for_run_raises_run_not_found_when_mantid_finds_nothing(monkeypatch):
    def _find_runs(hint):
        raise RuntimeError(f"Unable to find file: search object '{hint}'")

    monkeypatch.setattr("lr_reduction.io.run_loader.FileFinder", SimpleNamespace(findRuns=_find_runs))

    with pytest.raises(RunNotFoundError, match="12345"):
        RunLoader().get_filepath_for_run(12345)


@pytest.mark.datarepo
def test_get_filepath_for_run_finds_a_run_in_the_data_repository(nexus_dir):
    with amend_config(data_dir=nexus_dir):
        path = RunLoader().get_filepath_for_run(198409)

    assert path == Path(nexus_dir) / "REF_L_198409.nxs.h5"


def _ads_names() -> set[str]:
    return set(mtd.getObjectNames())


def test_load_single_workspace_raises_run_not_found_for_a_missing_file(tmp_path):
    before = _ads_names()

    with pytest.raises(RunNotFoundError, match="REF_L_12345"):
        RunLoader()._load_single_workspace(tmp_path / "REF_L_12345.nxs.h5")

    assert _ads_names() == before


@pytest.mark.datarepo
def test_load_single_workspace_loads_events_and_rejected_events(nexus_dir):
    before = _ads_names()

    name, error_events_name = RunLoader()._load_single_workspace(Path(nexus_dir) / "REF_L_198409.nxs.h5")
    try:
        assert re.fullmatch(r"REF_L_198409__[0-9a-f]{8}", name)
        assert error_events_name == f"{name}_errors"
        assert _ads_names() - before == {name, error_events_name}
        assert mtd[name].getNumberEvents() > 0
        assert mtd[error_events_name].getNumberEvents() > 0
        assert mtd[name].blocksize() == 1
    finally:
        DeleteWorkspace(name)
        DeleteWorkspace(error_events_name)


@pytest.mark.datarepo
def test_load_single_workspace_draws_a_new_name_per_load(nexus_dir):
    path = Path(nexus_dir) / "REF_L_198409.nxs.h5"
    loader = RunLoader()

    first = loader._load_single_workspace(path)
    second = loader._load_single_workspace(path)
    try:
        assert set(first).isdisjoint(second)
    finally:
        for name in (*first, *second):
            DeleteWorkspace(name)


@pytest.mark.datarepo
def test_load_single_workspace_without_rejected_events(nexus_dir, tmp_path):
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "REF_L_201285.nxs.h5"
    shutil.copy(Path(nexus_dir) / path.name, path)
    with h5py.File(path, "a") as nexus:
        del nexus["entry/bank_error_events"]
    before = _ads_names()

    name, error_events_name = RunLoader()._load_single_workspace(path)
    try:
        assert error_events_name is None
        assert _ads_names() - before == {name}
    finally:
        DeleteWorkspace(name)


@pytest.mark.datarepo
def test_load_single_workspace_deletes_the_events_when_rejected_events_fail(nexus_dir, monkeypatch):
    def _fail(**_kwargs):
        raise RuntimeError("LoadErrorEventsNexus-v1: unreadable")

    monkeypatch.setattr("lr_reduction.io.run_loader.LoadErrorEventsNexus", _fail)
    before = _ads_names()

    with pytest.raises(RuntimeError, match="unreadable"):
        RunLoader()._load_single_workspace(Path(nexus_dir) / "REF_L_198409.nxs.h5")

    assert _ads_names() == before


def _assert_loaded_run_198409(run, path, before):
    """The RunData of run 198409 loaded from *path*, adding only its two workspaces to the ADS."""
    assert run.run_numbers == (198409,)
    assert isinstance(run.run_numbers[0], int)
    assert run.source_paths == (path,)
    assert _ads_names() - before == {run.workspace, run.error_events_workspace}
    assert mtd[run.workspace].getNumberEvents() > 0
    assert mtd[run.error_events_workspace].getNumberEvents() > 0


@pytest.mark.datarepo
def test_load_reads_a_run_from_the_data_repository_by_number(nexus_dir):
    before = _ads_names()

    with amend_config(data_dir=nexus_dir):
        run = RunLoader().load(198409)
    try:
        _assert_loaded_run_198409(run, Path(nexus_dir) / "REF_L_198409.nxs.h5", before)
    finally:
        DeleteWorkspace(run.workspace)
        DeleteWorkspace(run.error_events_workspace)


@pytest.mark.datarepo
def test_load_from_path_reads_a_run_from_the_data_repository(nexus_dir):
    path = Path(nexus_dir) / "REF_L_198409.nxs.h5"
    before = _ads_names()

    run = RunLoader().load_from_path(path)
    try:
        _assert_loaded_run_198409(run, path, before)
    finally:
        DeleteWorkspace(run.workspace)
        DeleteWorkspace(run.error_events_workspace)
