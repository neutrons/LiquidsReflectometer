import re
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest
from mantid.kernel import amend_config
from mantid.simpleapi import DeleteWorkspace, mtd

from lr_reduction.exceptions import RunNotFoundError
from lr_reduction.io.run_loader import RunLoader
from lr_reduction.models.run_data import RunData


def test_load_returns_run_data():
    assert isinstance(RunLoader().load(12345), RunData)


def test_load_from_path_returns_run_data(tmp_path):
    nexus_file = tmp_path / "REF_L_12345.nxs.h5"
    nexus_file.touch()

    assert isinstance(RunLoader().load_from_path(nexus_file), RunData)


def test_resolve_path_returns_the_file_mantid_finds(monkeypatch, tmp_path):
    found = tmp_path / "REF_L_12345.nxs.h5"
    searched = []

    def _find_runs(hint):
        searched.append(hint)
        return [str(found)]

    monkeypatch.setattr("lr_reduction.io.run_loader.FileFinder", SimpleNamespace(findRuns=_find_runs))

    assert RunLoader().resolve_path(12345) == found
    assert searched == ["REF_L_12345"]


def test_resolve_path_raises_run_not_found_when_mantid_finds_nothing(monkeypatch):
    def _find_runs(hint):
        raise RuntimeError(f"Unable to find file: search object '{hint}'")

    monkeypatch.setattr("lr_reduction.io.run_loader.FileFinder", SimpleNamespace(findRuns=_find_runs))

    with pytest.raises(RunNotFoundError, match="12345"):
        RunLoader().resolve_path(12345)


@pytest.mark.datarepo
def test_resolve_path_finds_a_run_in_the_data_repository(nexus_dir):
    with amend_config(data_dir=nexus_dir):
        path = RunLoader().resolve_path(198409)

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
