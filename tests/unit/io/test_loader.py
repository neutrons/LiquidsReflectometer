from pathlib import Path
from types import SimpleNamespace

import pytest
from mantid.kernel import amend_config

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
