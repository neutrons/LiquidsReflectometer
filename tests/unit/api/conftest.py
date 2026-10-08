from pathlib import Path

import pytest
from mantid.simpleapi import CreateSampleWorkspace, mtd

from lr_reduction.models.run_data import RunData
from lr_reduction.utils.sample_logs import SampleLogs


def _synthetic_run(run_number: int, source_paths: tuple[Path, ...] = ()) -> RunData:
    """A one-event run recording sequence_number 1, standing in for a loaded NeXus file."""
    name = mtd.unique_hidden_name()
    CreateSampleWorkspace(WorkspaceType="Event", NumBanks=1, BankPixelWidth=1, NumEvents=1, OutputWorkspace=name)
    SampleLogs(name).insert("sequence_number", 1)
    return RunData(workspace=name, run_numbers=(run_number,), source_paths=source_paths)


@pytest.fixture(autouse=True)
def _synthetic_run_loader(monkeypatch):
    """Entrypoint tests reduce made-up runs, so RunLoader hands back synthetic runs instead of reading files."""
    monkeypatch.setattr(
        "lr_reduction.io.run_loader.RunLoader.load",
        lambda _self, run_number: _synthetic_run(run_number),
    )
    monkeypatch.setattr(
        "lr_reduction.io.run_loader.RunLoader.load_from_path",
        lambda _self, nexus_file_path: _synthetic_run(0, (Path(nexus_file_path),)),
    )
