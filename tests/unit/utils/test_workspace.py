import pytest
from mantid.simpleapi import CreateWorkspace, DeleteWorkspace, mtd

from lr_reduction.utils.workspace import workspace_handle


@pytest.fixture
def workspace():
    ws = CreateWorkspace(DataX=[0, 1], DataY=[0, 10], OutputWorkspace=mtd.unique_hidden_name())
    yield ws
    DeleteWorkspace(ws)


def test_object_passes_through(workspace):
    assert workspace_handle(workspace) is workspace


def test_name_resolves_to_the_same_workspace(workspace):
    assert workspace_handle(str(workspace)).name() == workspace.name()


def test_unknown_name_raises():
    with pytest.raises(KeyError):
        workspace_handle("no_such_workspace")
