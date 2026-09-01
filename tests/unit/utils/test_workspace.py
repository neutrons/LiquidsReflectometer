import pytest
from mantid.simpleapi import CreateWorkspace, DeleteWorkspace, mtd

from lr_reduction.exceptions import LrReductionError, NotFoundError, WorkspaceNotFoundError
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
    with pytest.raises(WorkspaceNotFoundError, match="no_such_workspace"):
        workspace_handle("no_such_workspace")


def test_unknown_name_raise_stays_in_the_package_family():
    """Mantid signals a missing workspace with a bare KeyError, which sits outside the
    exception family. It reaches callers by surprising routes — `name in SampleLogs(name)`
    resolves the name first — so a membership test must not be able to raise KeyError."""
    with pytest.raises(LrReductionError):
        workspace_handle("no_such_workspace")
    with pytest.raises(NotFoundError):
        workspace_handle("no_such_workspace")


def test_a_deleted_workspace_is_reported_as_missing():
    """The failure that actually bites: a name that resolved earlier stops resolving once
    the workspace is deleted. Self-contained rather than using the `workspace` fixture,
    whose teardown cannot delete a workspace this test has already deleted."""
    name = mtd.unique_hidden_name()
    CreateWorkspace(DataX=[0, 1], DataY=[0, 10], OutputWorkspace=name)
    assert workspace_handle(name).name() == name

    DeleteWorkspace(name)

    with pytest.raises(WorkspaceNotFoundError):
        workspace_handle(name)
