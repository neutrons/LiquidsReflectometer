"""Mantid workspace helpers shared by the reduction workflow.

Named `workspace` rather than `mantid` so no module here shares a name with the
third-party `mantid` package.
"""

from mantid.api import AnalysisDataService, Workspace

from lr_reduction.exceptions import WorkspaceNotFoundError
from lr_reduction.types import MantidWorkspace, MantidWorkspaceName


def workspace_handle(workspace: MantidWorkspace) -> Workspace:
    r"""Resolve a workspace name or object to a workspace handle.

    Parameters
    ----------
    workspace
        Name of the workspace in the analysis data service, or the workspace object.

    Returns
    -------
        The Workspace instance.

    Raises
    ------
    WorkspaceNotFoundError
        If a name is given that the analysis data service does not hold — the workspace
        was never created, or has since been deleted or replaced.
    """
    if isinstance(workspace, str):
        try:
            return AnalysisDataService[workspace]
        except KeyError as exc:
            # Mantid signals this with a bare KeyError, which sits outside this package's
            # exception family. It reaches callers by surprising routes — a membership
            # test like `"sequence_number" in SampleLogs(name)` resolves the name first —
            # so translate it rather than letting a KeyError escape a bool-returning call.
            raise WorkspaceNotFoundError(f"No workspace named {workspace!r} in the analysis data service") from exc
    return workspace


def workspace_exists(name: MantidWorkspaceName) -> bool:
    """Whether the analysis data service holds a workspace of this name.

    The membership half of `workspace_handle`, for a caller that wants to check a name
    rather than resolve it. Kept here so this module stays the one place that touches the
    analysis data service directly.
    """
    return AnalysisDataService.doesExist(name)
