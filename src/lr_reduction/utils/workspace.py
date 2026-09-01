"""Mantid workspace helpers shared by the reduction workflow.

Named `workspace` rather than `mantid` so no module here shares a name with the
third-party `mantid` package.
"""

from mantid.api import AnalysisDataService, Workspace

from lr_reduction.exceptions import WorkspaceNotFoundError
from lr_reduction.types import MantidWorkspace


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
