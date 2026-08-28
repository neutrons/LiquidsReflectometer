"""Mantid workspace helpers shared by the reduction workflow.

Named `workspace` rather than `mantid` so no module here shares a name with the
third-party `mantid` package.
"""

from mantid.api import AnalysisDataService, Workspace

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
    KeyError
        If a name is given that the analysis data service does not hold.
    """
    if isinstance(workspace, str):
        return AnalysisDataService[workspace]
    return workspace
