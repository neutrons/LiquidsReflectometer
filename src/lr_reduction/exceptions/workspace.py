"""Exceptions raised while resolving a Mantid workspace."""

from lr_reduction.exceptions.base import LrReductionError, NotFoundError


class WorkspaceError(LrReductionError):
    """Base for any failure resolving a workspace handle."""


class WorkspaceNotFoundError(WorkspaceError, NotFoundError):
    """The analysis data service holds no workspace of the requested name.

    Mantid signals this with a bare `KeyError`, which sits outside this package's
    exception family; `workspace_handle` translates it so a caller can catch a deleted or
    replaced workspace alongside every other `LrReductionError`.
    """
