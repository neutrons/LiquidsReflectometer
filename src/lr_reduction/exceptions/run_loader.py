"""Exceptions raised while resolving or loading a run."""

from lr_reduction.exceptions.base import LrReductionError, NotFoundError


class LoaderError(LrReductionError):
    """Base for any failure resolving or loading a run."""


class RunNotFoundError(LoaderError, NotFoundError):
    """The run number or path does not resolve to an existing NeXus file."""
