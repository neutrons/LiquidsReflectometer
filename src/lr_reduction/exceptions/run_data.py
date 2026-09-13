"""Exceptions raised while constructing or validating a RunData."""

from lr_reduction.exceptions.base import LrReductionError, LrValidationError


class RunDataError(LrReductionError):
    """Base for any failure constructing or validating a RunData.

    Deliberately does not gain a log-lookup member: the metadata a `RunData` exposes is
    reached through `SampleLogs`, whose own `SampleLogsError` family already covers a
    missing or ambiguous log. Those propagate from `RunData`'s properties unwrapped.
    """


class IncompleteRunDataError(RunDataError, LrValidationError):
    """RunData was constructed without a required field (workspace, run_numbers)."""
