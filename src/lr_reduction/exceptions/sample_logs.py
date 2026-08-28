"""Exceptions raised while reading or writing a workspace's sample logs."""

from lr_reduction.exceptions.base import (
    LrReductionError,
    LrValidationError,
    NotFoundError,
    UnsupportedFormatError,
)


class SampleLogsError(LrReductionError):
    """Base for any failure resolving a sample log through `SampleLogs`.

    Deliberately its own family rather than a member of the `RunData` family: `SampleLogs`
    is a general workspace accessor and is not tied to any one model.
    """


class LogNotFoundError(SampleLogsError, NotFoundError):
    """The requested sample log does not exist on the wrapped workspace."""


class AmbiguousLogError(SampleLogsError, LrValidationError):
    """The log cannot be reduced to a single value without the caller choosing how.

    Raised by `SampleLogs.__getitem__` for a time series whose entries are not all equal,
    in which case the caller must state the reduction explicitly via
    `single_value(name, operation=...)` or `time_average(name)`; and by any accessor for a
    series that records no values at all, where no reduction would help.
    """


class LogUnitError(SampleLogsError, LrValidationError):
    """The log's recorded unit is not one the caller was willing to accept."""


class LogTypeError(SampleLogsError, UnsupportedFormatError):
    """A log's contents and the requested operation are incompatible.

    Either a numeric operation was requested on a log that does not hold numbers, or a
    value was written that Mantid cannot record as a log at all.
    """
