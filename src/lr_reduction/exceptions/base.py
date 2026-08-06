"""Generic, package-wide exceptions.

Domain-specific modules (e.g. `exceptions.config`) extend these rather than raising them
directly, so callers can catch failures either by domain (`ConfigError`) or by kind
(`ParseError`) regardless of domain.
"""


class LrReductionError(Exception):
    """Root exception for all lr_reduction-raised errors."""


class NotFoundError(LrReductionError):
    """A required resource (file, run, record, ...) could not be located."""


class ParseError(LrReductionError):
    """Input could not be parsed: malformed, wrong type, or otherwise unreadable."""


class LrValidationError(LrReductionError):
    """Parsed input failed schema, range, or referential-integrity validation."""


class UnsupportedFormatError(LrReductionError):
    """The input's format/type is not one this package supports."""
