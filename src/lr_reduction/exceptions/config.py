"""Exceptions raised while loading or validating a reduction configuration."""

from __future__ import annotations

from lr_reduction.exceptions.base import (
    LrReductionError,
    NotFoundError,
    ParseError,
    UnsupportedFormatError,
    ValidationError,
)


class ConfigError(LrReductionError):
    """Base for any failure loading, validating, or querying a reduction configuration.

    Callers outside this package (notably the GUI) should catch this family instead of
    exceptions native to a specific format's loader (`yaml.YAMLError`, `pydantic.ValidationError`,
    …), so they don't need to know which loader handled a given file.
    """

    def __init__(self, message: str, *, filepath: str | None = None):
        self.filepath = filepath
        super().__init__(message)


class ConfigNotFoundError(ConfigError, NotFoundError):
    """The configuration file does not exist."""


class ConfigFileTypeError(ConfigError, UnsupportedFormatError):
    """The file's extension is not a recognized configuration format."""


class ConfigParseError(ConfigError, ParseError):
    """The file could not be parsed: invalid YAML, not a mapping, or otherwise malformed."""


class ConfigValidationError(ConfigError, ValidationError):
    """The parsed configuration failed Pydantic schema/range/referential validation.

    `errors` carries the raw Pydantic error dicts (`pydantic.ValidationError.errors()`) so a
    caller that wants structured per-field detail (e.g. a GUI form) doesn't have to re-parse
    the message string.
    """

    def __init__(self, message: str, *, filepath: str | None = None, errors: list[dict] | None = None):
        self.errors = errors or []
        super().__init__(message, filepath=filepath)
