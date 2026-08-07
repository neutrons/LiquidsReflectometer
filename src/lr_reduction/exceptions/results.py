"""Exceptions raised while constructing or validating reduction results."""

from lr_reduction.exceptions.base import LrReductionError, LrValidationError


class ResultError(LrReductionError):
    """Base for any failure constructing or validating a reduction result."""


class IncompleteDataError(ResultError, LrValidationError):
    """Only part of a column group that must be supplied together was provided."""


class MalformedDataError(ResultError, LrValidationError):
    """A data column has the wrong dimensionality or an inconsistent length."""
