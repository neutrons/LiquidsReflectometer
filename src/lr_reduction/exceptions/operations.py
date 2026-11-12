"""Exceptions raised by the reduction operations."""

from lr_reduction.exceptions.base import LrReductionError, LrValidationError


class OperationError(LrReductionError):
    """Base for any failure inside a reduction operation."""


class OperationInputError(OperationError, LrValidationError):
    """An operation's data or configuration cannot be processed; the message names each run and cause."""
