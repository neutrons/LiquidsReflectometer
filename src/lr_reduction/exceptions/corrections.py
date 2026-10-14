"""Exceptions raised while preparing or applying a correction to a workspace."""

from lr_reduction.exceptions.base import LrReductionError, LrValidationError


class CorrectionError(LrReductionError):
    """Base for any failure preparing or applying a correction to a workspace."""


class EmissionTimeCorrectionError(CorrectionError, LrValidationError):
    """The emission-time coefficients could not be resolved for a run: the configuration
    asks for the correction, supplies no coefficients, and the run's chopper logs do not
    say what they are."""
