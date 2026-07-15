"""Helpers shared by entrypoint leaves.

Both groups here stand in for logic that doesn't exist yet (§6.4, §1.2): locating a
run's standard configuration, and fabricating a placeholder result until Op 1/Op 2/Op 3
are implemented. Centralizing them keeps the leaves thin wrappers (§6.4.3) instead of
duplicating stub logic across manual.py/live.py/autoreduce.py.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import mantid
import numpy as np
from orsopy.fileio import Software

from lr_reduction import __version__
from lr_reduction.config.model import ReductionConfig
from lr_reduction.types import SequenceResult, SingleRunResult
from lr_reduction.utils import get_logger

logger = get_logger(__name__)

_LR_REDUCTION_SOFTWARE = Software(name="lr_reduction", version=__version__)
_MANTID_SOFTWARE = Software(name="mantid", version=mantid.__version__)


def locate_standard_configuration(run_number: int) -> Path:
    """Resolve the standard configuration path for *run_number* (§6.4.5.1).

    Placeholder: the beamline's standard directory-layout convention is not yet defined.
    """
    logger.info(f"Resolving standard configuration for run number {run_number}")
    return Path(f"run_{run_number}.yaml")


def locate_configuration_relative_to(output_directory: Path) -> Path:
    """Resolve the configuration path relative to *output_directory* (§6.4.1).

    Placeholder: the beamline's standard directory-layout convention is not yet defined.
    """
    logger.info(f"Resolving configuration relative to output directory {output_directory}")
    return Path(output_directory) / "config_0.yaml"


def parse_numeric_identifier(identifier: str | int, *, field_name: str) -> int:
    """Convert an identifier to int and raise a field-specific error when conversion fails."""
    try:
        return int(identifier)
    except ValueError as error:
        raise ValueError(f"{field_name} must be an integer, got {identifier!r}") from error


def placeholder_single_run_result(config: ReductionConfig, sequence_number: int = 0) -> SingleRunResult:
    """Fabricate a SingleRunResult until Op 1 (direct beam) / Op 2 (single-run reduce) exist (§1.2)."""
    return SingleRunResult(
        reduction_output=np.empty((0, 4)),
        sequence_id=config.sequence_id,
        sequence_number=sequence_number,
        lr_reduction_info=_LR_REDUCTION_SOFTWARE,
        mantid_info=_MANTID_SOFTWARE,
        reduction_timestamp=datetime.now(),
        experiment_id="PLACEHOLDER EXPERIMENT ID",
        configuration_yaml={},
    )


def placeholder_sequence_result(config: ReductionConfig, reflected_runs: list[str]) -> SequenceResult:
    """Fabricate a SequenceResult until Op 3 (assembly) exists (§1.2)."""
    return SequenceResult(
        sequence_id=config.sequence_id,
        reflected_runs=reflected_runs,
        assembly_result=np.empty((0, 4)),
        stitching_scale_factors=[1.0] * len(reflected_runs),
        html_report="",
    )
