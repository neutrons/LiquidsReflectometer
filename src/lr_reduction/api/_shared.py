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
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult
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


def placeholder_single_run_result(config: ReductionConfig) -> ReductionResult:
    """Fabricate a ReductionResult until Op 1 (direct beam) / Op 2 (single-run reduce) exist (§1.2)."""
    return ReductionResult(
        reflectivity=np.empty((0, 4)),
        reduction_config=config,
        html_report=None,
        lr_reduction_info=_LR_REDUCTION_SOFTWARE,
        mantid_info=_MANTID_SOFTWARE,
        reduction_timestamp=datetime.now(),
    )


def placeholder_sequence_result(config: ReductionConfig) -> CombinedReductionResult:
    """Fabricate a CombinedReductionResult until Op 3 (assembly) exists (§1.2)."""
    return CombinedReductionResult(
        reflectivity=np.empty((0, 4)),
        reduction_config=config,
        html_report=None,
        lr_reduction_info=_LR_REDUCTION_SOFTWARE,
        mantid_info=_MANTID_SOFTWARE,
        reduction_timestamp=datetime.now(),
        assembly_config=config.assembly,
    )
