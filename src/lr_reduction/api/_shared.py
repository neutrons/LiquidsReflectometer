"""Helpers shared by entrypoint leaves.

Both groups here stand in for logic that doesn't exist yet (§6.4, §1.2): locating a
run's standard configuration, and fabricating a placeholder result until Op 1/Op 2/Op 3
are implemented. Centralizing them keeps the leaves thin wrappers (§6.4.3) instead of
duplicating stub logic across manual.py/live.py/autoreduce.py.
"""

from __future__ import annotations

from pathlib import Path

from lr_reduction.io import RunLoader
from lr_reduction.models.config import DirectBeamConfig, ReductionConfig
from lr_reduction.models.run_data import RunData
from lr_reduction.types import ID
from lr_reduction.utils import get_logger

logger = get_logger(__name__)


def get_direct_beam_config(sequence_number: ID, config: ReductionConfig) -> DirectBeamConfig:
    """Get the direct beam configuration for a given run, according to the configuration.

    `DirectBeamConfig.direct_beam_run_numbers` gives the run numbers to load; the object itself is
    what `DirectBeamCompositionOperation` needs, so callers can carry a single resolved
    value through `load_data` -> `call_operations` instead of re-deriving it from `config`.
    """
    db_name = config.runs[sequence_number].direct_beam
    return config.direct_beams[db_name]




def get_direct_beams(run_loader: RunLoader, db_config: DirectBeamConfig) -> list[RunData]:
    """Load the direct beam runs for a given run, according to the configuration."""
    return [run_loader.load(run_number) for run_number in db_config.direct_beam_run_numbers]


def locate_standard_configuration(run_number: ID) -> Path:
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
