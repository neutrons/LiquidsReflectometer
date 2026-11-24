"""Helpers shared by entrypoint leaves.

Both groups here stand in for logic that doesn't exist yet (§6.4, §1.2): locating a
run's standard configuration, and fabricating a placeholder result until Op 1/Op 2/Op 3
are implemented. Centralizing them keeps the leaves thin wrappers (§6.4.3) instead of
duplicating stub logic across manual.py/live.py/autoreduce.py.
"""

from __future__ import annotations

import logging
from pathlib import Path

from lr_reduction.io import RunLoader
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.run_data import RunData
from lr_reduction.types import ID

logger = logging.getLogger(__name__)


def get_direct_beam_name(sequence_number: ID, config: ReductionConfig) -> str:
    """The name of the composite direct beam the reflected run `sequence_number` is normalized by.

    The name selects the entry of `config.direct_beams` to load and is what
    `DirectBeamCompositionOperation` takes, so callers carry it through `load_data` ->
    `call_operations` instead of re-deriving it from `config`.
    """
    return config.runs[sequence_number].direct_beam


def get_direct_beams(run_loader: RunLoader, config: ReductionConfig, direct_beam_name: str) -> list[RunData]:
    """Load the runs of the composite direct beam `direct_beam_name`."""
    return [run_loader.load(run_number) for run_number in config.direct_beams[direct_beam_name].direct_beam_run_numbers]


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
