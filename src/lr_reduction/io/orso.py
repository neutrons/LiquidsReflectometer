"""ORSO format I/O module."""

from pathlib import Path

from orsopy import fileio

from lr_reduction.config.model import ReductionConfig
from lr_reduction.types import SequenceResult, SingleRunResult
from lr_reduction.utils.logging import get_logger

logger = get_logger(__name__)

header = fileio.orso.Orso.empty()


def read(filepath: str) -> ReductionConfig:
    """Read an ORSO file and return a ReductionConfig populated with data from the ORSO file."""
    ...


def read_partials(partial_dir: str, sequence_id: int) -> list[SingleRunResult]:
    """Discover and load partial ORSO files for a given sequence ID."""
    ...


def write(results: SingleRunResult | SequenceResult, output_dir: str | Path = "."): ...
