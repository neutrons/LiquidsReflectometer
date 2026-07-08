"""ORSO format I/O module."""

import datetime
from pathlib import Path

import numpy as np
from orsopy import fileio
from orsopy.fileio import (
    DataSource,
    Experiment,
    InstrumentSettings,
    Measurement,
    Orso,
    Person,
    Reduction,
    Sample,
    Software,
    Value,
)

from lr_reduction import __version__ as lr_reduction_version
from lr_reduction.config.model import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
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
