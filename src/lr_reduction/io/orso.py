"""ORSO format I/O module."""

from orsopy import fileio

from lr_reduction.config.model import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.types import SequenceResult, SingleRunResult
from lr_reduction.utils import get_logger, get_sequence_id_from_path

logger = get_logger(__name__)

header = fileio.orso.Orso.empty()


def read_config(filepath: str) -> ReductionConfig:
    """Read an ORSO file and return a ReductionConfig populated with metadata from the ORSO file."""
    sequence_id = get_sequence_id_from_path(filepath)
    dbs = DirectBeamConfig(name="PLACEHOLDER NAME", db_runs=["PLACEHOLDER DB RUNS"])
    refs = ReflectedRunConfig(run_id="PLACEHOLDER REF RUN", direct_beam=dbs.name)
    return ReductionConfig(
        sequence_id=sequence_id,
        direct_beams=[dbs],
        reflected_runs=[refs],
    )


def read_single_run(filepath: str) -> SingleRunResult:
    """Read an ORSO file and return a SingleRunResult populated with data from the ORSO file."""
    ...


def read_partials(partial_dir: str, sequence_id: int) -> list[SingleRunResult]:
    """Discover and load partial ORSO files for a given sequence ID."""
    ...


def write(result: SingleRunResult | SequenceResult, output_dir: str = ".") -> str:
    """Persist *result* as an ORSO file under *output_dir* and return the path written."""
    logger.info(f"Writing ORSO reduced data for {type(result).__name__} to {output_dir}")
    ...
