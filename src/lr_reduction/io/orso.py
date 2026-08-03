"""ORSO format I/O module."""

from orsopy import fileio

from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult
from lr_reduction.utils import get_logger
from lr_reduction.utils.files import get_sequence_id_from_path

logger = get_logger(__name__)

header = fileio.orso.Orso.empty()


def read_config(filepath: str) -> ReductionConfig:
    """Read an ORSO file and return a ReductionConfig populated with metadata from the ORSO file.

    Placeholder: real metadata extraction from the ORSO header's recorded configuration
    (§3.3.8.1) is not yet implemented; this fabricates a minimal valid single-run config from
    the run number embedded in the filename.
    """
    run_number = int(get_sequence_id_from_path(filepath))
    direct_beam_name = "PLACEHOLDER_DB"
    return ReductionConfig(
        direct_beams={direct_beam_name: DirectBeamConfig(db_runs=[run_number])},
        runs=[ReflectedRunConfig(sequence_number=1, direct_beam=direct_beam_name, run_number=run_number)],
    )


def read_single_run(filepath: str) -> ReductionResult:
    """Read an ORSO file and return a ReductionResult populated with data from the ORSO file."""
    ...


def read_partials(partial_dir: str, sequence_id: int) -> list[ReductionResult]:
    """Discover and load partial ORSO files for a given sequence ID."""
    ...


def write(result: ReductionResult | CombinedReductionResult, output_dir: str = ".") -> str:
    """Persist *result* as an ORSO file under *output_dir* and return the path written."""
    logger.info(f"Writing ORSO reduced data for {type(result).__name__} to {output_dir}")
    ...
