"""ORSO format I/O module."""

from orsopy import fileio

from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult
from lr_reduction.utils import get_logger

logger = get_logger(__name__)

# the custom Orso header attribute (via Orso.__init__'s **user_data) that carries the complete
# configuration used for the reduction, for reproducibility (§4.2.3); mirrors
# ReductionResult.reduction_config
CONFIG_HEADER_KEY = "configuration_yaml"

header = fileio.orso.Orso.empty()


def read_config(filepath: str) -> ReductionConfig:
    """Read a prior ORSO output file and return the ReductionConfig embedded in its header (§3.3.8.1(b), §4.2.3)."""
    # Placeholder implementation; read the CONFIG_HEADER_KEY block from the ORSO header via
    # orsopy.fileio.load_orso(filepath)[...].info.user_data and construct a ReductionConfig
    ...


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
