"""Legacy XML I/O module. Reads legacy RefRed XML templates as configuration input."""

from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.utils import get_logger
from lr_reduction.utils.files import get_sequence_id_from_path

logger = get_logger(__name__)


def read_config(filepath: str) -> ReductionConfig:
    """Read a legacy XML template and convert it to a ReductionConfig.

    Placeholder: real XML parsing (§7, tools/xml_to_yaml.py) is not yet implemented; this
    fabricates a minimal valid single-run config from the run number embedded in the filename.
    """
    logger.info(f"Reading legacy XML configuration from {filepath}")
    run_number = int(get_sequence_id_from_path(filepath))
    direct_beam_name = "PLACEHOLDER_DB"
    return ReductionConfig(
        direct_beams={direct_beam_name: DirectBeamConfig(db_runs=[run_number])},
        runs=[ReflectedRunConfig(sequence_number=1, direct_beam=direct_beam_name, run_number=run_number)],
    )
