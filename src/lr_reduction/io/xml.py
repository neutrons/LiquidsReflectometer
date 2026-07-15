"""Legacy XML I/O module. Reads legacy RefRed XML templates as configuration input."""

from lr_reduction.models.config import ReductionConfig
from lr_reduction.utils import get_logger

logger = get_logger(__name__)


def read_config(filepath: str) -> ReductionConfig:
    """Read a legacy XML template and convert it to a ReductionConfig."""
    logger.info(f"Reading legacy XML configuration from {filepath}")
    sequence_id = get_sequence_id_from_path(filepath)
    dbs = DirectBeamConfig(name="PLACEHOLDER NAME", db_runs=["PLACEHOLDER DB RUNS"])
    refs = ReflectedRunConfig(run_id=sequence_id, direct_beam=dbs.name)
    return ReductionConfig(
        sequence_id=sequence_id,
        direct_beams=[dbs],
        reflected_runs=[refs],
    )
