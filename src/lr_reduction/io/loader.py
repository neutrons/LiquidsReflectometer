from pathlib import Path

from lr_reduction.config.model import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.utils.deprecated import deprecated
from lr_reduction.utils.logging import get_logger

logger = get_logger(__name__)


def load(config: str) -> ReductionConfig:
    """Load a ReductionConfig from the specified path."""
    path = Path(config)
    logger.info(f"Loading configuration from {path}")
    match path.suffix:
        case ".ort":
            return load_orso(path)
        case ".yaml" | ".yml":
            return load_config(path)
        case ".xml":
            return load_xml(path)
        case _:
            raise ValueError(f"Unsupported configuration file format: {path.suffix}")


def load_config(path: Path) -> ReductionConfig:
    """Load a ReductionConfig from a configuration file."""
    sequence_id = path.stem.split("_")[1]
    dbs = DirectBeamConfig(name="PLACEHOLDER NAME", db_runs=["PLACEHOLDER DB RUNS"])
    refs = ReflectedRunConfig(run_id="PLACEHOLDER REF RUN", direct_beam=dbs.name)
    return ReductionConfig(
        sequence_id=sequence_id,
        direct_beams=[dbs],
        reflected_runs=[refs],
    )


def load_orso(path: Path) -> ReductionConfig:
    """Load a ReductionConfig from an ORSO file."""
    sequence_id = path.stem.split("_")[1]
    dbs = DirectBeamConfig(name="PLACEHOLDER NAME", db_runs=["PLACEHOLDER DB RUNS"])
    refs = ReflectedRunConfig(run_id="PLACEHOLDER REF RUN", direct_beam=dbs.name)
    return ReductionConfig(
        sequence_id=sequence_id,
        direct_beams=[dbs],
        reflected_runs=[refs],
    )


@deprecated("XML config files are deprecated. Please use YAML or ORSO format instead.")
def load_xml(path: Path) -> ReductionConfig:
    """Legacy method to load a ReductionConfig from an XML file, for backward compatibility."""
    sequence_id = path.stem.split("_")[1]
    dbs = DirectBeamConfig(name="PLACEHOLDER NAME", db_runs=["PLACEHOLDER DB RUNS"])
    refs = ReflectedRunConfig(run_id="PLACEHOLDER REF RUN", direct_beam=dbs.name)
    return ReductionConfig(
        sequence_id=sequence_id,
        direct_beams=[dbs],
        reflected_runs=[refs],
    )
