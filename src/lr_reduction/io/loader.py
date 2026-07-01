from pathlib import Path

from lr_reduction.utils.deprecated import deprecated
from lr_reduction.utils.logging import get_logger

logger = get_logger(__name__)


def load_config(path: Path):
    """Load a ReductionConfig from a configuration file."""
    pass


def load_orso(path: Path):
    """Load a ReductionConfig from an ORSO file."""
    pass


def load(config: str):
    """Load a ReductionConfig from the specified path."""
    path = Path(config)
    logger.info(f"Loading configuration from {path}")
    match path.suffix:
        case ".orso":
            return load_orso(path)
        case ".yaml" | ".yml":
            return load_config(path)
        case ".xml":
            return load_xml(path)


def load_partials(partial_dir: Path, sequence_id: int):
    """Load partial results from the specified directory for the given sequence ID."""
    pass


@deprecated("XML config files are deprecated. Please use YAML or ORSO format instead.")
def load_xml(path: Path):
    """Legacy method to load a ReductionConfig from an XML file, for backward compatibility."""
    pass
