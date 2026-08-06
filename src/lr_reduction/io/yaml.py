"""YAML I/O module. YAML is the native configuration format."""

import yaml

from lr_reduction.exceptions.config import ConfigParseError
from lr_reduction.models.config import ReductionConfig
from lr_reduction.utils import get_logger

logger = get_logger(__name__)


def read_config(filepath: str) -> ReductionConfig:
    """Read a YAML configuration file and return a validated ReductionConfig.

    Uses `yaml.BaseLoader` so every scalar stays a string and Pydantic -- not PyYAML's
    YAML-1.1 heuristics -- decides its type (avoids e.g. `012345` silently becoming octal 5349).

    Raises
    ------
    OSError
        if `filepath` does not exist or cannot be opened (e.g. `FileNotFoundError`).
    ConfigParseError
        if the file's contents are not valid YAML, or do not parse to a mapping.
    pydantic.ValidationError
        if the parsed data fails `ReductionConfig`'s schema/range/referential validation.
    """
    logger.info(f"Reading YAML configuration from {filepath}")
    with open(filepath) as f:
        try:
            data = yaml.load(f, Loader=yaml.BaseLoader)
        except yaml.YAMLError as exc:
            logger.error(f"Failed to parse YAML in {filepath}: {exc}")
            raise ConfigParseError(f"{filepath}: invalid YAML ({exc})", filepath=filepath) from exc
    if not isinstance(data, dict):
        raise ConfigParseError(
            f"{filepath}: expected a YAML mapping, got {type(data).__name__}", filepath=filepath
        )
    return ReductionConfig(**data)
