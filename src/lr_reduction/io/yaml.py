"""YAML I/O module. YAML is the native configuration format."""

import yaml

from lr_reduction.models.config import ReductionConfig
from lr_reduction.utils import get_logger

logger = get_logger(__name__)


def read_config(filepath: str) -> ReductionConfig:
    """Read a YAML configuration file and return a validated ReductionConfig.

    Uses `yaml.BaseLoader` so every scalar stays a string and Pydantic -- not PyYAML's
    YAML-1.1 heuristics -- decides its type (avoids e.g. `012345` silently becoming octal 5349).
    """
    logger.info(f"Reading YAML configuration from {filepath}")
    with open(filepath) as f:
        data = yaml.load(f, Loader=yaml.BaseLoader)
    if not isinstance(data, dict):
        raise ValueError(f"{filepath}: expected a YAML mapping, got {type(data).__name__}")
    return ReductionConfig(**data)
