from pathlib import Path

from lr_reduction.io.interfaces import ConfigLoaderInterface
from lr_reduction.io.yaml import read_config as read_yaml_config
from lr_reduction.models.config import ReductionConfig
from lr_reduction.utils import get_logger

logger = get_logger(__name__)


class ConfigLoader(ConfigLoaderInterface):
    """A loader for reduction configurations from various file formats."""

    def __init__(self):
        self._loaders = {
            ".yaml": self._load_config,
            ".yml": self._load_config,
        }

    def load(self, path: str) -> ReductionConfig:
        """Load a ReductionConfig from the specified path."""
        logger.info(f"Loading configuration from {path}")
        fp = Path(path)
        loader = self._loaders.get(fp.suffix)
        if not loader:
            raise ValueError(f"Unsupported configuration file format: {fp.suffix}")

        return loader(str(fp))

    def _load_config(self, path: str) -> ReductionConfig:
        """Load a ReductionConfig from a YAML configuration file."""
        return read_yaml_config(path)
