from pathlib import Path

from lr_reduction.config.model import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.io.interfaces import ConfigLoaderInterface
from lr_reduction.io.orso import read_config as read_orso_config
from lr_reduction.utils import deprecated, get_logger, get_sequence_id_from_path

logger = get_logger(__name__)


class ConfigLoader(ConfigLoaderInterface):
    """A loader for reduction configurations from various file formats."""

    def __init__(self):
        self._loaders = {
            ".ort": self._load_orso,
            ".yaml": self._load_config,
            ".yml": self._load_config,
            ".xml": self._load_xml,
        }

    def load(self, path: str) -> ReductionConfig:
        """Load a ReductionConfig from the specified path."""
        logger.info(f"Loading configuration from {path}")
        fp = Path(path)
        loader = self._loaders.get(fp.suffix)
        if not loader:
            raise ValueError(f"Unsupported configuration file format: {fp.suffix}")

        sequence_id = get_sequence_id_from_path(fp)
        return loader(sequence_id)

    def _load_config(self, path: str) -> ReductionConfig:
        """Load a ReductionConfig from a configuration file."""
        sequence_id = get_sequence_id_from_path(path)
        dbs = DirectBeamConfig(name="PLACEHOLDER NAME", db_runs=["PLACEHOLDER DB RUNS"])
        refs = ReflectedRunConfig(run_id="PLACEHOLDER REF RUN", direct_beam=dbs.name)
        return ReductionConfig(
            sequence_id=sequence_id,
            direct_beams=[dbs],
            reflected_runs=[refs],
        )

    def _load_orso(self, path: str) -> ReductionConfig:
        """Load a ReductionConfig from an ORSO file."""
        return read_orso_config(path)

    @deprecated("XML config files are deprecated. Please use YAML or ORSO format instead.")
    def _load_xml(self, path: str) -> ReductionConfig:
        """Legacy method to load a ReductionConfig from an XML file, for backward compatibility."""
        sequence_id = get_sequence_id_from_path(path)
        dbs = DirectBeamConfig(name="PLACEHOLDER NAME", db_runs=["PLACEHOLDER DB RUNS"])
        refs = ReflectedRunConfig(run_id="PLACEHOLDER REF RUN", direct_beam=dbs.name)
        return ReductionConfig(
            sequence_id=sequence_id,
            direct_beams=[dbs],
            reflected_runs=[refs],
        )
