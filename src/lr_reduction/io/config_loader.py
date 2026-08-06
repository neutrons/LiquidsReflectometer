from pathlib import Path

from pydantic import ValidationError

from lr_reduction.exceptions.config import (
    ConfigFileTypeError,
    ConfigNotFoundError,
    ConfigValidationError,
)
from lr_reduction.io.interfaces import ConfigLoaderInterface
from lr_reduction.io.orso import read_config as read_orso_config
from lr_reduction.io.xml import read_config as read_xml_config
from lr_reduction.io.yaml import read_config as read_yaml_config
from lr_reduction.models.config import ReductionConfig
from lr_reduction.utils import deprecated, get_logger

logger = get_logger(__name__)


def _validation_error_fields(exc: ValidationError) -> str:
    """A human-readable, comma-separated list of the dotted field paths that failed."""
    locations = [".".join(str(part) for part in error["loc"]) for error in exc.errors()]
    return ", ".join(locations) if locations else "<root>"


class ConfigLoader(ConfigLoaderInterface):
    """A loader for reduction configurations from various file formats.

    Dispatches to a format-specific reader (YAML, ORSO, legacy XML) based on the file
    extension, and normalizes whatever that reader raises into a `ConfigError` subclass, so
    callers can branch on failure kind (missing file vs. bad syntax vs. schema violation)
    without depending on the internals of any one loader.
    """

    def __init__(self):
        self._loaders = {
            ".ort": self._load_orso,
            ".yaml": self._load_config,
            ".yml": self._load_config,
            ".xml": self._load_xml,
        }

    def load(self, path: str) -> ReductionConfig:
        """Load a ReductionConfig from the specified path.

        Raises
        ------
        ConfigNotFoundError
            if the file does not exist.
        ConfigFileTypeError
            if the file's extension is not a recognized configuration format.
        ConfigParseError
            if the file exists but cannot be parsed (invalid YAML, not a mapping, …).
        ConfigValidationError
            if the parsed data fails schema, range, or referential-integrity validation.
        """
        logger.info(f"Loading configuration from {path}")
        fp = Path(path)
        loader = self._loaders.get(fp.suffix)
        if not loader:
            raise ConfigFileTypeError(f"Unsupported configuration file format: {fp.suffix!r}", filepath=path)
        if not fp.is_file():
            raise ConfigNotFoundError(f"Configuration file not found: {path}", filepath=path)

        try:
            return loader(str(fp))
        except ValidationError as exc:
            fields = _validation_error_fields(exc)
            logger.error(f"Validation failed for {path}, field(s): {fields}")
            raise ConfigValidationError(
                f"{path}: configuration validation failed for field(s): {fields}",
                filepath=path,
                errors=exc.errors(),
            ) from exc
        except FileNotFoundError as exc:
            logger.error(f"Configuration file not found: {path}")
            raise ConfigNotFoundError(f"Configuration file not found: {path}", filepath=path) from exc

    def _load_config(self, path: str) -> ReductionConfig:
        """Load a ReductionConfig from a YAML configuration file."""
        return read_yaml_config(path)

    def _load_orso(self, path: str) -> ReductionConfig:
        """Load a ReductionConfig from an ORSO file."""
        return read_orso_config(path)

    @deprecated("XML config files are deprecated. Please use YAML or ORSO format instead.")
    def _load_xml(self, path: str) -> ReductionConfig:
        """Legacy method to load a ReductionConfig from an XML file, for backward compatibility."""
        return read_xml_config(path)
