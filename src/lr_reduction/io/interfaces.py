from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

from lr_reduction.config.model import ReductionConfig
from lr_reduction.io.run_data import RunData


class Loader(ABC):
    """Base class for all loaders."""


class RunLoaderInterface(Loader):
    """
    Abstract loader for a single experimental run.

    Accepts a run number and returns a RunData instance containing
    neutron event data. Concrete implementations are responsible for
    locating and reading data.
    """

    @abstractmethod
    def load(self, run_number: int) -> RunData:
        """Load raw event data for *run_number* and return it as RunData."""


class ConfigLoaderInterface(Loader):
    """
    Abstract loader for a reduction configuration file.

    Accepts a file-system path and returns a validated ReductionConfig.
    Concrete implementations may support different file formats (YAML, JSON,
    XML legacy template, …) while sharing this interface.
    """

    @abstractmethod
    def load(self, path: Path) -> ReductionConfig:
        """Parse the configuration file at *path* and return a ReductionConfig."""
