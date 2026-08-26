from __future__ import annotations

from pathlib import Path

from lr_reduction.io.interfaces import RunLoaderInterface
from lr_reduction.models.run_data import RunData
from lr_reduction.types import ID
from lr_reduction.utils.logging import get_logger

logger = get_logger(__name__)


class RunLoader(RunLoaderInterface):
    """Loader for single experimental run."""

    def load(self, run_number: ID) -> RunData:
        """Load raw event data for *run_number* and return it as RunData."""
        logger.info(f"Loading run data for run number {run_number}")
        # Placeholder implementation; replace with actual data loading logic
        return RunData()

    def load_from_path(self, nexus_file_path: str | Path) -> RunData:
        """Load raw event data directly from a NeXus file path and return it as RunData."""
        logger.info(f"Loading run data from path {nexus_file_path}")
        # Placeholder implementation; replace with actual data loading logic
        return RunData()
