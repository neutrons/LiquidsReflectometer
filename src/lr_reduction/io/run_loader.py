from lr_reduction.io.interfaces import RunLoaderInterface
from lr_reduction.io.run_data import RunData
from lr_reduction.utils.logging import get_logger

logger = get_logger(__name__)


class RunLoader(RunLoaderInterface):
    """Loader for single experimental run."""

    def load(self, run_number: int) -> RunData:
        """Load raw event data for *run_number* and return it as RunData."""
        logger.info(f"Loading run data for run number {run_number}")
        # Placeholder implementation; replace with actual data loading logic
        return RunData()
