from __future__ import annotations

from pathlib import Path

from mantid.simpleapi import CreateSampleWorkspace, mtd

from lr_reduction.io.interfaces import RunLoaderInterface
from lr_reduction.models.run_data import RunData
from lr_reduction.types import ID, MantidWorkspaceName
from lr_reduction.utils.logging import get_logger
from lr_reduction.utils.sample_logs import SampleLogs

logger = get_logger(__name__)

#: Stands in for the run's recorded sequence_number until the real load reads it.
_PLACEHOLDER_SEQUENCE_NUMBER = 1


def _placeholder_workspace() -> MantidWorkspaceName:
    """Name of an empty event workspace carrying a fabricated sequence_number log.

    TODO: replaced by the real NeXus-backed workspace with RunLoader's own implementation.
    """
    name = mtd.unique_hidden_name()
    CreateSampleWorkspace(
        WorkspaceType="Event",
        NumBanks=1,
        BankPixelWidth=1,
        NumEvents=1,
        OutputWorkspace=name,
    )
    SampleLogs(name).insert("sequence_number", _PLACEHOLDER_SEQUENCE_NUMBER)
    return name


class RunLoader(RunLoaderInterface):
    """Loader for single experimental run."""

    def load(self, run_number: ID) -> RunData:
        """Load raw event data for *run_number* and return it as RunData."""
        logger.info(f"Loading run data for run number {run_number}")
        # Placeholder implementation; replace with actual data loading logic
        return RunData.from_workspace(_placeholder_workspace(), run_numbers=(run_number,))

    def load_from_path(self, nexus_file_path: str | Path) -> RunData:
        """Load raw event data directly from a NeXus file path and return it as RunData."""
        logger.info(f"Loading run data from path {nexus_file_path}")
        # Placeholder implementation; replace with actual data loading logic. The run
        # number is not known until the file is read, so it stands in as 0.
        return RunData.from_workspace(
            _placeholder_workspace(),
            run_numbers=(0,),
            source_paths=(Path(nexus_file_path),),
        )
