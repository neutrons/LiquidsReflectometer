from __future__ import annotations

from os import PathLike
from pathlib import Path
from secrets import token_hex

from mantid.api import FileFinder
from mantid.simpleapi import CreateSampleWorkspace, DeleteWorkspace, LoadErrorEventsNexus, LoadEventNexus, mtd

from lr_reduction.exceptions import RunNotFoundError
from lr_reduction.io.interfaces import RunLoaderInterface
from lr_reduction.models.run_data import RunData
from lr_reduction.types import ID, MantidWorkspaceName
from lr_reduction.utils.logging import get_logger
from lr_reduction.utils.sample_logs import SampleLogs

logger = get_logger(__name__)

#: Stands in for the run's recorded sequence_number until the real load reads it.
_PLACEHOLDER_SEQUENCE_NUMBER = 1

#: The NeXus entry holding rejected events; `LoadErrorEventsNexus` names it when a file has none.
_ERROR_EVENTS_ENTRY = "bank_error_events"


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

    def resolve_path(self, run_number: ID) -> Path:
        """The NeXus file to load for *run_number*.

        Found by Mantid's own facility/archive search, the same search `LoadEventNexus`
        runs on a bare `REF_L_<n>` token, so the path the loader reads is known before
        loading.

        Raises
        ------
        RunNotFoundError
            Mantid's search finds no file for *run_number*.
        """
        try:
            found = FileFinder.findRuns(f"REF_L_{run_number}")
        except RuntimeError as error:
            raise RunNotFoundError(f"No NeXus file found for run {run_number}: {error}") from error
        if not found:
            raise RunNotFoundError(f"No NeXus file found for run {run_number}")
        return Path(found[0])

    def _load_single_workspace(self, path: Path) -> tuple[MantidWorkspaceName, MantidWorkspaceName | None]:
        """Load the events in the NeXus file at *path* and, if it has them, its rejected events.

        The names start with the file's run token and a random token drawn once per load,
        e.g. `REF_L_212345__3fa9c1d2` and `REF_L_212345__3fa9c1d2_errors`, so concurrent
        reductions cannot collide in the analysis data service (§11.1.6).

        Returns
        -------
            The name of the events workspace, and the name of the rejected-events workspace,
            or None when the file records no rejected events.

        Raises
        ------
        RunNotFoundError
            No file exists at *path*.
        """
        if not path.is_file():
            raise RunNotFoundError(f"No NeXus file at {path}")
        # TODO: we will probably revisit this unique naming scheme.
        name = f"{path.name.split('.', 1)[0]}__{token_hex(4)}"
        error_events_name = f"{name}_errors"
        LoadEventNexus(Filename=str(path), OutputWorkspace=name)
        try:
            LoadErrorEventsNexus(Filename=str(path), OutputWorkspace=error_events_name)
        except BaseException as error:
            if isinstance(error, RuntimeError) and f"{_ERROR_EVENTS_ENTRY} does not exist" in str(error):
                logger.info(f"{path} records no rejected events")
                return name, None
            DeleteWorkspace(name)
            raise
        return name, error_events_name

    def load(self, run_number: ID) -> RunData:
        """Load raw event data for *run_number* and return it as RunData."""
        logger.info(f"Loading run data for run number {run_number}")
        # Placeholder implementation; replace with actual data loading logic
        return RunData(workspace=_placeholder_workspace(), run_numbers=(run_number,))

    def load_from_path(self, nexus_file_path: str | PathLike[str]) -> RunData:
        """Load raw event data directly from a NeXus file path and return it as RunData."""
        logger.info(f"Loading run data from path {nexus_file_path}")
        # Placeholder implementation; replace with actual data loading logic. The run
        # number is not known until the file is read, so it stands in as 0.
        return RunData(
            workspace=_placeholder_workspace(),
            run_numbers=(0,),
            source_paths=(Path(nexus_file_path),),
        )
