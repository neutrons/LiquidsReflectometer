from __future__ import annotations

from pathlib import Path
from secrets import token_hex

from mantid.api import FileFinder
from mantid.simpleapi import DeleteWorkspace, LoadErrorEventsNexus, LoadEventNexus

from lr_reduction.exceptions import RunNotFoundError
from lr_reduction.io.interfaces import RunLoaderInterface
from lr_reduction.models.run_data import RunData
from lr_reduction.types import ID, MantidWorkspaceName
from lr_reduction.utils.logging import get_logger
from lr_reduction.utils.workspace import workspace_handle

logger = get_logger(__name__)

#: The NeXus entry holding rejected events; `LoadErrorEventsNexus` names it when a file has none.
_ERROR_EVENTS_ENTRY = "bank_error_events"


class RunLoader(RunLoaderInterface):
    """Loader for single experimental run."""

    def get_filepath_for_run(self, run_number: ID) -> Path:
        """Return the path to the NeXus file to load for *run_number*.

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
        except RuntimeError as error:
            if f"{_ERROR_EVENTS_ENTRY} does not exist" in str(error):
                logger.info(f"{path} records no rejected events")
                return name, None
            DeleteWorkspace(name)
            raise
        return name, error_events_name

    def load(self, run_number: ID) -> RunData:
        """Load raw event data for *run_number* and return it as RunData.

        Raises
        ------
        RunNotFoundError
            Mantid's search finds no file for *run_number*.
        """
        logger.info(f"Loading run data for run number {run_number}")
        return self.load_from_path(self.get_filepath_for_run(run_number))

    def load_from_path(self, nexus_file_path: str | Path) -> RunData:
        """Load raw event data directly from a NeXus file path and return it as RunData.

        The run number is read from the loaded workspace's `run_number` log.

        Raises
        ------
        RunNotFoundError
            No file exists at *nexus_file_path*.
        """
        path = Path(nexus_file_path)
        logger.info(f"Loading run data from path {path}")
        workspace, error_events_workspace = self._load_single_workspace(path)
        return RunData(
            workspace=workspace,
            run_numbers=(workspace_handle(workspace).getRunNumber(),),
            error_events_workspace=error_events_workspace,
            source_paths=(path,),
        )
