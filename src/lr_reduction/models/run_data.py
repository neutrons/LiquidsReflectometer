"""The unit of "a loaded run" handed from the loader to the operation layer."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from lr_reduction.exceptions import IncompleteRunDataError, WorkspaceNotFoundError
from lr_reduction.models.config import RunFilter
from lr_reduction.types import ID, MantidWorkspaceName
from lr_reduction.utils.sample_logs import SampleLogs
from lr_reduction.utils.workspace import workspace_exists


@dataclass
class RunData:
    """Raw neutron event data for a single run, plus the provenance a merged or filtered
    workspace's own sample logs can no longer reliably carry.

    Returned by `RunLoaderInterface`. All other run metadata — sequence identity, run
    title, experiment identifier, incidence and detector angles, slit settings, wavelength
    request — is read from `logs`, not duplicated here. That is deliberate: a copy taken at
    load time could silently drift from the workspace once a merge or a filter touches the
    logs, so only what the logs *cannot* reconstruct afterwards is stored.

    Attributes
    ----------
    workspace
        Name of the primary event workspace in the analysis data service. Per-pixel-
        resolved and never histogrammed here; this is the medium the corrections and
        operations exchange. A name rather than a workspace object because §11.1.6 requires
        workspaces to be passed by name, and because a name is re-resolved on every read:
        an algorithm that writes its output back under the same name replaces the entry,
        and a held object would go on reading the workspace it replaced. A caller holding a
        workspace object passes `workspace.name()`; one that was never registered has no
        name to pass, and is rejected rather than silently wrapped.
    run_numbers
        The constituent run number(s). Tracked explicitly rather than read back off the
        workspace's own `run_number` log, because Mantid's `MergeRuns` keeps a single
        value, so after the loader sums several source runs that log no longer identifies
        every run in the sum. A tuple because it is part of this object's identity and must
        not be mutated in place; any sequence may be passed and is stored as one, so a
        caller's list cannot go on aliasing it. Note `ReductionResult.run_numbers` is a
        `list`, so there is a deliberate conversion at that seam.
    error_events_workspace
        Name of the paired rejected-event workspace (Mantid `LoadErrorEventsNexus`), an input to
        the dead-time correction. Optional: not every file has one, and not every load
        needs one. When several source runs are summed, the loader merges their
        rejected-event companions in step, so this stays the rejected-event population of
        the whole merged run.
    source_paths
        The NeXus file path(s) actually read. A provenance and debugging aid, independent
        of whether the run was addressed by number or by path. Stored as a tuple, on the
        same terms as `run_numbers`.
    applied_filter
        The time or log-value filter the loader applied before producing this RunData, if
        any. Lets downstream and diagnostic code see what happened to this workspace
        without re-deriving it from configuration.

    Raises
    ------
    IncompleteRunDataError
        `workspace` is empty, or `run_numbers` is.
    WorkspaceNotFoundError
        The analysis data service holds no workspace of a given name.
    """

    workspace: MantidWorkspaceName
    run_numbers: tuple[ID, ...]
    error_events_workspace: MantidWorkspaceName | None = None
    source_paths: tuple[Path, ...] = ()
    applied_filter: RunFilter | None = None

    def __post_init__(self):
        # An empty name is rejected along with None: that is what `.name()` returns for a
        # workspace which was never registered, and it is the one way a nameless workspace
        # could reach here from a caller holding an object.
        if not self.workspace:
            raise IncompleteRunDataError("RunData requires a workspace name")
        if not self.run_numbers:
            raise IncompleteRunDataError("RunData requires at least one run number")
        # Normalized after the guards above, not before: `tuple(None)` raises a bare
        # TypeError, and `None` from an unset optional is the likely bad input -- the same
        # trap `SampleLogs._insertable_sequence` documents having fallen into.
        self.run_numbers = tuple(self.run_numbers)
        self.source_paths = tuple(self.source_paths)
        # Fail at the seam that introduced a name the analysis data service cannot resolve,
        # rather than much later at the first log read. `WorkspaceNotFoundError` rather than
        # this module's own family: it is the same failure `workspace_handle` reports, and
        # it is about the workspace, not about RunData's own fields.
        self._require_registered(self.workspace)
        if self.error_events_workspace is not None:
            self._require_registered(self.error_events_workspace)

    @staticmethod
    def _require_registered(name: MantidWorkspaceName) -> None:
        """Raise unless the analysis data service holds a workspace of this name."""
        if not workspace_exists(name):
            raise WorkspaceNotFoundError(f"No workspace named {name!r} in the analysis data service")

    @property
    def logs(self) -> SampleLogs:
        """This run's sample logs — the single source of truth for its metadata.

        Built fresh on each access rather than cached: `SampleLogs` re-resolves the name on
        every read by design, and constructing one is a single assignment.
        """
        return SampleLogs(self.workspace)

    @property
    def label(self) -> str:
        """Short human-readable identifier for logging and diagnostics.

        ``"12345"`` for a single run, ``"12345+12346"`` for a sum. A rendering of
        `run_numbers`, not state of its own.
        """
        return "+".join(str(run_number) for run_number in self.run_numbers)

    @property
    def sequence_id(self) -> ID:
        """sequence_id as recorded in this run's NeXus logs (§3.1.2)."""
        return int(self.logs["sequence_id"])

    @property
    def sequence_number(self) -> ID:
        """sequence_number as RECORDED IN THIS RUN'S LOGS (§3.1.2).

        Not necessarily the value this run is being reduced against: per §3.1.2.1 an
        entrypoint may override it for a given invocation, in which case the effective
        value is `ReflectedRunConfig.sequence_number`. Anything selecting a configuration
        entry, or populating `ReductionResult.sequence_number` for output, reads it from
        the config; this property is for diagnostics, logging, and the default
        (no-override) case.
        """
        # The DAS may record the log as a double; a float would key `ReductionConfig.runs`
        # by 1.0 and reach ORSO output as "REF_L_<run>_1.0". Spelled `int` and not `ID`:
        # `ID` is a PEP 695 type alias, which is a typing construct and not callable.
        return int(self.logs["sequence_number"])
