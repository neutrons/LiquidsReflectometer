"""The unit of "a loaded run" handed from the loader to the operation layer."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from lr_reduction.exceptions import IncompleteRunDataError
from lr_reduction.models.config import RunFilter
from lr_reduction.types import ID, MantidWorkspace
from lr_reduction.utils.sample_logs import SampleLogs


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
        The primary event workspace. Per-pixel-resolved and never histogrammed here; this
        is the medium the corrections and operations exchange.
    run_numbers
        The constituent run number(s). Tracked explicitly rather than read back off the
        workspace's own `run_number` log, because Mantid's `MergeRuns` keeps a single
        value, so after the loader sums several source runs that log no longer identifies
        every run in the sum. A tuple because it is part of this object's identity and must
        not be mutated in place; note `ReductionResult.run_numbers` is a `list`, so there
        is a deliberate conversion at that seam.
    error_events_workspace
        The paired rejected-event workspace (Mantid `LoadErrorEventsNexus`), an input to
        the dead-time correction. Optional: not every file has one, and not every load
        needs one. When several source runs are summed, the loader merges their
        rejected-event companions in step, so this stays the rejected-event population of
        the whole merged run.
    source_paths
        The NeXus file path(s) actually read. A provenance and debugging aid, independent
        of whether the run was addressed by number or by path.
    applied_filter
        The time or log-value filter the loader applied before producing this RunData, if
        any. Lets downstream and diagnostic code see what happened to this workspace
        without re-deriving it from configuration.
    """

    workspace: MantidWorkspace
    run_numbers: tuple[ID, ...]
    error_events_workspace: MantidWorkspace | None = None
    source_paths: tuple[Path, ...] = ()
    applied_filter: RunFilter | None = None

    def __post_init__(self):
        if self.workspace is None:
            raise IncompleteRunDataError("RunData requires a workspace")
        if not self.run_numbers:
            raise IncompleteRunDataError("RunData requires at least one run number")

    @classmethod
    def from_workspace(
        cls,
        workspace: MantidWorkspace,
        *,
        run_numbers: Sequence[ID],
        source_paths: Sequence[Path] = (),
        error_events_workspace: MantidWorkspace | None = None,
        applied_filter: RunFilter | None = None,
    ) -> RunData:
        """Build a RunData around an already-loaded workspace.

        The validated, ergonomic entry point, and the one the loader always constructs
        through: it accepts any sequence for the two tuple fields and normalizes them, so a
        caller holding a list need not convert.
        """
        return cls(
            workspace=workspace,
            run_numbers=tuple(run_numbers),
            error_events_workspace=error_events_workspace,
            source_paths=tuple(source_paths),
            applied_filter=applied_filter,
        )

    @property
    def logs(self) -> SampleLogs:
        """This run's sample logs — the single source of truth for its metadata.

        Built fresh on each access rather than cached: `SampleLogs` re-resolves its
        workspace on every read by design, and constructing one is a single assignment.
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
        return ID(self.logs["sequence_id"])

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
        # `ID` is an int and the DAS may record the log as a double; a float would key
        # `ReductionConfig.runs` by 1.0 and reach ORSO output as "REF_L_<run>_1.0".
        return ID(self.logs["sequence_number"])
