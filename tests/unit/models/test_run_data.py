from pathlib import Path

import pytest
from mantid.kernel import FloatTimeSeriesProperty, Int32TimeSeriesProperty
from mantid.simpleapi import CreateSampleWorkspace, DeleteWorkspace, mtd

from lr_reduction.exceptions import IncompleteRunDataError, LogNotFoundError, WorkspaceNotFoundError
from lr_reduction.models.config import RunFilter
from lr_reduction.models.run_data import RunData
from lr_reduction.utils.sample_logs import SampleLogs


def _event_workspace() -> str:
    """A small event workspace, registered under a unique name.

    Returns the name, not the object: `RunData` addresses its workspace by name (§11.1.6).
    """
    name = mtd.unique_hidden_name()
    CreateSampleWorkspace(
        WorkspaceType="Event",
        NumBanks=1,
        BankPixelWidth=1,
        NumEvents=1,
        OutputWorkspace=name,
    )
    return name


@pytest.fixture
def workspace():
    """Name of a run workspace carrying the two logs `RunData`'s passthrough properties read.

    `sequence_number` is a constant time series, as the DAS records it; `sequence_id` is a
    scalar. Function-scoped so no test can contaminate another's logs.
    """
    name = _event_workspace()
    run = mtd[name].getRun()

    recorded = Int32TimeSeriesProperty("sequence_number")
    recorded.addValue("2020-01-01T00:00:00", 3)
    recorded.addValue("2020-01-01T00:00:01", 3)
    run.addProperty("sequence_number", recorded, True)
    run.addProperty("sequence_id", 778, "", True)

    yield name
    DeleteWorkspace(name)


# --- construction and validation -------------------------------------------------------


def test_a_workspace_name_is_required():
    with pytest.raises(IncompleteRunDataError, match="workspace name"):
        RunData(workspace=None, run_numbers=(12345,))


def test_an_empty_workspace_name_is_rejected():
    """`.name()` is the empty string for a workspace that was never registered, so this is
    how a nameless workspace reaches construction from a caller holding an object."""
    with pytest.raises(IncompleteRunDataError, match="workspace name"):
        RunData(workspace="", run_numbers=(12345,))


def test_at_least_one_run_number_is_required(workspace):
    """The one piece of identity the loader must track itself: a merged workspace's own
    `run_number` log no longer names every constituent run (§3.1.3)."""
    with pytest.raises(IncompleteRunDataError, match="run number"):
        RunData(workspace=workspace, run_numbers=())


def test_a_name_the_analysis_data_service_does_not_hold_is_rejected():
    with pytest.raises(WorkspaceNotFoundError, match="no_such_workspace"):
        RunData(workspace="no_such_workspace", run_numbers=(12345,))


def test_a_deleted_workspace_is_rejected_at_construction():
    """Fails at the seam that supplied the stale name, rather than later at the first log
    read, which is where the failure would otherwise surface."""
    name = _event_workspace()
    DeleteWorkspace(name)

    with pytest.raises(WorkspaceNotFoundError, match=name):
        RunData(workspace=name, run_numbers=(12345,))


def test_the_error_events_workspace_name_is_validated_too(workspace):
    with pytest.raises(WorkspaceNotFoundError, match="no_such_workspace"):
        RunData(workspace=workspace, run_numbers=(12345,), error_events_workspace="no_such_workspace")


def test_constructor_defaults_the_optional_provenance(workspace):
    run = RunData(workspace=workspace, run_numbers=(12345,))

    assert run.error_events_workspace is None
    assert run.source_paths == ()
    assert run.applied_filter is None


def test_from_workspace_normalizes_sequences_to_tuples(workspace):
    """The ergonomic entry point: a caller holding lists need not convert them, and the
    stored identity is still immutable."""
    run = RunData.from_workspace(
        workspace,
        run_numbers=[12345, 12346],
        source_paths=[Path("REF_L_12345.nxs.h5"), Path("REF_L_12346.nxs.h5")],
    )

    assert run.run_numbers == (12345, 12346)
    assert run.source_paths == (Path("REF_L_12345.nxs.h5"), Path("REF_L_12346.nxs.h5"))


def test_from_workspace_carries_every_field(workspace):
    error_events = _event_workspace()
    run_filter = RunFilter(start_time=0.0, stop_time=10.0)

    run = RunData.from_workspace(
        workspace,
        run_numbers=(12345,),
        source_paths=(Path("REF_L_12345.nxs.h5"),),
        error_events_workspace=error_events,
        applied_filter=run_filter,
    )

    assert run.workspace == workspace
    assert run.error_events_workspace == error_events
    assert run.applied_filter is run_filter
    DeleteWorkspace(error_events)


def test_from_workspace_validates(workspace):
    with pytest.raises(IncompleteRunDataError, match="run number"):
        RunData.from_workspace(workspace, run_numbers=())


# --- label -----------------------------------------------------------------------------


def test_label_of_a_single_run(workspace):
    assert RunData.from_workspace(workspace, run_numbers=(12345,)).label == "12345"


def test_label_of_a_summed_run(workspace):
    assert RunData.from_workspace(workspace, run_numbers=(12345, 12346)).label == "12345+12346"


# --- logs ------------------------------------------------------------------------------


def test_logs_reads_through_to_the_workspace(workspace):
    run = RunData.from_workspace(workspace, run_numbers=(12345,))

    assert isinstance(run.logs, SampleLogs)
    assert run.logs["sequence_id"] == 778


def test_logs_is_built_fresh_on_each_access(workspace):
    """Not cached: `SampleLogs` re-resolves the name on every read so a workspace an
    algorithm rewrote under that name is picked up, rather than the one it replaced."""
    run = RunData.from_workspace(workspace, run_numbers=(12345,))

    assert run.logs is not run.logs


def test_logs_follows_the_name_when_the_workspace_is_replaced(workspace):
    """The reason the field holds a name and not an object: an algorithm writing its output
    back under the same name replaces the analysis data service entry, and the RunData must
    read the replacement."""
    run = RunData.from_workspace(workspace, run_numbers=(12345,))
    assert run.sequence_id == 778

    CreateSampleWorkspace(WorkspaceType="Event", NumBanks=1, BankPixelWidth=1, NumEvents=1, OutputWorkspace=workspace)
    mtd[workspace].getRun().addProperty("sequence_id", 779, "", True)

    assert run.sequence_id == 779


# --- sequence identity -----------------------------------------------------------------


def test_sequence_id_reads_the_recorded_log(workspace):
    run = RunData.from_workspace(workspace, run_numbers=(12345,))

    assert run.sequence_id == 778
    assert isinstance(run.sequence_id, int)


def test_sequence_number_reads_the_recorded_log(workspace):
    run = RunData.from_workspace(workspace, run_numbers=(12345,))

    assert run.sequence_number == 3
    assert isinstance(run.sequence_number, int)


def test_sequence_number_is_coerced_from_a_double_log():
    """The DAS may record it as a double; a float would key `ReductionConfig.runs` by 1.0
    and reach ORSO output as "REF_L_<run>_1.0"."""
    name = _event_workspace()
    recorded = FloatTimeSeriesProperty("sequence_number")
    recorded.addValue("2020-01-01T00:00:00", 2.0)
    mtd[name].getRun().addProperty("sequence_number", recorded, True)

    run = RunData.from_workspace(name, run_numbers=(12345,))

    assert run.sequence_number == 2
    assert isinstance(run.sequence_number, int)
    DeleteWorkspace(name)


def test_a_missing_log_raises_the_sample_logs_error_unwrapped():
    """`SampleLogs` owns log lookup failures; `RunData` does not re-wrap them into its own
    family, so a caller catching `LogNotFoundError` still sees one here."""
    name = _event_workspace()
    run = RunData.from_workspace(name, run_numbers=(12345,))

    with pytest.raises(LogNotFoundError, match="sequence_number"):
        _ = run.sequence_number

    DeleteWorkspace(name)
