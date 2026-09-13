from pathlib import Path

import pytest
from mantid.kernel import FloatTimeSeriesProperty, Int32TimeSeriesProperty
from mantid.simpleapi import CreateSampleWorkspace, DeleteWorkspace, mtd

from lr_reduction.exceptions import IncompleteRunDataError, LogNotFoundError
from lr_reduction.models.config import RunFilter
from lr_reduction.models.run_data import RunData
from lr_reduction.utils.sample_logs import SampleLogs


def _event_workspace():
    """A small event workspace, named so it is reachable through the analysis data service."""
    return CreateSampleWorkspace(
        WorkspaceType="Event",
        NumBanks=1,
        BankPixelWidth=1,
        NumEvents=1,
        OutputWorkspace=mtd.unique_hidden_name(),
    )


@pytest.fixture
def workspace():
    """A run workspace carrying the two logs `RunData`'s passthrough properties read.

    `sequence_number` is a constant time series, as the DAS records it; `sequence_id` is a
    scalar. Function-scoped so no test can contaminate another's logs.
    """
    ws = _event_workspace()
    run = ws.getRun()

    recorded = Int32TimeSeriesProperty("sequence_number")
    recorded.addValue("2020-01-01T00:00:00", 3)
    recorded.addValue("2020-01-01T00:00:01", 3)
    run.addProperty("sequence_number", recorded, True)
    run.addProperty("sequence_id", 778, "", True)

    yield ws
    DeleteWorkspace(ws)


# --- construction and validation -------------------------------------------------------


def test_workspace_is_required():
    with pytest.raises(IncompleteRunDataError, match="workspace"):
        RunData(workspace=None, run_numbers=(12345,))


def test_at_least_one_run_number_is_required(workspace):
    """The one piece of identity the loader must track itself: a merged workspace's own
    `run_number` log no longer names every constituent run (§3.1.3)."""
    with pytest.raises(IncompleteRunDataError, match="run number"):
        RunData(workspace=workspace, run_numbers=())


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

    assert run.workspace is workspace
    assert run.error_events_workspace is error_events
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
    """Not cached: `SampleLogs` re-resolves its workspace on every read so a name whose
    analysis-data-service entry was replaced does not read stale logs."""
    run = RunData.from_workspace(workspace, run_numbers=(12345,))

    assert run.logs is not run.logs


def test_logs_resolves_a_workspace_given_by_name(workspace):
    """`MantidWorkspace` is a name or an object, and both must reach the same logs."""
    run = RunData.from_workspace(workspace.name(), run_numbers=(12345,))

    assert run.sequence_id == 778


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
    ws = _event_workspace()
    recorded = FloatTimeSeriesProperty("sequence_number")
    recorded.addValue("2020-01-01T00:00:00", 2.0)
    ws.getRun().addProperty("sequence_number", recorded, True)

    run = RunData.from_workspace(ws, run_numbers=(12345,))

    assert run.sequence_number == 2
    assert isinstance(run.sequence_number, int)
    DeleteWorkspace(ws)


def test_a_missing_log_raises_the_sample_logs_error_unwrapped():
    """`SampleLogs` owns log lookup failures; `RunData` does not re-wrap them into its own
    family, so a caller catching `LogNotFoundError` still sees one here."""
    ws = _event_workspace()
    run = RunData.from_workspace(ws, run_numbers=(12345,))

    with pytest.raises(LogNotFoundError, match="sequence_number"):
        _ = run.sequence_number

    DeleteWorkspace(ws)
