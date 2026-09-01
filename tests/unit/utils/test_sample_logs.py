import datetime

import numpy as np
import pytest
from mantid.kernel import (
    BoolTimeSeriesProperty,
    FloatTimeSeriesProperty,
    Int32TimeSeriesProperty,
    StringTimeSeriesProperty,
)
from mantid.simpleapi import AddSampleLog, CreateWorkspace, DeleteWorkspace, mtd

from lr_reduction.exceptions import (
    AmbiguousLogError,
    LogNotFoundError,
    LogTypeError,
    LogUnitError,
    WorkspaceNotFoundError,
)
from lr_reduction.utils.sample_logs import SampleLogs


def _float_series(name, entries):
    """A float time series, `entries` being (ISO timestamp, value) pairs."""
    prop = FloatTimeSeriesProperty(name)
    for timestamp, value in entries:
        prop.addValue(timestamp, value)
    return prop


def _uneven_series():
    """0.0 recorded at t=0, then 10.0 at t=100 s and again at t=101 s.

    The arithmetic mean is 20/3. The time-weighted mean is 20/102 = 0.196: Mantid weights
    each entry by its own interval and repeats the final interval's duration, giving
    weights of 100 s, 1 s and 1 s. The third entry is therefore load-bearing — with only
    the first two, both intervals are 100 s and the time-weighted mean is 5.0.

    The two disagreeing is the whole point of `mean` vs `time_average` being separate
    methods.
    """
    return _float_series(
        "uneven",
        [
            ("2020-01-01T00:00:00", 0.0),
            ("2020-01-01T00:01:40", 10.0),
            ("2020-01-01T00:01:41", 10.0),
        ],
    )


@pytest.fixture
def workspace():
    """A workspace carrying one log of every shape this class must handle.

    Function-scoped on purpose: `insert` mutates the workspace, so a shared workspace
    would let the write tests contaminate the read tests.
    """
    ws = CreateWorkspace(DataX=[0, 1], DataY=[0, 10], OutputWorkspace=mtd.unique_hidden_name())
    AddSampleLog(ws, LogName="scalar_int", LogText="42", LogType="Number", NumberType="Int")
    AddSampleLog(ws, LogName="scalar_float", LogText="3.14", LogType="Number", NumberType="Double")
    AddSampleLog(ws, LogName="scalar_str", LogText="hello", LogType="String")
    AddSampleLog(ws, LogName="with_units", LogText="1.5", LogType="Number", NumberType="Double", LogUnit="degree")

    run = ws.getRun()

    constant_ints = Int32TimeSeriesProperty("constant_int_series")
    constant_ints.addValue("2020-01-01T00:00:00", 7)
    constant_ints.addValue("2020-01-01T00:00:01", 7)
    run.addProperty("constant_int_series", constant_ints, True)

    constant_strings = StringTimeSeriesProperty("constant_str_series")
    constant_strings.addValue("2020-01-01T00:00:00", "steady")
    constant_strings.addValue("2020-01-01T00:00:01", "steady")
    run.addProperty("constant_str_series", constant_strings, True)

    varying_strings = StringTimeSeriesProperty("varying_str_series")
    varying_strings.addValue("2020-01-01T00:00:00", "first")
    varying_strings.addValue("2020-01-01T00:00:01", "second")
    run.addProperty("varying_str_series", varying_strings, True)

    run.addProperty("uneven", _uneven_series(), True)
    run.addProperty("empty", FloatTimeSeriesProperty("empty"), True)
    run.addProperty("all_nan", _float_series("all_nan", [("2020-01-01T00:00:00", float("nan"))]), True)

    # Vector-valued logs: a plain sequence with no timestamps, which is what `insert`
    # writes and what a composite run's per-run values are recorded as. Mantid gives
    # these a `Vector*PropertyWithValue`, not a time series.
    run.addProperty("vector", [1.0, 2.0, 6.0], "", True)
    run.addProperty("vector_str", ["first", "second"], "", True)

    # A series with one NaN among good entries. Distinct from `all_nan`: this one has a
    # real arithmetic mean if the NaN is skipped, so it separates "reject NaN" from
    # "reject a log with nothing readable in it".
    run.addProperty(
        "some_nan",
        _float_series(
            "some_nan",
            [
                ("2020-01-01T00:00:00", 1.0),
                ("2020-01-01T00:00:01", float("nan")),
                ("2020-01-01T00:00:02", 3.0),
            ],
        ),
        True,
    )

    # Boolean logs, which `insert` writes and `[...]` reads. `np.bool_` is not a subtype
    # of `np.number`, so these are the case a dtype check silently excludes.
    run.addProperty("flag_scalar", True, "", True)
    flag_series = BoolTimeSeriesProperty("flag_series")
    flag_series.addValue("2020-01-01T00:00:00", True)
    flag_series.addValue("2020-01-01T00:00:01", False)
    run.addProperty("flag_series", flag_series, True)

    # An empty *string* series. `np.asarray([])` is float64 whatever the log holds, so
    # this is the log a dtype-first check mistakes for an empty numeric one.
    run.addProperty("empty_str", StringTimeSeriesProperty("empty_str"), True)

    yield ws
    DeleteWorkspace(ws)


@pytest.fixture
def logs(workspace):
    return SampleLogs(workspace)


# --- construction -----------------------------------------------------------------


def test_accepts_a_workspace_object(workspace):
    assert SampleLogs(workspace)["scalar_int"] == 42


def test_accepts_a_workspace_name(workspace):
    assert SampleLogs(str(workspace))["scalar_int"] == 42


def test_a_named_workspace_is_resolved_on_every_read():
    """Regression guard: an algorithm writing its output back to the same name replaces
    the entry in the analysis data service with a new workspace object. A wrapper that
    resolved the name once would keep reading the old object's logs — or report a log
    that is plainly present as missing."""
    name = mtd.unique_hidden_name()
    try:
        before = CreateWorkspace(DataX=[0, 1], DataY=[0, 10], OutputWorkspace=name)
        AddSampleLog(before, LogName="tag", LogText="before", LogType="String")
        logs = SampleLogs(name)
        assert logs["tag"] == "before"

        after = CreateWorkspace(DataX=[0, 1], DataY=[0, 10], OutputWorkspace=name)
        AddSampleLog(after, LogName="tag", LogText="after", LogType="String")

        assert logs["tag"] == "after"
    finally:
        DeleteWorkspace(name)


# --- __contains__ -----------------------------------------------------------------


def test_contains_reports_presence(logs):
    assert "scalar_int" in logs
    assert "no_such_log" not in logs


def test_presence_does_not_imply_readable(logs):
    """`in` is about existence; `[...]` additionally requires an unambiguous value."""
    assert "uneven" in logs
    with pytest.raises(AmbiguousLogError):
        logs["uneven"]


# --- __getitem__ ------------------------------------------------------------------


def test_getitem_returns_scalars(logs):
    assert logs["scalar_int"] == 42
    assert logs["scalar_float"] == 3.14
    assert logs["scalar_str"] == "hello"


def test_getitem_returns_constant_series_value(logs):
    assert logs["constant_int_series"] == 7


def test_getitem_returns_constant_string_series_value(logs):
    """Regression guard: a string series is a Python list, so the constancy check must
    go through np.asarray or `value == value[0]` is False for every such log."""
    assert logs["constant_str_series"] == "steady"


def test_getitem_returns_python_natives_not_numpy_scalars(logs):
    """Regression guard for numpy leakage: np.int32 is not an int and is not JSON
    serializable, and these values flow into ReductionResult and ORSO output."""
    value = logs["constant_int_series"]
    assert isinstance(value, int)
    assert not isinstance(value, np.generic)


def test_getitem_rejects_a_varying_series(logs):
    with pytest.raises(AmbiguousLogError, match="single_value"):
        logs["uneven"]


def test_getitem_rejects_an_empty_series(logs):
    with pytest.raises(AmbiguousLogError, match="no values"):
        logs["empty"]


def test_getitem_reports_nan_entries_distinctly(logs):
    """An all-NaN log is unreadable because NaN != NaN, not because it varies."""
    with pytest.raises(AmbiguousLogError, match="NaN"):
        logs["all_nan"]


def test_getitem_rejects_a_missing_log(logs):
    with pytest.raises(LogNotFoundError, match="no_such_log"):
        logs["no_such_log"]


# --- property ---------------------------------------------------------------------


def test_property_returns_the_raw_object(logs):
    prop = logs.property("with_units")
    assert prop.value == 1.5
    assert prop.units == "degree"


def test_property_rejects_a_missing_log(logs):
    with pytest.raises(LogNotFoundError):
        logs.property("no_such_log")


# --- mean and time_average --------------------------------------------------------


def test_mean_is_arithmetic_and_time_average_is_weighted(logs):
    """Pins the distinction the design document got wrong: getStatistics().mean is the
    unweighted arithmetic mean, not the time-weighted one."""
    assert logs.mean("uneven") == pytest.approx(20.0 / 3.0)
    assert logs.time_average("uneven") == pytest.approx(0.19607843, rel=1e-6)
    assert logs.mean("uneven") != pytest.approx(logs.time_average("uneven"))


def test_mean_rejects_a_string_log(logs):
    """Mantid returns nan here rather than raising, which would be silently meaningless."""
    with pytest.raises(LogTypeError):
        logs.mean("scalar_str")


def test_time_average_rejects_a_string_log(logs):
    with pytest.raises(LogTypeError):
        logs.time_average("scalar_str")


def test_mean_rejects_an_empty_series(logs):
    """Mantid returns nan for the statistics of an empty series rather than raising,
    which would propagate a silently meaningless number."""
    with pytest.raises(AmbiguousLogError, match="no values"):
        logs.mean("empty")


def test_time_average_rejects_an_empty_series(logs):
    with pytest.raises(AmbiguousLogError, match="no values"):
        logs.time_average("empty")


def test_mean_rejects_a_missing_log(logs):
    with pytest.raises(LogNotFoundError):
        logs.mean("no_such_log")


def test_mean_averages_a_vector_log(logs):
    """Regression guard: Mantid's getStatistics reports nan for a vector-valued log
    rather than averaging it, so `mean` reduces the values itself."""
    assert logs.mean("vector") == pytest.approx(3.0)


def test_mean_of_a_vector_log_returns_a_python_native(logs):
    assert not isinstance(logs.mean("vector"), np.generic)


def test_mean_rejects_a_vector_of_strings(logs):
    with pytest.raises(LogTypeError):
        logs.mean("vector_str")


def test_time_average_rejects_a_vector_log(logs):
    """A vector log records no times, so there is nothing to weight by. Mantid returns
    nan for it, which would be a silently meaningless number."""
    with pytest.raises(LogTypeError, match="no times"):
        logs.time_average("vector")


def test_time_average_still_weights_a_time_series(logs):
    """Guards the vector rejection against over-reach: a series must still be weighted,
    and so must a scalar, whose statistics Mantid reports correctly."""
    assert logs.time_average("uneven") == pytest.approx(0.19607843, rel=1e-6)
    assert logs.time_average("scalar_float") == pytest.approx(3.14)


# --- NaN handling in the fixed-reduction accessors ---------------------------------


def test_mean_rejects_a_series_containing_nan(logs):
    """`np.mean` propagates NaN, so this used to return nan — a silently meaningless
    number reaching Q, the same failure class as the vector-log statistics."""
    with pytest.raises(AmbiguousLogError, match="NaN"):
        logs.mean("some_nan")


def test_time_average_rejects_a_series_containing_nan(logs):
    with pytest.raises(AmbiguousLogError, match="NaN"):
        logs.time_average("some_nan")


def test_nan_rejection_names_the_escape_hatch(logs):
    with pytest.raises(AmbiguousLogError, match="nanmean"):
        logs.mean("some_nan")


def test_single_value_still_honors_a_nan_aware_operation(logs):
    """The reduction is the caller's to state, so `single_value` does not reject NaN --
    it is the escape hatch the other accessors point at."""
    assert logs.single_value("some_nan", operation=np.nanmean) == pytest.approx(2.0)


def test_single_value_propagates_nan_under_the_default_operation(logs):
    """Documented rather than guarded: the default `np.mean` propagates NaN."""
    assert np.isnan(logs.single_value("some_nan"))


def test_find_log_with_units_rejects_nan_even_with_a_nan_aware_operation(logs):
    """The cost of putting the check in the shared numeric layer, pinned so it is a
    decision rather than a surprise: the unit assertion runs first, so a NaN-aware
    operation does not get the chance to skip them."""
    with pytest.raises(AmbiguousLogError, match="NaN"):
        logs.find_log_with_units("some_nan", operation=np.nanmean)


# --- boolean logs ------------------------------------------------------------------


def test_mean_averages_a_boolean_series(logs):
    """`np.bool_` is not a subtype of `np.number`, so a dtype check alone rejects a log
    `insert` itself writes. The mean of a flag is the fraction of entries it held."""
    assert logs.mean("flag_series") == pytest.approx(0.5)


def test_time_average_gives_a_boolean_series_duty_cycle(logs):
    assert logs.time_average("flag_series") == pytest.approx(0.5)


def test_find_log_with_units_reads_a_boolean_log(logs):
    assert logs.find_log_with_units("flag_scalar") is True


def test_getitem_still_reads_boolean_logs(logs):
    assert logs["flag_scalar"] is True


def test_a_written_boolean_log_round_trips_through_every_accessor(logs):
    """Regression guard for the inconsistency itself: the class must not reject a log it
    just wrote.

    Note the two shapes of answer, both intended: the accessors that read a value give the
    flag back as `True`, while `mean` reduces it and so gives the number 1.0.
    """
    logs.insert("dead_time_applied", True)

    assert logs["dead_time_applied"] is True
    assert logs.single_value("dead_time_applied") is True
    assert logs.find_log_with_units("dead_time_applied") is True
    assert logs.mean("dead_time_applied") == pytest.approx(1.0)


# --- the empty-string-series mislabelling -----------------------------------------


def test_mean_reports_an_empty_string_series_as_a_type_error(logs):
    """`np.asarray([])` is float64 whatever the log holds, so a dtype-first check passed
    this log as numeric and blamed its emptiness. The type is the real problem, and a
    caller that branches on the exception kind needs to be told which."""
    with pytest.raises(LogTypeError, match="strings"):
        logs.mean("empty_str")


def test_mean_still_reports_an_empty_numeric_series_as_empty(logs):
    """The other half of the distinction: an empty float series really is just empty."""
    with pytest.raises(AmbiguousLogError, match="no values"):
        logs.mean("empty")


def test_time_average_reports_an_empty_string_series_as_a_type_error(logs):
    with pytest.raises(LogTypeError, match="strings"):
        logs.time_average("empty_str")


# --- single_value -----------------------------------------------------------------


def test_single_value_defaults_to_mean(logs):
    assert logs.single_value("uneven") == pytest.approx(20.0 / 3.0)


def test_single_value_honors_the_operation(logs):
    assert logs.single_value("uneven", operation=np.max) == 10.0
    assert logs.single_value("uneven", operation=np.min) == 0.0


def test_single_value_returns_python_natives(logs):
    assert not isinstance(logs.single_value("uneven", operation=np.max), np.generic)


def test_single_value_passes_scalars_through(logs):
    assert logs.single_value("scalar_int") == 42
    assert logs.single_value("scalar_str") == "hello"


def test_single_value_reads_a_constant_string_series(logs):
    assert logs.single_value("constant_str_series") == "steady"


def test_single_value_rejects_a_varying_string_series(logs):
    """A reduction function is meaningless over strings, so ambiguity still stands."""
    with pytest.raises(AmbiguousLogError):
        logs.single_value("varying_str_series")


def test_single_value_rejects_a_missing_log(logs):
    with pytest.raises(LogNotFoundError):
        logs.single_value("no_such_log")


# --- find_log_with_units ----------------------------------------------------------


def test_find_log_with_units_accepts_a_matching_unit(logs):
    assert logs.find_log_with_units("with_units", "degree") == 1.5


def test_find_log_with_units_accepts_any_of_several_spellings(logs):
    """Real files spell the same unit inconsistently: ths is 'degree', thi is 'deg'."""
    assert logs.find_log_with_units("with_units", ("deg", "degree")) == 1.5


def test_find_log_with_units_does_not_split_a_string_into_characters(logs):
    with pytest.raises(LogUnitError):
        logs.find_log_with_units("with_units", "d")


def test_find_log_with_units_rejects_a_mismatched_unit(logs):
    with pytest.raises(LogUnitError, match="degree"):
        logs.find_log_with_units("with_units", "mm")


def test_find_log_with_units_skips_the_check_when_no_unit_given(logs):
    assert logs.find_log_with_units("with_units") == 1.5


def test_find_log_with_units_requires_an_unambiguous_value_by_default(logs):
    with pytest.raises(AmbiguousLogError):
        logs.find_log_with_units("uneven")


def test_find_log_with_units_reduces_when_given_an_operation(logs):
    assert logs.find_log_with_units("uneven", operation=np.max) == 10.0


def test_find_log_with_units_rejects_a_non_numeric_log(logs):
    with pytest.raises(LogTypeError):
        logs.find_log_with_units("scalar_str")


def test_find_log_with_units_rejects_a_missing_log(logs):
    with pytest.raises(LogNotFoundError):
        logs.find_log_with_units("no_such_log")


# --- insert -----------------------------------------------------------------------


def test_insert_writes_a_float_with_units(logs):
    logs.insert("db_pixel", 2.5, unit="mm")
    assert logs["db_pixel"] == 2.5
    assert logs.property("db_pixel").units == "mm"


def test_insert_writes_an_int(logs):
    logs.insert("db_run", 12345)
    assert logs["db_run"] == 12345


def test_insert_writes_a_string(logs):
    logs.insert("comment", "composite")
    assert logs["comment"] == "composite"


def test_insert_writes_a_boolean(logs):
    logs.insert("dead_time_applied", True)
    assert logs["dead_time_applied"] is True


def test_insert_replaces_an_existing_log(logs):
    logs.insert("scalar_int", 99)
    assert logs["scalar_int"] == 99


def test_insert_writes_a_sequence_as_a_readable_series(logs):
    logs.insert("db_runs", [11111, 22222], unit="")
    assert list(logs.property("db_runs").value) == [11111, 22222]


def test_insert_writes_a_constant_sequence_readable_through_getitem(logs):
    logs.insert("tthd_values", [1.5, 1.5])
    assert logs["tthd_values"] == 1.5


def test_insert_writes_to_a_workspace_outside_the_analysis_data_service():
    """Regression guard: the intermediate workspaces inside an algorithm are not in the
    ADS, and `AddSampleLog` cannot write to those at all — it fails when storing its
    output, after having already added the property."""
    unmanaged = CreateWorkspace(DataX=[0, 1], DataY=[0, 10], StoreInADS=False)
    logs = SampleLogs(unmanaged)

    logs.insert("db_pixel", 2.5, unit="mm")

    assert logs["db_pixel"] == 2.5
    assert logs.property("db_pixel").units == "mm"


@pytest.mark.parametrize(
    "value, expected",
    [(np.True_, True), (np.int32(7), 7), (np.float64(1.5), 1.5)],
    ids=["bool", "int", "float"],
)
def test_insert_converts_numpy_scalars(logs, value, expected):
    """Regression guard: Mantid registers no converter for numpy scalars and rejects
    them outright, so the write side must normalize just as the read side does."""
    logs.insert("converted", value)

    assert logs["converted"] == expected
    assert not isinstance(logs["converted"], np.generic)


def test_insert_converts_numpy_array_elements(logs):
    logs.insert("db_runs", np.array([11111, 22222]))
    assert list(logs.property("db_runs").value) == [11111, 22222]


def test_insert_writes_a_constant_string_sequence_readable_through_getitem(logs):
    """Round-trips the write side against the read side's string-series trap: a string
    vector's value comes back as a Python list, not an array."""
    logs.insert("db_configuration", ["standard", "standard"])
    assert logs["db_configuration"] == "standard"


def test_insert_rejects_an_empty_sequence(logs):
    """Mantid raises a bare RuntimeError from its mapping layer for a zero-length
    sequence log; this class translates Mantid's errors rather than leaking them."""
    with pytest.raises(LogTypeError, match="empty sequence"):
        logs.insert("nothing", [])


@pytest.mark.parametrize(
    "value",
    [None, datetime.datetime(2020, 1, 1), object(), 3 + 4j],
    ids=["none", "datetime", "object", "complex"],
)
def test_insert_rejects_an_unsupported_value_with_a_package_error(logs, value):
    """Iteration used to decide this, so an unsupported value failed as a bare TypeError
    ("'NoneType' object is not iterable") rather than the LogTypeError the docstring
    promises. `None` from an unset optional is the one a caller will actually hit."""
    with pytest.raises(LogTypeError):
        logs.insert("unsupported", value)


def test_insert_rejects_a_mapping_rather_than_writing_its_keys(logs):
    """Worse than an error before this: iterating a mapping yields its keys, so
    `insert("x", {"a": 1})` silently recorded ["a"] and discarded the values."""
    with pytest.raises(LogTypeError, match="mapping"):
        logs.insert("from_mapping", {"a": 1, "b": 2})

    assert "from_mapping" not in logs


def test_insert_rejects_a_zero_dimensional_array(logs):
    """`np.array(5.0)` is not a sequence — iterating it raises "iteration over a 0-d
    array" — but it is the shape a caller lands on after an accidental reduction."""
    with pytest.raises(LogTypeError, match="dimensional"):
        logs.insert("zero_dim", np.array(5.0))


def test_insert_rejects_a_nested_sequence(logs):
    """Mantid rejects this one itself, with a bare ValueError from its mapping layer."""
    with pytest.raises(LogTypeError):
        logs.insert("nested", [[1, 2], [3, 4]])


def test_insert_still_accepts_a_scalar_after_the_shape_checks(logs):
    """Guards the new validation against over-reach."""
    logs.insert("still_scalar", 2.5, unit="mm")
    logs.insert("still_sequence", [1.5, 2.5])

    assert logs["still_scalar"] == 2.5
    assert list(logs.property("still_sequence").value) == [1.5, 2.5]


# --- workspace resolution ----------------------------------------------------------


def test_contains_on_a_deleted_workspace_raises_a_package_error():
    """`AnalysisDataService[name]` raises a bare KeyError, which a membership test has no
    business raising and which sits outside the package's exception family. A name that
    resolved earlier stops resolving once the workspace is deleted or replaced."""
    name = mtd.unique_hidden_name()
    CreateWorkspace(DataX=[0, 1], DataY=[0, 10], OutputWorkspace=name)
    logs = SampleLogs(name)
    assert "run_title" in logs

    DeleteWorkspace(name)

    with pytest.raises(WorkspaceNotFoundError):
        "run_title" in logs  # noqa: B015 (the raise is the assertion)


def test_reading_a_deleted_workspace_raises_a_package_error():
    name = mtd.unique_hidden_name()
    CreateWorkspace(DataX=[0, 1], DataY=[0, 10], OutputWorkspace=name)
    logs = SampleLogs(name)

    DeleteWorkspace(name)

    with pytest.raises(WorkspaceNotFoundError):
        logs["run_title"]
