import numpy as np
import pytest
from mantid.kernel import FloatTimeSeriesProperty, Int32TimeSeriesProperty, StringTimeSeriesProperty
from mantid.simpleapi import AddSampleLog, CreateWorkspace, DeleteWorkspace, mtd

from lr_reduction.exceptions import AmbiguousLogError, LogNotFoundError, LogTypeError, LogUnitError
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
