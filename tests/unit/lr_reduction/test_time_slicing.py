"""Time slicing (upstream add-time-slicing, intake time-slicing-reconcile): the events, error events and proton
charge a time window selects; how the window travels from reduce_from_file to the file read; how slices are made,
named, reduced and written.

The inputs are **synthetic**. Each test writes a tiny NeXus file with ten pulses, one per second from 0 s, and three
events per pulse whose ids encode their pulse, so a test can name the exact events a window must select. The pulses
sit exactly on the window edges, so the ±1 conventions at an edge are visible (plan L14). The physics of a slice
(NR_Reduction._reduce_single_run) is a stub, as in test_prior_combination.py: its R(Q) is made up, not a measurement.
The window selection, the slice bookkeeping, the header writing and the file and path handling are the real code.

No test may reach /SNS: every path is a tmp_path, and the config's /SNS default points into it (`no_facility_paths`).
"""

import json
import logging
import re
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

import lr_reduction.binary_processing as BP
import lr_reduction.new_reduction_from_file as nrff
import lr_reduction.save_reduced_data as save_fn
from lr_reduction.nr_reduction_calc import NR_Reduction
from lr_reduction.nr_reduction_config import NRReductionConfig

PULSES = np.arange(10.0)  # pulse times (event_time_zero), seconds from the run's start; the last at the run's end
PER_PULSE = 3  # events per pulse
CHARGE = (np.arange(10.0) + 1.0) * 10.0  # per-pulse proton charge (DASlogs/proton_charge/value): each pulse its own
TOTAL = CHARGE.sum() + 0.5  # entry/proton_charge, deliberately not the per-pulse sum: a test can tell which a path reads
CLOSED_AT_9 = float(np.nextafter(9.0, np.inf))  # a final window's end, moved just past the default run's last pulse
RUN = 230001
SEQ_ID = 230001
EXPERIMENT = "IPTS-00001"
LOWRES = [0, 255]


def ids_of(pulses):
    """The event ids of the given pulses, in order: pulse p holds 100 p, 100 p + 1, 100 p + 2."""
    return [100 * p + k for p in pulses for k in range(PER_PULSE)]


def offsets_of(pulses):
    """The events' time-of-flight offsets (µs): 1000 + id."""
    return [1000.0 + i for i in ids_of(pulses)]


def error_offsets_of(pulses):
    """The error bank holds one event per pulse, at 500 + pulse µs."""
    return [500.0 + p for p in pulses]


def write_run(path, events=True, pulses=PULSES, error_pulses=None, charge_times=None, charge=None, seq=1, counts=None,
              duration=None):
    """A run's NeXus file: the logs reduce_from_file groups by, its duration and charge, and (unless `events` is
    False) the event and error-event banks with their pulse times and per-pulse indices. By default every bank has
    the same pulses, one error event per pulse and one charge entry per pulse, as on the real runs measured; a test
    can give the error bank or the charge log pulses of their own, each pulse's number of events (`counts`, 0 for an
    empty pulse), and entry/duration (float32, as on the real files). Pulse i's events have ids 100 i + k."""
    n = len(pulses)
    counts = [PER_PULSE] * n if counts is None else list(counts)
    error_pulses = pulses if error_pulses is None else error_pulses
    charge_times = pulses if charge_times is None else charge_times
    charge = (np.arange(n) + 1.0) * 10.0 if charge is None else charge
    duration = (pulses[-1] if n else 0.0) if duration is None else duration
    with h5py.File(path, "w") as f:
        f["entry/title"] = [b"slicing test"]
        f["entry/duration"] = np.array([duration], dtype=np.float32)  # float32 on the real files (REF_L_184981)
        f["entry/proton_charge"] = np.array([TOTAL])
        logs = f.create_group("entry/DASlogs")
        logs["BL4B:CS:Autoreduce:Sequence:Num/value"] = np.array([seq])
        logs["BL4B:CS:Autoreduce:Sequence:Id/value"] = np.array([SEQ_ID])
        logs["BL4B:Mot:ths.RBV/value"] = np.array([0.6])
        logs["BL4B:Mot:thi.RBV/value"] = np.array([0.6])
        logs["proton_charge/value"] = charge
        logs["proton_charge/time"] = charge_times
        if events:
            ids = np.array([100 * i + k for i in range(n) for k in range(counts[i])], dtype=np.int64)
            f["entry/bank1_events/event_id"] = ids
            f["entry/bank1_events/event_time_offset"] = (1000.0 + ids).astype(np.float32)
            f["entry/bank1_events/event_time_zero"] = np.asarray(pulses, dtype=float)
            f["entry/bank1_events/event_index"] = np.concatenate([[0], np.cumsum(counts)[:-1]]).astype(np.int64) if n else []
            f["entry/bank_error_events/event_time_offset"] = (500.0 + np.asarray(error_pulses)).astype(np.float32)
            f["entry/bank_error_events/event_time_zero"] = np.asarray(error_pulses, dtype=float)
            f["entry/bank_error_events/event_index"] = np.arange(len(error_pulses))


#: Runs shaped as real files are (I-62), each (pulses, events per pulse or None for three, entry/duration). The default run
#: (integer pulses in order, the duration equal to the last pulse) has none of these irregularities.
LAST_184981 = 76.780626  # REF_L_184981's last pulse; its float32 entry/duration, 76.78062439, rounds below it
PAST_184981 = LAST_184981 + 1e-4  # a duration 100 µs past that pulse: on 27 of the 63 real runs it passes the latest pulse
REAL_SHAPED = {
    "float32 duration below the last pulse": (np.append(np.arange(0.0, 77.0), LAST_184981), None, LAST_184981),
    "float32 duration just past the last pulse": (np.append(np.arange(0.0, 77.0), LAST_184981), None, PAST_184981),
    "unsorted tail": (np.concatenate([np.arange(0.0, 17.0), [17.10, 10.72, 10.73]]), None, 17.2),
    "unsorted tail, the latest pulse past the duration": (
        np.concatenate([np.arange(0.0, 70.0), [LAST_184981, 70.5, 70.6]]), None, LAST_184981),
    "empty pulses": (PULSES, [3, 0, 3, 0, 0, 3, 3, 0, 3, 0], PULSES[-1]),
    "beam-off tail": (PULSES, None, 14.5),
    "beam-off tail, a quarter on the last pulse": (PULSES, None, 18.0),
    "last pulse at the duration": (PULSES, None, PULSES[-1]),
    "charge log past the last event pulse": (PULSES, None, 9.5),
}
#: The banks other than the events, where a shape gives them pulses of their own: here the charge log has an entry after
#: the last event pulse, so the run's end is the charge log's (read_run_end takes the latest over the three banks).
SHAPED_BANKS = {"charge log past the last event pulse": {"charge_times": np.append(PULSES, 9.5),
                                                         "charge": (np.arange(11) + 1.0) * 10.0}}


def write_shaped(path, shape):
    pulses, counts, duration = REAL_SHAPED[shape]
    write_run(path, pulses=pulses, counts=counts, duration=duration, **SHAPED_BANKS.get(shape, {}))
    return pulses


@pytest.fixture(autouse=True)
def no_facility_paths(tmp_path, monkeypatch):
    """No test reads or writes the facility tree.
    - The config's /SNS/REF_L default (NRReductionConfig.base_path) points into tmp_path. A test that wrote under it
      fails at teardown: every folder a test uses is given explicitly. Only computing the default is harmless
      (reduce_from_file's hasattr() evaluates the path property).
    - A path spelled under /SNS in the code (the contribution's literals, a mutant restoring one) is refused before
      the file is opened or written: h5py.File and numpy.savetxt raise for it. The attempt is recorded and read at
      teardown, so a handler that catches the raise cannot hide it."""
    default = tmp_path / "facility-default"
    monkeypatch.setattr(NRReductionConfig, "base_path", property(lambda self: default / str(self.experiment_id)))
    attempts, real_file, real_savetxt = [], h5py.File, np.savetxt

    def refuse(path):
        if str(path).startswith("/SNS"):
            attempts.append(str(path))
            raise AssertionError(f"a test reached the facility tree: {path}")

    def guarded_file(name, *args, **kwargs):
        refuse(name)
        return real_file(name, *args, **kwargs)

    def guarded_savetxt(fname, *args, **kwargs):
        refuse(fname)
        return real_savetxt(fname, *args, **kwargs)

    monkeypatch.setattr(h5py, "File", guarded_file)
    monkeypatch.setattr(np, "savetxt", guarded_savetxt)
    yield
    assert attempts == [], f"the facility tree was reached: {attempts}"
    assert not default.exists(), sorted(str(path.relative_to(default)) for path in default.rglob("*"))


@pytest.fixture
def run_file(tmp_path, monkeypatch):
    """One run with events. The per-run log reader is not under test: it returns what convert_to_binary needs."""
    path = tmp_path / f"REF_L_{RUN}.nxs.h5"
    write_run(path)
    monkeypatch.setattr(BP, "get_log_values", lambda _fname: {"start_time": "2025-03-04T11:22:33-05:00"})
    return path


# ---------------------------------------------------------------------------
# T1-T4: the window, at the file read (binary_processing.load_and_extract / convert_to_binary)


def test_no_windows_is_the_whole_run(run_file):
    """T1a: without windows every event, every error event and every pulse's charge is read, as today; the run's
    charge is the file's entry/proton_charge, not a sum recomputed from the per-pulse log (N1)."""
    whole = BP.load_and_extract(run_file)
    explicit = BP.load_and_extract(run_file, None, None)
    for read, again in zip(whole[:5], explicit[:5]):
        np.testing.assert_array_equal(read, again)
    e_offset, event_id, error_offset, pcharge, cpc, _ = whole
    assert event_id.tolist() == ids_of(range(10))
    assert e_offset.tolist() == offsets_of(range(10))
    assert error_offset.tolist() == error_offsets_of(range(10))
    assert np.asarray(pcharge).tolist() == [TOTAL]
    assert cpc.tolist() == CHARGE.tolist()


def test_windows_select_by_pulse_time(run_file):
    """T1b: a window [start, end) selects the events of the pulses with start <= t < end, exactly: [2, 5) is pulses
    2, 3 and 4. The error events and the proton charge are selected by the same predicate."""
    e_offset, event_id, error_offset, pcharge, cpc, _ = BP.load_and_extract(run_file, [2.0], [5.0])
    assert event_id.tolist() == ids_of([2, 3, 4])
    assert e_offset.tolist() == offsets_of([2, 3, 4])
    assert error_offset.tolist() == error_offsets_of([2, 3, 4])
    assert cpc.tolist() == CHARGE[2:5].tolist()
    assert np.asarray(pcharge).tolist() == [CHARGE[2:5].sum()]


def test_start_zero_is_not_special(run_file):
    """T1c: a window from 0 starts at the first pulse like any other start: [0, 3) is pulses 0, 1 and 2."""
    _, event_id, error_offset, _, cpc, _ = BP.load_and_extract(run_file, [0.0], [3.0])
    assert event_id.tolist() == ids_of([0, 1, 2])
    assert error_offset.tolist() == error_offsets_of([0, 1, 2])
    assert cpc.tolist() == CHARGE[0:3].tolist()


def test_contiguous_slices_partition_the_run(run_file):
    """T2a: contiguous windows over the run, the final one closed as the slicing functions close it (close_final_window:
    its end is moved past the last pulse, which sits at the run's end here), select disjoint event sets whose union is
    every event, in order; no event is lost or duplicated at an edge. The slices' charges add up to the run's
    per-pulse total."""
    windows = BP.close_final_window([(0.0, 3.0), (3.0, 6.0), (6.0, 9.0)], *BP.read_run_end(run_file))
    assert windows == [(0.0, 3.0), (3.0, 6.0), (6.0, CLOSED_AT_9)]
    slices = [BP.load_and_extract(run_file, [a], [b]) for a, b in windows]
    assert sum((s[1].tolist() for s in slices), []) == ids_of(range(10))
    assert sum((s[2].tolist() for s in slices), []) == error_offsets_of(range(10))
    assert sum(float(np.asarray(s[3])[0]) for s in slices) == CHARGE.sum()
    assert sum((s[4].tolist() for s in slices), []) == CHARGE.tolist()


@pytest.mark.parametrize("end, pulses", [(9.0, [6, 7, 8]), (12.0, [6, 7, 8, 9])],
                         ids=["end at the last pulse", "end after the run"])
def test_the_file_read_is_half_open_at_every_end(run_file, end, pulses):
    """T1'/T2' (B-3): at the file read every window is half-open, the last pulse's time included: [6, 9) leaves out
    the pulse at 9 s, and an end past it takes it. A window that ends on the last pulse cannot tell there whether it
    is the run's final window or the one before [9, ...); the slicing functions, which know every window, close the
    final one (test_contiguous_slices_partition_every_real_shaped_run)."""
    _, event_id, error_offset, _, cpc, _ = BP.load_and_extract(run_file, [6.0], [end])
    assert event_id.tolist() == ids_of(pulses)
    assert error_offset.tolist() == error_offsets_of(pulses)
    assert cpc.tolist() == CHARGE[pulses].tolist()


def test_pulses_are_selected_by_their_own_time(tmp_path, monkeypatch):
    """T1'a (B-2): a pulse belongs to a window by its own time, wherever it sits in the file. Here the last two pulses
    go back in time (17.10 s, then 10.72 and 10.73 s, as on REF_L_198410): [10, 17) takes the pulses at 10 to 16 s and
    the two at 10.72 and 10.73 s, in file order, and not the one at 17.10 s; the error bank and the charge alike."""
    path = tmp_path / f"REF_L_{RUN}.nxs.h5"
    pulses = write_shaped(path, "unsorted tail")
    monkeypatch.setattr(BP, "get_log_values", lambda _fname: {})
    inside = [i for i, t in enumerate(pulses) if 10.0 <= t < 17.0]
    assert inside == [10, 11, 12, 13, 14, 15, 16, 18, 19]
    _, event_id, error_offset, _, cpc, _ = BP.load_and_extract(path, [10.0], [17.0])
    assert event_id.tolist() == ids_of(inside)
    assert error_offset.tolist() == np.float32(500.0 + pulses[inside]).tolist()
    assert cpc.tolist() == ((np.array(inside) + 1.0) * 10.0).tolist()


@pytest.mark.parametrize("starts, ends, pulses", [
    ([0.0, 5.0], [2.0, 7.0], [0, 1, 5, 6]),
    ([5.0, 0.0], [7.0, 2.0], [5, 6, 0, 1]),
    ([0.0, 2.0], [3.0, 4.0], [0, 1, 2, 2, 3]),
], ids=["two disjoint", "given out of order", "overlapping"])
def test_several_windows_concatenate_in_order(run_file, starts, ends, pulses):
    """T3: several windows in one call concatenate in the order given. Overlapping windows are allowed (§10 A2), and a
    pulse in two of them is read twice: nothing is deduplicated."""
    _, event_id, error_offset, pcharge, cpc, _ = BP.load_and_extract(run_file, starts, ends)
    assert event_id.tolist() == ids_of(pulses)
    assert error_offset.tolist() == error_offsets_of(pulses)
    assert cpc.tolist() == CHARGE[pulses].tolist()
    assert np.asarray(pcharge).tolist() == [CHARGE[pulses].sum()]


@pytest.mark.parametrize("starts, ends", [
    ([1.0], None), (None, [2.0]), ([1.0, 2.0], [3.0]), ([3.0], [3.0]), ([4.0], [3.0]), ([], []),
], ids=["no end", "no start", "unequal lengths", "start at end", "start after end", "no windows"])
def test_invalid_windows_raise_before_reading_events(tmp_path, monkeypatch, starts, ends):
    """T4a: an invalid window list is a ValueError before any event is read: one of the lists missing, unequal
    lengths, a start not before its end, or two empty lists (§10 A1: None is the spelling of "no filter"). The file
    here has no event banks, so a reader that read them first would fail otherwise."""
    path = tmp_path / f"REF_L_{RUN}.nxs.h5"
    write_run(path, events=False)
    monkeypatch.setattr(BP, "get_log_values", lambda _fname: {})
    with pytest.raises(ValueError):
        BP.load_and_extract(path, starts, ends)


@pytest.mark.parametrize("one", [float, np.array], ids=["numbers", "0-d arrays"])
def test_scalars_are_one_window(run_file, one):
    """T4b: a window given as two numbers is the one-window list. A 0-d array is a number: reduce_time_slices'
    arithmetic on entry/duration makes them."""
    scalars = BP.load_and_extract(run_file, one(2.0), one(5.0))
    lists = BP.load_and_extract(run_file, [2.0], [5.0])
    for a, b in zip(scalars[:5], lists[:5]):
        np.testing.assert_array_equal(a, b)
    assert scalars[1].tolist() == ids_of([2, 3, 4])


@pytest.mark.parametrize("start, end", [(20.0, 30.0), (2.25, 2.75)], ids=["after the run", "between two pulses"])
def test_a_window_without_charge_is_a_named_error(run_file, start, end):
    """T4c (the read): a window that selects no pulse with charge cannot be reduced. convert_to_binary says so, naming
    the file and the window, instead of returning None for every caller to unpack (N2)."""
    with pytest.raises(ValueError) as raised:
        BP.convert_to_binary(run_file, LOWRES, start_times=[start], end_times=[end])
    message = str(raised.value)
    assert run_file.name in message and str(start) in message and str(end) in message, message


def test_a_window_after_the_run_selects_nothing(run_file):
    """T4c (the selection): a window after the last pulse, or between two pulses, reads no event, no error event and
    no charge; convert_to_binary then says so (above)."""
    for start, end in [(20.0, 30.0), (2.25, 2.75)]:
        e_offset, event_id, error_offset, pcharge, cpc, _ = BP.load_and_extract(run_file, [start], [end])
        assert (len(e_offset), len(event_id), len(error_offset), len(cpc)) == (0, 0, 0, 0), (start, end)
        assert np.asarray(pcharge).tolist() == [0.0]


def test_a_run_without_pulses_reads_nothing(tmp_path, monkeypatch):
    """A run that recorded no pulse: a window reads nothing (no IndexError on the empty pulse list), and the reduction
    says there is no charge."""
    path = tmp_path / f"REF_L_{RUN}.nxs.h5"
    write_run(path, pulses=np.array([]), charge=np.array([]))
    monkeypatch.setattr(BP, "get_log_values", lambda _fname: {})
    _, event_id, error_offset, _, cpc, _ = BP.load_and_extract(path, [0.0], [5.0])
    assert (len(event_id), len(error_offset), len(cpc)) == (0, 0, 0)
    with pytest.raises(ValueError, match="No proton charge"):
        BP.convert_to_binary(path, LOWRES, start_times=[0.0], end_times=[5.0])


def test_each_bank_is_selected_by_its_own_pulse_times(tmp_path, monkeypatch):
    """T1: the error events and the charge are selected by the window's predicate on their own pulse times, not by
    the event pulses' positions. Here the error bank has a pulse every other second and the charge log one entry more
    at the start (at -1 s): [2, 5) is error pulses 2 and 4, and the charge logged at 2, 3 and 4 s."""
    path = tmp_path / f"REF_L_{RUN}.nxs.h5"
    charge_times = np.concatenate([[-1.0], PULSES])
    charge = np.concatenate([[1000.0], CHARGE])
    write_run(path, error_pulses=np.arange(0.0, 10.0, 2.0), charge_times=charge_times, charge=charge)
    monkeypatch.setattr(BP, "get_log_values", lambda _fname: {})
    _, event_id, error_offset, pcharge, cpc, _ = BP.load_and_extract(path, [2.0], [5.0])
    assert event_id.tolist() == ids_of([2, 3, 4])
    assert error_offset.tolist() == [502.0, 504.0]
    assert cpc.tolist() == CHARGE[2:5].tolist()
    assert np.asarray(pcharge).tolist() == [CHARGE[2:5].sum()]


# ---------------------------------------------------------------------------
# The window's path: reduce_from_file -> reduce -> _reduce_single_run -> _load_and_extract_lambda
# -> _make_binary_files -> convert_to_binary -> load_and_extract (F5). Each link is asked separately: the next
# step records what it received, then stops the call.


class StopError(Exception):
    """Raised by a recording stand-in to end the call once it has what it came for."""


WINDOW = {"start_times": [1.0, 4.0], "end_times": [2.0, 5.0]}


def settings_for(tmp_path, **extra):
    """A one-position reduction settings file (written with save_config_json's keys), as in test_prior_combination."""
    settings = {
        "method_per_run": ["meanTheta"],
        "DBname": ["DB_A1.dat"],
        "RB_Ymin": [146],
        "RB_Ymax": [155],
        "Sname": "blank",
        "experiment_id": EXPERIMENT,
        "AutoScale": False,
        "LambdaMin": [2.7],
        "LambdaMax": [9.5],
        "ThetaShift": [0],
        "ScaleFactor": [1],
        "LambdaMinUse": 2.7,
        "LambdaMaxUse": 9.5,
        **extra,
    }
    path = tmp_path / "reduce_settings.json"
    path.write_text(json.dumps(settings))
    return path


@pytest.fixture
def reducer(tmp_path):
    """An NR_Reduction for one run whose folders are all under tmp_path."""
    config = nrff.json_to_config(nrff.load_from_file(settings_for(tmp_path))["config"])
    config.RBnum = [RUN]
    config.experiment_id = EXPERIMENT
    config.NEXUSpathRB = tmp_path
    config.Spath = tmp_path
    config.DBpath = tmp_path
    config.BINpath = tmp_path
    return NR_Reduction(config)


def recorder(record):
    def stand_in(*args, **kwargs):
        record.append((args, kwargs))
        raise StopError
    return stand_in


def received_window(record):
    (args, kwargs), = record
    return {name: kwargs.get(name) for name in ("start_times", "end_times")}, args


@pytest.mark.usefixtures("run_file")  # the grouping reads the run
def test_the_window_reaches_reduce_from_reduce_from_file(tmp_path, monkeypatch):
    record = []
    monkeypatch.setattr(NR_Reduction, "reduce", recorder(record))
    with pytest.raises(StopError):
        nrff.reduce_from_file([RUN], settings_for(tmp_path), EXPERIMENT, datapath=tmp_path, plot=False,
                              override_params={"Spath": tmp_path}, **WINDOW)
    assert received_window(record)[0] == WINDOW


def test_the_window_reaches_each_run_from_reduce(monkeypatch, reducer):
    record = []
    monkeypatch.setattr(NR_Reduction, "_reduce_single_run", recorder(record))
    with pytest.raises(StopError):
        reducer.reduce(plot=False, **WINDOW)
    window, args = received_window(record)
    assert window == WINDOW and args[1:] == (0, RUN)


def test_the_window_reaches_the_lambda_load_from_a_run(monkeypatch, reducer):
    record = []
    monkeypatch.setattr(NR_Reduction, "_load_and_extract_lambda", recorder(record))
    with pytest.raises(StopError):
        reducer._reduce_single_run(0, RUN, **WINDOW)
    window, args = received_window(record)
    assert window == WINDOW and args[1:] == (0, RUN)


def test_the_window_reaches_the_binary_files_from_the_lambda_load(monkeypatch, reducer):
    record = []
    monkeypatch.setattr(NR_Reduction, "_make_binary_files", recorder(record))
    with pytest.raises(StopError):
        reducer._load_and_extract_lambda(0, RUN, **WINDOW)
    window, args = received_window(record)
    assert window == WINDOW and args[1] == RUN


def test_the_window_reaches_convert_to_binary_from_the_binary_files(monkeypatch, reducer, run_file):
    record = []
    monkeypatch.setattr(BP, "convert_to_binary", recorder(record))
    with pytest.raises((StopError, RuntimeError)):  # _make_binary_files wraps what its call raises
        reducer._make_binary_files(RUN, 0, 50000, **WINDOW)
    window, args = received_window(record)
    assert window == WINDOW and Path(args[0]) == run_file


def test_the_window_reaches_the_file_read_from_convert_to_binary(monkeypatch, run_file):
    record = []
    monkeypatch.setattr(BP, "load_and_extract", recorder(record))
    with pytest.raises(StopError):
        BP.convert_to_binary(run_file, LOWRES, **WINDOW)
    window, args = received_window(record)
    assert window == WINDOW and Path(args[0]) == run_file


# ---------------------------------------------------------------------------
# T4c (the slice), T5, T6: slices made, named, reduced and written (new_reduction_time_resolved, the header writer)


@pytest.fixture
def slicing(tmp_path, monkeypatch):
    """A NeXus folder holding the run, an output folder, a settings file, and a physics stub that records the
    window, subname and folders each reduction received."""
    import lr_reduction.new_reduction_time_resolved as nrtr

    nexus, out = tmp_path / "nexus", tmp_path / "out"
    nexus.mkdir()
    out.mkdir()
    write_run(nexus / f"REF_L_{RUN}.nxs.h5")
    calls, fail = [], []  # fail: (window, exception) pairs; a window may hold lists, so not a dict key

    def fake_reduce_single_run(self, i, rb_num, save=True, start_times=None, end_times=None):  # noqa: ARG001 -- the real method's signature; save is not used
        calls.append({"window": (start_times, end_times), "subname": self.config.subname,
                      "Spath": Path(self.config.Spath), "nexus": Path(self.config.NEXUSpathRB)})
        if start_times is not None and end_times is not None:  # what the real file read selects for this window
            _, ids, errors, _, cpc, _ = BP.load_and_extract(
                Path(self.config.NEXUSpathRB) / f"REF_L_{rb_num}.nxs.h5", start_times, end_times)
            calls[-1]["selected"] = (ids.tolist(), errors.tolist(), cpc.tolist())
        for window, error in fail:
            if window == (start_times, end_times):
                raise error
        q = np.geomspace(0.01, 0.1, 20)
        r = 1e-6 * q**-4 * (1.0 + len(calls))
        zeros = np.zeros_like(q)
        self.log_values = {"title": "slicing test", "ths": 0.6, "thi": 0.6, "ThCen": 0.6}
        self.config.LambdaMinUse = self.config.LambdaMin[i]
        self.config.LambdaMaxUse = self.config.LambdaMax[i]
        result = {"q": q, "r": r, "dr": 0.05 * r, "dq": 0.02 * q, "t": zeros, "l": zeros, "dt": zeros, "dl": zeros}
        return result, self.config, self.log_values

    monkeypatch.setattr(NR_Reduction, "_reduce_single_run", fake_reduce_single_run)
    monkeypatch.setattr(BP, "get_log_values", lambda _fname: {})
    return SimpleNamespace(nrtr=nrtr, nexus=nexus, out=out, calls=calls, fail=fail, settings=settings_for(tmp_path))


def header_lines(path):
    return [line[2:] for line in Path(path).read_text().splitlines() if line.startswith("# ")]


def assert_one_window_line_and_the_marker(out):
    """The time-sliced outputs' header invariant (plan §4): each file has exactly one "Time resolved:" line and one
    format marker, the window line before the marker."""
    files = sorted(Path(out).glob("*.dat"))
    assert files
    for path in files:
        lines = header_lines(path)
        windows = [i for i, line in enumerate(lines) if line.startswith("Time resolved: ")]
        markers = [i for i, line in enumerate(lines) if line.startswith("Header format: ")]
        assert len(windows) == 1 and len(markers) == 1 and windows[0] < markers[0], (path.name, lines)


def test_reduce_time_slices_makes_equal_windows_named_i_of_n(slicing, monkeypatch):
    """T5a: n slices are n equal windows over entry/duration, the last closed at the run's end. Each is reduced
    through reduce_time_list under the name slice_<i>of<n>, and the kinetic plot gets each window's mid-point. One
    flat result pack per slice comes back."""
    shown = []
    monkeypatch.setattr(slicing.nrtr, "plot_kinetic", lambda outputs, run, times, show=True, figure=None: shown.append(
        (len(outputs), run, list(times), show, figure)) or "figure")
    outputs, plots = slicing.nrtr.reduce_time_slices(RUN, slicing.settings, EXPERIMENT, 3, savepath=slicing.out,
                                                     datapath=slicing.nexus, show_plots=False)
    # the last window is the run's final one: its end is moved just past the last pulse, at 9 s (T2')
    assert [call["window"] for call in slicing.calls] == [(0.0, 3.0), (3.0, 6.0), (6.0, CLOSED_AT_9)]
    assert [call["subname"] for call in slicing.calls] == [
        "_slice_1of3_slice_0_3", "_slice_2of3_slice_3_6", "_slice_3of3_slice_6_9"]
    assert len(outputs) == 3 and all(isinstance(pack, list) and pack and isinstance(pack[0], dict) for pack in outputs)
    assert plots == "figure" and shown == [(3, RUN, [1.5, 4.5, 7.5], False, None)]
    assert_one_window_line_and_the_marker(slicing.out)


def test_the_figure_given_reaches_the_plot(slicing, monkeypatch):
    """T8': the figure a caller gives (the tab's own) reaches plot_kinetic through reduce_time_slices and through
    reduce_time_list."""
    drawn = []
    monkeypatch.setattr(slicing.nrtr, "plot_kinetic",
                        lambda *_args, figure=None, **_kwargs: drawn.append(figure) or figure)
    own = object()
    slicing.nrtr.reduce_time_slices(RUN, slicing.settings, EXPERIMENT, 2, savepath=slicing.out, datapath=slicing.nexus,
                                    show_plots=False, figure=own)
    slicing.nrtr.reduce_time_list(RUN, slicing.settings, EXPERIMENT, [0.0], [3.0], savepath=slicing.out,
                                  datapath=slicing.nexus, show_plots=False, figure=own)
    assert drawn == [own, own]


@pytest.mark.parametrize("n", [0, -1])
def test_reduce_time_slices_needs_at_least_one_slice(slicing, n):
    """T5a: zero or a negative number of slices is a ValueError, before anything is reduced."""
    with pytest.raises(ValueError):
        slicing.nrtr.reduce_time_slices(RUN, slicing.settings, EXPERIMENT, n, savepath=slicing.out,
                                        datapath=slicing.nexus, plot_time=False)
    assert slicing.calls == []


def test_reduce_time_list_returns_a_pack_per_window_in_order(slicing):
    """T5a (reduce_time_list): one flat result pack per window, in the order given, each reduced with its own
    window and named by it."""
    outputs, plots = slicing.nrtr.reduce_time_list(RUN, slicing.settings, EXPERIMENT, [4.0, 0.0], [6.0, 2.0],
                                                   savepath=slicing.out, datapath=slicing.nexus, plot_time=False)
    assert [call["window"] for call in slicing.calls] == [(4.0, 6.0), (0.0, 2.0)]
    assert [call["subname"] for call in slicing.calls] == ["_slice_4_6", "_slice_0_2"]
    assert len(outputs) == 2 and plots is None


@pytest.mark.parametrize("starts, ends", [([0.0, 5.0], [2.0, 7.0]), ([5.0, 0.0], [7.0, 2.0])],
                         ids=["in order", "given out of order"])
def test_nested_windows_are_one_slice_named_by_its_span(slicing, starts, ends):
    """T3 through reduce_time_list: a window entry may itself be a list of windows (the function's own comment), one
    output slice concatenating them, named by its span: the earliest start and the latest end, whatever the order the
    windows are given in. int() of the list raised (N4)."""
    slicing.nrtr.reduce_time_list(RUN, slicing.settings, EXPERIMENT, [starts], [ends],
                                  savepath=slicing.out, datapath=slicing.nexus, plot_time=False)
    assert [call["window"] for call in slicing.calls] == [(starts, ends)]
    assert [call["subname"] for call in slicing.calls] == ["_slice_0_7"]


def whole_run(path):
    """Every event id, error event and pulse charge of the run, by the file read without a window."""
    _, ids, errors, _, cpc, _ = BP.load_and_extract(path)
    return ids.tolist(), errors.tolist(), cpc.tolist()


def assert_partition(slices, whole):
    """The slices' selections partition the run in every bank: each event, each error event and each pulse's charge
    exactly once (each pulse's charge is its own multiple of 10, so a charge names its pulse)."""
    ids, errors, cpc = whole
    assert Counter(i for s in slices for i in s[0]) == Counter(ids)
    assert Counter(e for s in slices for e in s[1]) == Counter(errors)
    assert Counter(c for s in slices for c in s[2]) == Counter(cpc)


@pytest.mark.parametrize("n", [1, 4, 10])
@pytest.mark.parametrize("shape", list(REAL_SHAPED))
def test_contiguous_slices_partition_every_real_shaped_run(slicing, shape, n):
    """T2'a (B-1, B-2): reduce_time_slices' n windows over entry/duration partition the run in every bank, on every
    real-shaped run: a float32 duration below the last pulse or just past it, a tail that goes back in time, empty
    pulses, a beam-off tail, the last pulse at the duration. Each slice's selection is the real file read's, for the
    window it got."""
    path = slicing.nexus / f"REF_L_{RUN}.nxs.h5"
    write_shaped(path, shape)
    slicing.nrtr.reduce_time_slices(RUN, slicing.settings, EXPERIMENT, n, savepath=slicing.out,
                                    datapath=slicing.nexus, plot_time=False)
    assert len(slicing.calls) == n
    assert_partition([call["selected"] for call in slicing.calls], whole_run(path))


#: For each of the three banks the selection reads, a run in which that bank alone holds the latest pulse, at 9.5 s,
#: half a second after the other two banks' last (9 s), with entry/duration ending there.
LATEST_IN_ONE_BANK = {
    "error events": {"error_pulses": np.append(PULSES, 9.5)},
    "events": {"pulses": np.append(PULSES, 9.5), "error_pulses": PULSES, "charge_times": PULSES, "charge": CHARGE},
    "charge log": SHAPED_BANKS["charge log past the last event pulse"],
}


@pytest.mark.parametrize("n", [1, 4])
@pytest.mark.parametrize("bank", list(LATEST_IN_ONE_BANK))
def test_the_run_ends_at_the_latest_pulse_of_each_bank(slicing, bank, n):
    """T2'a, bank by bank (plan v3 P1, P2): the run's end is the latest pulse time over the three banks, whichever bank
    holds it. Here one bank alone holds it, where entry/duration ends: reduce_time_slices' n slices still partition the
    run in every bank, the final slice taking that bank's latest pulse (I-64 B-1: an error pulse later than every event
    pulse and charge time; B-2: an event pulse later than the charge log's last)."""
    path = slicing.nexus / f"REF_L_{RUN}.nxs.h5"
    write_run(path, duration=9.5, **LATEST_IN_ONE_BANK[bank])
    slicing.nrtr.reduce_time_slices(RUN, slicing.settings, EXPERIMENT, n, savepath=slicing.out,
                                    datapath=slicing.nexus, plot_time=False)
    assert len(slicing.calls) == n
    assert_partition([call["selected"] for call in slicing.calls], whole_run(path))


@pytest.mark.parametrize("boundary", ["the latest pulse", "the last entry"])
@pytest.mark.parametrize("shape", ["beam-off tail", "unsorted tail"])
def test_a_boundary_on_the_last_pulse_does_not_duplicate_it(slicing, shape, boundary):
    """T2'b (B-3): user windows [0, t) and [t, duration], with t the latest pulse's time (or the last entry's, which
    differ on an unsorted tail), partition the run: the pulse at t is in the second window only. Only the final window
    is closed; the first stays half-open although it ends on the last pulse."""
    path = slicing.nexus / f"REF_L_{RUN}.nxs.h5"
    pulses = write_shaped(path, shape)
    t = float(max(pulses)) if boundary == "the latest pulse" else float(pulses[-1])
    duration = float(np.float32(REAL_SHAPED[shape][2]))
    slicing.nrtr.reduce_time_list(RUN, slicing.settings, EXPERIMENT, [0.0, t], [t, duration], savepath=slicing.out,
                                  datapath=slicing.nexus, plot_time=False)
    assert len(slicing.calls) == 2
    assert_partition([call["selected"] for call in slicing.calls], whole_run(path))
    # both ends are as given: the first is not final, and the final one's end passes the last pulse already
    assert [call["window"] for call in slicing.calls] == [(0.0, t), (t, duration)]


def test_overlapping_windows_that_end_last_both_reach_the_run_end(slicing):
    """T2'/T3: two overlapping windows that end together at the run's last pulse are both final: each takes the pulses
    in its own range, the last one included, so the overlap is counted twice as T3 says, and nothing is lost."""
    slicing.nrtr.reduce_time_list(RUN, slicing.settings, EXPERIMENT, [0.0, 5.0], [9.0, 9.0], savepath=slicing.out,
                                  datapath=slicing.nexus, plot_time=False)
    assert [call["selected"][0] for call in slicing.calls] == [ids_of(range(10)), ids_of(range(5, 10))]


def test_reduce_time_list_closes_its_final_window(slicing):
    """T2' through reduce_time_list: user windows ending at the run's recorded duration, which rounds below the last
    pulse (as on REF_L_184981), partition the run: the final one takes the last pulse."""
    path = slicing.nexus / f"REF_L_{RUN}.nxs.h5"
    write_shaped(path, "float32 duration below the last pulse")
    duration = float(np.float32(LAST_184981))
    slicing.nrtr.reduce_time_list(RUN, slicing.settings, EXPERIMENT, [0.0, 40.0], [40.0, duration],
                                  savepath=slicing.out, datapath=slicing.nexus, plot_time=False)
    assert_partition([call["selected"] for call in slicing.calls], whole_run(path))


@pytest.mark.parametrize("end", ["the latest pulse", "between the latest pulse and the duration"])
def test_a_window_ending_on_the_latest_pulse_is_final_when_the_duration_passes_it(slicing, end):
    """T2', the first arm (plan v3 P3, I-64 B-3): a window reaches the run's end when its end is at or past the earlier of
    the latest pulse and entry/duration. On a run whose float32 duration passes its latest pulse (27 of the 63 real
    runs), reduce_time_list's [0, max(t)] is therefore the final window, closed just past the latest pulse: it takes
    that pulse's events, its error event and its charge, as [0, e] does for an e between the two, which is already
    past the pulse and stays as given. Each partitions the run."""
    path = slicing.nexus / f"REF_L_{RUN}.nxs.h5"
    write_shaped(path, "float32 duration just past the last pulse")
    duration = float(np.float32(PAST_184981))
    assert LAST_184981 < duration
    e = LAST_184981 if end == "the latest pulse" else (LAST_184981 + duration) / 2
    slicing.nrtr.reduce_time_list(RUN, slicing.settings, EXPERIMENT, [0.0], [e], savepath=slicing.out,
                                  datapath=slicing.nexus, plot_time=False)
    closed = float(np.nextafter(LAST_184981, np.inf)) if end == "the latest pulse" else e
    assert [call["window"] for call in slicing.calls] == [(0.0, closed)]
    assert_partition([call["selected"] for call in slicing.calls], whole_run(path))


def test_a_nested_entrys_final_window_is_closed(slicing):
    """T2'/T3: when a slice is a list of windows, the run's final window may be one of them: it is closed there too, and
    the entry stays a list."""
    slicing.nrtr.reduce_time_list(RUN, slicing.settings, EXPERIMENT, [[0.0, 5.0]], [[2.0, 9.0]], savepath=slicing.out,
                                  datapath=slicing.nexus, plot_time=False)
    assert [call["window"] for call in slicing.calls] == [([0.0, 5.0], [2.0, CLOSED_AT_9])]
    assert slicing.calls[0]["selected"][0] == ids_of([0, 1, 5, 6, 7, 8, 9])


def test_an_invalid_window_is_reported_with_the_others(slicing):
    """T4/T4c through reduce_time_list: a window that is not valid ([5, 5)) is left out of the choice of the final
    window and reported after the others have been reduced, with its window; it does not end the call."""
    with pytest.raises(ValueError) as raised:
        slicing.nrtr.reduce_time_list(RUN, slicing.settings, EXPERIMENT, [0.0, 5.0, 3.0], [3.0, 5.0, 6.0],
                                      savepath=slicing.out, datapath=slicing.nexus, plot_time=False)
    assert [call["window"] for call in slicing.calls] == [(0.0, 3.0), (3.0, 6.0)]
    assert "[5.0, 5.0)" in str(raised.value), str(raised.value)


@pytest.mark.parametrize("n", [1, 4])
def test_the_final_window_takes_every_remaining_pulse_despite_float32_duration(slicing, n):
    """T2'c (B-1): entry/duration is float32 and here rounds below the last pulse (76.78062439 against 76.780626, as on
    REF_L_184981). The final window still takes the last pulse's events, its error event and its charge."""
    path = slicing.nexus / f"REF_L_{RUN}.nxs.h5"
    pulses = write_shaped(path, "float32 duration below the last pulse")
    assert float(np.float32(LAST_184981)) < LAST_184981
    slicing.nrtr.reduce_time_slices(RUN, slicing.settings, EXPERIMENT, n, savepath=slicing.out,
                                    datapath=slicing.nexus, plot_time=False)
    last = len(pulses) - 1
    ids, errors, cpc = slicing.calls[-1]["selected"]
    assert set(ids_of([last])) <= set(ids)
    assert float(np.float32(500.0 + LAST_184981)) in errors
    assert (last + 1.0) * 10.0 in cpc


def test_a_window_that_cannot_be_reduced_is_reported_and_the_others_are_reduced(slicing):
    """T4c (the slice): a window whose reduction fails (here as the real chain fails for a window with no charge:
    _make_binary_files' RuntimeError) is reported after every other window has been reduced and written; the error
    names the window."""
    slicing.fail.append(((20.0, 30.0), RuntimeError(f"Failed to compute binary data for run {RUN}: no proton charge")))
    with pytest.raises(ValueError) as raised:
        slicing.nrtr.reduce_time_list(RUN, slicing.settings, EXPERIMENT, [0.0, 20.0, 3.0], [3.0, 30.0, 6.0],
                                      savepath=slicing.out, datapath=slicing.nexus, plot_time=False)
    assert [call["window"] for call in slicing.calls] == [(0.0, 3.0), (20.0, 30.0), (3.0, 6.0)]
    assert "[20.0, 30.0)" in str(raised.value) and "no proton charge" in str(raised.value), str(raised.value)
    written = sorted(path.name for path in slicing.out.glob("*.dat"))
    assert any("_slice_0_3" in name for name in written) and any("_slice_3_6" in name for name in written), written
    assert not any("_slice_20_30" in name for name in written), written


def test_a_slice_that_cannot_be_reduced_is_reported_after_the_others(slicing):
    """T4c through reduce_time_slices: a slice that fails does not end the call. The slices after it are reduced and
    written, and the failure is reported at the end with its window."""
    slicing.fail.append(((3.0, 6.0), RuntimeError(f"Failed to compute binary data for run {RUN}: no proton charge")))
    with pytest.raises(ValueError) as raised:
        slicing.nrtr.reduce_time_slices(RUN, slicing.settings, EXPERIMENT, 3, savepath=slicing.out,
                                        datapath=slicing.nexus, plot_time=False)
    assert [call["window"] for call in slicing.calls] == [(0.0, 3.0), (3.0, 6.0), (6.0, CLOSED_AT_9)]
    assert "[3.0, 6.0)" in str(raised.value) and "no proton charge" in str(raised.value), str(raised.value)
    written = sorted(path.name for path in slicing.out.glob("*.dat"))
    assert any("_slice_3of3_" in name for name in written) and not any("_slice_2of3_" in name for name in written)


def test_flatten_reduced_results_takes_each_shape(slicing):
    """The frame: reduce_from_file returns one dict, a list of dicts, or lists of them nested (with priors); each is
    one flat list of the dicts, and anything else is dropped."""
    flatten = slicing.nrtr.flatten_reduced_results
    a, b, c = {"Q": 1}, {"Q": 2}, {"Q": 3}
    assert flatten(a) == [a]
    assert flatten([a, b]) == [a, b]
    assert flatten([a, [b, (c,)], "not a result"]) == [a, b, c]
    assert flatten("not a result") == []


def test_a_window_with_no_reduced_data_is_reported(slicing, monkeypatch):
    """T4c: a window whose reduction returns no data is reported with its window, as one that raises is."""
    monkeypatch.setattr(slicing.nrtr.reduction, "reduce_from_file", lambda *_args, **_kwargs: ([], [], [], None))
    with pytest.raises(ValueError) as raised:
        slicing.nrtr.reduce_time_list(RUN, slicing.settings, EXPERIMENT, [0.0], [3.0], savepath=slicing.out,
                                      datapath=slicing.nexus, plot_time=False)
    assert "[0.0, 3.0)" in str(raised.value) and "no reduced result data" in str(raised.value), str(raised.value)


def test_the_priors_re_saves_name_the_window(tmp_path, slicing):
    """T6: a windowed reduction that combines with the files already in its folder (check_for_prior, as autoreduction
    does) rewrites them; every file it writes names the window. The sequence has two positions, reduced one after the
    other, so the second reduction finds the first's file and re-saves both and the combined file."""
    second = RUN + 1
    write_run(slicing.nexus / f"REF_L_{second}.nxs.h5", seq=2)
    two = {key: value * 2 for key, value in json.loads(slicing.settings.read_text()).items() if isinstance(value, list)}
    settings = settings_for(tmp_path, **two)
    for run in (RUN, second):
        nrff.reduce_from_file([run], settings, EXPERIMENT, datapath=slicing.nexus, plot=False,
                              override_params={"Spath": slicing.out}, check_for_prior=True,
                              start_times=[0.0], end_times=[3.0])
    assert len(slicing.calls) == 2 and sorted(path.name for path in slicing.out.glob("*.dat")) == [
        f"REFL_{SEQ_ID}_1_{RUN}.dat", f"REFL_{SEQ_ID}_2_{second}.dat", f"REFL_{SEQ_ID}_combined.dat"]
    assert_one_window_line_and_the_marker(slicing.out)


def test_the_kinetic_map_shows_r_and_says_so(slicing):
    """T9'a (B-5): plot_kinetic's colour map is the slices' R, row by row in slice order, under a colour bar labelled
    "R" (the offset panel beside it plots R). The contribution mapped dR there under the same label."""
    q = np.geomspace(0.01, 0.1, 5)
    outputs = [[{"Q": q, "R": (k + 1) * 1e-3 * np.ones(5), "dR": (k + 1) * 1e-5 * np.ones(5), "dQ": 0.02 * q}]
               for k in range(3)]
    pyplot = __import__("matplotlib.pyplot").pyplot
    registered = len(pyplot.get_fignums())
    figure = slicing.nrtr.plot_kinetic(outputs, RUN, times=[1.0, 2.0, 3.0], show=False)
    assert len(pyplot.get_fignums()) == registered, "a plot that is not shown leaves no pyplot figure (I-62 A-1)"
    try:
        image = next(axes.images[0] for axes in figure.axes if axes.images)
        np.testing.assert_array_equal(np.asarray(image.get_array()), [pack[0]["R"] for pack in outputs])
        bar = image.colorbar
        assert bar is not None and bar.ax.get_ylabel() == "R"
    finally:
        __import__("matplotlib.pyplot").pyplot.close(figure)


def test_the_run_file_is_resolved_the_campaigns_way(tmp_path, slicing):
    """T5b (F9): the run's NeXus file and the output folder come from the reduction's own config, never a /SNS/REF_L
    literal: an explicit datapath, else the settings' NEXUSpathRB; savepath, else the settings' Spath. The module
    holds no facility path (and no_facility_paths fails any test that reaches the config's default)."""
    settings = settings_for(tmp_path, NEXUSpathRB=str(slicing.nexus), Spath=str(slicing.out))
    slicing.nrtr.reduce_time_slices(RUN, settings, EXPERIMENT, 2, plot_time=False)
    assert {(call["nexus"], call["Spath"]) for call in slicing.calls} == {(slicing.nexus, slicing.out)}
    assert sorted(slicing.out.glob("*.dat"))
    source = Path(slicing.nrtr.__file__).read_text()
    assert "/SNS" not in source


def test_the_file_is_closed(slicing, monkeypatch):
    """T5c: reduce_time_slices reads the run's duration and closes the file. Only the module's own h5py is replaced:
    the grouping in reduce_from_file opens the run too, and is not this slug's."""
    opened = []

    def recording_file(*args, **kwargs):
        handle = h5py.File(*args, **kwargs)
        opened.append(handle)
        return handle

    monkeypatch.setattr(slicing.nrtr, "h5py", SimpleNamespace(File=recording_file))
    slicing.nrtr.reduce_time_slices(RUN, slicing.settings, EXPERIMENT, 2, savepath=slicing.out,
                                    datapath=slicing.nexus, plot_time=False)
    assert opened and not any(handle.id.valid for handle in opened)


@pytest.mark.parametrize("full", [True, False], ids=["full header", "short header"])
def test_the_header_names_the_window_beside_the_marker(reducer, full):
    """T6a: with a window, the header (both writers) has one line naming it, the contribution's text, before the
    format marker; without one there is no such line (§10 A3), so today's headers are unchanged."""
    log_values = {"title": ["t"], "ths": [0.6], "thi": [0.6], "ThCen": [0.6]}
    windowed = save_fn._build_header(reducer.config, log_values, full=full, time_window=([600.0], [1200.0]))
    lines = windowed.splitlines()
    assert [line for line in lines if line.startswith("Time resolved")] == ["Time resolved: [600.0], [1200.0]"]
    window_at = lines.index("Time resolved: [600.0], [1200.0]")
    marker_at = next(i for i, line in enumerate(lines) if line.startswith("Header format: "))
    assert window_at < marker_at
    plain = save_fn._build_header(reducer.config, log_values, full=full)
    assert "Time resolved" not in plain
    assert plain.splitlines() == [line for line in lines if not line.startswith("Time resolved")]


def test_reduce_does_not_write_windows_into_the_callers_config(slicing, reducer):
    """T6b: reduce() with a window leaves the caller's config without one (M2's R4/R5: what a call used is a record
    for its header, not a write into the shared config)."""
    reducer.config.NEXUSpathRB = slicing.nexus  # the run the stub reads the window's selection from
    reducer.reduce(plot=False, start_times=[0.0], end_times=[3.0])
    assert getattr(reducer.config, "start_times", None) is None
    assert getattr(reducer.config, "end_times", None) is None
    assert slicing.calls[-1]["window"] == ([0.0], [3.0])


def test_a_run_without_windows_writes_todays_header(slicing):
    """§10 A3 and acceptance 2: reduced without windows, a run's header has no "Time resolved" line and its Config
    line has no window keys (the window is not a config field, N3)."""
    nrff.reduce_from_file([RUN], slicing.settings, EXPERIMENT, datapath=slicing.nexus, plot=False,
                          override_params={"Spath": slicing.out})
    for path in sorted(slicing.out.glob("*.dat")):
        lines = header_lines(path)
        assert not any(line.startswith("Time resolved") for line in lines), path.name
        config = json.loads(next(line for line in lines if line.startswith("Config: "))[len("Config: "):])
        assert "start_times" not in config and "end_times" not in config, sorted(config)


# ---------------------------------------------------------------------------
# T9: progress goes to logging, not print


@pytest.mark.parametrize("module", ["new_reduction_time_resolved", "binary_processing", "direct_beam_maker"])
def test_no_print_in_library_code(module):
    """T9: the library modules the intake touches report through logging; a print() writes to whatever stdout the
    caller has (the launcher's console, an autoreduction log) and cannot be filtered."""
    path = Path(BP.__file__).with_name(f"{module}.py")
    prints = [n for n, line in enumerate(path.read_text().splitlines(), 1) if re.search(r"\bprint\(", line)]
    assert prints == [], f"{path.name}: print() at lines {prints}"


def test_the_window_filter_logs_at_info(run_file, caplog):
    """T9: the selection says what it did through the module's logger."""
    with caplog.at_level(logging.INFO, logger=BP.__name__):
        BP.load_and_extract(run_file, [2.0], [5.0])
    assert any(record.name == BP.__name__ and record.levelno == logging.INFO for record in caplog.records)
