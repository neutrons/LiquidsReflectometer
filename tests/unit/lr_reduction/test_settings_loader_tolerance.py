"""The settings loader reads a file written by any version of the code.

A settings JSON, and the ``# Config:`` line of a reduced ``.dat``, are written from ``config.__dict__``: every
attribute the config had when the file was written is a key of the file. A version that has since dropped an
attribute (or has not gained it yet) used to refuse the whole file with ``AttributeError: <key> is not a valid config
parameter``. A key this version does not have is now skipped, with one warning that names the key, its value and the
file. ``override_params`` stays strict: those are typed by the caller in code, and a typo there should fail at once.
"""

import json
import logging

import pytest

import lr_reduction.new_reduction_from_file as nrff
from lr_reduction.nr_reduction_calc import NR_Reduction
from lr_reduction.nr_reduction_config import NRReductionConfig

#: A key no version of NRReductionConfig has had.
UNKNOWN = "a_setting_of_another_version"

LOADER_LOGGER = "lr_reduction.new_reduction_from_file"


def settings(**extra):
    """One sequence position, enough for NR_Reduction to accept the config."""
    return {"Sname": "from_file", "qmin": 0.008, "DBname": ["DB_A1.dat"], "RB_Ymin": [146], "RB_Ymax": [155], **extra}


def write_json(tmp_path, values):
    path = tmp_path / "REFL_1_settings.json"
    path.write_text(json.dumps(values, indent=2))
    return path


def write_dat(tmp_path, values):
    path = tmp_path / "REFL_1_combined.dat"
    path.write_text(f"# NR_runs = [1]\n# Config: {json.dumps(values)}\n# columns = Q, R, dR, dQ (sigma)\n"
                    "0.01 1.0 0.1 0.001\n0.02 0.5 0.05 0.002\n")
    return path


def warnings_of(caplog):
    return [record for record in caplog.records if record.name == LOADER_LOGGER and record.levelno == logging.WARNING]


@pytest.mark.parametrize("write", [write_json, write_dat], ids=["settings json", "dat header"])
def test_a_file_with_a_key_this_version_lacks_loads(tmp_path, caplog, write):
    path = write(tmp_path, settings(**{UNKNOWN: [1, 2]}))
    with caplog.at_level(logging.WARNING, logger=LOADER_LOGGER):
        config = nrff.json_to_config(nrff.load_from_file(path)["config"], source=path)
    assert isinstance(config, NRReductionConfig)
    assert (config.Sname, config.qmin, config.DBname) == ("from_file", 0.008, ["DB_A1.dat"])
    assert not hasattr(config, UNKNOWN)
    (warning,) = warnings_of(caplog)
    message = warning.getMessage()
    assert UNKNOWN in message and "[1, 2]" in message and str(path) in message


def test_each_skipped_key_is_named_once(caplog):
    with caplog.at_level(logging.WARNING, logger=LOADER_LOGGER):
        nrff.json_to_config(settings(first_gone=1, second_gone="two"))
    messages = [record.getMessage() for record in warnings_of(caplog)]
    assert len(messages) == 2
    assert sum("first_gone" in message for message in messages) == 1
    assert sum("second_gone" in message for message in messages) == 1


def test_a_file_of_this_version_loads_without_a_word(caplog):
    with caplog.at_level(logging.WARNING, logger=LOADER_LOGGER):
        config = nrff.json_to_config(dict(NRReductionConfig().__dict__))
    assert config.__dict__ == NRReductionConfig().__dict__
    assert warnings_of(caplog) == []


def test_a_key_that_names_a_read_only_property_is_skipped_too(caplog):
    """``base_path`` passes ``hasattr`` and has no setter: it is a key the config cannot take, like any other."""
    with caplog.at_level(logging.WARNING, logger=LOADER_LOGGER):
        config = nrff.json_to_config(settings(base_path="/elsewhere"))
    assert str(config.base_path) != "/elsewhere"
    (warning,) = warnings_of(caplog)
    assert "base_path" in warning.getMessage()


class StopError(Exception):
    """Raised by the stand-in reducer: the call got as far as reducing."""


@pytest.fixture
def one_group(monkeypatch):
    """reduce_from_file up to the reduction itself, without a NeXus file: one run at sequence position 1."""
    monkeypatch.setattr(nrff, "group_runs", lambda *_args: {"run_nums": [[1]], "seq_nums": [[1]], "seq_ids": [[1]]})
    received = []

    def stand_in(self, *_args, **_kwargs):
        received.append(self.config)
        raise StopError

    monkeypatch.setattr(NR_Reduction, "reduce", stand_in)
    return received


@pytest.mark.parametrize("write", [write_json, write_dat], ids=["settings json", "dat header"])
def test_a_reduction_from_such_a_file_runs_and_names_the_file(tmp_path, caplog, one_group, write):
    path = write(tmp_path, settings(**{UNKNOWN: 3}))
    with caplog.at_level(logging.WARNING, logger=LOADER_LOGGER), pytest.raises(StopError):
        nrff.reduce_from_file([1], path, "IPTS-00001", datapath=tmp_path, plot=False)
    (config,) = one_group
    assert config.qmin == 0.008
    (warning,) = warnings_of(caplog)
    assert UNKNOWN in warning.getMessage() and str(path) in warning.getMessage()


def test_the_time_resolved_folders_resolve_from_such_a_file(tmp_path):
    import lr_reduction.new_reduction_time_resolved as nrtr

    path = write_json(tmp_path, settings(**{UNKNOWN: 3}))
    nexus, out = nrtr.run_folders(path, "IPTS-00001", datapath=tmp_path / "nexus", savepath=tmp_path / "out")
    assert (nexus, out) == (tmp_path / "nexus", tmp_path / "out")


def test_an_override_the_config_lacks_still_raises(tmp_path, one_group):
    """Typed by the caller, in code: a typo fails at once, before anything is reduced."""
    path = write_json(tmp_path, settings())
    with pytest.raises(AttributeError, match=f"{UNKNOWN} is not a valid config parameter"):
        nrff.reduce_from_file([1], path, "IPTS-00001", datapath=tmp_path, plot=False, override_params={UNKNOWN: 3})
    assert one_group == []


def test_the_loader_hands_back_what_it_did_not_use():
    """A caller that shows the settings to someone (the settings editor) has to say what the file held and the
    config does not: the keys skipped, with their values, in the dict it passes as ``unused``."""
    unused = {}
    nrff.json_to_config(settings(first_gone=1, base_path="/elsewhere"), unused=unused)
    assert unused == {"first_gone": 1, "base_path": "/elsewhere"}
    nothing = {}
    nrff.json_to_config(settings(), unused=nothing)
    assert nothing == {}


# ---------------------------------------------------------------------------
# The time window is not a setting (PR #223 review, D1 option A). The add-time-slicing branch (upstream PR #205) held
# it in the config and wrote it into every file: null for a run reduced whole, the window for a slice.

WINDOW_KEYS = ("start_times", "end_times")


def test_the_config_has_no_window_fields():
    """The window is an argument of a reduction call and a line of the header it writes, never a config attribute:
    an attribute is a key of every file written, and a key of a file is read back as a setting."""
    assert not set(WINDOW_KEYS) & set(NRReductionConfig().__dict__)


@pytest.mark.parametrize("write", [write_json, write_dat], ids=["settings json", "dat header"])
def test_a_file_that_recorded_no_window_loads_as_if_the_keys_were_absent(tmp_path, caplog, write):
    """A whole-run file of that branch holds both keys as null: nothing to say, nothing ignored."""
    path = write(tmp_path, settings(start_times=None, end_times=None))
    with caplog.at_level(logging.DEBUG, logger=LOADER_LOGGER):
        config = nrff.json_to_config(nrff.load_from_file(path)["config"], source=path)
    assert config.qmin == 0.008 and not set(WINDOW_KEYS) & set(config.__dict__)
    assert warnings_of(caplog) == []


@pytest.mark.parametrize("window", [(120, 136), ([0.0], [3.0]), ([[0, 5], [2, 7]], None)],
                         ids=["numbers", "lists", "one of the two"])
def test_a_file_that_recorded_a_window_loads_and_says_the_window_is_not_applied(tmp_path, caplog, window):
    """A slice's file of that branch holds the window it was reduced with. Read as settings it reduces the whole
    run, as it did on that branch (nothing there read the keys back), and the one warning says so and names the
    arguments that do reduce a window."""
    start, end = window
    path = write_dat(tmp_path, settings(start_times=start, end_times=end))
    with caplog.at_level(logging.WARNING, logger=LOADER_LOGGER):
        config = nrff.json_to_config(nrff.load_from_file(path)["config"], source=path)
    assert not set(WINDOW_KEYS) & set(config.__dict__)
    (warning,) = warnings_of(caplog)
    message = warning.getMessage()
    assert str(path) in message and repr(start) in message and repr(end) in message
    assert "whole run" in message and "start_times=" in message and "end_times=" in message


def test_a_reduction_from_a_slices_file_is_handed_no_window(tmp_path, caplog, monkeypatch):
    """End to end: the recorded window of the settings file reaches neither the config nor the reduction call."""
    monkeypatch.setattr(nrff, "group_runs", lambda *_args: {"run_nums": [[1]], "seq_nums": [[1]], "seq_ids": [[1]]})
    received = []

    def stand_in(self, *_args, **kwargs):
        received.append((dict(self.config.__dict__), kwargs))
        raise StopError

    monkeypatch.setattr(NR_Reduction, "reduce", stand_in)
    path = write_dat(tmp_path, settings(start_times=120, end_times=136))
    with caplog.at_level(logging.WARNING, logger=LOADER_LOGGER), pytest.raises(StopError):
        nrff.reduce_from_file([1], path, "IPTS-00001", datapath=tmp_path, plot=False)
    ((config, kwargs),) = received
    assert not set(WINDOW_KEYS) & set(config)
    assert kwargs["start_times"] is None and kwargs["end_times"] is None
    assert len(warnings_of(caplog)) == 1


@pytest.mark.parametrize("key", WINDOW_KEYS)
def test_a_window_among_the_overrides_raises(tmp_path, one_group, key):
    """override_params={"start_times": ...} is a call that means to slice and would not: it fails at once."""
    path = write_json(tmp_path, settings())
    with pytest.raises(AttributeError, match=f"{key} is not a valid config parameter"):
        nrff.reduce_from_file([1], path, "IPTS-00001", datapath=tmp_path, plot=False, override_params={key: [0.0]})
    assert one_group == []


def test_a_recorded_window_is_handed_back_as_not_used():
    """The settings editor has to say that the file names a window and that it is not applied."""
    unused = {}
    nrff.json_to_config(settings(start_times=120, end_times=136), unused=unused)
    assert unused == {"start_times": 120, "end_times": 136}
    nothing = {}
    nrff.json_to_config(settings(start_times=None, end_times=None), unused=nothing)
    assert nothing == {}


# ---------------------------------------------------------------------------
# From the test review of this change: the cases its first tests let through.


@pytest.mark.parametrize("values, handed_back", [
    ({"start_times": 0, "end_times": None}, {"start_times": 0, "end_times": None}),
    ({"start_times": 0}, {"start_times": 0}),
    ({"end_times": [3.0]}, {"end_times": [3.0]}),
], ids=["zero and null", "start only", "end only"])
def test_half_a_recorded_window_is_still_a_window(caplog, values, handed_back):
    """A start of 0 is a window (not "nothing"), and one key without the other is still named: one warning, and the
    caller is handed back what the file held, the null half included."""
    unused = {}
    with caplog.at_level(logging.WARNING, logger=LOADER_LOGGER):
        config = nrff.json_to_config(settings(**values), unused=unused)
    assert unused == handed_back
    assert len(warnings_of(caplog)) == 1
    assert not set(WINDOW_KEYS) & set(config.__dict__)


def test_the_warning_is_logged_whether_or_not_the_caller_collects(caplog):
    with caplog.at_level(logging.WARNING, logger=LOADER_LOGGER):
        nrff.json_to_config(settings(first_gone=1), unused={})
    assert len(warnings_of(caplog)) == 1


def test_the_time_resolved_folders_name_the_file_they_read(tmp_path, caplog):
    import lr_reduction.new_reduction_time_resolved as nrtr

    path = write_json(tmp_path, settings(**{UNKNOWN: 3}))
    with caplog.at_level(logging.WARNING, logger=LOADER_LOGGER):
        nrtr.run_folders(path, "IPTS-00001", datapath=tmp_path)
    (warning,) = warnings_of(caplog)
    assert str(path) in warning.getMessage()


@pytest.mark.parametrize("not_settings", [None, ["Sname", "x"], "Sname"], ids=["none", "a list", "text"])
def test_what_is_not_a_mapping_of_settings_is_refused_by_name(tmp_path, not_settings):
    """A .dat without a "# Config:" line gives None, a JSON file may hold a list: neither is a settings file of any
    version, and the error says so and names the file (it was a TypeError or a bare AttributeError)."""
    with pytest.raises(ValueError, match="holds no reduction settings") as raised:
        nrff.json_to_config(not_settings, source=tmp_path / "REFL_1.dat")
    assert "REFL_1.dat" in str(raised.value)


# load_from_file reads to the end looking for the line, then hands loadtxt the exhausted file: its warning, not this change's
@pytest.mark.filterwarnings("ignore:loadtxt. input contained no data:UserWarning")
def test_a_dat_without_a_config_line_is_refused_by_name(tmp_path, one_group):
    path = tmp_path / "bare.dat"
    path.write_text("# NR_runs = [1]\n0.01 1.0 0.1 0.001\n0.02 0.5 0.05 0.002\n")
    with pytest.raises(ValueError, match="holds no reduction settings"):
        nrff.reduce_from_file([1], path, "IPTS-00001", datapath=tmp_path, plot=False)
    assert one_group == []
