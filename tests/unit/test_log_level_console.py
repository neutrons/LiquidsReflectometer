"""LR_REDUCTION_LOG_LEVEL: the package's own log records on stdout, on by default at INFO; NOTSET switches it off.

Unset (or blank), one stdout handler on the ``lr_reduction`` logger prints every ``lr_reduction.*`` record at INFO and
above, so diagnostics reach the console and nothing reaches stderr (a process that treats any stderr output as a
failure stays green). A level name or a positive integer sets the level. ``NOTSET`` or ``0`` switches it off, which is
the old behaviour: Python's lastResort handler prints WARNING and above on stderr and drops INFO/DEBUG. Any other value
uses INFO and says so once; no value can make ``import lr_reduction`` fail.
"""

import logging
import os
import subprocess
import sys

import pytest

import lr_reduction
from lr_reduction import configure_console_logging

HANDLER_NAME = "lr_reduction.console"
FORMAT = "%(levelname)s %(name)s: %(message)s"


def _console_handlers(logger_name="lr_reduction"):
    return [h for h in logging.getLogger(logger_name).handlers if h.get_name() == HANDLER_NAME]


@pytest.fixture(autouse=True)
def _reset():
    """After each test: no handler of this name on the package logger or on root, and the package level NOTSET.

    Root is cleared too, so a handler attached to the wrong logger cannot leak from one test into the next.
    """
    yield
    for name in ("lr_reduction", None):
        logger = logging.getLogger(name)
        for handler in [h for h in logger.handlers if h.get_name() == HANDLER_NAME]:
            logger.removeHandler(handler)
    logging.getLogger("lr_reduction").setLevel(logging.NOTSET)


@pytest.mark.parametrize("value", [None, "", "  "], ids=["none", "empty", "blank"])
def test_unset_or_blank_turns_on_at_info_on_stdout(value):
    handler = configure_console_logging(value)
    assert isinstance(handler, logging.StreamHandler)
    assert handler.get_name() == HANDLER_NAME
    assert handler.stream is sys.stdout
    assert handler.formatter._fmt == FORMAT
    assert _console_handlers() == [handler]
    assert handler.level == logging.INFO
    assert logging.getLogger("lr_reduction").level == logging.INFO


@pytest.mark.parametrize("value", ["NOTSET", "notset", "0"])
def test_notset_or_zero_switches_off(value):
    assert configure_console_logging(value) is None
    assert _console_handlers() == []
    assert logging.getLogger("lr_reduction").level == logging.NOTSET
    configure_console_logging("INFO")
    assert len(_console_handlers()) == 1
    assert configure_console_logging(value) is None
    assert _console_handlers() == []
    assert logging.getLogger("lr_reduction").level == logging.NOTSET


def test_child_info_reaches_stdout_only(capsys):
    configure_console_logging("INFO")
    logging.getLogger("lr_reduction.unit_probe").info("ping")
    out, err = capsys.readouterr()
    assert out.splitlines() == ["INFO lr_reduction.unit_probe: ping"]
    assert err == ""


@pytest.mark.parametrize(
    "value, level",
    [("debug", logging.DEBUG), ("WARNING", logging.WARNING), ("WARN", logging.WARNING), ("error", logging.ERROR),
     ("10", 10), ("１０", 10)],
)
def test_level_names_and_integers(value, level):
    handler = configure_console_logging(value)
    assert handler.level == level
    assert logging.getLogger("lr_reduction").level == level


@pytest.mark.parametrize("value", ["²", "①", "-5", "1.5", "verbose", " 7 x"])
def test_invalid_values_use_info_and_warn_once_and_never_raise(value, capsys):
    handler = configure_console_logging(value)
    assert handler.level == logging.INFO
    assert logging.getLogger("lr_reduction").level == logging.INFO
    out, err = capsys.readouterr()
    lines = out.splitlines()
    assert len(lines) == 1, lines
    assert lines[0].startswith("WARNING lr_reduction: ")
    assert repr(value) in lines[0]
    for accepted in ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "positive integer", "NOTSET"):
        assert accepted in lines[0]
    assert err == ""


def test_second_call_updates_both_levels_without_a_second_handler():
    configure_console_logging("INFO")
    configure_console_logging("DEBUG")
    handlers = _console_handlers()
    assert len(handlers) == 1
    assert handlers[0].level == logging.DEBUG
    assert logging.getLogger("lr_reduction").level == logging.DEBUG


def test_a_warning_is_printed_once_on_stdout(capsys):
    configure_console_logging("INFO")
    logging.getLogger("lr_reduction.unit_probe").warning("w")
    out, err = capsys.readouterr()
    assert out.splitlines().count("WARNING lr_reduction.unit_probe: w") == 1
    assert err == ""


def _fresh_process(level):
    env = {k: v for k, v in os.environ.items() if k != "LR_REDUCTION_LOG_LEVEL"}
    if level is not None:
        env["LR_REDUCTION_LOG_LEVEL"] = level
    # the child imports this checkout's package, whatever an installed copy points at
    package_root = os.path.dirname(os.path.dirname(os.path.abspath(lr_reduction.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [package_root, env.get("PYTHONPATH")]))
    program = ("import logging, lr_reduction; L = logging.getLogger('lr_reduction.probe'); "
               "L.info('ping'); L.warning('w')")
    return subprocess.run([sys.executable, "-c", program], env=env, capture_output=True, text=True, timeout=120)


def test_fresh_process_unset_prints_on_stdout_and_nothing_on_stderr():
    result = _fresh_process(None)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["INFO lr_reduction.probe: ping", "WARNING lr_reduction.probe: w"]
    assert result.stderr == ""


def test_fresh_process_notset_is_the_old_behaviour():
    result = _fresh_process("NOTSET")
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == "w\n"


def test_fresh_process_invalid_value_imports_and_warns():
    result = _fresh_process("²")
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert len(lines) == 3, lines
    assert lines[0].startswith("WARNING lr_reduction: ") and "'²'" in lines[0]
    assert lines[1:] == ["INFO lr_reduction.probe: ping", "WARNING lr_reduction.probe: w"]
    assert result.stderr == ""


def test_other_packages_are_untouched(capsys):
    configure_console_logging("INFO")
    logging.getLogger("not_lr_reduction").info("x")
    out, err = capsys.readouterr()
    assert out == "" and err == ""
