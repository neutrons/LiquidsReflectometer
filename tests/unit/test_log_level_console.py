"""LR_REDUCTION_LOG_LEVEL: the package's records below WARNING on stdout; WARNING and above on stderr, as Python prints them.

Unset (or blank), the ``lr_reduction`` logger has two handlers. One prints its records below WARNING on stdout (the
informational messages that were ``print()`` calls), from INFO by default. The other prints WARNING and above on stderr
exactly as Python's last-resort handler does when nothing is configured (the bare message), so a process that treats
stderr output as a failure sees the package's warnings, errors and logged exceptions as before; uncaught exceptions are
untouched. A level name or a positive integer sets the level. ``NOTSET`` or ``0`` switches both off: Python's default,
WARNING and above on stderr and nothing below. Any other value uses INFO and says so once, at INFO; no value can make
``import lr_reduction`` fail.
"""

import contextlib
import io
import logging
import os
import subprocess
import sys

import pytest

import lr_reduction
from lr_reduction import configure_console_logging

HANDLER_NAME = "lr_reduction.console"
STDERR_NAME = "lr_reduction.stderr"
FORMAT = "%(levelname)s %(name)s: %(message)s"


def _console_handlers(logger_name="lr_reduction"):
    return [h for h in logging.getLogger(logger_name).handlers if h.get_name() == HANDLER_NAME]


def _stderr_handlers(logger_name="lr_reduction"):
    return [h for h in logging.getLogger(logger_name).handlers if h.get_name() == STDERR_NAME]


def _clear():
    for name in ("lr_reduction", None):
        logger = logging.getLogger(name)
        for handler in [h for h in logger.handlers if h.get_name() in (HANDLER_NAME, STDERR_NAME)]:
            logger.removeHandler(handler)
    logging.getLogger("lr_reduction").setLevel(logging.NOTSET)


def _package_state():
    """What a test may change on the package logger and this module gives back: its handlers, level, propagate and
    disabled. Not saved: its filters, and handlers on other loggers (root's handler of this name is only cleared)."""
    package = logging.getLogger("lr_reduction")
    return list(package.handlers), package.level, package.propagate, package.disabled


@contextlib.contextmanager
def _isolated_package_logger():
    """Inside: no handler of this name on the package logger or on root, and the package level NOTSET. On exit: the
    package logger's state (`_package_state`) as it was on entry."""
    package = logging.getLogger("lr_reduction")
    handlers, level, propagate, disabled = _package_state()
    _clear()
    try:
        yield
    finally:
        _clear()
        for handler in list(package.handlers):
            package.removeHandler(handler)
        for handler in handlers:
            package.addHandler(handler)
        package.setLevel(level)
        package.propagate, package.disabled = propagate, disabled


@pytest.fixture(autouse=True, scope="module")
def _module_gives_the_package_logger_back():
    """The module leaves the package logger as it found it, whatever the order of its tests: if the per-test fixture
    below stops restoring, the module ends in an error here, not in a later module that finds package logging off."""
    before = _package_state()
    yield
    assert _package_state() == before


@pytest.fixture(autouse=True)
def _reset():
    """Each test starts with no handler of this name on the package logger or on root and the package level NOTSET,
    and the package logger is given back as it was.

    Cleared before, because the import-time call has already attached a handler bound to the stdout of the moment, so a
    test run alone, or first, must not reuse it; cleared after, so nothing leaks into the next test (root too, so a
    handler attached to the wrong logger cannot leak); restored after, so the test modules that run later in the
    session see the import-time configuration, not package logging switched off.
    """
    with _isolated_package_logger():
        yield


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
    [stderr] = _stderr_handlers()
    assert stderr.stream is sys.stderr
    assert stderr.level == logging.WARNING
    assert stderr.formatter._fmt == "%(message)s"


@pytest.mark.parametrize("value", ["NOTSET", "notset", "0"])
def test_notset_or_zero_switches_off(value):
    assert configure_console_logging(value) is None
    assert _console_handlers() == [] and _stderr_handlers() == []
    assert logging.getLogger("lr_reduction").level == logging.NOTSET
    configure_console_logging("INFO")
    assert len(_console_handlers()) == 1 and len(_stderr_handlers()) == 1
    assert configure_console_logging(value) is None
    assert _console_handlers() == [] and _stderr_handlers() == []
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


@pytest.mark.parametrize("value", ["²", "①", "-5", "1.5", "verbose", " 7 x", pytest.param("9" * 5000, id="5000-digits")])
def test_invalid_values_use_info_and_say_so_once_at_info_and_never_raise(value, capsys):
    handler = configure_console_logging(value)
    assert handler.level == logging.INFO
    assert logging.getLogger("lr_reduction").level == logging.INFO
    out, err = capsys.readouterr()
    lines = out.splitlines()
    assert len(lines) == 1, lines
    assert lines[0].startswith("INFO lr_reduction: ")
    assert repr(value) in lines[0]
    for accepted in ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "positive integer", "NOTSET"):
        assert accepted in lines[0]
    assert err == ""


def test_surrounding_whitespace_is_ignored(capsys):
    handler = configure_console_logging(" debug\n")
    assert handler.level == logging.DEBUG
    assert logging.getLogger("lr_reduction").level == logging.DEBUG
    assert capsys.readouterr() == ("", "")


@pytest.mark.parametrize("value", [None, "DEBUG", "NOTSET", "verbose"])
def test_records_still_propagate(value):
    configure_console_logging(value)
    assert logging.getLogger("lr_reduction").propagate is True


def test_a_second_call_replaces_the_handler():
    first = configure_console_logging("INFO")
    [first_stderr] = _stderr_handlers()
    second = configure_console_logging("DEBUG")
    assert second is not first
    assert _console_handlers() == [second]
    assert second.level == logging.DEBUG
    [second_stderr] = _stderr_handlers()
    assert second_stderr is not first_stderr


def test_a_second_call_prints_on_the_stdout_of_the_moment(monkeypatch):
    configure_console_logging("INFO")
    now = io.StringIO()
    monkeypatch.setattr(sys, "stdout", now)
    configure_console_logging("INFO")
    logging.getLogger("lr_reduction.unit_probe").info("ping")
    assert now.getvalue().splitlines() == ["INFO lr_reduction.unit_probe: ping"]


def test_the_package_logger_is_given_back_after_each_test():
    handler = configure_console_logging("INFO")
    with _isolated_package_logger():
        assert _console_handlers() == []
        assert configure_console_logging("NOTSET") is None
    assert _console_handlers() == [handler]
    assert logging.getLogger("lr_reduction").level == logging.INFO


def test_second_call_updates_both_levels_without_a_second_handler():
    configure_console_logging("INFO")
    configure_console_logging("DEBUG")
    handlers = _console_handlers()
    assert len(handlers) == 1
    assert handlers[0].level == logging.DEBUG
    assert logging.getLogger("lr_reduction").level == logging.DEBUG


def test_a_warning_goes_to_stderr_bare_and_not_to_stdout(capsys):
    configure_console_logging("DEBUG")
    logging.getLogger("lr_reduction.unit_probe").warning("w")
    logging.getLogger("lr_reduction.unit_probe").error("e")
    out, err = capsys.readouterr()
    assert out == ""
    assert err == "w\ne\n"


def test_the_stderr_handler_writes_to_the_stderr_of_the_moment(monkeypatch):
    configure_console_logging("INFO")
    now = io.StringIO()
    monkeypatch.setattr(sys, "stderr", now)
    logging.getLogger("lr_reduction.unit_probe").warning("w")
    assert now.getvalue() == "w\n"


#: Every kind of record a reduction makes: informational, a warning, an error, a caught exception logged.
_EVERY_KIND = ("L.debug('dbg'); L.info('ping'); L.warning('w'); L.error('e')\n"
               "try:\n    1 / 0\nexcept ZeroDivisionError:\n    L.exception('x')\n")


def _fresh_process(level, body="L.info('ping'); L.warning('w')", import_package=True):
    env = {k: v for k, v in os.environ.items() if k != "LR_REDUCTION_LOG_LEVEL"}
    if level is not None:
        env["LR_REDUCTION_LOG_LEVEL"] = level
    # the child imports this checkout's package, whatever an installed copy points at
    package_root = os.path.dirname(os.path.dirname(os.path.abspath(lr_reduction.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [package_root, env.get("PYTHONPATH")]))
    imports = "import logging, lr_reduction" if import_package else "import logging"
    program = f"{imports}\nL = logging.getLogger('lr_reduction.probe')\n{body}"
    return subprocess.run([sys.executable, "-c", program], env=env, capture_output=True, text=True, timeout=120)


def _without_the_package():
    """What Python prints for the same records with no handler configured anywhere (upstream's state)."""
    return _fresh_process(None, _EVERY_KIND, import_package=False)


def test_fresh_process_unset_prints_info_on_stdout_and_the_rest_on_stderr_as_python_does():
    result = _fresh_process(None, _EVERY_KIND)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["INFO lr_reduction.probe: ping"]
    assert result.stderr == _without_the_package().stderr
    assert result.stderr.startswith("w\ne\nx\nTraceback")


def test_fresh_process_an_uncaught_exception_is_untouched():
    result = _fresh_process(None, "L.info('ping')\n1 / 0")
    assert result.returncode == 1
    assert result.stdout.splitlines() == ["INFO lr_reduction.probe: ping"]
    assert result.stderr.startswith("Traceback") and result.stderr.rstrip().endswith("ZeroDivisionError: division by zero")


@pytest.mark.parametrize("level", ["NOTSET", "0"])
def test_fresh_process_notset_is_python_s_default_on_both_streams(level):
    result = _fresh_process(level, _EVERY_KIND)
    baseline = _without_the_package()
    assert result.returncode == 0, result.stderr
    assert (result.stdout, result.stderr) == (baseline.stdout, baseline.stderr)


def test_fresh_process_at_warning_stdout_is_empty_and_stderr_is_python_s():
    result = _fresh_process("WARNING", _EVERY_KIND)
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == _without_the_package().stderr


def test_fresh_process_at_debug_debug_and_info_reach_stdout():
    result = _fresh_process("DEBUG", _EVERY_KIND)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["DEBUG lr_reduction.probe: dbg", "INFO lr_reduction.probe: ping"]
    assert result.stderr == _without_the_package().stderr


def test_fresh_process_invalid_value_imports_and_says_so_on_stdout_only():
    result = _fresh_process("²", "L.info('ping')")
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert len(lines) == 2, lines
    assert lines[0].startswith("INFO lr_reduction: ") and "'²'" in lines[0]
    assert lines[1] == "INFO lr_reduction.probe: ping"
    assert result.stderr == ""


def test_fresh_process_very_long_number_imports_and_says_so_on_stdout_only():
    result = _fresh_process("9" * 5000, "L.info('ping')")
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert len(lines) == 2, [line[:80] for line in lines]
    assert lines[0].startswith("INFO lr_reduction: ") and "'9999" in lines[0]
    assert lines[1] == "INFO lr_reduction.probe: ping"
    assert result.stderr == ""


def test_other_packages_are_untouched(capsys):
    configure_console_logging("INFO")
    logging.getLogger("not_lr_reduction").info("x")
    out, err = capsys.readouterr()
    assert out == "" and err == ""
