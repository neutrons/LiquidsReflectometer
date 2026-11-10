import logging
import logging.handlers
import os
import re

import pytest

from lr_reduction.utils.logging import (
    LOG_FMT,
    add_log_fh,
    log_to_file,
    remove_log_fh,
    setup_logging,
)

_LEVELS = (logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL)


@pytest.fixture(autouse=True)
def isolate_root_logger(monkeypatch, tmp_path):
    """Give each test a clean root logger and restore global logging state afterwards."""
    root = logging.getLogger()
    saved_handlers = root.handlers[:]
    saved_level = root.level
    saved_level_names = {lvl: logging.getLevelName(lvl) for lvl in _LEVELS}

    # Keep setup_logging(rotate=True) from writing to the real ~/lr_reduction.log
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))

    for handler in saved_handlers:
        root.removeHandler(handler)

    yield

    for handler in root.handlers[:]:
        root.removeHandler(handler)
        handler.close()
    for handler in saved_handlers:
        root.addHandler(handler)
    root.setLevel(saved_level)
    for lvl, name in saved_level_names.items():
        logging.addLevelName(lvl, name)


def _file_handlers(path) -> list[logging.FileHandler]:
    resolved = os.path.abspath(path)
    return [
        h
        for h in logging.getLogger().handlers
        if isinstance(h, logging.FileHandler) and os.path.abspath(h.baseFilename) == resolved
    ]


def _stream_handlers() -> list[logging.StreamHandler]:
    # Exact type match: FileHandler and pytest's LogCaptureHandler subclass StreamHandler.
    return [h for h in logging.getLogger().handlers if type(h) is logging.StreamHandler]


def _record(name: str, pathname: str = "/pkg/some_module.py", msg: str = "hello") -> logging.LogRecord:
    return logging.LogRecord(name, logging.INFO, pathname, 1, msg, None, None)


# --------------------------------------------------------------------------- #
# setup_logging
# --------------------------------------------------------------------------- #


def test_setup_logging_reset_handlers_is_idempotent():
    setup_logging(rotate=False)
    setup_logging(rotate=False)

    assert len(_stream_handlers()) == 1


def test_setup_logging_without_reset_keeps_existing_handlers():
    existing = logging.NullHandler()
    logging.getLogger().addHandler(existing)

    setup_logging(rotate=False, reset_handlers=False)

    assert existing in logging.getLogger().handlers
    assert len(_stream_handlers()) == 1


def test_setup_logging_rotate_adds_timed_rotating_handler(tmp_path):
    setup_logging(rotate=True, num_backups=3)

    rotating = [h for h in logging.getLogger().handlers if isinstance(h, logging.handlers.TimedRotatingFileHandler)]
    assert len(rotating) == 1
    handler = rotating[0]
    assert handler.backupCount == 3
    assert handler.when == "MIDNIGHT"
    assert os.path.abspath(handler.baseFilename) == os.path.abspath(tmp_path / "lr_reduction.log")


def test_setup_logging_level_filters_output(capsys):
    setup_logging(level="WARNING", rotate=False)

    root = logging.getLogger()
    assert root.level == logging.WARNING
    assert all(h.level == logging.WARNING for h in _stream_handlers())

    root.info("hidden info message")
    root.warning("visible warning message")

    out = capsys.readouterr().out
    assert "hidden info message" not in out
    assert "visible warning message" in out


# --------------------------------------------------------------------------- #
# OriginFormatter / format
# --------------------------------------------------------------------------- #


def test_origin_formatter_uses_module_for_root_logger():
    formatted = LOG_FMT.format(_record("root", pathname="/pkg/some_module.py"))

    assert "| some_module " in formatted


def test_file_output_matches_log_format(tmp_path):
    setup_logging(rotate=False)
    logfile = tmp_path / "format.log"
    add_log_fh(str(logfile))

    logging.getLogger("my.logger").info("formatted message")

    line = logfile.read_text(encoding="utf-8").strip()
    assert re.fullmatch(
        r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \| INFO \| my\.logger\s+ \| formatted message",
        line,
    ), line


def test_add_log_fh_works_without_prior_setup(tmp_path):
    root = logging.getLogger()
    # pytest attaches its capture handlers after fixtures run; drop them to get a bare root.
    for handler in root.handlers[:]:
        root.removeHandler(handler)
    root.setLevel(logging.INFO)
    logfile = tmp_path / "no_setup.log"

    handler = add_log_fh(str(logfile))
    root.info("no setup message")

    assert handler.formatter is LOG_FMT
    assert "no setup message" in logfile.read_text(encoding="utf-8")


# --------------------------------------------------------------------------- #
# add_log_fh / remove_log_fh
# --------------------------------------------------------------------------- #


def test_add_log_fh_replaces_existing_handler_for_same_file(tmp_path):
    logfile = tmp_path / "replace.log"

    first = add_log_fh(str(logfile))
    second = add_log_fh(str(logfile))

    assert _file_handlers(logfile) == [second]
    assert first is not second
    assert first.stream is None  # closed


def test_add_log_fh_normalizes_relative_paths(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    add_log_fh("relative.log")
    add_log_fh(str(tmp_path / "relative.log"))

    assert len(_file_handlers(tmp_path / "relative.log")) == 1


@pytest.mark.parametrize(("mode", "keeps_first"), [("w", False), ("a", True)])
def test_add_log_fh_mode(tmp_path, mode, keeps_first):
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    logfile = tmp_path / f"mode_{mode}.log"

    add_log_fh(str(logfile), mode="w")
    root.info("first message")
    add_log_fh(str(logfile), mode=mode)
    root.info("second message")
    remove_log_fh(str(logfile))

    content = logfile.read_text(encoding="utf-8")
    assert ("first message" in content) is keeps_first
    assert "second message" in content


def test_remove_log_fh_only_removes_matching_handler(tmp_path):
    setup_logging(rotate=False)
    keep = add_log_fh(str(tmp_path / "keep.log"))
    drop = add_log_fh(str(tmp_path / "drop.log"))

    remove_log_fh(str(tmp_path / "drop.log"))

    handlers = logging.getLogger().handlers
    assert keep in handlers
    assert drop not in handlers
    assert drop.stream is None
    assert len(_stream_handlers()) == 1


# --------------------------------------------------------------------------- #
# log_to_file
# --------------------------------------------------------------------------- #


def test_log_to_file_adds_and_removes_file_handler(tmp_path):
    setup_logging(rotate=False)
    logfile = tmp_path / "logging_alt_context.log"

    with log_to_file(str(logfile)):
        logging.getLogger().info("message written in context")
        assert len(_file_handlers(logfile)) == 1

    assert _file_handlers(logfile) == []
    assert "message written in context" in logfile.read_text(encoding="utf-8")


def test_log_to_file_removes_handler_on_exception(tmp_path):
    setup_logging(rotate=False)
    logfile = tmp_path / "exception.log"

    with pytest.raises(RuntimeError, match="boom"), log_to_file(str(logfile)):
        handler = _file_handlers(logfile)[0]
        raise RuntimeError("boom")

    assert _file_handlers(logfile) == []
    assert handler.stream is None
