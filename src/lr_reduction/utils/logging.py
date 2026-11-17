"""
Custom logging setup for lr_reduction with optional log rotation and enhanced formatting.
Configures the root logger to log to both console and a rotating file,
with a custom format that includes timestamps, log levels, and source modules.

Usage
-----

Call ``setup_logging`` once, at program start (the CLI entry points already do this):

.. code-block:: python

    from lr_reduction.utils.logging import setup_logging

    setup_logging(level="INFO", rotate=True)

Then in any module, use the standard logging module to log messages:

.. code-block:: python

    import logging

    logger = logging.getLogger(__name__)
    logger.info("This is an info message.")
"""

import logging
import logging.handlers
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Literal

from lr_reduction import app_config


class OriginFormatter(logging.Formatter):
    """Display the source module for root logger records."""

    def format(self, record):
        record.origin = record.module if record.name == "root" else record.name
        return super().format(record)


def _detach_handler(handler: logging.Handler) -> None:
    """Remove a handler from the root logger and release its resources."""
    logging.getLogger().removeHandler(handler)
    handler.close()


_FMT = "%(asctime)s | %(levelname)-4s | %(origin)-12s | %(message)s"
_DATE_FMT = "%Y-%m-%d %H:%M:%S"
LOG_FMT = OriginFormatter(fmt=_FMT, datefmt=_DATE_FMT)


def setup_logging(
    level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] | None = None,
    rotate: bool | None = None,
    num_backups: int | None = None,
    log_file: str | Path | None = None,
    reset_handlers: bool = True,
):
    """Configure logging with a custom format and optional timed rotating file handler.

    Arguments
    ---------
    level: str
        The logging level to set for the root logger. Options are "DEBUG", "INFO", "WARNING", "ERROR", and "CRITICAL".
    rotate: bool
        If True, add a TimedRotatingFileHandler that rotates logs at midnight and keeps N backups.
    num_backups: int
        The number of backup files to keep. Default is 15. Only relevant if `rotate` is True.
    log_file: str | Path
        The path to the log file. Default is "~/lr_reduction.log". Only relevant if `rotate` is True.
    reset_handlers: bool
        If True, remove existing handlers from the root logger before adding new ones.
        This prevents duplicate logs if setup_logging is called multiple times.
    """

    level = app_config["logging.level"] if level is None else level
    rotate = app_config["logging.rotating_log_files"] if rotate is None else rotate
    num_backups = app_config["logging.num_backups"] if num_backups is None else num_backups
    log_file = app_config["logging.log_file"] if log_file is None else log_file

    logging.getLogger().setLevel(level)

    # Shorten level names to 4 characters for better log formatting
    logging.addLevelName(logging.DEBUG, "DBUG")
    logging.addLevelName(logging.INFO, "INFO")
    logging.addLevelName(logging.WARNING, "WARN")
    logging.addLevelName(logging.ERROR, "ERR")
    logging.addLevelName(logging.CRITICAL, "CRIT")

    if reset_handlers:
        # Remove existing handlers to prevent duplicate logs
        for handler in logging.getLogger().handlers[:]:
            _detach_handler(handler)

    # Setup stream handler
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(level)
    sh.setFormatter(LOG_FMT)
    logging.getLogger().addHandler(sh)

    if rotate:
        # Setup timed rotating file handler. It opens the file immediately, so an unwritable
        # location fails here; logging then continues to stdout only (warned about below).
        try:
            fh = logging.handlers.TimedRotatingFileHandler(
                log_file,
                when="midnight",
                backupCount=num_backups,
            )
        except OSError as e:  # PermissionError, FileNotFoundError, ... are all OSError
            logging.getLogger(__name__).warning(f"Cannot write log file {log_file}: {e};\nLogging to stdout only.")
        else:
            fh.setLevel(level)
            fh.setFormatter(LOG_FMT)
            logging.getLogger().addHandler(fh)

    # Configure Mantid logging to integrate with the Python logging system,
    # ensuring that logs from Mantid are captured and formatted consistently with the rest of the application.
    from mantid.utils.logging import log_to_python

    log_to_python(
        level=app_config["logging.mantid_log_level"].lower(),
        pattern="%s | %t",
    )


def add_log_fh(logfile: str | Path, mode: Literal["w", "a"] = "w") -> logging.Handler:
    """Add a file handler to the root logger, replacing any existing one for the same file.

    `mode` follows FileHandler semantics, typically "w" (overwrite) or "a" (append).
    """

    resolved_logfile = Path(logfile).resolve()
    resolved_logfile.parent.mkdir(parents=True, exist_ok=True)

    # Keep only one open handler per logfile so callers can toggle logging on/off.
    remove_log_fh(resolved_logfile)

    file_handler = logging.FileHandler(resolved_logfile, mode=mode)
    file_handler.setFormatter(LOG_FMT)
    logging.getLogger().addHandler(file_handler)

    return file_handler


def remove_log_fh(logfile: str | Path):
    """Remove a file handler from the root logger that matches the given logfile."""

    resolved_logfile = Path(logfile).resolve()

    for handler in logging.getLogger().handlers[:]:
        if not isinstance(handler, logging.FileHandler):
            continue

        if Path(handler.baseFilename).resolve() != resolved_logfile:
            continue

        _detach_handler(handler)


@contextmanager
def log_to_file(logfile: str | Path, mode: Literal["w", "a"] = "w") -> Iterator[None]:
    """Context manager to temporarily add a file handler to the root logger.

    The file handler is removed when the context exits, even if an exception occurs.
    """

    handler = add_log_fh(logfile, mode=mode)
    try:
        yield
    finally:
        _detach_handler(handler)
