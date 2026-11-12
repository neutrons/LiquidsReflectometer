"""
Custom logging setup for lr_reduction with optional log rotation and enhanced formatting.
Configures the root logger to log to both console and a rotating file,
with a custom format that includes timestamps, log levels, and source modules.

Usage:

In __init__.py or main.py
```python
from python_project.utils.logging_alt import setup_logging
setup_logging(level="DEBUG", rotate=True)
```

Then in any module, use the standard logging module to log messages:
```python
import logging
logger = logging.getLogger(__name__)
logger.info("This is an info message.")
```
"""

import logging
import logging.handlers
import os
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from typing import Literal


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
    level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO",
    rotate: bool = True,
    num_backups: int = 15,
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
    reset_handlers: bool
        If True, remove existing handlers from the root logger before adding new ones.
        This prevents duplicate logs if setup_logging is called multiple times.
    """
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

    if rotate:
        # Setup timed rotating file handler
        fh = logging.handlers.TimedRotatingFileHandler(
            os.path.join(os.path.expanduser("~"), "lr_reduction.log"),
            when="midnight",
            backupCount=num_backups,
        )
        fh.setLevel(level)
        fh.setFormatter(LOG_FMT)
        logging.getLogger().addHandler(fh)

    # Setup stream handler
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(level)
    sh.setFormatter(LOG_FMT)
    logging.getLogger().addHandler(sh)


def add_log_fh(logfile: str, mode: Literal["w", "a"] = "w") -> logging.Handler:
    """Add a file handler to the root logger, replacing any existing one for the same file.

    `mode` follows FileHandler semantics, typically "w" (overwrite) or "a" (append).
    """

    resolved_logfile = os.path.abspath(logfile)
    os.makedirs(os.path.dirname(resolved_logfile), exist_ok=True)

    # Keep only one open handler per logfile so callers can toggle logging on/off.
    remove_log_fh(resolved_logfile)

    file_handler = logging.FileHandler(resolved_logfile, mode=mode)
    file_handler.setFormatter(LOG_FMT)
    logging.getLogger().addHandler(file_handler)

    return file_handler


def remove_log_fh(logfile: str):
    """Remove a file handler from the root logger that matches the given logfile."""

    resolved_logfile = os.path.abspath(logfile)

    for handler in logging.getLogger().handlers[:]:
        if not isinstance(handler, logging.FileHandler):
            continue

        if os.path.abspath(handler.baseFilename) != resolved_logfile:
            continue

        _detach_handler(handler)


@contextmanager
def log_to_file(logfile: str, mode: Literal["w", "a"] = "w") -> Iterator[None]:
    """Context manager to temporarily add a file handler to the root logger.

    The file handler is removed when the context exits, even if an exception occurs.
    """

    handler = add_log_fh(logfile, mode=mode)
    try:
        yield
    finally:
        _detach_handler(handler)
