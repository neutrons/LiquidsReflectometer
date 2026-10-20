import logging
import os
import sys

try:
    from ._version import __version__
except ImportError:
    __version__ = "unknown"

_HANDLER_NAME = "lr_reduction.console"
_FORMAT = "%(levelname)s %(name)s: %(message)s"
_ACCEPTED = "DEBUG, INFO, WARNING, ERROR, CRITICAL, a positive integer, or NOTSET to switch off"


def _parse_level(text: str) -> tuple[int | None, bool]:
    """``(level, invalid)`` for a stripped, non-empty ``LR_REDUCTION_LOG_LEVEL``; ``level`` is None for "off".

    Every ``str`` lands in exactly one branch, so nothing here raises: what ``int()`` accepts is a number (0 is
    NOTSET's number: off; a negative is invalid); otherwise ``NOTSET`` is off and a name ``logging`` knows is that
    level; anything else is invalid.
    """
    try:
        number = int(text)
    except ValueError:
        number = None
    if number is not None:
        if number == 0:
            return None, False
        return (number, False) if number > 0 else (logging.INFO, True)
    name = text.upper()
    if name == "NOTSET":
        return None, False
    level = logging.getLevelName(name)
    if isinstance(level, int) and level > 0:
        return level, False
    return logging.INFO, True


def configure_console_logging(value: str | None) -> logging.Handler | None:
    """Print the package's own log records on stdout, at INFO unless ``value`` says otherwise.

    ``value`` is what ``LR_REDUCTION_LOG_LEVEL`` holds. Unset (``None``), empty or blank: on at INFO. A level name (any
    case) or a positive integer: on at that level. ``NOTSET`` or ``0``: off, which is Python's default (its lastResort
    handler prints WARNING and above on stderr); a handler attached earlier is removed. Any other value: on at INFO,
    and one WARNING says so. One handler, named ``lr_reduction.console``, on the ``lr_reduction`` logger that every
    ``lr_reduction.*`` logger propagates to; no other logger is touched. A second call updates its level.
    """
    package = logging.getLogger(__name__)
    handler = next((h for h in package.handlers if h.get_name() == _HANDLER_NAME), None)
    text = "" if value is None else value.strip()
    level, invalid = (logging.INFO, False) if not text else _parse_level(text)
    if level is None:
        if handler is not None:
            package.removeHandler(handler)
        package.setLevel(logging.NOTSET)
        return None
    if handler is None:
        handler = logging.StreamHandler(sys.stdout)
        handler.set_name(_HANDLER_NAME)
        handler.setFormatter(logging.Formatter(_FORMAT))
        package.addHandler(handler)
    handler.setLevel(level)
    package.setLevel(level)
    if invalid:
        package.warning("LR_REDUCTION_LOG_LEVEL=%r is not a level name (%s); using INFO", value, _ACCEPTED)
    return handler


configure_console_logging(os.environ.get("LR_REDUCTION_LOG_LEVEL"))
