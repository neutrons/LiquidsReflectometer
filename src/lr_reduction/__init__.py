import logging
import os
import sys

try:
    from ._version import __version__
except ImportError:
    __version__ = "unknown"

_HANDLER_NAME = "lr_reduction.console"
_STDERR_NAME = "lr_reduction.stderr"
_ACCEPTED = "DEBUG, INFO, WARNING, ERROR, CRITICAL, a positive integer, or NOTSET to switch off"


class _StderrHandler(logging.StreamHandler):
    """Writes to the ``sys.stderr`` of the moment, as Python's last-resort handler does."""

    def __init__(self, level=logging.NOTSET):
        logging.Handler.__init__(self, level)

    @property
    def stream(self):
        return sys.stderr


def configure_console_logging(value: str | None) -> logging.Handler | None:
    """Print the package's records below WARNING on stdout; leave WARNING and above on stderr, as Python prints them.

    The package's informational records (what were ``print()`` calls) go to stdout, which autoreduction does not read
    as a failure. Its warnings, errors and logged exceptions go to stderr exactly as Python's last-resort handler
    prints them when nothing is configured, the bare message; uncaught exceptions are not logging's and are unchanged.
    A handler that filters out WARNING still counts as a handler, so the last-resort handler would no longer run: the
    second handler, on stderr, is what keeps those records there.

    ``value`` is ``LR_REDUCTION_LOG_LEVEL``: unset or blank is INFO; a level name or a number sets the level;
    ``NOTSET`` or ``0`` switches both handlers off (Python's default: WARNING and above on stderr, nothing below).
    Never raises. Returns the stdout handler, or None when switched off.
    """
    package = logging.getLogger(__name__)
    for old in [h for h in package.handlers if h.get_name() in (_HANDLER_NAME, _STDERR_NAME)]:
        package.removeHandler(old)
    text = (value or "").strip().upper() or "INFO"
    # A number is read only when it is all decimal digits, every one of which int() parses, and shorter than ten
    # (no level needs more), far inside int()'s digit limit (4300 by default, never set below 640): neither call can
    # raise. Signed or underscored forms ("+5", "1_0"), which int() would accept, take the invalid path.
    level = int(text) if text.isdecimal() and len(text) < 10 else logging.getLevelNamesMapping().get(text)
    invalid = level is None
    level = logging.INFO if invalid else level
    package.setLevel(level)
    if level == logging.NOTSET:
        return None
    handler = logging.StreamHandler(sys.stdout)
    handler.set_name(_HANDLER_NAME)
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter("%(levelname)s %(name)s: %(message)s"))
    handler.addFilter(lambda record: record.levelno < logging.WARNING)
    package.addHandler(handler)
    stderr = _StderrHandler(logging.WARNING)
    stderr.set_name(_STDERR_NAME)
    stderr.setFormatter(logging.Formatter("%(message)s"))
    package.addHandler(stderr)
    if invalid:
        # A notice about the configuration, not a failure of the reduction: INFO, so it never puts a run on stderr.
        package.info("LR_REDUCTION_LOG_LEVEL=%r is not a level name (%s); using INFO", value, _ACCEPTED)
    return handler


configure_console_logging(os.environ.get("LR_REDUCTION_LOG_LEVEL"))
