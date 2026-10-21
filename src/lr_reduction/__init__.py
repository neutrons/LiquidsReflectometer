import logging
import os
import sys

try:
    from ._version import __version__
except ImportError:
    __version__ = "unknown"

_HANDLER_NAME = "lr_reduction.console"
_ACCEPTED = "DEBUG, INFO, WARNING, ERROR, CRITICAL, a positive integer, or NOTSET to switch off"


def configure_console_logging(value: str | None) -> logging.Handler | None:
    """Print the package's own log records on stdout (autoreduction treats anything on stderr as a failure).

    ``value`` is ``LR_REDUCTION_LOG_LEVEL``: unset or blank is INFO; a level name or a number sets the level;
    ``NOTSET`` or ``0`` switches it off (Python's default: WARNING and above on stderr). Never raises.
    """
    package = logging.getLogger(__name__)
    for old in [h for h in package.handlers if h.get_name() == _HANDLER_NAME]:
        package.removeHandler(old)
    text = (value or "").strip().upper() or "INFO"
    # str.isdecimal() is exactly what int() accepts as an unsigned number, and int() refuses more than 4300 digits,
    # so a number is read only when it is short (no level needs ten digits): neither call can raise.
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
    package.addHandler(handler)
    if invalid:
        package.warning("LR_REDUCTION_LOG_LEVEL=%r is not a level name (%s); using INFO", value, _ACCEPTED)
    return handler


configure_console_logging(os.environ.get("LR_REDUCTION_LOG_LEVEL"))
