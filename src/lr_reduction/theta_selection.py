"""
Shared theta log selection rules for Liquids Reflectometer reduction.
"""

_MISSING = object()


def _has_log(log_source, log_name: str) -> bool:
    return log_name in log_source


def _read_log(log_source, log_name: str, default=_MISSING):
    """
    Read a sample log value from either a Mantid run object or a mapping-like log source.
    """
    if hasattr(log_source, "getProperty"):
        if log_name in log_source:
            prop = log_source.getProperty(log_name)
            if hasattr(prop, "size"):
                return prop.value[-1]
            return prop.value
    elif log_name in log_source:
        return log_source[log_name]

    if default is not _MISSING:
        return default
    raise KeyError(log_name)


def uses_incident_theta(log_source) -> bool:
    """
    Return True when theta should be derived from `thi` instead of `ths`.
    """
    if _has_log(log_source, "BL4B:CS:Mode:Coordinates"):
        coordinates_mode = _read_log(log_source, "BL4B:CS:Mode:Coordinates")
        try:
            return int(coordinates_mode) == 0
        except (TypeError, ValueError):
            return False

    return _read_log(log_source, "BL4B:CS:ExpPl:OperatingMode", default="") == "Free Liquid"


def theta_log_name(log_source) -> str:
    """
    Return the sample log name that should be used for theta.
    """
    return "thi" if uses_incident_theta(log_source) else "ths"
