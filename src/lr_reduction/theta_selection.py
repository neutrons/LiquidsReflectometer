"""
Shared theta log selection rules for Liquids Reflectometer reduction.
"""

from typing import Any, Protocol

from lr_reduction.mantid_utils import SampleLogValues
from lr_reduction.types import MantidWorkspace


class SampleLogSource(Protocol):
    """Value-oriented sample log accessor."""

    def __contains__(self, log_name: str) -> bool: ...

    def __getitem__(self, log_name: str) -> Any: ...


def _sample_logs(log_source: MantidWorkspace) -> SampleLogSource:
    """
    Normalize workspace and mapping inputs to a value-based sample-log accessor.

    Parameters
    ----------
    log_source : MantidWorkspace
        Workspace or mapping-like object containing sample log values.

    Returns
    -------
    SampleLogSource
        Sample-log accessor that returns values instead of Mantid property objects.
    """
    if hasattr(log_source, "getRun") or isinstance(log_source, str):
        return SampleLogValues(log_source)
    return log_source


def is_earth_centered_geometry(log_source: MantidWorkspace) -> bool:
    """
    Determine whether the run used earth-centered geometry.

    Parameters
    ----------
    log_source : MantidWorkspace | SampleLogSource
        Workspace or mapping-like object containing the coordinate-mode and
        operating-mode logs used to determine the beamline geometry.

    Returns
    -------
    bool
        True for earth-centered geometry, False for beam-centered geometry.
        When the coordinate-mode log is unavailable, falls back to the legacy
        operating-mode log and treats ``"Free Liquid"`` as earth-centered.
    """
    sample_logs = _sample_logs(log_source)

    if "BL4B:CS:Mode:Coordinates" in sample_logs:
        coordinates_mode = sample_logs["BL4B:CS:Mode:Coordinates"]
        try:
            return int(coordinates_mode) == 0
        except (TypeError, ValueError):
            return False

    return "BL4B:CS:ExpPl:OperatingMode" in sample_logs and sample_logs["BL4B:CS:ExpPl:OperatingMode"] == "Free Liquid"


def theta_log_name(log_source: MantidWorkspace | SampleLogSource) -> str:
    """
    Return the sample log name that should be used for theta.

    Parameters
    ----------
    log_source : MantidWorkspace | SampleLogSource
        Workspace or mapping-like object containing the relevant sample log values.

    Returns
    -------
    str
        ``"thi"`` for earth-centered geometry and ``"ths"`` otherwise.
    """
    return "thi" if is_earth_centered_geometry(log_source) else "ths"
