from enum import IntEnum

import numpy as np

from lr_reduction.typing import MantidWorkspace
from lr_reduction.utils import workspace_handle


def _log_value(run, log_name: str, default=None):
    """
    Extract a value from a Mantid workspace run log.

    Parameters
    ----------
    run : mantid.api.Run
        The Mantid workspace run object containing log properties.
    log_name : str
        The name of the log property to retrieve from the workspace run.
    default : optional
        Default value to return if the log is not found. If None and the log
        is not found, a RuntimeError is raised.

    Returns
    -------
    The value from the log property. For log series (int series, float series,
    string series), returns the first value. For single values, returns the
    value directly.

    Raises
    ------
    RuntimeError
        If the log is not found and no default value is provided.
    """
    if log_name in run:
        p = run.getProperty(log_name)
        if hasattr(p, "size"):  # a log series (int series, float series, string series)
            # logs can have some pre-run data in them at the start
            # and the last value is typically more reliable as a reference for the run
            return p.value[-1]
        else:
            return p.value
    elif default is not None:
        return default
    else:
        raise RuntimeError(f"Log {log_name} not found")


class GravityDirection(IntEnum):

    DOWN = -1
    OFF = 0
    UP = 1


    @classmethod
    def from_value(cls, value):
        """Create GravityDirection from value, returning None if value is None"""
        return None if value is None else cls(value)

    @classmethod
    def find_direction(cls, workspace: MantidWorkspace) -> 'GravityDirection':
        """(docstring here)"""

        try:
            ws = workspace_handle(workspace)
            run = ws.getRun()
            assert ws.getInstrument().getName() == "REF_L", "Gravity direction can only be determined for REF_L"
            # Determine whether reflect up or down based on the sign of sample angle (ths) as don't have other flag for this.
            # For theta-theta geometry ths=0.0 and is reflect up so include 0.0 as if positive value of ths.
            # This is a workaround without an alternative flag.
            ths_val_RB = _log_value(run, "BL4B:Mot:ths.RBV")
            if ths_val_RB < -0.001:
                return cls.DOWN
            else:
                return cls.UP
        except Exception:  # noqa: BLE001
            return cls.OFF


def _theta_in(workspace: MantidWorkspace) -> float:
    """
    Calculate the incident angle (theta_in) for a given Mantid workspace.

    Parameters
    ----------
    workspace : MantidWorkspace
        The Mantid workspace from which to extract metadata for the calculation.

    Returns
    -------
    float
        The calculated incident angle (theta_in) in degrees.
    """
    ws = workspace_handle(workspace)
    run = ws.getRun()

    # Angle calculated from thi and a flag on earth-centered vs beam-centered
    thi = _log_value(run, "BL4B:Mot:thi.RBV")
    if "BL4B:CS:Mode:Coordinates" in run:
        if _log_value(run, "BL4B:CS:Mode:Coordinates") == 0: # Earth-centered=0
            theta_in = thi
        else:
            theta_in = thi - 4.0  # Beamline optics gives -4 deg incline. In future will have PV.
    elif _log_value(run, "BL4B:CS:ExpPl:OperatingMode", default="") == "Free Liquid":
        theta_in = thi
    else:
        theta_in = thi - 4.0
    return abs(theta_in)


def _theta_sample(workspace: MantidWorkspace, wavelengths: np.ndarray, theta_in: float) -> np.ndarray:
    """
    Calculate the sample angle corrected for gravity effects.

    This function computes the gravity-corrected angle at the sample position for
    neutron events with given wavelengths. The calculation accounts for the neutron
    beam path through the instrument slits and the gravitational deflection based on
    the ILL paper methodology for inclined beams.

    Parameters
    ----------
    workspace : MantidWorkspace
        The Mantid workspace containing instrument geometry and run metadata.
    wavelengths : np.ndarray
        Array of neutron wavelengths in Angstroms for each event.
    theta_in : float
        The incident angle in degrees.

    Returns
    -------
    np.ndarray
        Array of gravity-corrected sample angles in degrees, one for each wavelength.

    Notes
    -----
    The calculation uses:
    - Instrument parameters for slit positions and distances
    - Neutron velocity calculated from de Broglie wavelength
    - Parabolic trajectory correction under gravitational acceleration
    - The sample position is defined as the coordinate origin (x=0, y=0)

    The algorithm determines the parabolic path of neutrons through the instrument
    slits and calculates the corrected angle at the sample position using the
    derivative of the parabolic trajectory.
    """
    ws = workspace_handle(workspace)
    instrument = ws.getInstrument()
    run = ws.getRun()

    def instrument_parameter(parameter_name: str, default):
        if instrument.hasParameter(parameter_name):
            return instrument.getNumberParameter(parameter_name)[0]
        else:
            return default

    # Xi reference would be the position of xi if the si slit were to be positioned
    # at the sample. The distance from the sample to si is then xi_reference - xi.
    # Default reference value based on instrument setup. xi would read 445 at the sample point.
    xi_reference = instrument_parameter("xi-reference", default=445)  # in mm
    # Distance between slit s1 and the sample
    s1_sample_distance = instrument_parameter("s1-sample-distance", default=1.485) * 1000  # in mm
    # Slit i position given by xi motor. Reference value is xi motor readback at the sample position.
    xi = _log_value(run, "BL4B:Mot:xi.RBV", default=310) # in mm
    sample_si_distance = xi_reference - xi
    slit_distance = s1_sample_distance - sample_si_distance

    # Calculation from the ILL paper. This works for inclined beams.
    # Calculated theta is the angle on the sample

    g = 9.8067  # m/s^2
    h = 6.6260715e-34  # Js=kg m^2/s
    mn = 1.67492749804e-27  # kg

    v = h / (mn * wavelengths * 1e-10)
    k = g / (2 * v**2)

    # Define the sample position as x=0, y=0. increasing x is towards moderator
    xs = 0

    # positions of slits, in meters
    x1 = sample_si_distance / 1000
    x2 = (sample_si_distance + slit_distance) / 1000

    # height of slits determined by incident theta, y=0 is the sample height
    y1 = x1 * np.tan(theta_in * np.pi / 180)
    y2 = x2 * np.tan(theta_in * np.pi / 180)

    # This is the location of the top of the parabola
    x0 = (y1 - y2 + k * (x1**2 - x2**2)) / (2 * k * (x1 - x2))
    y0 = y2 + k * (x2 - x0)**2

    # Shift in x for where neutron hits the sample plane
    xs = x0 - np.sqrt(y0 / k)

    # Angle is arctan(dy/dx) at sample
    return np.arctan(2 * k * (x0 - xs)) * 180 / np.pi


def gravity_correction(workspace: MantidWorkspace,
                       wavelengths: np.ndarray,
                       theta_in: float = None,
                       gravity_direction: GravityDirection = None) -> np.ndarray:
    """
    Gravity correction for each event

    Parameters
    ----------
    workspace : mantid.api.Workspace
        Mantid workspace to extract correction meta-data from.
    wavelengths : np.ndarray
        Array of wavelengths for each event.
    theta_in : float
        theta angle.. If not specified, it will be calculated with the metadata in the workspace.
    gravity_direction : GravityDirection
        Specify gravity direction. If not specified, it will be determined from the workspace metadata.


    Returns
    -------
    Array of gravity-corrected theta values for each event, in radians.
    """
    if theta_in is None:
        theta_in = _theta_in(workspace)
    theta_sample = _theta_sample(workspace, wavelengths, theta_in)
    if gravity_direction is None:
        gravity_direction = GravityDirection.find_direction(workspace)
    return int(gravity_direction)  * (theta_sample - theta_in) * np.pi / 180.0
