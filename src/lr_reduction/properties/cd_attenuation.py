"""Cd-foil transmission correction for direct-beam runs.

A direct-beam run is measured through a stack of Cd foils chosen to keep the detector on-scale.
Dividing each event's weight by the stack's transmission, following the Beer-Lambert law
T(λ) = exp(-μ(λ)·d), recovers the unattenuated incident flux the sample would see.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
from mantid.api import IEventWorkspace, MatrixWorkspace
from mantid.simpleapi import ConvertToHistogram, CreateWorkspace, Divide, FilterByXValue

from lr_reduction.exceptions import LogNotFoundError, LrValidationError, SampleLogsError
from lr_reduction.types import MantidWorkspace, MantidWorkspaceName
from lr_reduction.utils.sample_logs import SampleLogs
from lr_reduction.utils.workspace import workspace_exists, workspace_handle

#######################################################################################
### Constants

# Path to the Cd attenuation table, in wavelength (Å) vs. linear attenuation coefficient μ (cm^-1).
CD_ATTENUATION_FILE = Path(__file__).parent / "Cd-attenuation.csv"

# Cumulative Cd thickness in microns that each attenuator stage puts in the beam, in log order.
CD_FOILS = (57.5, 126.5, 126.5 + 123.0, 2 * (126.5 + 123.0))

# Sample log names for the attenuator flags, in order of the stages in `CD_FOILS`. The Menu log is preferred,
MENU_LOG = "BL4B:Actuator:Menu"
STAGE_LOGS = [f"BL4B:Actuator:{stage}MRb" for stage in ("50", "100", "200", "400")]
THICKNESS_LOG = "cd_thickness"

# Maximum (mu*d) before the transmission correction overflows to inf,
# since Mantid stores event weights and squared errors as float32
_float32_exponent_limit = 0.5 * np.log(np.finfo(np.float32).max)  # ≈ 44.36
_headroom = np.log(1e3)  # input weights >1, later sums
_MAX_ATTENUATION_EXPONENT = _float32_exponent_limit - _headroom  # ≈ 37.4

#######################################################################################
### Private helpers


@lru_cache(maxsize=1)
def _load_cd_attenuation_data() -> tuple[np.ndarray, np.ndarray]:
    """Read the Cd attenuation table.

    Cached, so the file is read once per process; the arrays are read-only because every caller
    shares them. A failed read raises and is not cached, so it is retried on the next call.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The table's wavelengths in Angstrom and the Cd linear attenuation coefficient μ at each,
        in cm^-1.
    """
    wavelength, mu = np.loadtxt(CD_ATTENUATION_FILE, delimiter="\t", skiprows=1, unpack=True)
    wavelength.setflags(write=False)
    mu.setflags(write=False)
    return wavelength, mu


def _is_transmission_for(name: MantidWorkspaceName, thickness: float) -> bool:
    """Whether the workspace of this name is the transmission this module built for `thickness`.

    The name is shared, process-global state: a workspace written over it — an `apply_correction`
    output, say — or a rounded name shared by two thicknesses must be rebuilt, not reused.
    """
    if not workspace_exists(name):
        return False
    ws = workspace_handle(name)
    if not isinstance(ws, MatrixWorkspace) or not ws.isHistogramData():
        return False
    logs = SampleLogs(ws)
    return THICKNESS_LOG in logs and logs[THICKNESS_LOG] == thickness


def _overflow_cutoff_wavelength(ws: IEventWorkspace, cd_thickness: float, bin_edges: np.ndarray) -> float | None:
    """Check that a workspace's events are within the Cd attenuation table and won't overflow the float32 weight limit.

    Mantid stores weighted-event weights and squared errors as float32.

    Returns
    -------
    cutoff: float | None
        The wavelength in Angstrom just below the first bin where the exponent exceeds the safe limit,
        or None if all are safe. Events above it must be dropped before the correction.

    Raises
    ------
    LrValidationError
        If all events are outside the Cd attenuation table, or if every event's wavelength would cause the
        transmission correction to overflow the float32 event weights.
    """
    # Despite the name, these return the events' X values in the workspace's current unit;
    # ConvertUnits rewrites each event's "tof" field in place, so here they are wavelengths.
    wavelength_min, wavelength_max = ws.getTofMin(), ws.getTofMax()
    if wavelength_min < bin_edges[0] or wavelength_max >= bin_edges[-1]:
        raise LrValidationError(
            f"Events span {wavelength_min:.3f}-{wavelength_max:.3f} Angstrom, beyond the Cd attenuation "
            f"table's {bin_edges[0]:.3f}-{bin_edges[-1]:.3f} Angstrom"
        )

    # Enumerate the bins that contain events, and compute the exponent μ·d for each.
    _, mu = _load_cd_attenuation_data()
    first_bin, last_bin = np.searchsorted(bin_edges, [wavelength_min, wavelength_max], side="right") - 1
    exponent_list = mu[first_bin : last_bin + 1] * cd_thickness

    # If all events are in bins that would overflow, raise an error
    if np.all(exponent_list > _MAX_ATTENUATION_EXPONENT):
        raise LrValidationError(
            f"All events would overflow the float32 event weights after Cd attenuation correction "
            f"(μ·d = {exponent_list.min():.3f}-{exponent_list.max():.3f} > {_MAX_ATTENUATION_EXPONENT:.3f})"
        )
    # If all events are in bins that are safe, return None
    elif np.all(exponent_list <= _MAX_ATTENUATION_EXPONENT):
        return None
    # Find the lower edge of the first bin where the exponent exceeds the safe limit
    else:
        first_overflow = first_bin + np.argmax(exponent_list > _MAX_ATTENUATION_EXPONENT)
        # FilterByXValue keeps x == XMax; step just below the edge so the overflowing bin is fully dropped.
        return float(np.nextafter(bin_edges[first_overflow], -np.inf))


#######################################################################################
### Public functions


def get_attenuator_log(ws) -> list[int]:
    """Get the attenuator log from a direct-beam workspace.

    Returns a list of 0/1 flags, one per attenuator stage, in the order of `CD_FOILS`.
    """

    def _to_bits(n: int) -> list[int]:
        """Convert an integer 0-15 to a list of attenuator flags (LSB-first binary)."""
        if not 0 <= n <= 15:
            raise SampleLogsError(f"Attenuator code must be 0-15, got {n}")
        return [(n >> i) & 1 for i in range(4)]

    try:
        atten = _to_bits(int(SampleLogs(ws)["BL4B:Actuator:Menu"]))
        return atten
    except LogNotFoundError:
        # Older data files lack the Actuator:Menu log, and used a different name for the attenuator flags.
        logs = [
            SampleLogs(ws).single_value("BL4B:Actuator:50MRb", operation=np.max),
            SampleLogs(ws).single_value("BL4B:Actuator:100MRb", operation=np.max),
            SampleLogs(ws).single_value("BL4B:Actuator:200MRb", operation=np.max),
            SampleLogs(ws).single_value("BL4B:Actuator:400MRb", operation=np.max),
        ]
        atten = [int(i) for i in logs]
        return atten


def get_transmission_workspace(cd_thickness: float) -> MatrixWorkspace:
    """The Cd transmission T(λ) = exp(-μ(λ)·d) through a given thickness, as a histogram.

    The table's wavelengths become bin centres, so μ is piecewise constant across each bin
    rather than interpolated. The workspace is kept in the analysis data service under a name
    derived from the thickness and reused by later calls.

    Parameters
    ----------
    cd_thickness
        The total thickness of Cadmium in the beam, in centimeters.

    Returns
    -------
    MatrixWorkspace
        One spectrum of transmission against wavelength, carrying no uncertainty.

    Raises
    ------
    LrValidationError
        `cd_thickness` is not a number, is negative, or is not finite.
    """
    try:
        thickness = float(cd_thickness)
    except (TypeError, ValueError) as exc:
        raise LrValidationError(f"Cd thickness must be a number of centimeters, got {cd_thickness!r}") from exc
    if not np.isfinite(thickness) or thickness < 0:
        raise LrValidationError(f"Cd thickness must be a non-negative number of centimeters, got {thickness}")
    # The name is for readability only; reuse is decided by the exact thickness logged on the workspace.
    name = f"transmission_{thickness * 1e4:g}um"
    if _is_transmission_for(name, thickness):
        return workspace_handle(name)
    wavelength, mu = _load_cd_attenuation_data()
    CreateWorkspace(OutputWorkspace=name, DataX=wavelength, DataY=np.exp(-mu * thickness), UnitX="Wavelength")
    transmission = ConvertToHistogram(InputWorkspace=name, OutputWorkspace=name)
    SampleLogs(transmission).insert(THICKNESS_LOG, thickness, unit="cm")
    return transmission


def compute_cd_thickness(atten: list[int], flip_atten: bool = False) -> float:
    """Total Cd thickness in the beam, from a run's attenuator sample logs.

    Parameters
    ----------
    atten
        One flag per attenuator stage, ordered as `CD_FOILS`: 1 in the beam, 0 out.
    flip_atten
        Whether the log was recorded with inverted polarity, 0 meaning in the beam.

    Returns
    -------
    float
        The summed thickness of the stages in the beam, in centimeters.

    Raises
    ------
    SampleLogsError
        Issue with list of attenuator flags: is not a flat sequence of numbers,
        holds fewer entries than there are stages, or holds a value other than 0 or 1.
    """
    try:
        states = np.asarray(atten, dtype=float)
    except (TypeError, ValueError) as exc:
        raise SampleLogsError(
            f"`atten` list must hold one numeric flag per attenuator stage, got {type(atten).__name__}"
        ) from exc
    if states.ndim != 1 or states.size != len(CD_FOILS):
        raise SampleLogsError(
            f"`atten` list must be a 1D array of length {len(CD_FOILS)}, "
            f"got a {states.ndim}-d array with length {states.size}."
        )
    states = states[: len(CD_FOILS)]
    if not np.isin(states, (0.0, 1.0)).all():
        raise SampleLogsError(f"`atten` list values must each be 0 or 1, got {states.tolist()}")
    if flip_atten:
        states = 1.0 - states
    return float(np.dot(CD_FOILS, states)) * 1e-4  # microns -> cm


def apply_correction(
    workspace: MantidWorkspace, cd_thickness: float, output_workspace: MantidWorkspaceName
) -> MatrixWorkspace:
    """Correct a workspace for the Cd attenuation of its beam, according to the Beer-Lambert law.

        I_corrected = I_measured / T(λ, d)

    On an event workspace this is applied per event: each event's weight is divided by the
    transmission at its wavelength, and its squared error by the transmission squared.

    Parameters
    ----------
    workspace
        The workspace to correct, in units of wavelength.
    cd_thickness
        The total thickness of Cadmium in the beam, in centimeters.
    output_workspace
        Name to register the corrected workspace under; may be the input's own name.

    Returns
    -------
    MatrixWorkspace
        The corrected workspace. An event workspace stays one.

    Raises
    ------
    LrValidationError
        `workspace` is not in wavelength; its events reach beyond the attenuation table; or the
        correction would overflow the float32 event weights.
    ValueError
        Raised by Mantid for a histogram workspace binned differently from the transmission.
    """
    ws = workspace_handle(workspace)
    unit = ws.getAxis(0).getUnit().unitID()
    if unit.lower() != "wavelength":
        raise LrValidationError(f"Cd attenuation correction needs a workspace in Wavelength, got {unit!r}")
    transmission = get_transmission_workspace(cd_thickness)
    bin_edges = transmission.readX(0)
    if (isinstance(ws, IEventWorkspace) and ws.getNumberEvents() > 0) and (
        (cutoff := _overflow_cutoff_wavelength(ws, cd_thickness, bin_edges)) is not None
    ):
        ws = FilterByXValue(InputWorkspace=ws, XMax=cutoff, OutputWorkspace=output_workspace)
    return Divide(LHSWorkspace=ws, RHSWorkspace=transmission, OutputWorkspace=output_workspace)
