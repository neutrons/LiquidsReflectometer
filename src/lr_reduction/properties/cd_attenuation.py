"""Cd-foil transmission correction for direct-beam runs.

A direct-beam run is measured through a stack of Cd foils chosen to keep the detector on-scale.
Dividing each event's weight by the stack's transmission, following the Beer-Lambert law
T(λ) = exp(-μ(λ)·d), recovers the unattenuated incident flux the sample would see.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
from mantid.api import IEventWorkspace, MatrixWorkspace
from mantid.simpleapi import ConvertToHistogram, CreateWorkspace, Divide

from lr_reduction.exceptions import LrValidationError, SampleLogsError
from lr_reduction.types import MantidWorkspace, MantidWorkspaceName
from lr_reduction.utils.workspace import workspace_exists, workspace_handle

cd_attenuation_file = Path(__file__).parent / "Cd-attenuation.csv"

#: Cumulative Cd thickness in microns that each attenuator stage puts in the beam, in log order.
CD_FOILS = (57.5, 126.5, 126.5 + 123.0, 2 * (126.5 + 123.0))

# Mantid stores weighted-event weights and squared errors as float32. Dividing by T = exp(-μ·d)
# scales a squared error by exp(2·μ·d), which overflows to inf beyond this exponent.
_MAX_ATTENUATION_EXPONENT = 0.5 * np.log(np.finfo(np.float32).max)


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
    wavelength, mu = np.loadtxt(cd_attenuation_file, delimiter="\t", skiprows=1, unpack=True)
    wavelength.setflags(write=False)
    mu.setflags(write=False)
    return wavelength, mu


def compute_cd_thickness(atten, flip_atten: bool = False) -> float:
    """Total Cd thickness in the beam, from a run's ``Atten`` log value.

    Parameters
    ----------
    atten
        One flag per attenuator stage, ordered as `CD_FOILS`: 1 in the beam, 0 out. Entries
        beyond ``len(CD_FOILS)`` are ignored.
    flip_atten
        Whether the log was recorded with inverted polarity, 0 meaning in the beam.

    Returns
    -------
    float
        The summed thickness of the stages in the beam, in centimeters.

    Raises
    ------
    SampleLogsError
        `atten` is not a flat sequence of numbers, holds fewer entries than there are stages,
        or holds a value other than 0 or 1.
    """
    try:
        states = np.asarray(atten, dtype=float)
    except (TypeError, ValueError) as exc:
        raise SampleLogsError(
            f"Atten log must hold one numeric flag per attenuator stage, got {type(atten).__name__}"
        ) from exc
    if states.ndim != 1:
        raise SampleLogsError(f"Atten log must hold one flag per attenuator stage, got a {states.ndim}-d value")
    if states.size < len(CD_FOILS):
        raise SampleLogsError(f"Atten log holds {states.size} values, fewer than the {len(CD_FOILS)} attenuator stages")
    states = states[: len(CD_FOILS)]
    if not np.isin(states, (0.0, 1.0)).all():
        raise SampleLogsError(f"Atten log values must each be 0 or 1, got {states.tolist()}")
    if flip_atten:
        states = 1.0 - states
    return float(np.dot(CD_FOILS, states)) * 1e-4  # microns -> cm


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
        `cd_thickness` is negative or not finite.
    """
    if not np.isfinite(cd_thickness) or cd_thickness < 0:
        raise LrValidationError(f"Cd thickness must be a non-negative number of centimeters, got {cd_thickness}")
    # Keyed to 1e-8 cm: a coarser key would hand one thickness the cached transmission of a nearby one.
    name = f"transmission_{cd_thickness:.8f}"
    if workspace_exists(name):
        return workspace_handle(name)
    wavelength, mu = _load_cd_attenuation_data()
    CreateWorkspace(OutputWorkspace=name, DataX=wavelength, DataY=np.exp(-mu * cd_thickness), UnitX="Wavelength")
    return ConvertToHistogram(InputWorkspace=name, OutputWorkspace=name)


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
    if unit != "Wavelength":
        raise LrValidationError(f"Cd attenuation correction needs a workspace in Wavelength, got {unit!r}")
    transmission = get_transmission_workspace(cd_thickness)
    if isinstance(ws, IEventWorkspace) and ws.getNumberEvents() > 0:
        _check_event_wavelengths(ws, cd_thickness, transmission.readX(0))
    return Divide(LHSWorkspace=ws, RHSWorkspace=transmission, OutputWorkspace=output_workspace)


def _check_event_wavelengths(ws: IEventWorkspace, cd_thickness: float, bin_edges: np.ndarray) -> None:
    """Raise unless the correction reaches every event and keeps its weight finite.

    Both would otherwise fail silently: Mantid leaves an event outside the transmission's bins
    unscaled, and a weight or squared error past the float32 range becomes inf.
    """
    wavelength_min, wavelength_max = ws.getTofMin(), ws.getTofMax()
    if wavelength_min < bin_edges[0] or wavelength_max >= bin_edges[-1]:
        raise LrValidationError(
            f"Events span {wavelength_min:.3f}-{wavelength_max:.3f} Angstrom, beyond the Cd attenuation "
            f"table's {bin_edges[0]:.3f}-{bin_edges[-1]:.3f} Angstrom"
        )
    _, mu = _load_cd_attenuation_data()
    first_bin, last_bin = np.searchsorted(bin_edges, [wavelength_min, wavelength_max], side="right") - 1
    exponent = mu[first_bin : last_bin + 1].max() * cd_thickness
    if exponent > _MAX_ATTENUATION_EXPONENT:
        raise LrValidationError(
            f"Cd transmission falls to exp(-{exponent:.1f}) within the events' wavelength range; "
            f"correcting it would overflow the float32 event weights"
        )
