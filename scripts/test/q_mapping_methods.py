"""
Pixel-to-Q mapping methods for a reflected-beam R(y, lambda) histogram, on one common input.

The three candidate methods for Requirements §2.3.9, and the event-mode reference they approximate:

* **A, point assignment** (`map_point`) — Mantid `MagnetismReflectometryReduction` with
  `ConstantQBinning=True`: each cell's center Q goes into exactly one bin, the bin value is
  `N_rows · ΣR / n_cells`, and bins with fewer than `trim · max(n_cells)` cells are dropped.
* **B, prototype spread** (`map_prototype`) — `NR_Reduction._reduce_single_run` on LR-exp `exp`: each
  cell's Q interval is shared equally among the bins whose centers it contains, with weight
  `|J| · Δλ / n_overlapped`, where `|J| = √((∂Q/∂θ)² + (∂Q/∂λ)²)` with θ in degrees; the bin value is
  `r / (Jsum / N_rows)`, multiplied by the Q-line coverage fraction.
* **C, fractional overlap** (`map_overlap`) — the Cubitt et al. (2015) constant-Q projection used by
  Mantid `ReflectometrySumInQ`, written for an input that is already divided by the direct beam:
  each cell carries `R · w_λ`, where `w_λ` is its Q width at the row angle, spread uniformly over
  its Q footprint and shared by overlap length; a bin's value is the sum of its shares over `ΔQ`.
* **E, event reference** (`map_events`) — the per-event estimator of legacy
  `EventReflectivity.specular_weighted(q_summing=True)`: each reflected event is weighted by
  `1 / I_DB(λ) · |∂Q/∂λ|`, histogrammed in Q and divided by `ΔQ`. Here `I_DB` is the same
  direct-beam histogram the cells are divided by, so E differs from A, B and C only in mapping each
  event rather than each cell.

* **T, constant TOF** (`map_constant_tof`) — Requirements `constantTOF` and legacy
  `specular(q_summing=False)`: the rows are summed at each λ first, and the sum is mapped to Q at
  the single angle `θ_central`. It is not a constant-Q method; the beam's angular spread stays in the
  result as Q resolution.

C is the histogram limit of E: for a row at a fixed angle, the events of a cell carry
`R · Δλ · |∂Q/∂λ| = R · w_λ` in total, so sharing that over the cell's Q interval by overlap is what E
does when the events are uniform within the cell.

The per-row scattering angle is a separate input (`row_angles`, or `mean_theta_row_angles` for the
prototype's slit-weighted mean angle), so the binning method and the angle model can be varied
independently. Gravity enters as an optional per-wavelength angle offset added to
every row; for E it is a function of the event wavelength.

Units: angles in radians, wavelength in Å, Q in 1/Å, unless a name says otherwise.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import gaussian_filter1d, uniform_filter1d

FOUR_PI = 4.0 * np.pi

ANGLE_CONVENTIONS = ("half", "full")
"""`half`: δ = arctan(Δy·p/L)/2 (legacy, and the Mantid magnetism code to first order).
`full`: δ = arctan(Δy·p/L) (the prototype's `constantq`, and Requirements §2.3.6).
The third angle model, the prototype's `meantheta`, needs the slits as well (`mean_theta_row_angles`)."""


@dataclass(frozen=True)
class Events:
    """One beam's events after TOF→λ: detector row, wavelength and event weight."""

    row: np.ndarray
    wavelength: np.ndarray
    weight: np.ndarray

    def __post_init__(self) -> None:
        if not (len(self.row) == len(self.wavelength) == len(self.weight)):
            raise ValueError("row, wavelength and weight must have the same length")


@dataclass(frozen=True)
class RowAngles:
    """Scattering angle of each detector row: at the row center and at its two edges (radians)."""

    rows: np.ndarray
    center: np.ndarray
    lo: np.ndarray
    hi: np.ndarray


@dataclass(frozen=True)
class Cells:
    """The reflected-beam histogram R(y, λ), already divided by the direct beam.

    `r[i, l]` is the reflected counts of row `rows[i]` in wavelength bin `l`, divided by the
    direct-beam counts in the same bin, each normalized by its own `norm` (proton charge, or incident
    neutrons for synthetic data). Summing `r` over rows gives the reflectivity at that wavelength for a
    single-angle measurement. `db` is the direct-beam counts per norm in each bin; cells where it is not
    positive have no defined ratio and are excluded by `valid`.
    """

    rows: np.ndarray
    lam_edges: np.ndarray
    r: np.ndarray
    var: np.ndarray
    db: np.ndarray

    @property
    def lam_centers(self) -> np.ndarray:
        return 0.5 * (self.lam_edges[1:] + self.lam_edges[:-1])

    @property
    def lam_widths(self) -> np.ndarray:
        return np.diff(self.lam_edges)

    @property
    def valid(self) -> np.ndarray:
        return np.broadcast_to(self.db > 0, self.r.shape)


@dataclass
class QCurve:
    """R(Q) on the bins `q_edges`; bins with no estimate are NaN."""

    method: str
    q_edges: np.ndarray
    r: np.ndarray
    dr: np.ndarray
    diagnostics: dict[str, np.ndarray] = field(default_factory=dict)

    @property
    def q_centers(self) -> np.ndarray:
        return np.sqrt(self.q_edges[1:] * self.q_edges[:-1])


# ---------------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------------


def angle_offset(dy_pixels: np.ndarray, pixel_width: float, distance: float, convention: str) -> np.ndarray:
    """Angle offset of a detector position `dy_pixels` from the reference pixel, in radians.

    `pixel_width` and `distance` must share a unit.
    """
    if convention not in ANGLE_CONVENTIONS:
        raise ValueError(f"convention must be one of {ANGLE_CONVENTIONS}, got {convention!r}")
    full = np.arctan(np.asarray(dy_pixels, dtype=float) * pixel_width / distance)
    return full / 2.0 if convention == "half" else full


def row_angles(
    rows: np.ndarray,
    theta_central: float,
    y_ref: float,
    pixel_width: float,
    distance: float,
    convention: str = "half",
    sign: float = 1.0,
) -> RowAngles:
    """Per-row scattering angles: `θ_central + sign · offset(y − y_ref)`, row `y` spanning `y ± 0.5`.

    `sign` is legacy's `sign(ths)`, which flips the pixel direction for reflection down.
    """
    rows = np.asarray(rows)
    dy = rows.astype(float) - y_ref

    def at(offset_pixels: float) -> np.ndarray:
        return theta_central + sign * angle_offset(dy + offset_pixels, pixel_width, distance, convention)

    lower, upper = at(-0.5), at(0.5)
    return RowAngles(rows=rows, center=at(0.0), lo=np.minimum(lower, upper), hi=np.maximum(lower, upper))


@dataclass(frozen=True)
class SlitGeometry:
    """The incident slits, the distances and the detector response that `mean_theta_offset` uses.

    `si` is the gap of the slit nearer the sample and `s1` that of the slit upstream of it. All lengths
    share one unit, the unit of the pixel width they are used with. `det_sigma` is the σ of the detector
    response, `det_res_fn` its shape: `rectangular` (the prototype default, a box of width σ·√12),
    `gaussian`, or `none`.
    """

    s1: float
    si: float
    interslit: float
    si_sample: float
    sample_detector: float
    det_sigma: float = 0.0
    det_res_fn: str = "none"


def mean_theta_offset(dy: np.ndarray, geometry: SlitGeometry) -> tuple[np.ndarray, np.ndarray]:
    """Mean and σ of the incidence-angle offset of the neutrons measured at detector offset `dy` (radians).

    A vectorized port of the prototype's `meantheta` angle model (LR-exp `calc_mean_theta`,
    `_calc_detector_convolution` and `nr_tools.calc_beam_geometry_from_slits`). The neutrons that pass
    both slits fill a parallelogram in (angle at Si, position on the detector), with the detector a
    straight-line distance `si_sample + sample_detector` from Si, the reflection unfolded. The
    parallelogram is sampled on a 1000 × 1000 grid, each angle's column is convolved with the detector
    response, and at each `dy` the nearest grid row gives the mean and σ of the angle. Offsets beyond the
    convolved beam take the edge row's values.

    The slit gaps go in as the signature names them; the prototype's `meantheta` call passes them
    swapped. `det_res_fn = "none"`, which the prototype does not offer, skips the convolution.
    """
    g = geometry
    if g.det_res_fn not in ("rectangular", "gaussian", "none"):
        raise ValueError(f"det_res_fn must be 'rectangular', 'gaussian' or 'none', got {g.det_res_fn!r}")
    # The parallelogram's corners, angles in degrees as in the prototype.
    a_plus = np.arctan((g.s1 + g.si) / (2 * g.interslit))
    a_minus = np.arctan((g.s1 - g.si) / (2 * g.interslit))
    a_rad = np.array([a_plus, a_minus, -a_plus, -a_minus])
    length = g.si_sample + g.sample_detector
    y = np.array([g.si / 2, -g.si / 2, -g.si / 2, g.si / 2]) + length * np.tan(a_rad)
    a = np.degrees(a_rad)
    slope = (y - np.roll(y, -1)) / (a - np.roll(a, -1))
    icept = y - slope * a
    nxt = np.roll(np.arange(4), -1)
    t_corner = (icept[nxt] - icept) / (slope - slope[nxt])
    y_corner = slope * t_corner + icept

    if g.det_res_fn == "rectangular":
        pad = width = g.det_sigma * np.sqrt(12)
    elif g.det_res_fn == "gaussian":
        pad = g.det_sigma * 4.0
    else:
        pad = 0.0
    n = 1000
    t_step = (t_corner.max() - t_corner.min()) / n
    y_step = ((y_corner.max() + pad) - (y_corner.min() - pad)) / n
    t_vec = np.linspace(t_corner.min(), t_corner.max() + t_step, n)
    y_vec = np.linspace(y_corner.min() - pad, y_corner.max() + pad + y_step, n)
    hi = np.minimum(t_vec * slope[3] + icept[3], t_vec * slope[2] + icept[2])
    lo = np.maximum(t_vec * slope[0] + icept[0], t_vec * slope[1] + icept[1])
    smear = ((y_vec[:, None] >= lo[None, :]) & (y_vec[:, None] <= hi[None, :])).astype(float)
    if g.det_res_fn == "rectangular":
        smear = uniform_filter1d(smear, size=int(width / y_step), axis=0)
    elif g.det_res_fn == "gaussian":
        smear = gaussian_filter1d(smear, sigma=g.det_sigma / y_step, axis=0)
    keep = np.any(smear != 0, axis=1)
    smear, y_vec = smear[keep], y_vec[keep]

    dy = np.atleast_1d(np.asarray(dy, dtype=float))
    nearest = np.argmin(np.abs(dy[:, None] - y_vec[None, :]), axis=1)
    weights = smear[nearest]
    total = weights.sum(axis=1)
    mean = weights @ t_vec / total
    var = np.sum(weights * (t_vec[None, :] - mean[:, None]) ** 2, axis=1) / total
    return np.radians(mean), np.radians(np.sqrt(var))


def mean_theta_row_angles(
    rows: np.ndarray,
    theta_central: float,
    y_ref: float,
    pixel_width: float,
    geometry: SlitGeometry,
    sign: float = 1.0,
) -> RowAngles:
    """Per-row angles under the prototype's `meantheta` model: `θ_central + sign · mean offset(y − y_ref)`.

    `pixel_width` shares the unit of `geometry`. The row edges take the model at `y ± 0.5`; where the
    model is flat, they coincide with the center.
    """
    rows = np.asarray(rows)
    dy = rows.astype(float) - y_ref

    def at(offset_pixels: float) -> np.ndarray:
        return theta_central + sign * mean_theta_offset((dy + offset_pixels) * pixel_width, geometry)[0]

    lower, upper = at(-0.5), at(0.5)
    return RowAngles(rows=rows, center=at(0.0), lo=np.minimum(lower, upper), hi=np.maximum(lower, upper))


def histogram_cells(
    reflected: Events,
    direct: Events,
    rows: np.ndarray,
    lam_edges: np.ndarray,
    reflected_norm: float = 1.0,
    direct_norm: float = 1.0,
    include_direct_beam_error: bool = False,
) -> Cells:
    """Histogram both beams on `lam_edges` and divide: the single Rebin and the divide of §2.3.5.

    The direct beam is summed over all its rows. By default the cell variance carries only the
    reflected-beam counting error, so every method and the event reference E are compared on the same
    error model; `include_direct_beam_error` adds the direct-beam term to each cell.
    """
    rows = np.asarray(rows)
    lam_edges = np.asarray(lam_edges, dtype=float)
    n_lam = len(lam_edges) - 1

    row_index = _row_index(reflected.row, rows)
    lam_index = np.searchsorted(lam_edges, reflected.wavelength, side="right") - 1
    keep = (row_index >= 0) & (lam_index >= 0) & (lam_index < n_lam)
    flat = row_index[keep] * n_lam + lam_index[keep]
    size = len(rows) * n_lam
    counts = np.bincount(flat, weights=reflected.weight[keep], minlength=size).reshape(len(rows), n_lam)
    counts_var = np.bincount(flat, weights=reflected.weight[keep] ** 2, minlength=size).reshape(len(rows), n_lam)

    db_index = np.searchsorted(lam_edges, direct.wavelength, side="right") - 1
    db_keep = (db_index >= 0) & (db_index < n_lam)
    db = np.bincount(db_index[db_keep], weights=direct.weight[db_keep], minlength=n_lam) / direct_norm
    db_var = np.bincount(db_index[db_keep], weights=direct.weight[db_keep] ** 2, minlength=n_lam) / direct_norm**2

    with np.errstate(divide="ignore", invalid="ignore"):
        positive = db > 0
        r = np.where(positive, counts / reflected_norm / db, 0.0)
        var = np.where(positive, counts_var / reflected_norm**2 / db**2, 0.0)
        if include_direct_beam_error:
            var = var + np.where(positive, r**2 * db_var / db**2, 0.0)
    return Cells(rows=rows, lam_edges=lam_edges, r=r, var=var, db=db)


def log_q_edges(q_min: float, q_max: float, dq_over_q: float) -> np.ndarray:
    """Logarithmic bin edges with constant `ΔQ/Q`, the same rule as legacy `get_q_binning(q_step<0)`."""
    step = 1.0 + abs(dq_over_q)
    n = int(np.log(q_max / q_min) / np.log(step)) + 1
    return q_min * step ** np.arange(n)


def qline_coverage(q_centers: np.ndarray, angles: RowAngles, lam_min: float, lam_max: float) -> np.ndarray:
    """Fraction of each constant-Q line's λ span (across the rows) inside `[lam_min, lam_max]`, §2.3.10.

    The prototype's `Qline_fraction`: `(min(λmax, λ(θmax)) − max(λmin, λ(θmin))) / (λ(θmax) − λ(θmin))`.
    """
    theta_min = np.min(np.abs(angles.center))
    theta_max = np.max(np.abs(angles.center))
    lam_lo = FOUR_PI * np.sin(theta_min) / q_centers
    lam_hi = FOUR_PI * np.sin(theta_max) / q_centers
    span = lam_hi - lam_lo
    with np.errstate(divide="ignore", invalid="ignore"):
        fraction = (np.minimum(lam_max, lam_hi) - np.maximum(lam_min, lam_lo)) / span
    # A single row, or rows at one angle: the line is a point, covered when it lies in the window.
    point = span <= 0
    fraction[point] = ((lam_lo[point] >= lam_min) & (lam_lo[point] <= lam_max)).astype(float)
    return fraction


# ---------------------------------------------------------------------------------------------------
# Method A: point assignment (Mantid MagnetismReflectometryReduction, ConstantQBinning=True)
# ---------------------------------------------------------------------------------------------------


def map_point(
    cells: Cells,
    angles: RowAngles,
    q_edges: np.ndarray,
    gravity: np.ndarray | None = None,
    trim: float = 0.5,
) -> QCurve:
    """Method A. `trim` is Mantid's `ConstQTrim`; 0 keeps every bin with at least one cell.

    Ported from `constant_q` in Mantid 6.16 `MagnetismReflectometryReduction.py`, with the angle model
    taken from `angles` instead of the code's own `θ + θ_f` and `cos(θ_f)` factor.
    """
    theta = _cell_angles(angles.center, cells, gravity)
    q = FOUR_PI * np.sin(theta) / cells.lam_centers[None, :]
    n_q = len(q_edges) - 1
    k = np.digitize(q, q_edges) - 1
    use = cells.valid & (k >= 0) & (k < n_q)

    total = np.bincount(k[use], weights=cells.r[use], minlength=n_q)
    total_var = np.bincount(k[use], weights=cells.var[use], minlength=n_q)
    n = np.bincount(k[use], minlength=n_q).astype(float)

    n_rows = cells.r.shape[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        r = n_rows * total / n
        dr = n_rows * np.sqrt(total_var) / n
    dropped = (n == 0) | (n < trim * n.max(initial=0.0))
    r[dropped] = np.nan
    dr[dropped] = np.nan
    return QCurve("A point", np.asarray(q_edges), r, dr, {"n_cells": n})


# ---------------------------------------------------------------------------------------------------
# Method B: the prototype's Jacobian-weighted equal-share spread
# ---------------------------------------------------------------------------------------------------


def prototype_jacobian(lam: np.ndarray, theta_degrees: np.ndarray, include_dqdtheta: bool = True) -> np.ndarray:
    """`NR_Reduction._calculate_jacobian`: ∂Q/∂θ per degree and ∂Q/∂λ per Å, combined in quadrature."""
    dqdtheta = (FOUR_PI / lam) * np.cos(np.radians(theta_degrees)) * np.pi / 180
    dqdlambda = -(FOUR_PI / lam**2) * np.sin(np.radians(theta_degrees))
    if include_dqdtheta:
        return np.abs(np.sqrt(dqdtheta**2 + dqdlambda**2))
    return np.abs(dqdlambda)


def map_prototype(
    cells: Cells,
    angles: RowAngles,
    q_edges: np.ndarray,
    gravity: np.ndarray | None = None,
    qline_threshold: float = 1.0,
    lam_window: tuple[float, float] | None = None,
) -> QCurve:
    """Method B, a vectorized port of the `constantq` / `meantheta` branch of the prototype's Q loop.

    The prototype maps onto Q *points* (`self.q`, from `log_qvector`); here those points are the
    geometric centers of `q_edges`, so all methods share one output grid. As in the prototype:

    * a cell's interval is `[4π sin(θ − Δθ/2)/(λ + Δλ/2), 4π sin(θ + Δθ/2)/(λ − Δλ/2)]` with
      `Δθ = |θ(y+1) − θ(y)|` (the last row repeats the one before) and `Δλ` the bin width;
    * the points are first restricted to those whose Q-line coverage (§2.3.10, over `lam_window`,
      the λ grid's range when omitted) is at least `qline_threshold`;
    * a cell is used only if its whole interval lies within the range of the remaining points;
    * it is shared equally among the points inside its interval, or given to the nearest point
      when none is;
    * the result is multiplied by the coverage fraction.

    Not ported: the prototype sets cells with no direct beam to 0 and includes them; here they are
    excluded, as in every method.
    """
    q_all = np.sqrt(q_edges[1:] * q_edges[:-1])
    lam_min, lam_max = lam_window if lam_window is not None else (cells.lam_edges[0], cells.lam_edges[-1])
    fraction_all = qline_coverage(q_all, angles, lam_min, lam_max)
    selected = np.flatnonzero(fraction_all >= qline_threshold)
    r_all = np.full(len(q_all), np.nan)
    dr_all = np.full(len(q_all), np.nan)
    jsum_all = np.zeros(len(q_all))
    if len(selected) == 0:
        return QCurve(
            "B prototype", np.asarray(q_edges), r_all, dr_all, {"jsum": jsum_all, "qline_fraction": fraction_all}
        )
    q_points = q_all[selected]
    fraction = fraction_all[selected]
    n_q = len(q_points)
    n_rows = cells.r.shape[0]

    theta = _cell_angles(angles.center, cells, gravity)
    theta_deg = np.abs(np.degrees(theta))
    theta_bin = np.abs(np.diff(np.degrees(angles.center)))
    theta_bin = np.concatenate([theta_bin, theta_bin[-1:]]) if len(theta_bin) else np.zeros(1)
    half_theta = (theta_bin / 2.0)[:, None]

    lam = cells.lam_centers[None, :]
    dlam = cells.lam_widths[None, :]
    jac = prototype_jacobian(lam, theta_deg, include_dqdtheta=True)
    qcen = FOUR_PI * np.sin(np.radians(theta_deg)) / lam
    qlo = FOUR_PI * np.sin(np.radians(theta_deg - half_theta)) / (lam + dlam / 2)
    qhi = FOUR_PI * np.sin(np.radians(theta_deg + half_theta)) / (lam - dlam / 2)

    use = cells.valid & (qhi <= q_points[-1]) & (qlo >= q_points[0])
    first = np.searchsorted(q_points, qlo, side="left")
    stop = np.searchsorted(q_points, qhi, side="right")
    empty = stop <= first
    nearest = _nearest_index(q_points, qcen)
    first = np.where(empty, nearest, first)
    stop = np.where(empty, nearest + 1, stop)
    share = jac / (stop - first) * dlam

    signal = _range_add(first[use], stop[use], (cells.r * share)[use], n_q)
    err2 = _range_add(first[use], stop[use], (np.sqrt(cells.var) * share)[use] ** 2, n_q)
    jsum = _range_add(first[use], stop[use], np.broadcast_to(share, cells.r.shape)[use], n_q)
    # The difference array leaves round-off where no cell deposits; count deposits exactly instead.
    n_cells = np.rint(_range_add(first[use], stop[use], np.ones(int(use.sum())), n_q))

    factor = jsum / n_rows
    with np.errstate(divide="ignore", invalid="ignore"):
        r = signal / factor * fraction
        dr = np.sqrt(err2) / factor * fraction
    dropped = ~np.isfinite(r) | (n_cells == 0)
    r[dropped] = np.nan
    dr[dropped] = np.nan
    r_all[selected], dr_all[selected], jsum_all[selected] = r, dr, jsum
    return QCurve("B prototype", np.asarray(q_edges), r_all, dr_all, {"jsum": jsum_all, "qline_fraction": fraction_all})


# ---------------------------------------------------------------------------------------------------
# Method C: fractional overlap of each cell's Q footprint
# ---------------------------------------------------------------------------------------------------


def map_overlap(
    cells: Cells,
    angles: RowAngles,
    q_edges: np.ndarray,
    gravity: np.ndarray | None = None,
    include_angular_width: bool = True,
) -> QCurve:
    """Method C.

    A cell's footprint is `[4π sin θ_lo/λ_hi, 4π sin θ_hi/λ_lo]` with `include_angular_width` (the
    Cubitt / `ReflectometrySumInQ` corner projection), or `[4π sin θ/λ_hi, 4π sin θ/λ_lo]` at the row
    angle without it. Either way the cell carries `R · w_λ`, with `w_λ` the second interval's width,
    so including the angular width spreads the same amount over a wider range. Bin `k` gets
    `Σ R · w_λ · overlap_k / (footprint · ΔQ_k)`; parts of a footprint outside the grid are dropped.

    `row_coverage` is, per bin, the mean over rows of the fraction of the bin the row's cells cover:
    1 where every row spans the whole bin.

    With `include_angular_width`, the uniform box over the corner footprint is wider than the spacing
    of the cells, so even a uniform input comes out with a ripple of order `w_λ / footprint` before
    the Q bins average it; the true distribution of Q in a cell is the λ box convolved with the θ
    box, a trapezoid. The total is conserved either way.
    """
    q_edges = np.asarray(q_edges, dtype=float)
    n_q = len(q_edges) - 1
    widths_q = np.diff(q_edges)
    n_rows = cells.r.shape[0]

    lam_lo = cells.lam_edges[None, :-1]
    lam_hi = cells.lam_edges[None, 1:]
    theta = _cell_angles(angles.center, cells, gravity)
    q_lam_lo = FOUR_PI * np.sin(theta) / lam_hi
    q_lam_hi = FOUR_PI * np.sin(theta) / lam_lo
    mass_width = q_lam_hi - q_lam_lo
    if include_angular_width:
        foot_lo = FOUR_PI * np.sin(_cell_angles(angles.lo, cells, gravity)) / lam_hi
        foot_hi = FOUR_PI * np.sin(_cell_angles(angles.hi, cells, gravity)) / lam_lo
    else:
        foot_lo, foot_hi = q_lam_lo, q_lam_hi
    foot = foot_hi - foot_lo

    use = cells.valid & (foot > 0) & (foot_hi > q_edges[0]) & (foot_lo < q_edges[-1])
    lo, hi = foot_lo[use], foot_hi[use]
    density = (cells.r * mass_width / foot)[use]
    density_sigma = (np.sqrt(cells.var) * mass_width / foot)[use]
    coverage_density = (mass_width / foot)[use]

    k_first = np.clip(np.searchsorted(q_edges, lo, side="right") - 1, 0, n_q - 1)
    k_last = np.clip(np.searchsorted(q_edges, hi, side="left") - 1, 0, n_q - 1)
    signal = np.zeros(n_q)
    err2 = np.zeros(n_q)
    covered = np.zeros(n_q)
    for step in range(int((k_last - k_first).max(initial=0)) + 1):
        k = k_first + step
        inside = k <= k_last
        kk = k[inside]
        overlap = np.minimum(hi[inside], q_edges[kk + 1]) - np.maximum(lo[inside], q_edges[kk])
        overlap = np.clip(overlap, 0.0, None)
        signal += np.bincount(kk, weights=density[inside] * overlap, minlength=n_q)
        err2 += np.bincount(kk, weights=(density_sigma[inside] * overlap) ** 2, minlength=n_q)
        covered += np.bincount(kk, weights=coverage_density[inside] * overlap, minlength=n_q)

    r = signal / widths_q
    dr = np.sqrt(err2) / widths_q
    row_coverage = covered / widths_q / n_rows
    empty = covered == 0
    r[empty] = np.nan
    dr[empty] = np.nan
    name = "C overlap" if include_angular_width else "C overlap (no angular width)"
    return QCurve(name, q_edges, r, dr, {"row_coverage": row_coverage})


# ---------------------------------------------------------------------------------------------------
# T: constant TOF
# ---------------------------------------------------------------------------------------------------


def map_constant_tof(
    cells: Cells, theta_central: float, q_edges: np.ndarray, gravity: np.ndarray | None = None
) -> QCurve:
    """T: sum the rows at each λ, then map each λ bin to Q at the one angle `θ_central`.

    Requirements §2.3.5/§2.3.6 `constantTOF`, and legacy `specular(q_summing=False)` in its histogram
    limit. `R(λ) = Σ_y r[y, λ]` is spread over `[4π sin θ_c/λ_hi, 4π sin θ_c/λ_lo]` by overlap, which is
    method C0 applied to a single row at `θ_central`. Every row of a λ bin goes to the same Q, so the
    result is R averaged over the beam's angular spread at fixed λ; the spread becomes Q resolution.
    """
    collapsed = Cells(
        rows=np.array([0]),
        lam_edges=cells.lam_edges,
        r=cells.r.sum(axis=0, keepdims=True),
        var=cells.var.sum(axis=0, keepdims=True),
        db=cells.db,
    )
    theta = np.array([theta_central], dtype=float)
    one_row = RowAngles(rows=np.array([0]), center=theta, lo=theta, hi=theta)
    curve = map_overlap(collapsed, one_row, q_edges, gravity=gravity, include_angular_width=False)
    curve.method = "T constant TOF"
    return curve


# ---------------------------------------------------------------------------------------------------
# E: event-mode reference
# ---------------------------------------------------------------------------------------------------


def map_events(
    reflected: Events,
    cells: Cells,
    angles: RowAngles,
    q_edges: np.ndarray,
    reflected_norm: float = 1.0,
    gravity: Callable[[np.ndarray], np.ndarray] | None = None,
) -> QCurve:
    """E: each reflected event at its own Q, weighted by `1/I_DB(λ) · |∂Q/∂λ| / ΔQ`.

    `I_DB(λ)` is `cells.db` per unit wavelength in the event's bin, so the direct beam is the same step
    function A, B and C divide by. `gravity`, when given, maps event wavelengths to an angle offset.
    The error is the reflected-beam counting error only.
    """
    q_edges = np.asarray(q_edges, dtype=float)
    n_q = len(q_edges) - 1
    lam = reflected.wavelength
    row_index = _row_index(reflected.row, angles.rows)
    lam_index = np.searchsorted(cells.lam_edges, lam, side="right") - 1
    n_lam = len(cells.lam_edges) - 1
    keep = (row_index >= 0) & (lam_index >= 0) & (lam_index < n_lam)
    keep[keep] &= cells.db[lam_index[keep]] > 0

    lam = lam[keep]
    theta = angles.center[row_index[keep]]
    if gravity is not None:
        theta = theta + gravity(lam)
    q = FOUR_PI * np.sin(theta) / lam
    density = cells.db[lam_index[keep]] / cells.lam_widths[lam_index[keep]]
    weight = reflected.weight[keep] / reflected_norm / density * (q / lam)

    k = np.digitize(q, q_edges) - 1
    inside = (k >= 0) & (k < n_q)
    widths_q = np.diff(q_edges)
    signal = np.bincount(k[inside], weights=weight[inside], minlength=n_q) / widths_q
    err2 = np.bincount(k[inside], weights=weight[inside] ** 2, minlength=n_q) / widths_q**2
    n_events = np.bincount(k[inside], minlength=n_q).astype(float)
    r = signal
    dr = np.sqrt(err2)
    r[n_events == 0] = np.nan
    dr[n_events == 0] = np.nan
    return QCurve("E events", q_edges, r, dr, {"n_events": n_events})


# ---------------------------------------------------------------------------------------------------
# Angle-model check: do the rows agree in Q?
# ---------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class RowAlignment:
    """Agreement of the per-row R(Q) curves, each scaled to the common mean (see `row_alignment`)."""

    chi2_per_dof: float
    rms_relative: float
    rows: np.ndarray
    n_bins: int


def row_alignment(
    cells: Cells, angles: RowAngles, q_edges: np.ndarray, min_row_share: float = 0.03, iterations: int = 5
) -> RowAlignment:
    """Map each row on its own (method C without angular width) and measure how well the rows agree.

    Each row sees `p_y · R(Q)`, so with the right angle model every row's curve is the same shape up to
    a scale. The rows holding at least `min_row_share` of the signal are scaled to a common weighted
    mean on the bins every one of them fully covers, and the scatter about it is reported as
    `χ²/dof` (1 when the rows agree within counting statistics) and as an unweighted relative RMS.
    A wrong angle model moves each row's sharp features to a different Q and raises both.
    """
    share = cells.r.sum(axis=1) / cells.r.sum()
    use_rows = np.flatnonzero(share >= min_row_share)
    curves, sigmas, full = [], [], []
    for i in use_rows:
        one = Cells(cells.rows[i : i + 1], cells.lam_edges, cells.r[i : i + 1], cells.var[i : i + 1], cells.db)
        one_angles = RowAngles(
            angles.rows[i : i + 1], angles.center[i : i + 1], angles.lo[i : i + 1], angles.hi[i : i + 1]
        )
        curve = map_overlap(one, one_angles, q_edges, include_angular_width=False)
        curves.append(curve.r)
        sigmas.append(curve.dr)
        full.append(curve.diagnostics["row_coverage"] > 1 - 1e-9)
    r = np.array(curves)
    dr = np.array(sigmas)
    common = np.all(np.array(full), axis=0) & np.all(np.isfinite(r) & (dr > 0) & (r > 0), axis=0)
    r, dr = r[:, common], dr[:, common]
    n_rows, n_bins = r.shape
    if n_rows < 2 or n_bins < 3:
        return RowAlignment(np.nan, np.nan, cells.rows[use_rows], n_bins)

    scale = r.mean(axis=1) / r.mean()
    for _ in range(iterations):
        w = (scale[:, None] / dr) ** 2
        mean = np.sum(w * r / scale[:, None], axis=0) / np.sum(w, axis=0)
        scale = np.sum(r * mean / dr**2, axis=1) / np.sum(mean**2 / dr**2, axis=1)
    resid = (r / scale[:, None] - mean) / (dr / scale[:, None])
    dof = n_rows * n_bins - n_bins - (n_rows - 1)
    rms = float(np.sqrt(np.mean((r / scale[:, None] / mean - 1.0) ** 2)))
    return RowAlignment(float(np.sum(resid**2) / dof), rms, cells.rows[use_rows], n_bins)


# ---------------------------------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------------------------------


def _row_index(event_rows: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """Position of each event's row in `rows`, or −1 when the row is not in it."""
    rows = np.asarray(rows)
    order = np.argsort(rows)
    sorted_rows = rows[order]
    pos = np.searchsorted(sorted_rows, event_rows)
    pos_clipped = np.clip(pos, 0, len(rows) - 1)
    found = sorted_rows[pos_clipped] == event_rows
    return np.where(found, order[pos_clipped], -1)


def _cell_angles(row_values: np.ndarray, cells: Cells, gravity: np.ndarray | None) -> np.ndarray:
    """Broadcast a per-row angle over the λ axis and add the per-λ gravity offset, if any."""
    theta = np.broadcast_to(np.asarray(row_values, dtype=float)[:, None], cells.r.shape)
    if gravity is None:
        return theta
    gravity = np.asarray(gravity, dtype=float)
    if gravity.shape != cells.lam_centers.shape:
        raise ValueError("gravity must hold one angle offset per wavelength bin")
    return theta + gravity[None, :]


def _nearest_index(points: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Index of the nearest of the sorted `points` to each value (first one on a tie, like argmin)."""
    right = np.clip(np.searchsorted(points, values, side="left"), 0, len(points) - 1)
    left = np.clip(right - 1, 0, len(points) - 1)
    closer_left = np.abs(values - points[left]) <= np.abs(points[right] - values)
    return np.where(closer_left, left, right)


def _range_add(first: np.ndarray, stop: np.ndarray, values: np.ndarray, size: int) -> np.ndarray:
    """Add `values[i]` to every index in `[first[i], stop[i])`, via a difference array."""
    diff = np.zeros(size + 1)
    np.add.at(diff, first, values)
    np.add.at(diff, stop, -values)
    return np.cumsum(diff[:-1])
