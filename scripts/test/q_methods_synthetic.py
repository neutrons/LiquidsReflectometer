"""
Synthetic single-angle reflectivity measurement with a known R(Q), for `compare_q_methods.py`.

Two beam geometries (`SyntheticSetup.geometry`):

* `point` — a ray's detector row is set by its angle alone: a ray incident at `θ_c + δ` lands at
  `y_ref + tan(2δ)·L/p` pixels, the `half` angle convention. This is the relation for a rotation of the
  sample by δ; a point sample in a divergent beam puts the ray at `tan(δ)·L`.
* `slits` — rays are traced through the two slits: a ray crossing S1 at height `h1` and Si at `hi`
  has `δ = (hi − h1)/D` and lands at `y_ref + (hi + (d_Si,sample + L)·tan δ)/p`, the reflection unfolded.
  A row then holds a spread of angles, and the mean angle per row is what the prototype's `meantheta`
  model computes.

In both:

* optional Gaussian detector blur moves a ray off its row;
* the incident angles follow the trapezoid that two slits of gaps `s1`, `si` a distance `D` apart
  produce, the sum of two uniform distributions;
* the flux spectrum is a smooth moderator-like shape, and wavelengths are exact (no TOF resolution);
* there is no gravity and no background.

Events are returned as `q_mapping_methods.Events`. In `noise_free` mode they are a dense deterministic
grid of weighted pseudo-events, so histogram cells and the event reference E see the same expected
counts; under `slits` the pseudo-events of one row and wavelength are summed into one, which every
method sees the same way. In Poisson mode they are individually sampled unit-weight events.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from q_mapping_methods import FOUR_PI, Events, SlitGeometry
from scipy.special import ndtr

GEOMETRIES = ("point", "slits")


@dataclass(frozen=True)
class Layer:
    sld: float  # 1/Å², e.g. 2.07e-6 for Si
    thickness: float  # Å
    roughness: float  # Å, of the interface on top of this layer


def parratt(q: np.ndarray, layers: list[Layer], substrate: Layer, ambient_sld: float = 0.0) -> np.ndarray:
    """Specular reflectivity by the Parratt recursion with Névot–Croce roughness."""
    q = np.asarray(q, dtype=float)
    sld = np.array([ambient_sld] + [layer.sld for layer in layers] + [substrate.sld])
    thickness = np.array([0.0] + [layer.thickness for layer in layers] + [0.0])
    sigma = np.array([layer.roughness for layer in layers] + [substrate.roughness])
    kz = np.sqrt((q[:, None] / 2.0) ** 2 - FOUR_PI * (sld - sld[0]) + 0j)
    x = np.zeros(len(q), dtype=complex)
    for j in range(len(sld) - 2, -1, -1):
        rj = (
            (kz[:, j] - kz[:, j + 1])
            / (kz[:, j] + kz[:, j + 1])
            * np.exp(-2.0 * kz[:, j] * kz[:, j + 1] * sigma[j] ** 2)
        )
        phase = np.exp(2j * kz[:, j + 1] * thickness[j + 1])
        x = (rj + x * phase) / (1.0 + rj * x * phase)
    return np.abs(x) ** 2


DEFAULT_FILM = [Layer(sld=4.0e-6, thickness=400.0, roughness=4.0)]
DEFAULT_SUBSTRATE = Layer(sld=2.07e-6, thickness=0.0, roughness=3.0)


@dataclass(frozen=True)
class SyntheticSetup:
    theta_central_deg: float = 0.45
    y_ref: float = 150.3  # fractional specular pixel
    pixel_width_mm: float = 0.70
    distance_mm: float = 1830.0  # sample to detector
    s1_mm: float = 2.0
    si_mm: float = 1.2
    interslit_mm: float = 1390.0
    si_sample_mm: float = 135.0  # Si to sample; enters only the `slits` geometry
    lam_min: float = 2.75
    lam_max: float = 9.5
    blur_pixels: float = 0.0  # Gaussian sigma of the detector response
    geometry: str = "point"

    def __post_init__(self) -> None:
        if self.geometry not in GEOMETRIES:
            raise ValueError(f"geometry must be one of {GEOMETRIES}, got {self.geometry!r}")

    @property
    def theta_central(self) -> float:
        return np.radians(self.theta_central_deg)

    def slit_geometry(self) -> SlitGeometry:
        """The slits and detector response in mm, as the `meantheta` angle model takes them."""
        return SlitGeometry(
            s1=self.s1_mm,
            si=self.si_mm,
            interslit=self.interslit_mm,
            si_sample=self.si_sample_mm,
            sample_detector=self.distance_mm,
            det_sigma=self.blur_pixels * self.pixel_width_mm,
            det_res_fn="gaussian" if self.blur_pixels > 0 else "none",
        )

    def delta_limits(self) -> tuple[float, float]:
        """Half-widths of the two uniform components of the incident-angle trapezoid (radians)."""
        return (self.s1_mm / 2.0 / self.interslit_mm, self.si_mm / 2.0 / self.interslit_mm)

    def reflectivity(self, q: np.ndarray) -> np.ndarray:
        return parratt(q, DEFAULT_FILM, DEFAULT_SUBSTRATE)


def flux(lam: np.ndarray) -> np.ndarray:
    """A smooth moderator-like spectrum, peaked near 3.5 Å."""
    lam_t = 3.0
    return (lam_t / lam) ** 5 * np.exp(-((lam_t / lam) ** 2)) * 50.0


def trapezoid_pdf(delta: np.ndarray, a: float, b: float) -> np.ndarray:
    """Density of `u + v` with `u ~ U(−a, a)`, `v ~ U(−b, b)`."""
    lo, hi = abs(a - b), a + b
    top = 1.0 / (2.0 * max(a, b))
    d = np.abs(delta)
    ramp = top * (hi - d) / (hi - lo) if hi > lo else np.zeros_like(d)
    return np.where(d <= lo, top, np.where(d < hi, ramp, 0.0))


def detector_y(setup: SyntheticSetup, delta: np.ndarray, h_si: np.ndarray | float = 0.0) -> np.ndarray:
    """Detector position in pixels of a ray at incidence offset `delta` crossing Si at height `h_si` (mm).

    `h_si` enters only the `slits` geometry.
    """
    if setup.geometry == "point":
        return setup.y_ref + np.tan(2.0 * delta) * setup.distance_mm / setup.pixel_width_mm
    length = setup.si_sample_mm + setup.distance_mm
    return setup.y_ref + (h_si + length * np.tan(delta)) / setup.pixel_width_mm


def rows_of(y: np.ndarray) -> np.ndarray:
    """Pixel row `j` spans `[j − 0.5, j + 0.5)`."""
    return np.floor(y + 0.5).astype(int)


def noise_free_events(setup: SyntheticSetup, n_delta: int = 240, lam_step: float = 0.0005) -> tuple[Events, Events]:
    """Dense weighted pseudo-events: a midpoint grid in (δ, λ), each weighted by its expected counts.

    The reflected weights are `p(δ)·dδ · φ(λ)·dλ · R(Q)` and the direct-beam weights `φ(λ)·dλ`, so both
    beams share one normalization of 1.
    """
    a, b = setup.delta_limits()
    edges_delta = np.linspace(-(a + b), a + b, n_delta + 1)
    delta = 0.5 * (edges_delta[1:] + edges_delta[:-1])
    p = trapezoid_pdf(delta, a, b) * np.diff(edges_delta)
    p /= p.sum()
    lam = np.arange(setup.lam_min + lam_step / 2, setup.lam_max, lam_step)
    phi = flux(lam) * lam_step

    q = FOUR_PI * np.sin(setup.theta_central + delta[:, None]) / lam[None, :]
    weight = p[:, None] * phi[None, :] * setup.reflectivity(q.ravel()).reshape(q.shape)
    direct = Events(np.zeros(len(lam), dtype=int), lam, phi)
    if setup.geometry == "slits":
        rows, fraction = _slit_row_fractions(setup, delta)
        summed = fraction.T @ weight  # (row, λ): the rows' shares of every δ, summed over δ
        return Events(np.repeat(rows, len(lam)), np.tile(lam, len(rows)), summed.ravel()), direct
    y = np.broadcast_to(detector_y(setup, delta)[:, None], q.shape)
    if setup.blur_pixels > 0:
        y, weight = _blur(y, weight, setup.blur_pixels)
    reflected = Events(rows_of(np.ravel(y)), np.broadcast_to(lam, y.shape).ravel().copy(), np.ravel(weight))
    return reflected, direct


def poisson_events(
    setup: SyntheticSetup,
    incident: float,
    incident_direct: float,
    rng: np.random.Generator,
    chunk: int = 10_000_000,
) -> tuple[Events, Events, float, float]:
    """Unit-weight events for `incident` reflected-beam and `incident_direct` direct-beam neutrons.

    Returns the two event sets and their normalizations (the incident counts), the role proton charge
    plays for measured data.
    """
    a, b = setup.delta_limits()
    lam_grid = np.linspace(setup.lam_min, setup.lam_max, 20001)
    cdf = np.concatenate([[0.0], np.cumsum(flux(0.5 * (lam_grid[1:] + lam_grid[:-1])) * np.diff(lam_grid))])
    cdf /= cdf[-1]

    def sample_lam(n: int) -> np.ndarray:
        return np.interp(rng.random(n), cdf, lam_grid)

    rows, lams = [], []
    remaining = int(rng.poisson(incident))
    while remaining > 0:  # in chunks, so memory stays bounded whatever `incident` is
        n = min(remaining, chunk)
        remaining -= n
        lam = sample_lam(n)
        if setup.geometry == "point":
            delta = rng.uniform(-a, a, n) + rng.uniform(-b, b, n)
            h_si = np.zeros(n)
        else:
            h_s1 = rng.uniform(-setup.s1_mm / 2, setup.s1_mm / 2, n)
            h_si = rng.uniform(-setup.si_mm / 2, setup.si_mm / 2, n)
            delta = (h_si - h_s1) / setup.interslit_mm  # the trapezoid of `delta_limits`
        q = FOUR_PI * np.sin(setup.theta_central + delta) / lam
        keep = rng.random(n) < setup.reflectivity(q)
        y = detector_y(setup, delta[keep], h_si[keep])
        if setup.blur_pixels > 0:
            y = y + rng.normal(0.0, setup.blur_pixels, len(y))
        rows.append(rows_of(y))
        lams.append(lam[keep])
    rb = Events(np.concatenate(rows), np.concatenate(lams), np.ones(sum(len(x) for x in lams)))

    n_db = rng.poisson(incident_direct)
    db = Events(np.zeros(n_db, dtype=int), sample_lam(n_db), np.ones(n_db))
    return rb, db, float(incident), float(incident_direct)


def bin_averaged_truth(setup: SyntheticSetup, q_edges: np.ndarray, samples_per_bin: int = 64) -> np.ndarray:
    """`(1/ΔQ) ∫ R(Q) dQ` over each bin: what an ideal constant-Q reduction returns."""
    q_edges = np.asarray(q_edges, dtype=float)
    t = (np.arange(samples_per_bin) + 0.5) / samples_per_bin
    q = q_edges[:-1, None] + np.diff(q_edges)[:, None] * t[None, :]
    return setup.reflectivity(q.ravel()).reshape(q.shape).mean(axis=1)


def divergence_averaged_truth(
    setup: SyntheticSetup, q_edges: np.ndarray, samples_per_bin: int = 64, n_delta: int = 240
) -> np.ndarray:
    """What constant TOF returns: R averaged over the incident-angle trapezoid at fixed λ, bin-averaged.

    A λ bin is labeled `Q = 4π sin θ_c/λ`, but the neutron at incidence offset δ has true Q
    `Q · sin(θ_c + δ)/sin θ_c`; constant TOF returns `∫ p(δ) R(that) dδ`.
    """
    q_edges = np.asarray(q_edges, dtype=float)
    a, b = setup.delta_limits()
    edges_delta = np.linspace(-(a + b), a + b, n_delta + 1)
    delta = 0.5 * (edges_delta[1:] + edges_delta[:-1])
    p = trapezoid_pdf(delta, a, b) * np.diff(edges_delta)
    p /= p.sum()
    t = (np.arange(samples_per_bin) + 0.5) / samples_per_bin
    q = (q_edges[:-1, None] + np.diff(q_edges)[:, None] * t[None, :]).ravel()
    scale = np.sin(setup.theta_central + delta) / np.sin(setup.theta_central)
    smeared = setup.reflectivity(np.outer(q, scale).ravel()).reshape(len(q), len(delta)) @ p
    return smeared.reshape(len(q_edges) - 1, samples_per_bin).mean(axis=1)


def _slit_row_fractions(setup: SyntheticSetup, delta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rows reached by the rays at each incidence offset `delta`, and each row's share, under `slits`.

    At fixed δ the rays cross Si uniformly over the heights both slits allow, so they land uniformly over
    an interval of the detector. Each row takes its overlap with that interval, convolved with the
    Gaussian blur when there is one. As in `_blur`, the blur reaches ±3σ, and the shares of each δ are
    renormalized over the rows within that reach. Returns the rows and a `(len(delta), len(rows))` array
    of shares.
    """
    shift = setup.interslit_mm * delta  # Si height minus S1 height
    h_lo = np.maximum(-setup.si_mm / 2, shift - setup.s1_mm / 2)
    h_hi = np.minimum(setup.si_mm / 2, shift + setup.s1_mm / 2)
    y_lo = detector_y(setup, delta, h_lo)
    y_hi = detector_y(setup, delta, h_hi)
    sigma = setup.blur_pixels
    rows = np.arange(rows_of(y_lo.min() - 3.0 * sigma), rows_of(y_hi.max() + 3.0 * sigma) + 1)
    lo_edge, hi_edge = rows[None, :] - 0.5, rows[None, :] + 0.5
    y_lo, y_hi = y_lo[:, None], y_hi[:, None]
    if sigma > 0:

        def ramp(x: np.ndarray) -> np.ndarray:  # ∫ Φ(x/σ) dx
            return x * ndtr(x / sigma) + sigma * np.exp(-0.5 * (x / sigma) ** 2) / np.sqrt(2 * np.pi)

        overlap = ramp(hi_edge - y_lo) - ramp(lo_edge - y_lo) - ramp(hi_edge - y_hi) + ramp(lo_edge - y_hi)
    else:
        overlap = np.clip(np.minimum(y_hi, hi_edge) - np.maximum(y_lo, lo_edge), 0.0, None)
    fraction = overlap / overlap.sum(axis=1, keepdims=True)
    used = fraction.max(axis=0) > 1e-12
    return rows[used], fraction[:, used]


def _blur(y: np.ndarray, weight: np.ndarray, sigma: float, n_points: int = 9) -> tuple[np.ndarray, np.ndarray]:
    """Spread each pseudo-event over Gaussian-weighted copies across ±3σ in y."""
    offsets = np.linspace(-3.0 * sigma, 3.0 * sigma, n_points)
    w = np.exp(-0.5 * (offsets / sigma) ** 2)
    w /= w.sum()
    return (y[None, ...] + offsets[:, None, None], weight[None, ...] * w[:, None, None])
