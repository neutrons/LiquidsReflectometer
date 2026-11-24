"""
Pixel-to-Q mapping methods on a noise-free synthetic measurement: a standalone, numpy-only version.

Builds the reflected-beam histogram r[row, λ] (already divided by the direct beam) for a sample with
R(Q) ≡ 1 and for a 400 Å film, maps it to Q with the methods E, A, B, C and C0, and prints the median
and maximum |deviation| from the truth on the fully covered Q bins.

    python q_mapping_minimal.py

Units: angles in radians, wavelength in Å, Q in 1/Å, lengths in mm.
"""

import numpy as np

FOUR_PI = 4.0 * np.pi

# Geometry of the synthetic measurement
THETA_C = np.radians(0.45)  # central incidence angle
Y_REF = 150.3  # pixel of the central ray
PIXEL, L = 0.70, 1830.0  # pixel width and sample-detector distance
S1, SI, D = 2.0, 1.2, 1390.0  # slit gaps and the distance between the slits
LAM_MIN, LAM_MAX = 2.75, 9.5


def flux(lam):
    return (3.0 / lam) ** 5 * np.exp(-((3.0 / lam) ** 2))


def film(q):
    """Parratt reflectivity of 400 Å of SLD 4e-6 on Si, Névot–Croce roughness 4 Å and 3 Å."""
    sld = np.array([0.0, 4.0e-6, 2.07e-6])
    thickness, sigma = np.array([0.0, 400.0, 0.0]), np.array([4.0, 3.0])
    kz = np.sqrt((q[:, None] / 2) ** 2 - FOUR_PI * sld + 0j)
    x = np.zeros(len(q), dtype=complex)
    for j in (1, 0):
        fresnel = (kz[:, j] - kz[:, j + 1]) / (kz[:, j] + kz[:, j + 1])
        rj = fresnel * np.exp(-2 * kz[:, j] * kz[:, j + 1] * sigma[j] ** 2)
        phase = np.exp(2j * kz[:, j + 1] * thickness[j + 1])
        x = (rj + x * phase) / (1 + rj * x * phase)
    return np.abs(x) ** 2


def pseudo_events(reflectivity, n_delta=240, lam_step=0.0005):
    """A dense grid in (incidence offset δ, λ), each point weighted by its expected counts.

    The incidence offsets follow the trapezoid that two uniformly lit slits give. A ray at offset δ
    lands at `Y_REF + L·tan(2δ)/PIXEL`, the `half` model, so the data and the reduction share one model.
    """
    a, b = S1 / 2 / D, SI / 2 / D
    edges = np.linspace(-(a + b), a + b, n_delta + 1)
    delta = 0.5 * (edges[1:] + edges[:-1])
    p = np.clip((a + b - np.abs(delta)) / (2 * min(a, b)), 0, 1)  # trapezoid, up to a constant
    p /= p.sum()
    lam = np.arange(LAM_MIN + lam_step / 2, LAM_MAX, lam_step)
    phi = flux(lam) * lam_step
    q = FOUR_PI * np.sin(THETA_C + delta[:, None]) / lam[None, :]
    weight = p[:, None] * phi[None, :] * reflectivity(q.ravel()).reshape(q.shape)
    row = np.floor(Y_REF + np.tan(2 * delta) * L / PIXEL + 0.5).astype(int)
    row = np.broadcast_to(row[:, None], q.shape).ravel()
    return row, np.broadcast_to(lam, q.shape).ravel(), weight.ravel(), lam, phi


def cells_from_events(row, lam, weight, db_lam, db_weight, rows, lam_edges):
    """r[i, l] = reflected counts in (rows[i], λ bin l) / direct-beam counts in λ bin l."""
    n_lam = len(lam_edges) - 1
    li = np.searchsorted(lam_edges, lam, side="right") - 1
    ok = (li >= 0) & (li < n_lam)
    counts = np.zeros((len(rows), n_lam))
    np.add.at(counts, (row[ok] - rows[0], li[ok]), weight[ok])
    db_i = np.searchsorted(lam_edges, db_lam, side="right") - 1
    ok = (db_i >= 0) & (db_i < n_lam)
    db = np.bincount(db_i[ok], weights=db_weight[ok], minlength=n_lam)
    return counts / db, db


def half_angles(rows, offset=0.0):
    return THETA_C + np.arctan((rows - Y_REF + offset) * PIXEL / L) / 2


# ----------------------------------------------------------------------------------------------------
# The methods. Each returns R on the bins q_edges, NaN where it has no value.
# ----------------------------------------------------------------------------------------------------


def map_events(row, lam, weight, rows, lam_edges, db, q_edges):
    """E: each event at its own Q, weighted by 1/I_DB(λ)·|∂Q/∂λ|, histogrammed and divided by ΔQ."""
    li = np.searchsorted(lam_edges, lam, side="right") - 1
    ok = (li >= 0) & (li < len(lam_edges) - 1)
    lam, li = lam[ok], li[ok]
    q = FOUR_PI * np.sin(half_angles(rows)[row[ok] - rows[0]]) / lam
    density = db[li] / np.diff(lam_edges)[li]  # direct beam per Å, the same step function as the cells
    w = weight[ok] / density * q / lam
    k = np.digitize(q, q_edges) - 1
    inside = (k >= 0) & (k < len(q_edges) - 1)
    total = np.bincount(k[inside], weights=w[inside], minlength=len(q_edges) - 1)
    n = np.bincount(k[inside], minlength=len(q_edges) - 1)
    return np.where(n > 0, total / np.diff(q_edges), np.nan)


def map_point(r, theta, lam_edges, q_edges):
    """A: each cell to the bin of its center Q; R = N_rows · (sum of cells) / (number of cells)."""
    lam_c = 0.5 * (lam_edges[1:] + lam_edges[:-1])
    k = np.digitize(FOUR_PI * np.sin(theta[:, None]) / lam_c[None, :], q_edges) - 1
    use = (k >= 0) & (k < len(q_edges) - 1)
    total = np.bincount(k[use], weights=r[use], minlength=len(q_edges) - 1)
    n = np.bincount(k[use], minlength=len(q_edges) - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, r.shape[0] * total / n, np.nan)


def map_prototype(r, theta, lam_edges, q_edges, covered):
    """B: equal shares to the Q points inside the cell's corner interval, weighted by |J|·Δλ.

    Points are the bin centers whose Q line is fully covered; θ and Δθ in degrees, as in the prototype.
    """
    points = np.flatnonzero(covered)
    qp = np.sqrt(q_edges[1:] * q_edges[:-1])[points]
    lam = 0.5 * (lam_edges[1:] + lam_edges[:-1])[None, :]
    dlam = np.diff(lam_edges)[None, :]
    th = np.degrees(theta)[:, None]
    half_dth = np.abs(np.diff(np.degrees(theta), append=2 * np.degrees(theta[-1]) - np.degrees(theta[-2])))
    half_dth = half_dth[:, None] / 2
    jac = np.hypot(FOUR_PI / lam * np.cos(np.radians(th)) * np.pi / 180, FOUR_PI / lam**2 * np.sin(np.radians(th)))
    qlo = FOUR_PI * np.sin(np.radians(th - half_dth)) / (lam + dlam / 2)
    qhi = FOUR_PI * np.sin(np.radians(th + half_dth)) / (lam - dlam / 2)
    qcen = FOUR_PI * np.sin(np.radians(th)) / lam
    first = np.searchsorted(qp, qlo, side="left")
    stop = np.searchsorted(qp, qhi, side="right")
    nearest = np.argmin(np.abs(qcen[..., None] - qp), axis=-1)
    empty = stop <= first
    first, stop = np.where(empty, nearest, first), np.where(empty, nearest + 1, stop)
    share = np.broadcast_to(jac * dlam / (stop - first), r.shape)
    use = (qhi <= qp[-1]) & (qlo >= qp[0])
    signal, jsum = np.zeros(len(qp) + 1), np.zeros(len(qp) + 1)
    for values, out in ((r * share, signal), (share, jsum)):
        np.add.at(out, first[use], values[use])
        np.add.at(out, stop[use], -values[use])
    signal, jsum = np.cumsum(signal[:-1]), np.cumsum(jsum[:-1])
    result = np.full(len(q_edges) - 1, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        result[points] = np.where(jsum > 1e-12 * jsum.max(), signal / (jsum / r.shape[0]), np.nan)
    return result


def map_overlap(r, theta, theta_lo, theta_hi, lam_edges, q_edges):
    """C (corner footprint) or C0 (theta_lo = theta_hi = theta): spread r·w_λ uniformly, share by overlap."""
    lam_lo, lam_hi = lam_edges[None, :-1], lam_edges[None, 1:]
    w_lam = FOUR_PI * np.sin(theta)[:, None] * (1 / lam_lo - 1 / lam_hi)  # cell's Q width at the row angle
    lo = FOUR_PI * np.sin(theta_lo)[:, None] / lam_hi
    hi = FOUR_PI * np.sin(theta_hi)[:, None] / lam_lo
    density = r * w_lam / (hi - lo)
    out = np.zeros(len(q_edges) - 1)
    for k in range(len(out)):
        overlap = np.clip(np.minimum(hi, q_edges[k + 1]) - np.maximum(lo, q_edges[k]), 0, None)
        out[k] = np.sum(density * overlap)
    return np.where(out > 0, out / np.diff(q_edges), np.nan)


# ----------------------------------------------------------------------------------------------------


def compare(name, reflectivity, truth_of, dq_over_q, lam_steps=(0.005, 0.02, 0.1)):
    row, lam, weight, db_lam, db_weight = pseudo_events(reflectivity)
    rows = np.arange(row.min(), row.max() + 1)
    theta, theta_lo, theta_hi = half_angles(rows), half_angles(rows, -0.5), half_angles(rows, 0.5)
    q_lo = FOUR_PI * np.sin(THETA_C) / LAM_MAX * 0.9
    n = int(np.log(FOUR_PI * np.sin(THETA_C) / LAM_MIN * 1.1 / q_lo) / np.log(1 + dq_over_q)) + 1
    q_edges = q_lo * (1 + dq_over_q) ** np.arange(n)
    q_c = np.sqrt(q_edges[1:] * q_edges[:-1])
    # A Q bin is fully covered when its constant-Q line stays inside the λ window at every row.
    covered = (FOUR_PI * np.sin(theta.min()) / q_c >= LAM_MIN) & (FOUR_PI * np.sin(theta.max()) / q_c <= LAM_MAX)
    truth = truth_of(q_edges)
    print(f"\n{name}, ΔQ/Q = {dq_over_q:.1%}: median |deviation| (max) on the fully covered bins")
    print(f"{'Δλ (Å)':>8} " + " ".join(f"{m:>16}" for m in ("E", "A", "B", "C", "C0")))
    for step in lam_steps:
        lam_edges = np.arange(LAM_MIN, LAM_MAX + 0.5 * step, step)
        lam_edges[-1] = min(lam_edges[-1], LAM_MAX)
        r, db = cells_from_events(row, lam, weight, db_lam, db_weight, rows, lam_edges)
        curves = {
            "E": map_events(row, lam, weight, rows, lam_edges, db, q_edges),
            "A": map_point(r, theta, lam_edges, q_edges),
            "B": map_prototype(r, theta, lam_edges, q_edges, covered),
            "C": map_overlap(r, theta, theta_lo, theta_hi, lam_edges, q_edges),
            "C0": map_overlap(r, theta, theta, theta, lam_edges, q_edges),
        }
        mask = covered & np.all([np.isfinite(c) for c in curves.values()], axis=0)
        cols = []
        for c in curves.values():
            dev = np.abs(c[mask] / truth[mask] - 1) * 100
            cols.append(f"{np.median(dev):6.3f}% ({dev.max():5.2f}%)")
        print(f"{step:>8} " + " ".join(f"{s:>16}" for s in cols) + f"   [{mask.sum()} bins]")


def bin_average(reflectivity, q_edges, samples=64):
    t = (np.arange(samples) + 0.5) / samples
    q = q_edges[:-1, None] + np.diff(q_edges)[:, None] * t[None, :]
    return reflectivity(q.ravel()).reshape(q.shape).mean(axis=1)


if __name__ == "__main__":
    flat = np.ones_like
    compare("R(Q) ≡ 1", flat, lambda e: np.ones(len(e) - 1), dq_over_q=0.005)
    compare("400 Å film", film, lambda e: bin_average(film, e), dq_over_q=0.01)
