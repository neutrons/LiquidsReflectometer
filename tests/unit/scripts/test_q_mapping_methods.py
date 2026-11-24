"""Tests for the pixel-to-Q comparison harness in scripts/test (q_mapping_methods, q_methods_synthetic)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts" / "test"))

import q_mapping_methods as m  # noqa: E402
import q_methods_synthetic as s  # noqa: E402

PIXEL, DISTANCE = 0.7, 1830.0


def _cells(r: np.ndarray, lam_edges: np.ndarray, rows: np.ndarray, var: np.ndarray | None = None) -> m.Cells:
    var = np.full_like(r, 1e-4) if var is None else var
    return m.Cells(rows=rows, lam_edges=lam_edges, r=r, var=var, db=np.ones(len(lam_edges) - 1))


# --------------------------------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------------------------------


def test_angle_offset_conventions():
    full = m.angle_offset(np.array([2.0]), PIXEL, DISTANCE, "full")
    half = m.angle_offset(np.array([2.0]), PIXEL, DISTANCE, "half")
    assert full[0] == pytest.approx(np.arctan(1.4 / 1830.0))
    assert half[0] == pytest.approx(full[0] / 2)
    with pytest.raises(ValueError, match="convention"):
        m.angle_offset(np.array([1.0]), PIXEL, DISTANCE, "double")


def test_row_angles_edges_bracket_center_for_either_sign():
    rows = np.arange(145, 156)
    for sign in (1.0, -1.0):
        angles = m.row_angles(rows, 0.01, 150.3, PIXEL, DISTANCE, "half", sign=sign)
        assert np.all(angles.lo < angles.center) and np.all(angles.center < angles.hi)
    up = m.row_angles(rows, 0.01, 150.3, PIXEL, DISTANCE, "half", sign=1.0)
    assert up.center[-1] > up.center[0]


def test_histogram_cells_divides_by_direct_beam_and_flags_empty_bins():
    rb = m.Events(row=np.array([1, 1, 2, 2]), wavelength=np.array([3.1, 3.2, 3.1, 4.5]), weight=np.ones(4))
    db = m.Events(row=np.zeros(5, dtype=int), wavelength=np.array([3.1, 3.3, 3.4, 3.6, 3.9]), weight=np.ones(5))
    cells = m.histogram_cells(rb, db, rows=np.array([1, 2]), lam_edges=np.array([3.0, 4.0, 5.0]), reflected_norm=2.0)
    # Row 1: 2 events / norm 2 over 5 direct-beam events; the second bin has no direct beam.
    assert cells.r[0, 0] == pytest.approx(0.2)
    assert cells.var[0, 0] == pytest.approx(2 / 4 / 25)
    assert cells.r[1, 0] == pytest.approx(0.1)
    assert not cells.valid[0, 1]
    assert cells.r[1, 1] == 0.0


def test_qline_coverage_inside_and_truncated():
    angles = m.row_angles(np.arange(148, 153), 0.01, 150.0, PIXEL, DISTANCE, "half")
    theta_min, theta_max = angles.center.min(), angles.center.max()
    q = np.array([4 * np.pi * np.sin(theta_max) / 5.0])
    # The line spans λ ∈ [5·sinθmin/sinθmax, 5]: fully inside [2, 9], and cut in half at its middle.
    assert m.qline_coverage(q, angles, 2.0, 9.0)[0] == pytest.approx(1.0)
    lam_lo = 5.0 * np.sin(theta_min) / np.sin(theta_max)
    middle = 0.5 * (lam_lo + 5.0)
    assert m.qline_coverage(q, angles, middle, 9.0)[0] == pytest.approx(0.5)


def test_log_q_edges_matches_legacy_rule():
    edges = m.log_q_edges(0.01, 0.1, 0.02)
    assert edges[0] == 0.01
    assert np.allclose(edges[1:] / edges[:-1], 1.02)
    assert edges[-1] <= 0.1 < edges[-1] * 1.02


# --------------------------------------------------------------------------------------------------
# Method A against a verbatim copy of Mantid MagnetismReflectometryReduction.constant_q
# --------------------------------------------------------------------------------------------------


def _mantid_constant_q_loop(signal, signal_err, qz, axis_z, trim_factor):
    """Lines 404-433 of Mantid 6.16 MagnetismReflectometryReduction.py; signal[ix][tof], qz[tof][ix]."""
    refl = np.zeros(len(axis_z) - 1)
    refl_err = np.zeros(len(axis_z) - 1)
    signal_n = np.zeros(len(axis_z) - 1)
    for tof in range(qz.shape[0]):
        z_inds = np.digitize(qz[tof], axis_z)
        z_inds -= 1
        for ix in range(signal.shape[0]):
            if z_inds[ix] < len(axis_z) - 1 and z_inds[ix] >= 0 and not np.isnan(signal[ix][tof]):
                refl[z_inds[ix]] += signal[ix][tof]
                refl_err[z_inds[ix]] += signal_err[ix][tof] * signal_err[ix][tof]
                signal_n[z_inds[ix]] += 1.0
    signal_n = np.where(signal_n > 0, signal_n, 1)
    refl = float(signal.shape[0]) * refl / signal_n
    refl_err = float(signal.shape[0]) * np.sqrt(refl_err) / signal_n
    for i in range(len(refl)):
        if signal_n[i] < np.max(signal_n) * trim_factor:
            refl[i] = 0.0
            refl_err[i] = 0.0
    return refl, refl_err


def test_map_point_matches_mantid_loop():
    rng = np.random.default_rng(1)
    rows = np.arange(145, 156)
    lam_edges = np.linspace(3.0, 9.0, 301)
    r = rng.uniform(0.0, 0.1, (len(rows), 300))
    err = rng.uniform(0.001, 0.01, r.shape)
    cells = _cells(r, lam_edges, rows, var=err**2)
    angles = m.row_angles(rows, np.radians(0.6), 150.4, PIXEL, DISTANCE, "half")
    q_edges = m.log_q_edges(0.014, 0.05, 0.01)

    got = m.map_point(cells, angles, q_edges, trim=0.5)

    qz = (4 * np.pi * np.sin(angles.center[:, None]) / cells.lam_centers[None, :]).T
    refl, refl_err = _mantid_constant_q_loop(r, err, qz, q_edges, 0.5)
    kept = np.isfinite(got.r)
    assert kept.sum() > 50
    np.testing.assert_allclose(got.r[kept], refl[kept], rtol=1e-12)
    np.testing.assert_allclose(got.dr[kept], refl_err[kept], rtol=1e-12)
    assert np.all(refl[~kept] == 0.0)


# --------------------------------------------------------------------------------------------------
# Method B against a verbatim copy of the prototype's Q loop (LR-exp nr_reduction_calc.py)
# --------------------------------------------------------------------------------------------------


def _prototype_loop(Rarr, REarr, Theta, LAMBDA, LambdaBinSize, q_all, lam_window, qline_threshold):  # noqa: N803
    """`_reduce_single_run`, constantq branch with useGravity False, from the Q-line mask to `r / FAC`.

    `self.q` → `q_all`, `self.config.LambdaMin/MaxUse` → `lam_window`, `tools.get_edges` → uniform
    edges; `_calculate_jacobian` and `_pixel_to_q` are inlined. dq is not compared.
    """
    ThetaBinSize = abs(np.diff(Theta))  # noqa: N806
    ThetaBinSize = np.concatenate([ThetaBinSize, [ThetaBinSize[-1]]])  # noqa: N806
    lmin = 4 * np.pi * np.sin(np.radians(min(abs(Theta)))) / q_all
    lmax = 4 * np.pi * np.sin(np.radians(max(abs(Theta)))) / q_all
    Qline_fraction = (np.minimum(lam_window[1], lmax) - np.maximum(lam_window[0], lmin)) / (lmax - lmin)  # noqa: N806
    mask = Qline_fraction >= qline_threshold
    q_vals = q_all[mask]
    Qline_fraction = Qline_fraction[mask]  # noqa: N806

    r = q_vals * 0
    dr = q_vals * 0
    Jsum = q_vals * 0  # noqa: N806
    dLambda_bins = np.full(len(LAMBDA), LambdaBinSize)  # noqa: N806
    ThetaGC = np.zeros(len(LAMBDA))  # noqa: N806
    for T in range(Rarr.shape[0]):  # noqa: N806
        Thv = abs(Theta[T] + ThetaGC)  # noqa: N806
        dqdtheta = (4 * np.pi / LAMBDA) * np.cos(np.radians(Thv)) * np.pi / 180
        dqdlambda = -(4 * np.pi / LAMBDA**2) * np.sin(np.radians(Thv))
        J = abs(np.sqrt(dqdtheta**2 + dqdlambda**2))  # noqa: N806
        theta_bin = ThetaBinSize[T]
        for L in range(Rarr.shape[1]):  # noqa: N806
            qcen = 4 * np.pi * np.sin(np.radians(Thv[L])) / LAMBDA[L]
            qlo = 4 * np.pi * np.sin(np.radians(Thv[L] - theta_bin / 2)) / (LAMBDA[L] + LambdaBinSize / 2)
            qhi = 4 * np.pi * np.sin(np.radians(Thv[L] + theta_bin / 2)) / (LAMBDA[L] - LambdaBinSize / 2)
            if (qhi <= max(q_vals)) & (qlo >= min(q_vals)):
                idx = np.flatnonzero((q_vals >= qlo) & (q_vals <= qhi))
                if idx.size == 0:
                    idx = np.array([np.argmin(np.abs(qcen - q_vals))])
                wt = J[L] / idx.size
                dl = dLambda_bins[L]
                r[idx] = r[idx] + Rarr[T, L] * wt * dl
                dr[idx] = dr[idx] + (REarr[T, L] * wt * dl) ** 2
                Jsum[idx] = Jsum[idx] + wt * dl
    FAC = Jsum / Rarr.shape[0]  # noqa: N806
    keep = np.isfinite(r) & (FAC != 0)
    out_r = np.full(len(q_vals), np.nan)
    out_dr = np.full(len(q_vals), np.nan)
    out_r[keep] = r[keep] / FAC[keep] * Qline_fraction[keep]
    out_dr[keep] = np.sqrt(dr[keep]) / FAC[keep] * Qline_fraction[keep]
    return mask, out_r, out_dr


@pytest.mark.parametrize("lam_step", [0.01, 0.05])
@pytest.mark.parametrize("qline_threshold", [1.0, 0.5])
def test_map_prototype_matches_prototype_loop(lam_step, qline_threshold):
    rng = np.random.default_rng(2)
    rows = np.arange(146, 155)
    lam_edges = np.arange(3.0, 9.0 + 1e-9, lam_step)
    n_lam = len(lam_edges) - 1
    r = rng.uniform(0.0, 0.1, (len(rows), n_lam))
    err = rng.uniform(0.001, 0.01, r.shape)
    cells = _cells(r, lam_edges, rows, var=err**2)
    angles = m.row_angles(rows, np.radians(0.6), 150.4, PIXEL, DISTANCE, "full")
    q_edges = m.log_q_edges(0.012, 0.06, 0.01)

    got = m.map_prototype(cells, angles, q_edges, qline_threshold=qline_threshold)

    mask, ref_r, ref_dr = _prototype_loop(
        r,
        err,
        np.degrees(angles.center),
        cells.lam_centers,
        lam_step,
        got.q_centers,
        (lam_edges[0], lam_edges[-1]),
        qline_threshold,
    )
    assert np.all(np.isnan(got.r[~mask]))
    assert np.isfinite(ref_r).sum() > 20
    np.testing.assert_allclose(got.r[mask], ref_r, rtol=1e-9, equal_nan=True)
    np.testing.assert_allclose(got.dr[mask], ref_dr, rtol=1e-9, equal_nan=True)


# --------------------------------------------------------------------------------------------------
# Method C and the event reference E
# --------------------------------------------------------------------------------------------------


def _uniform_rows():
    rows = np.arange(148, 153)
    lam_edges = np.linspace(3.0, 9.0, 241)
    row_fraction = np.array([0.1, 0.2, 0.4, 0.2, 0.1])
    cells = _cells(np.repeat(row_fraction[:, None], 240, axis=1) * 0.5, lam_edges, rows)
    angles = m.row_angles(rows, np.radians(0.6), 150.0, PIXEL, DISTANCE, "half")
    return cells, angles


def test_map_overlap_without_angular_width_is_exact_on_uniform_rows():
    cells, angles = _uniform_rows()
    got = m.map_overlap(cells, angles, m.log_q_edges(0.016, 0.04, 0.01), include_angular_width=False)
    full = got.diagnostics["row_coverage"] > 1 - 1e-9
    assert full.sum() > 30
    np.testing.assert_allclose(got.r[full], 0.5, rtol=1e-9)


@pytest.mark.parametrize("include_angular_width", [False, True])
def test_map_overlap_conserves_total(include_angular_width):
    # With the angular width each cell is a box over its corner footprint, wider than the cell
    # spacing, so a uniform input comes out with a ripple; the total is conserved either way.
    cells, angles = _uniform_rows()
    q_edges = m.log_q_edges(0.005, 0.1, 0.01)  # contains every footprint
    got = m.map_overlap(cells, angles, q_edges, include_angular_width=include_angular_width)
    q_width = (
        4 * np.pi * np.sin(angles.center[:, None]) * (1 / cells.lam_edges[None, :-1] - 1 / cells.lam_edges[None, 1:])
    )
    expected = np.sum(cells.r * q_width)
    assert np.nansum(got.r * np.diff(q_edges)) == pytest.approx(expected, rel=1e-12)


def test_map_overlap_drops_parts_outside_grid():
    rows = np.array([150])
    lam_edges = np.array([4.0, 5.0])
    cells = _cells(np.array([[1.0]]), lam_edges, rows)
    angles = m.row_angles(rows, np.radians(0.6), 150.0, PIXEL, DISTANCE, "half")
    q_lo = 4 * np.pi * np.sin(angles.center[0]) / 5.0
    q_hi = 4 * np.pi * np.sin(angles.center[0]) / 4.0
    q_edges = np.array([q_lo, 0.5 * (q_lo + q_hi)])  # only the cell's lower half is on the grid
    got = m.map_overlap(cells, angles, q_edges, include_angular_width=False)
    assert got.r[0] == pytest.approx(1.0)


def test_map_events_constant_reflectivity():
    # Evenly spaced events: E maps each event against a step-function direct beam, so it is exact
    # only when the events are uniform within each wavelength bin.
    n = 600_000
    lam = 3.0 + 6.0 * (np.arange(n) + 0.5) / n
    lam_edges = np.linspace(3.0, 9.0, 61)
    rows = np.array([150])
    direct = m.Events(np.zeros(len(lam), dtype=int), lam, np.ones(len(lam)))
    reflected = m.Events(np.full(len(lam), 150), lam, np.full(len(lam), 0.25))
    cells = m.histogram_cells(reflected, direct, rows, lam_edges)
    angles = m.row_angles(rows, np.radians(0.6), 150.0, PIXEL, DISTANCE, "half")
    q_edges = m.log_q_edges(0.016, 0.04, 0.02)
    got = m.map_events(reflected, cells, angles, q_edges)
    covered = np.isfinite(got.r)
    np.testing.assert_allclose(got.r[covered][2:-2], 0.25, rtol=1e-3)


# --------------------------------------------------------------------------------------------------
# Synthetic generator
# --------------------------------------------------------------------------------------------------


def test_parratt_single_interface_is_fresnel():
    sld = 2.07e-6
    q = np.linspace(0.002, 0.2, 400)
    got = s.parratt(q, [], s.Layer(sld, 0.0, 0.0))
    qc = 4 * np.sqrt(np.pi * sld)
    kz = q / 2
    ks = np.sqrt(kz**2 - (qc / 2) ** 2 + 0j)
    fresnel = np.abs((kz - ks) / (kz + ks)) ** 2
    np.testing.assert_allclose(got, fresnel, rtol=1e-12)
    assert np.all(got[q < qc] == pytest.approx(1.0))


def test_trapezoid_pdf_normalized():
    a, b = 7.0e-4, 4.0e-4
    d = np.linspace(-(a + b), a + b, 20001)
    assert np.trapezoid(s.trapezoid_pdf(d, a, b), d) == pytest.approx(1.0, rel=1e-6)


def test_event_reference_recovers_truth_on_noise_free_data():
    setup = s.SyntheticSetup()
    reflected, direct = s.noise_free_events(setup, n_delta=120, lam_step=0.002)
    rows = np.arange(reflected.row.min(), reflected.row.max() + 1)
    lam_edges = np.arange(setup.lam_min, setup.lam_max + 1e-9, 0.01)
    cells = m.histogram_cells(reflected, direct, rows, lam_edges)
    angles = m.row_angles(rows, setup.theta_central, setup.y_ref, setup.pixel_width_mm, setup.distance_mm, "half")
    q_edges = m.log_q_edges(0.009, 0.045, 0.01)
    covered = m.qline_coverage(np.sqrt(q_edges[1:] * q_edges[:-1]), angles, lam_edges[0], lam_edges[-1]) >= 1
    truth = s.bin_averaged_truth(setup, q_edges)

    events = m.map_events(reflected, cells, angles, q_edges)
    overlap = m.map_overlap(cells, angles, q_edges, include_angular_width=False)

    assert covered.sum() > 60
    assert np.median(np.abs(events.r[covered] / truth[covered] - 1)) < 0.005
    assert np.median(np.abs(overlap.r[covered] / events.r[covered] - 1)) < 0.005


def test_constant_tof_returns_the_divergence_averaged_reflectivity():
    setup = s.SyntheticSetup()
    reflected, direct = s.noise_free_events(setup, n_delta=120, lam_step=0.002)
    rows = np.arange(reflected.row.min(), reflected.row.max() + 1)
    lam_edges = np.arange(setup.lam_min, setup.lam_max + 1e-9, 0.01)
    cells = m.histogram_cells(reflected, direct, rows, lam_edges)
    q_edges = m.log_q_edges(0.009, 0.045, 0.01)
    q_c = np.sqrt(q_edges[1:] * q_edges[:-1])
    lam_c = 4 * np.pi * np.sin(setup.theta_central) / q_c
    inside = (lam_c > lam_edges[0] * 1.02) & (lam_c < lam_edges[-1] * 0.98)

    got = m.map_constant_tof(cells, setup.theta_central, q_edges)
    smeared = s.divergence_averaged_truth(setup, q_edges)
    truth = s.bin_averaged_truth(setup, q_edges)

    assert inside.sum() > 60
    assert np.median(np.abs(got.r[inside] / smeared[inside] - 1)) < 0.003
    # The film's fringes are washed out by the angular spread, so the unsmeared truth is missed.
    assert np.max(np.abs(got.r[inside] / truth[inside] - 1)) > 0.1


def test_row_alignment_prefers_the_convention_the_data_were_made_with():
    setup = s.SyntheticSetup()
    reflected, direct = s.noise_free_events(setup, n_delta=120, lam_step=0.002)
    rows = np.arange(reflected.row.min(), reflected.row.max() + 1)
    lam_edges = np.arange(setup.lam_min, setup.lam_max + 1e-9, 0.01)
    cells = m.histogram_cells(reflected, direct, rows, lam_edges)
    q_edges = m.log_q_edges(0.009, 0.045, 0.01)
    rms = {}
    for convention in m.ANGLE_CONVENTIONS:
        angles = m.row_angles(
            rows, setup.theta_central, setup.y_ref, setup.pixel_width_mm, setup.distance_mm, convention
        )
        rms[convention] = m.row_alignment(cells, angles, q_edges).rms_relative
    assert rms["half"] < 0.01
    assert rms["full"] > 10 * rms["half"]


# --------------------------------------------------------------------------------------------------
# The meantheta angle model against a verbatim copy of the prototype (LR-exp)
# --------------------------------------------------------------------------------------------------


def _prototype_beam_geometry(si_H, s1_H, d_s1_si, d_si_sam, d_sam_det):  # noqa: N803
    """`nr_tools.calc_beam_geometry_from_slits` with `radians=False`."""
    L = d_si_sam + d_sam_det  # noqa: N806
    a_plus = np.arctan((s1_H + si_H) / (2 * d_s1_si))
    a_minus = np.arctan((s1_H - si_H) / (2 * d_s1_si))
    a = np.array([a_plus, a_minus, -a_plus, -a_minus])
    y = np.array(
        [
            si_H / 2 + L * np.tan(a[0]),
            -si_H / 2 + L * np.tan(a[1]),
            -si_H / 2 + L * np.tan(a[2]),
            si_H / 2 + L * np.tan(a[3]),
        ]
    )
    a = np.degrees(a)
    m_ = (y - np.roll(y, -1)) / (a - np.roll(a, -1))
    b = y - m_ * a
    return a, y, m_, b


def _prototype_intersect(m1, b1, m2, b2):
    x = (b2 - b1) / (m1 - m2)
    return x, m1 * x + b1


def _prototype_detector_convolution(m_, b, ResFn, sigma):  # noqa: N803
    """`NR_Reduction._calc_detector_convolution` without the plot."""
    from scipy.ndimage import gaussian_filter1d, uniform_filter1d

    verts = np.array(
        [
            _prototype_intersect(m_[0], b[0], m_[1], b[1]),
            _prototype_intersect(m_[1], b[1], m_[2], b[2]),
            _prototype_intersect(m_[2], b[2], m_[3], b[3]),
            _prototype_intersect(m_[3], b[3], m_[0], b[0]),
        ]
    )
    if ResFn == "rectangular":
        pad = sigma * np.sqrt(12)
        width = sigma * np.sqrt(12)
    if ResFn == "gaussian":
        pad = sigma * 4.0
    nt = 1000
    ny = 1000
    tstep = ((max(verts[:, 0])) - (min(verts[:, 0]))) / nt
    ystep = ((max(verts[:, 1]) + pad) - (min(verts[:, 1]) - pad)) / ny
    tvec = np.linspace(min(verts[:, 0]), max(verts[:, 0]) + tstep, nt)
    yvec = np.linspace(min(verts[:, 1]) - pad, max(verts[:, 1]) + pad + ystep, ny)
    smear = np.zeros((ny, nt))
    for i in range(nt):
        hi = np.minimum(tvec[i] * m_[3] + b[3], tvec[i] * m_[2] + b[2])
        lo = np.maximum(tvec[i] * m_[0] + b[0], tvec[i] * m_[1] + b[1])
        p = np.where((yvec >= lo) & (yvec <= hi))
        d = yvec * 0
        d[p] = d[p] + 1
        if ResFn == "rectangular":
            smear[:, i] = uniform_filter1d(d, size=int(width / ystep))
        if ResFn == "gaussian":
            smear[:, i] = gaussian_filter1d(d, sigma=sigma / ystep)
    mask = np.any(smear != 0, axis=1)
    return yvec[mask], smear[mask], tvec


def _prototype_mean_theta(yvalues, si_H, s1_H, d_sam_det, d_si_sam, d_s1_si, sigma, ResFn):  # noqa: N803
    """`NR_Reduction.calc_mean_theta` without the figures; angles in degrees."""
    _, _, m_, b = _prototype_beam_geometry(si_H, s1_H, d_s1_si, d_si_sam, d_sam_det)
    yvec, smear, tvec = _prototype_detector_convolution(m_, b, ResFn, sigma)
    mean_theta, sigma_theta = [], []
    for Y in yvalues:  # noqa: N806
        pos = np.where(abs(Y - yvec) == min(abs(Y - yvec)))
        mu = np.sum(smear[pos, :] * tvec) / np.sum(smear[pos, :])
        var = np.sum(smear[pos, :] * (tvec - mu) ** 2) / np.sum(smear[pos, :])
        mean_theta.append(mu)
        sigma_theta.append(np.sqrt(var))
    return np.array(mean_theta), np.array(sigma_theta)


@pytest.mark.parametrize("det_res_fn", ["rectangular", "gaussian"])
@pytest.mark.parametrize(("si", "s1"), [(0.976, 1.5205), (0.25, 0.394), (2.0, 1.2)])
def test_mean_theta_offset_matches_prototype(det_res_fn, si, s1):
    geometry = m.SlitGeometry(
        s1=s1, si=si, interslit=1350.0, si_sample=135.0, sample_detector=DISTANCE, det_sigma=0.8, det_res_fn=det_res_fn
    )
    dy = (np.arange(-9, 10) + 0.26) * PIXEL
    mean, sigma = m.mean_theta_offset(dy, geometry)
    ref_mean, ref_sigma = _prototype_mean_theta(dy, si, s1, DISTANCE, 135.0, 1350.0, 0.8, det_res_fn)
    np.testing.assert_allclose(np.degrees(mean), ref_mean, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(np.degrees(sigma), ref_sigma, rtol=1e-9, atol=1e-12)


def test_mean_theta_offset_is_independent_of_the_length_unit():
    mm = m.SlitGeometry(1.52, 0.98, 1350.0, 135.0, DISTANCE, det_sigma=0.8, det_res_fn="rectangular")
    metre = m.SlitGeometry(1.52e-3, 0.98e-3, 1.35, 0.135, DISTANCE / 1000, det_sigma=0.8e-3, det_res_fn="rectangular")
    dy = (np.arange(-6, 7) + 0.4) * PIXEL
    np.testing.assert_allclose(m.mean_theta_offset(dy / 1000, metre)[0], m.mean_theta_offset(dy, mm)[0], rtol=1e-9)


def test_mean_theta_without_detector_response_lies_between_the_point_source_limits():
    # Every ray that reaches the detector at dy crossed Si and S1; for a narrow Si the mean slope tends
    # to 1/(d_Si,sample + L), for a narrow S1 to 1/(D + d_Si,sample + L).
    geometry = m.SlitGeometry(s1=1.52, si=0.98, interslit=1350.0, si_sample=135.0, sample_detector=DISTANCE)
    dy = np.linspace(-0.2, 0.2, 5)
    slope = np.polyfit(dy, m.mean_theta_offset(dy, geometry)[0], 1)[0]
    assert 1 / (1350.0 + 135.0 + DISTANCE) < slope < 1 / (135.0 + DISTANCE)


def test_mean_theta_is_flat_when_the_detector_response_is_wider_than_the_beam():
    geometry = m.SlitGeometry(0.394, 0.25, 1350.0, 135.0, DISTANCE, det_sigma=0.8, det_res_fn="rectangular")
    mean = m.mean_theta_offset(np.linspace(-0.3, 0.3, 7), geometry)[0]
    assert np.ptp(mean) < 1e-12


def test_mean_theta_row_angles_follow_the_sign():
    geometry = m.SlitGeometry(1.52, 0.98, 1350.0, 135.0, DISTANCE, det_sigma=0.8, det_res_fn="rectangular")
    rows = np.arange(146, 155)
    up = m.mean_theta_row_angles(rows, 0.04, 150.0, PIXEL, geometry, sign=1.0)
    down = m.mean_theta_row_angles(rows, 0.04, 150.0, PIXEL, geometry, sign=-1.0)
    np.testing.assert_allclose(up.center - 0.04, -(down.center - 0.04), atol=1e-15)
    # Odd about y_ref, up to the prototype's grid, whose angle axis runs one step past the last corner.
    np.testing.assert_allclose(up.center - 0.04, -(up.center[::-1] - 0.04), atol=5e-6)
    assert np.all(np.diff(up.center) >= 0)
    assert np.all(up.lo <= up.center) and np.all(up.center <= up.hi)


# --------------------------------------------------------------------------------------------------
# Slit-traced synthetic geometry
# --------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("blur", [0.0, 1.0])
def test_slit_geometry_keeps_the_point_geometry_total(blur):
    point = s.noise_free_events(s.SyntheticSetup(blur_pixels=blur), n_delta=120, lam_step=0.002)[0]
    slits = s.noise_free_events(s.SyntheticSetup(blur_pixels=blur, geometry="slits"), n_delta=120, lam_step=0.002)[0]
    assert slits.weight.sum() == pytest.approx(point.weight.sum(), rel=1e-12)


def test_slit_geometry_poisson_rows_match_noise_free():
    setup = s.SyntheticSetup(geometry="slits", blur_pixels=0.5)
    expected = s.noise_free_events(setup, n_delta=240, lam_step=0.002)[0]
    rng = np.random.default_rng(3)
    drawn = s.poisson_events(setup, 2e6, 1.0, rng)[0]
    rows = np.arange(expected.row.min(), expected.row.max() + 1)
    share = np.bincount(expected.row - rows[0], weights=expected.weight, minlength=len(rows))
    share /= share.sum()
    counts = np.bincount(np.clip(drawn.row - rows[0], 0, len(rows) - 1), minlength=len(rows))
    big = share > 0.01
    pull = (counts[big] - share[big] * counts.sum()) / np.sqrt(share[big] * counts.sum())
    assert np.max(np.abs(pull)) < 5


def test_row_alignment_prefers_mean_theta_on_slit_traced_data():
    setup = s.SyntheticSetup(geometry="slits")
    reflected, direct = s.noise_free_events(setup, n_delta=120, lam_step=0.002)
    rows = np.arange(reflected.row.min(), reflected.row.max() + 1)
    lam_edges = np.arange(setup.lam_min, setup.lam_max + 1e-9, 0.01)
    cells = m.histogram_cells(reflected, direct, rows, lam_edges)
    q_edges = m.log_q_edges(0.009, 0.045, 0.01)
    args = (rows, setup.theta_central, setup.y_ref, setup.pixel_width_mm)
    rms = {
        c: m.row_alignment(cells, m.row_angles(*args, setup.distance_mm, c), q_edges).rms_relative
        for c in ("half", "full")
    }
    rms["mean_theta"] = m.row_alignment(
        cells, m.mean_theta_row_angles(*args, setup.slit_geometry()), q_edges
    ).rms_relative
    assert rms["mean_theta"] < rms["half"] / 5
    assert rms["mean_theta"] < rms["full"] / 5
