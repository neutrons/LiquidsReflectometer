"""
Compare the pixel-to-Q mapping methods A, B and C (Requirements §2.3.9) on one common R(y, λ).

Methods (see `q_mapping_methods.py`):

    A   point assignment            Mantid MagnetismReflectometryReduction, ConstantQBinning=True
    B   prototype spread            LR-exp NR_Reduction, Jacobian-weighted equal share
    C   fractional overlap          Cubitt et al. 2015 / ReflectometrySumInQ, cell corner footprint
    C0  fractional overlap          C at the row angle only (the histogram limit of E)
    E   event reference             legacy specular_weighted(q_summing=True) estimator, same I_DB(λ)

Every method gets the same cells and the same per-row angles, so differences between them come from
the mapping alone. The angle model is a separate switch (`--convention half|full|mean_theta|both|all`):
`half` is legacy's arctan(Δy·p/L)/2, `full` the prototype's `constantq` arctan(Δy·p/L), and
`mean_theta` the prototype's `meantheta`, the mean incidence angle of the neutrons reaching each row
given the slits and the detector response. `both` is `half` and `full`; `all` adds `mean_theta`. The
per-row alignment check (`row_alignment`) reports, per convention, how well the rows agree in Q, which
real data can use to decide between them.

Two inputs:

* `synthetic` — a known Parratt R(Q) rendered as events; the reference is the bin-averaged truth.
  `--geometry point` (the default) places each ray by its angle alone under the `half` relation;
  `--geometry slits` traces rays through both slits, so each row holds a spread of angles.
  Noise-free by default; `--repeats N` draws N Poisson data sets instead and adds the bias and the
  pull width (std of (R − truth)/dR, 1 for honest error bars).
* `real` — a run reduced by legacy `EventReflectivity` with background, dead time and (by default)
  gravity off; the harness takes the same events through legacy's own TOF→λ, so geometry and
  wavelengths match legacy exactly. The reference is E; legacy's own curve is reported alongside.

Usage:
    pixi run python scripts/test/compare_q_methods.py synthetic --lam-step 0.005 0.02 0.1 --output out/
    pixi run python scripts/test/compare_q_methods.py synthetic --repeats 20 --incident 2e8
    pixi run python scripts/test/compare_q_methods.py real --run 198415 \
        --template tests/data/template.xml --nexus-dir tests/data/liquidsreflectometer-data/nexus \
        --lam-step 0.01 0.05 --output out/

Outputs (with --output): `<input>_summary.csv` (one row per method, convention and λ step),
`<input>_shared_bins_summary.csv` (the same on the bins every convention covers, when there are several),
`<input>_curves.csv` (every bin), `<input>_R_<convention>_<step>.png` (R(Q) and deviation from the
reference) and `<input>_sweep.png` (median deviation against λ step).
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import q_mapping_methods as m  # noqa: E402
import q_methods_synthetic as s  # noqa: E402

METHOD_KEYS = ("E", "A", "B", "C", "C0")

# Categorical slots 1-4 of the reference palette, validated for light mode; the two references are
# neutral ink. Each series also has its own marker and line style, so color never carries it alone.
STYLE = {
    "truth": dict(color="#2b2b2b", linestyle="-", marker=None, label="truth (bin average)"),
    "E": dict(color="#8a8a85", linestyle="--", marker=None, label="E events"),
    "A": dict(color="#2a78d6", linestyle="-", marker="o", label="A point"),
    "B": dict(color="#eb6834", linestyle="-", marker="s", label="B prototype"),
    "C": dict(color="#1baf7a", linestyle="-", marker="^", label="C overlap"),
    "C0": dict(color="#1baf7a", linestyle=":", marker="v", label="C0 overlap, row angle"),
    "legacy": dict(color="#eda100", linestyle="-.", marker="D", label="legacy"),
}


@dataclass
class Case:
    name: str
    reflected: m.Events
    direct: m.Events
    reflected_norm: float
    direct_norm: float
    rows: np.ndarray
    theta_central: float
    y_ref: float
    sign: float
    pixel_width: float
    distance: float
    lam_window: tuple[float, float]
    q_edges: np.ndarray
    truth: np.ndarray | None = None
    legacy: m.QCurve | None = None
    gravity: Callable[[np.ndarray], np.ndarray] | None = None
    slits: m.SlitGeometry | None = None  # in the unit of `pixel_width`; needed by `mean_theta`

    @property
    def reference_key(self) -> str:
        return "truth" if self.truth is not None else "E"


# ---------------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------------


def synthetic_case(setup: s.SyntheticSetup, dq: float, rng: np.random.Generator | None, incident: float) -> Case:
    if rng is None:
        reflected, direct = s.noise_free_events(setup)
        norms = (1.0, 1.0)
    else:
        reflected, direct, n_rb, n_db = s.poisson_events(setup, incident, incident, rng)
        norms = (n_rb, n_db)
    q_lo = 4 * np.pi * np.sin(setup.theta_central) / setup.lam_max * 0.9
    q_hi = 4 * np.pi * np.sin(setup.theta_central) / setup.lam_min * 1.1
    q_edges = m.log_q_edges(q_lo, q_hi, dq)
    rows = np.arange(reflected.row.min(), reflected.row.max() + 1)
    return Case(
        name="synthetic" if setup.geometry == "point" else f"synthetic-{setup.geometry}",
        reflected=reflected,
        direct=direct,
        reflected_norm=norms[0],
        direct_norm=norms[1],
        rows=rows,
        theta_central=setup.theta_central,
        y_ref=setup.y_ref,
        sign=1.0,
        pixel_width=setup.pixel_width_mm,
        distance=setup.distance_mm,
        lam_window=(setup.lam_min, setup.lam_max),
        q_edges=q_edges,
        truth=s.bin_averaged_truth(setup, q_edges),
        slits=setup.slit_geometry(),
    )


def real_case(
    run: int,
    template_path: Path,
    nexus_dir: Path,
    gravity: bool,
    det_sigma_mm: float = 0.8,
    det_res_fn: str = "rectangular",
) -> Case:
    """A run from a legacy template, with background and dead time off on both sides.

    The `mean_theta` slit geometry takes the slit gaps and `xi` from the run's logs, and `xi-reference`
    (445 mm) and `s1-sample-distance` (1485 mm) from the prototype's settings. `det_sigma_mm` and
    `det_res_fn` are the prototype's `DetSigma` and `DetResFn` defaults.
    """
    from mantid.simpleapi import LoadEventNexus

    from lr_reduction.legacy import event_reduction, peak_finding
    from lr_reduction.legacy.gravity_correction import GravityDirection, gravity_correction
    from lr_reduction.legacy.template import read_template
    from lr_reduction.legacy.theta_selection import is_earth_centered_geometry
    from lr_reduction.legacy.utils import amend_config

    nexus_dir = Path(nexus_dir).resolve()
    with amend_config(data_dir=str(nexus_dir)):
        ws_sc = LoadEventNexus(str(nexus_dir / f"REF_L_{run}.nxs.h5"), OutputWorkspace=f"q_methods_rb_{run}")
        sequence_number = 1
        if ws_sc.getRun().hasProperty("sequence_number"):
            sequence_number = ws_sc.getRun().getProperty("sequence_number").value[0]
        template = read_template(str(template_path), sequence_number)
        db_file = nexus_dir / f"REF_L_{template.norm_file}.nxs.h5"
        ws_db = LoadEventNexus(str(db_file), OutputWorkspace=f"q_methods_db_{template.norm_file}")
    thickness = event_reduction.get_attenuation_info(ws_db)
    if thickness > 0:
        ws_db = event_reduction.process_attenuation(ws_db, thickness=thickness)

    run_log = ws_sc.getRun()
    thi, ths = run_log["thi"].value[0], run_log["ths"].value[0]
    theta = np.radians((thi if is_earth_centered_geometry(ws_sc) else ths) + template.angle_offset)

    peak = template.data_peak_range
    _, x, y = peak_finding.process_data(ws_sc, summed=True, tof_step=200)
    start = int(np.argmax(y[peak[0] : peak[1]]) + peak[0])
    y_ref, _, _ = peak_finding.fit_signal_flat_bck(x, y, x_min=peak[0], x_max=peak[1], center=start, sigma=3.0)

    direction = None if gravity else GravityDirection.OFF
    refl = event_reduction.EventReflectivity(
        ws_sc,
        ws_db,
        signal_peak=peak,
        signal_bck=None,
        norm_peak=template.norm_peak_range,
        norm_bck=None,
        specular_pixel=y_ref,
        signal_low_res=template.data_x_range if template.data_x_range_flag else None,
        norm_low_res=template.norm_x_range if template.norm_x_range_flag else None,
        q_min=template.q_min,
        q_step=-template.q_step,
        q_max=None,
        tof_range=list(template.tof_range),
        theta=np.abs(theta),
        instrument=event_reduction.EventReflectivity.INSTRUMENT_4B,
        dead_time=False,
        use_emission_time=template.use_emission_time,
        gravity_direction=direction,
    )
    q_bins, r_legacy, dr_legacy = refl.specular_weighted(q_summing=True)
    legacy = m.QCurve("legacy", np.asarray(q_bins), np.where(r_legacy > 0, r_legacy, np.nan), dr_legacy)

    # The reflected events exactly as legacy `_reflectivity` collects them (4B pixel layout).
    rows, lams, weights = [], [], []
    ws = refl._ws_sc
    for i in range(refl.signal_low_res[0], int(refl.signal_low_res[1]) + 1):
        for j in range(peak[0], int(peak[1]) + 1):
            events = ws.getSpectrum(i * refl.n_y + j)
            if events.getNumberEvents() == 0:
                continue
            tofs = events.getTofs()
            if refl.use_emission_time:
                tofs = refl.emission_time_correction(ws, tofs=tofs)
            lams.append(tofs / refl.constant)
            weights.append(events.getWeights())
            rows.append(np.full(len(tofs), j))
    reflected = m.Events(np.concatenate(rows), np.concatenate(lams), np.concatenate(weights))
    db_lam, db_w = refl._get_events(refl._ws_db, refl.norm_peak, refl.norm_low_res)
    direct = m.Events(np.zeros(len(db_lam), dtype=int), db_lam, db_w)

    gravity_fn = None
    if gravity:

        def gravity_fn(lam: np.ndarray) -> np.ndarray:
            return gravity_correction(refl._ws_sc, np.asarray(lam, dtype=float), gravity_direction=direction)

    sign = float(np.sign(ws.getRun()["ths"].value[-1])) or 1.0
    to_m = 1e-3
    si_sample = 445.0 - float(run_log["BL4B:Mot:xi.RBV"].value[0])
    slits = m.SlitGeometry(
        s1=float(run_log["BL4B:Mot:s1:Y:Gap:Readback"].value[0]) * to_m,
        si=float(run_log["BL4B:Mot:si:Y:Gap:Readback"].value[0]) * to_m,
        interslit=(1485.0 - si_sample) * to_m,
        si_sample=si_sample * to_m,
        sample_detector=float(refl.sample_detector_distance),
        det_sigma=det_sigma_mm * to_m,
        det_res_fn=det_res_fn,
    )
    return Case(
        name=f"REF_L_{run}",
        reflected=reflected,
        direct=direct,
        reflected_norm=refl._ws_sc.getRun().getProtonCharge(),
        direct_norm=refl._ws_db.getRun().getProtonCharge(),
        rows=np.arange(peak[0], int(peak[1]) + 1),
        theta_central=float(refl.theta),
        y_ref=float(y_ref),
        sign=sign,
        pixel_width=refl.pixel_width,
        distance=refl.sample_detector_distance,
        lam_window=(float(refl.wl_range[0]), float(refl.wl_range[1])),
        q_edges=np.asarray(q_bins),
        legacy=legacy,
        gravity=gravity_fn,
        slits=slits,
    )


# ---------------------------------------------------------------------------------------------------
# Running and summarizing
# ---------------------------------------------------------------------------------------------------


CONVENTIONS = ("half", "full", "mean_theta")


def case_angles(case: Case, convention: str) -> m.RowAngles:
    """The rows' angles under one of `CONVENTIONS`."""
    if convention == "mean_theta":
        if case.slits is None:
            raise ValueError(f"{case.name} has no slit geometry for the mean_theta angle model")
        return m.mean_theta_row_angles(
            case.rows, case.theta_central, case.y_ref, case.pixel_width, case.slits, case.sign
        )
    return m.row_angles(
        case.rows, case.theta_central, case.y_ref, case.pixel_width, case.distance, convention, case.sign
    )


def run_methods(
    case: Case, convention: str, lam_step: float, trim: float = 0.0
) -> tuple[dict[str, m.QCurve], m.RowAlignment, np.ndarray]:
    """Run every method on one λ grid. `trim` is method A's `ConstQTrim` (Mantid's default is 0.5)."""
    angles = case_angles(case, convention)
    lam_edges = np.arange(case.lam_window[0], case.lam_window[1] + 0.5 * lam_step, lam_step)
    lam_edges[-1] = min(lam_edges[-1], case.lam_window[1])
    cells = m.histogram_cells(case.reflected, case.direct, case.rows, lam_edges, case.reflected_norm, case.direct_norm)
    gravity = case.gravity(cells.lam_centers) if case.gravity is not None else None
    curves = {
        "E": m.map_events(case.reflected, cells, angles, case.q_edges, case.reflected_norm, gravity=case.gravity),
        "A": m.map_point(cells, angles, case.q_edges, gravity=gravity, trim=trim),
        "B": m.map_prototype(cells, angles, case.q_edges, gravity=gravity, lam_window=case.lam_window),
        "C": m.map_overlap(cells, angles, case.q_edges, gravity=gravity, include_angular_width=True),
        "C0": m.map_overlap(cells, angles, case.q_edges, gravity=gravity, include_angular_width=False),
    }
    q_centers = np.sqrt(case.q_edges[1:] * case.q_edges[:-1])
    covered = m.qline_coverage(q_centers, angles, *case.lam_window) >= 1.0
    # Rows meet a given Q at different wavelengths, so direct-beam noise looks like row disagreement:
    # the alignment check carries it in the cell errors.
    cells_db_error = m.histogram_cells(
        case.reflected,
        case.direct,
        case.rows,
        lam_edges,
        case.reflected_norm,
        case.direct_norm,
        include_direct_beam_error=True,
    )
    return curves, m.row_alignment(cells_db_error, angles, case.q_edges), covered


def reference_values(case: Case, curves: dict[str, m.QCurve]) -> tuple[np.ndarray, np.ndarray | None]:
    if case.truth is not None:
        return case.truth, None
    return curves["E"].r, curves["E"].dr


def common_bins(case: Case, curves: dict[str, m.QCurve], covered: np.ndarray) -> np.ndarray:
    """Fully covered bins where every method, and legacy if present, has a value."""
    mask = covered.copy()
    for curve in curves.values():
        mask &= np.isfinite(curve.r)
    if case.legacy is not None:
        mask &= np.isfinite(case.legacy.r)
    return mask


def summarize(
    case: Case,
    convention: str,
    lam_step: float,
    curves: dict[str, m.QCurve],
    covered: np.ndarray,
    mask: np.ndarray | None = None,
) -> list[dict]:
    """Deviation of every method from the reference, on `mask` (by default `common_bins`)."""
    ref, ref_err = reference_values(case, curves)
    mask = common_bins(case, curves, covered) if mask is None else mask
    rows = []
    for key in METHOD_KEYS:
        if key == case.reference_key:
            continue
        curve = curves[key]
        rel = curve.r[mask] / ref[mask] - 1.0
        row = {
            "input": case.name,
            "convention": convention,
            "lam_step": lam_step,
            "method": key,
            "bins_finite": int(np.isfinite(curve.r).sum()),
            "bins_compared": int(mask.sum()),
            "median_abs_dev_pct": 100 * float(np.median(np.abs(rel))) if mask.any() else np.nan,
            "max_abs_dev_pct": 100 * float(np.max(np.abs(rel))) if mask.any() else np.nan,
        }
        if ref_err is not None:
            row["median_abs_dev_over_sigma"] = float(np.median(np.abs(curve.r[mask] - ref[mask]) / ref_err[mask]))
            row["median_dr_ratio_to_E"] = float(np.median(curve.dr[mask] / ref_err[mask]))
        if case.legacy is not None:
            row["median_abs_dev_vs_legacy_pct"] = 100 * float(
                np.median(np.abs(curve.r[mask] / case.legacy.r[mask] - 1))
            )
        rows.append(row)
    if case.legacy is not None:
        rel = case.legacy.r[mask] / ref[mask] - 1.0
        rows.append(
            {
                "input": case.name,
                "convention": convention,
                "lam_step": lam_step,
                "method": "legacy",
                "bins_finite": int(np.isfinite(case.legacy.r).sum()),
                "bins_compared": int(mask.sum()),
                "median_abs_dev_pct": 100 * float(np.median(np.abs(rel))),
                "max_abs_dev_pct": 100 * float(np.max(np.abs(rel))),
                "median_abs_dev_over_sigma": float(np.median(np.abs(case.legacy.r[mask] - ref[mask]) / ref_err[mask])),
                "median_dr_ratio_to_E": float(np.median(case.legacy.dr[mask] / ref_err[mask])),
            }
        )
    return rows


def poisson_summary(
    setup: s.SyntheticSetup, args: argparse.Namespace, conventions: tuple[str, ...], rng: np.random.Generator
) -> list[dict]:
    """Bias and pull width over `args.repeats` Poisson data sets, on bins valid in every repeat.

    Each data set is drawn once and reduced with every convention and λ step.
    """
    configs = [(c, step) for c in conventions for step in args.lam_step]
    values = {cfg: {k: [] for k in METHOD_KEYS} for cfg in configs}
    errors = {cfg: {k: [] for k in METHOD_KEYS} for cfg in configs}
    masks: dict[tuple, list[np.ndarray]] = {cfg: [] for cfg in configs}
    truth = None
    for repeat in range(args.repeats):
        case = synthetic_case(setup, args.dq, rng, args.incident)
        truth = case.truth
        for cfg in configs:
            curves, _, covered = run_methods(case, *cfg, args.trim)
            masks[cfg].append(common_bins(case, curves, covered))
            for key in METHOD_KEYS:
                values[cfg][key].append(curves[key].r)
                errors[cfg][key].append(curves[key].dr)
        print(f"data set {repeat + 1}/{args.repeats}: {len(case.reflected.row)} reflected events", flush=True)

    rows = []
    for cfg in configs:
        positive = [np.all(np.array(errors[cfg][k]) > 0, axis=0) for k in METHOD_KEYS]
        mask = np.all(masks[cfg], axis=0) & np.all(positive, axis=0)
        for key in METHOD_KEYS:
            r = np.array(values[cfg][key])[:, mask]
            dr = np.array(errors[cfg][key])[:, mask]
            rows.append(
                {
                    "input": f"{case.name}-poisson",
                    "convention": cfg[0],
                    "lam_step": cfg[1],
                    "method": key,
                    "repeats": args.repeats,
                    "bins_compared": int(mask.sum()),
                    "median_bias_pct": 100 * float(np.median(r.mean(axis=0) / truth[mask] - 1)),
                    "pull_std": float(np.std((r - truth[mask]) / dr)),
                    "median_scatter_over_mean_dr": float(np.median(r.std(axis=0) / dr.mean(axis=0))),
                }
            )
    return rows


# ---------------------------------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------------------------------


def print_table(rows: list[dict]) -> None:
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    widths = {k: max(len(k), *(len(_fmt(row.get(k, ""))) for row in rows)) for k in keys}
    print("  ".join(k.rjust(widths[k]) for k in keys))
    for row in rows:
        print("  ".join(_fmt(row.get(k, "")).rjust(widths[k]) for k in keys))


def _fmt(value) -> str:
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def write_csv(path: Path, rows: list[dict]) -> None:
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def curve_rows(
    case: Case, convention: str, lam_step: float, curves: dict[str, m.QCurve], covered: np.ndarray
) -> list[dict]:
    q = np.sqrt(case.q_edges[1:] * case.q_edges[:-1])
    out = []
    for k in range(len(q)):
        row = {
            "input": case.name,
            "convention": convention,
            "lam_step": lam_step,
            "q": q[k],
            "covered": bool(covered[k]),
        }
        if case.truth is not None:
            row["truth"] = case.truth[k]
        if case.legacy is not None:
            row["legacy"], row["legacy_dr"] = case.legacy.r[k], case.legacy.dr[k]
        for key in METHOD_KEYS:
            row[key], row[f"{key}_dr"] = curves[key].r[k], curves[key].dr[k]
        out.append(row)
    return out


def plot_curves(
    path: Path, case: Case, convention: str, lam_step: float, curves: dict[str, m.QCurve], covered: np.ndarray
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    q = np.sqrt(case.q_edges[1:] * case.q_edges[:-1])
    ref, _ = reference_values(case, curves)
    series = {"truth": case.truth} if case.truth is not None else {}
    series.update({key: curves[key].r for key in METHOD_KEYS})
    if case.legacy is not None:
        series["legacy"] = case.legacy.r

    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(8, 7), height_ratios=[3, 2])
    for key, r in series.items():
        style = STYLE[key]
        top.plot(q, r, linewidth=1.5, markersize=4, markevery=4, **style)
        if key != case.reference_key:
            bottom.plot(
                q[covered], 100 * (r[covered] / ref[covered] - 1), linewidth=1.5, markersize=4, markevery=4, **style
            )
    top.set_yscale("log")
    top.set_ylabel("R(Q)")
    top.set_title(f"{case.name}: {convention} angle convention, Δλ = {lam_step} Å")
    top.legend(frameon=False, fontsize=9)
    bottom.axhline(0.0, color="#8a8a85", linewidth=0.8)
    limit = 10.0
    bottom.set_ylim(-limit, limit)
    bottom.set_ylabel(
        f"deviation from {STYLE[case.reference_key]['label']} (%)\nfully covered bins, clipped at ±{limit:g}%"
    )
    bottom.set_xlabel("Q (1/Å)")
    for axis in (top, bottom):
        axis.grid(True, color="#e6e5e0", linewidth=0.6)
        axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_sweep(path: Path, name: str, rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    conventions = list(dict.fromkeys(row["convention"] for row in rows))
    fig, axes = plt.subplots(1, len(conventions), sharey=True, figsize=(5 * len(conventions), 4), squeeze=False)
    for axis, convention in zip(axes[0], conventions):
        for key in [k for k in (*METHOD_KEYS, "legacy") if any(r["method"] == k for r in rows)]:
            points = sorted(
                (r["lam_step"], r["median_abs_dev_pct"])
                for r in rows
                if r["method"] == key and r["convention"] == convention
            )
            steps, devs = zip(*points)
            axis.plot(steps, devs, linewidth=1.5, markersize=6, **STYLE[key])
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(f"{convention} angle convention")
        axis.set_xlabel("λ bin width (Å)")
        axis.grid(True, color="#e6e5e0", linewidth=0.6)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0][0].set_ylabel("median |deviation| from reference (%)")
    axes[0][-1].legend(frameon=False, fontsize=9)
    fig.suptitle(name)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest="input", required=True)
    for name in ("synthetic", "real"):
        p = sub.add_parser(name)
        p.add_argument("--convention", choices=(*CONVENTIONS, "both", "all"), default="both")
        p.add_argument("--lam-step", type=float, nargs="+", default=[0.005, 0.02, 0.1], help="λ bin widths (Å)")
        p.add_argument("--output", type=Path, help="directory for CSV and PNG output")
        p.add_argument(
            "--trim",
            type=float,
            default=0.0,
            help="method A's ConstQTrim (Mantid default 0.5; with a uniform λ grid and log Q bins it drops "
            "every bin above about twice the lowest covered Q, so it is off by default)",
        )
        if name == "synthetic":
            p.add_argument("--dq", type=float, default=0.01, help="ΔQ/Q of the log Q grid")
            p.add_argument("--theta", type=float, default=0.45, help="central angle (degrees)")
            p.add_argument("--blur", type=float, default=0.0, help="detector blur sigma (pixels)")
            p.add_argument("--geometry", choices=s.GEOMETRIES, default="point", help="synthetic beam geometry")
            p.add_argument("--s1", type=float, default=s.SyntheticSetup.s1_mm, help="S1 gap (mm)")
            p.add_argument("--si", type=float, default=s.SyntheticSetup.si_mm, help="Si gap (mm)")
            p.add_argument("--repeats", type=int, default=0, help="Poisson data sets (0: noise-free)")
            p.add_argument("--incident", type=float, default=5e7, help="incident neutrons per Poisson data set")
            p.add_argument("--seed", type=int, default=0)
        else:
            p.add_argument("--run", type=int, default=198415)
            p.add_argument("--template", type=Path, default=Path("tests/data/template.xml"))
            p.add_argument("--nexus-dir", type=Path, default=Path("tests/data/liquidsreflectometer-data/nexus"))
            p.add_argument("--gravity", action="store_true", help="apply legacy's gravity correction everywhere")
            p.add_argument("--det-sigma", type=float, default=0.8, help="mean_theta detector resolution σ (mm)")
            p.add_argument(
                "--det-res-fn",
                choices=("rectangular", "gaussian", "none"),
                default="rectangular",
                help="mean_theta detector response shape",
            )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    conventions = {"both": ("half", "full"), "all": CONVENTIONS}.get(args.convention, (args.convention,))
    if args.output:
        args.output.mkdir(parents=True, exist_ok=True)

    if args.input == "synthetic":
        setup = s.SyntheticSetup(
            theta_central_deg=args.theta, blur_pixels=args.blur, geometry=args.geometry, s1_mm=args.s1, si_mm=args.si
        )
        if args.repeats > 0:
            rng = np.random.default_rng(args.seed)
            rows = poisson_summary(setup, args, conventions, rng)
            print_table(rows)
            if args.output:
                write_csv(args.output / f"{rows[0]['input']}_summary.csv", rows)
            return
        case = synthetic_case(setup, args.dq, None, 0.0)
    else:
        case = real_case(args.run, args.template, args.nexus_dir, args.gravity, args.det_sigma, args.det_res_fn)

    summary, per_bin, alignment = [], [], []
    results = {}
    for convention in conventions:
        for step in args.lam_step:
            curves, aligned, covered = run_methods(case, convention, step, args.trim)
            results[convention, step] = (curves, covered)
            summary += summarize(case, convention, step, curves, covered)
            per_bin += curve_rows(case, convention, step, curves, covered)
            alignment.append(
                {
                    "convention": convention,
                    "lam_step": step,
                    "row_chi2_per_dof": aligned.chi2_per_dof,
                    "row_rms_pct": 100 * aligned.rms_relative,
                    "rows": f"{aligned.rows.min()}-{aligned.rows.max()}",
                    "bins": aligned.n_bins,
                }
            )
            if args.output:
                plot_curves(
                    args.output / f"{case.name}_R_{convention}_{step:g}.png", case, convention, step, curves, covered
                )
    print(f"\n{case.name}: deviation from {STYLE[case.reference_key]['label']} on fully covered bins")
    print_table(summary)
    print("\nPer-row alignment (lower is better; χ²/dof ≈ 1 when the rows agree within statistics)")
    print_table(alignment)
    shared = []
    if len(conventions) > 1:
        for step in args.lam_step:
            mask = np.logical_and.reduce([common_bins(case, *results[c, step]) for c in conventions])
            for convention in conventions:
                shared += summarize(case, convention, step, *results[convention, step], mask=mask)
        print(
            f"\n{case.name}: deviation from {STYLE[case.reference_key]['label']} on the bins every angle model covers"
        )
        print_table(shared)
    if args.output:
        write_csv(args.output / f"{case.name}_summary.csv", summary + alignment)
        if shared:
            write_csv(args.output / f"{case.name}_shared_bins_summary.csv", shared)
        write_csv(args.output / f"{case.name}_curves.csv", per_bin)
        plot_sweep(args.output / f"{case.name}_sweep.png", case.name, summary)
        print(f"\nwrote {args.output}/")


if __name__ == "__main__":
    main()
