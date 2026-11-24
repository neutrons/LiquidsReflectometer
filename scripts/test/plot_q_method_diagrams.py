"""
Explanatory figures for METHOD_Q_transform_comparison.md.

    python scripts/test/plot_q_method_diagrams.py [--output figures]

Writes:

* `q_methods_plane.png` — detector cells in the (λ, θ) plane, with constant-Q bin edges;
* `q_methods_plane_tof.png` — the region one Q bin collects under constant Q and under constant TOF;
* `q_methods_tof_film.png` — constant TOF next to C0 on the synthetic film;
* `q_methods_deposit.png` — how each method shares one cell among the Q bins;
* `q_methods_flat.png` — every method applied to a synthetic sample with R(Q) ≡ 1;
* `q_methods_row_weights.png` — the weight method A gives each row inside one Q bin;
* `q_methods_poisson.png` — bias and pull width of each method from the Poisson study CSV, when given;
* `q_methods_angle_models.png` — the row angles under `half`, `full` and the prototype's `meantheta`.

The deposits in `q_methods_deposit.png` follow the formulas in `q_mapping_methods.py`; the flat-R and
row-weight figures call those functions directly.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import q_mapping_methods as m
import q_methods_synthetic as s
from matplotlib.patches import Rectangle

FOUR_PI = 4.0 * np.pi
INK = "#2b2b2b"
MUTED = "#8a8a85"
GRID = "#e4e4e0"
COLORS = {"E": "#8a8a85", "A": "#2a78d6", "B": "#eb6834", "C": "#1baf7a", "C0": "#1baf7a"}
LABELS = {
    "E": "E events",
    "A": "A point",
    "B": "B prototype",
    "C": "C overlap, corner footprint",
    "C0": "C0 overlap, row angle",
}
MARKERS = {"E": None, "A": "o", "B": "s", "C": "^", "C0": "v"}
STYLES = {"E": "--", "A": "-", "B": "-", "C": "-", "C0": ":"}

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
    }
)


def q_of(theta: np.ndarray | float, lam: np.ndarray | float) -> np.ndarray:
    return FOUR_PI * np.sin(theta) / lam


# ---------------------------------------------------------------------------------------------------
# Figure 1: cells in the (λ, θ) plane
# ---------------------------------------------------------------------------------------------------


def plot_plane(path: Path) -> None:
    """A few rows × λ bins, the constant-Q bin edges through them, and one highlighted cell."""
    theta0 = np.radians(0.60)
    dtheta = np.radians(0.012)  # one row, exaggerated slightly for legibility
    rows = np.arange(-2, 3)
    lam_edges = np.arange(4.6, 5.41, 0.1)
    q_mid = q_of(theta0, 5.0)
    q_edges = q_mid * 1.01 ** np.arange(-12, 13)
    y_lo, y_hi = np.degrees(theta0 - 3.0 * dtheta), np.degrees(theta0 + 3.0 * dtheta)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for j in rows:
        for lo, hi in zip(lam_edges[:-1], lam_edges[1:]):
            ax.add_patch(
                Rectangle(
                    (lo, np.degrees(theta0 + (j - 0.5) * dtheta)),
                    hi - lo,
                    np.degrees(dtheta),
                    facecolor="none",
                    edgecolor=MUTED,
                    linewidth=0.8,
                )
            )
    lam_line = np.linspace(lam_edges[0], lam_edges[-1], 200)
    for q in q_edges:
        theta_line = np.degrees(np.arcsin(q * lam_line / FOUR_PI))
        ax.plot(lam_line, theta_line, color="#9fb6e3", linewidth=0.9, zorder=1)

    # Highlighted cell: row 0, λ bin 4.9–5.0.
    lo, hi = 4.9, 5.0
    t_lo, t_hi = theta0 - dtheta / 2, theta0 + dtheta / 2
    ax.add_patch(
        Rectangle(
            (lo, np.degrees(t_lo)),
            hi - lo,
            np.degrees(dtheta),
            facecolor="#1baf7a",
            alpha=0.18,
            edgecolor="#1baf7a",
            linewidth=1.6,
            zorder=3,
        )
    )
    (row_line,) = ax.plot([lo, hi], [np.degrees(theta0)] * 2, color="#1baf7a", linewidth=2.4, zorder=4)
    (center,) = ax.plot(0.5 * (lo + hi), np.degrees(theta0), "o", color="#2a78d6", markersize=8, zorder=5)
    (corners,) = ax.plot(
        [hi, lo], [np.degrees(t_lo), np.degrees(t_hi)], "^", color="#1baf7a", markersize=8, linestyle="none", zorder=5
    )
    (q_line,) = ax.plot([], [], color="#9fb6e3", linewidth=0.9)
    ax.legend(
        [q_line, center, corners, row_line],
        [
            "constant-Q bin edges, ΔQ/Q = 1%",
            "cell center: A puts the whole cell at this Q",
            "corners: the lowest and highest Q in the cell (B and C span these)",
            "row angle across the λ bin (C0 and E span this)",
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=2,
        fontsize=8,
        frameon=False,
    )
    ax.set_xlim(lam_edges[0] - 0.02, lam_edges[-1] + 0.02)
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel("wavelength λ (Å)")
    ax.set_ylabel("θ (deg), one band per row")
    ax.set_title("Detector cells (gray boxes): rows × λ bins, cut by constant-Q lines", color=INK, fontsize=10)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------
# Figure 1b: one Q bin under constant Q and under constant TOF
# ---------------------------------------------------------------------------------------------------


def plot_plane_tof(path: Path) -> None:
    """The cells of Figure 1, with the region one Q bin collects under constant Q and under constant TOF."""
    theta0 = np.radians(0.60)
    dtheta = np.radians(0.012)
    rows = np.arange(-2, 3)
    lam_edges = np.arange(4.6, 5.41, 0.1)
    q_lo = q_of(theta0, 5.0)
    q_hi = q_lo * 1.01
    t_min, t_max = theta0 - 2.5 * dtheta, theta0 + 2.5 * dtheta
    lam_line = np.linspace(lam_edges[0], lam_edges[-1], 400)
    theta_line = {q: np.degrees(np.arcsin(q * lam_line / FOUR_PI)) for q in (q_lo, q_hi)}
    band_color, tof_color = "#1baf7a", "#e87ba4"

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.9), sharey=True)
    for ax in axes:
        for j in rows:
            for lo, hi in zip(lam_edges[:-1], lam_edges[1:]):
                ax.add_patch(
                    Rectangle(
                        (lo, np.degrees(theta0 + (j - 0.5) * dtheta)),
                        hi - lo,
                        np.degrees(dtheta),
                        facecolor="none",
                        edgecolor=MUTED,
                        linewidth=0.8,
                        zorder=2,
                    )
                )
        ax.axhline(np.degrees(theta0), color=INK, linewidth=0.8, linestyle=(0, (1, 2)), zorder=2)
        ax.set_xlim(lam_edges[0] - 0.02, lam_edges[-1] + 0.02)
        ax.set_ylim(np.degrees(t_min) - 0.004, np.degrees(t_max) + 0.004)
        ax.set_xlabel("wavelength λ (Å)")
        ax.grid(False)

    # Constant Q: the bin is the band between its two constant-Q edges.
    ax = axes[0]
    ax.fill_between(lam_line, theta_line[q_lo], theta_line[q_hi], color=band_color, alpha=0.30, linewidth=0, zorder=1)
    for q in (q_lo, q_hi):
        ax.plot(lam_line, theta_line[q], color=band_color, linewidth=1.4, zorder=3)
    ax.set_ylabel("θ (deg), one band per row")
    ax.set_title("constant Q (A, B, C, C0, E): the bin follows constant-Q lines", color=INK, fontsize=10)
    ax.text(
        4.62,
        np.degrees(t_max) + 0.001,
        "each row meets the bin at its own λ",
        color=INK,
        fontsize=8,
        va="bottom",
    )

    # Constant TOF: the bin is a λ strip, the same for every row, labeled with Q at θ_central.
    ax = axes[1]
    lam_a = FOUR_PI * np.sin(theta0) / q_hi
    lam_b = FOUR_PI * np.sin(theta0) / q_lo
    ax.axvspan(lam_a, lam_b, color=tof_color, alpha=0.30, linewidth=0, zorder=1)
    for lam in (lam_a, lam_b):
        ax.axvline(lam, color=tof_color, linewidth=1.4, zorder=3)
    for q in (q_lo, q_hi):
        ax.plot(lam_line, theta_line[q], color=band_color, linewidth=1.0, linestyle="--", zorder=3)
    for j in rows:
        theta = theta0 + j * dtheta
        shift = 100 * (np.sin(theta) / np.sin(theta0) - 1)
        ax.text(
            lam_edges[-1] - 0.01,
            np.degrees(theta),
            f"true Q {shift:+.0f}%" if j else "Q as labeled",
            color=INK,
            fontsize=8,
            ha="right",
            va="center",
            zorder=4,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.0},
        )
    ax.text(lam_edges[0] + 0.02, np.degrees(theta0) + 0.0015, "θ_central", color=INK, fontsize=8, va="bottom")
    ax.set_title("constant TOF (T): the bin is a λ strip, Q taken at θ_central", color=INK, fontsize=10)

    handles = [
        plt.Line2D([], [], color=band_color, linewidth=1.4),
        plt.Line2D([], [], color=tof_color, linewidth=1.4),
        plt.Line2D([], [], color=band_color, linewidth=1.0, linestyle="--"),
    ]
    fig.legend(
        handles,
        [
            "constant-Q bin, ΔQ/Q = 1%",
            "constant-TOF bin with the same Q label",
            "that constant-Q bin, for comparison",
        ],
        loc="lower center",
        ncol=3,
        fontsize=8,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------
# Figure 2: one cell's deposit under each method
# ---------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class OneCell:
    theta: float
    dtheta: float
    lam: float
    dlam: float

    @property
    def q_center(self) -> float:
        return float(q_of(self.theta, self.lam))

    @property
    def lam_footprint(self) -> tuple[float, float]:
        return float(q_of(self.theta, self.lam + self.dlam / 2)), float(q_of(self.theta, self.lam - self.dlam / 2))

    @property
    def corner_footprint(self) -> tuple[float, float]:
        lo = q_of(self.theta - self.dtheta / 2, self.lam + self.dlam / 2)
        hi = q_of(self.theta + self.dtheta / 2, self.lam - self.dlam / 2)
        return float(lo), float(hi)

    def true_density(self, q: np.ndarray, n: int = 400) -> np.ndarray:
        """Density of Q for counts uniform in (θ, λ) over the cell, normalized to 1."""
        t = (np.arange(n) + 0.5) / n - 0.5
        theta = self.theta + self.dtheta * t
        lam = self.lam + self.dlam * t
        samples = q_of(theta[:, None], lam[None, :]).ravel()
        hist, edges = np.histogram(samples, bins=len(q) - 1, range=(q[0], q[-1]), density=True)
        return hist, edges


def box_shares(lo: float, hi: float, q_edges: np.ndarray) -> np.ndarray:
    overlap = np.clip(np.minimum(hi, q_edges[1:]) - np.maximum(lo, q_edges[:-1]), 0.0, None)
    return overlap / (hi - lo)


def plot_deposit(path: Path) -> None:
    cell = OneCell(theta=np.radians(0.60), dtheta=np.radians(0.011), lam=5.0, dlam=0.1)
    q_edges = cell.q_center * 1.005 ** np.arange(-7, 8) / 1.0025
    q_centers = np.sqrt(q_edges[1:] * q_edges[:-1])
    widths = np.diff(q_edges)

    # E: the Q distribution of the cell's own events, binned.
    fine = np.linspace(q_edges[0], q_edges[-1], 600)
    density, fine_edges = cell.true_density(fine)
    e_share = np.array(
        [
            density[(fine_edges[:-1] >= a) & (fine_edges[:-1] < b)].sum() * np.diff(fine_edges)[0]
            for a, b in zip(q_edges[:-1], q_edges[1:])
        ]
    )
    e_lam_only = box_shares(*cell.lam_footprint, q_edges)

    # A: all to the bin of the center.
    a_share = np.zeros(len(q_centers))
    a_share[np.digitize(cell.q_center, q_edges) - 1] = 1.0

    # B: equal shares to the bin centers inside the corner interval.
    lo, hi = cell.corner_footprint
    inside = (q_centers >= lo) & (q_centers <= hi)
    b_share = inside / inside.sum()

    c_share = box_shares(*cell.corner_footprint, q_edges)
    c0_share = box_shares(*cell.lam_footprint, q_edges)

    panels = [
        ("E", e_lam_only, "events at the row angle: uniform over the λ footprint"),
        ("A", a_share, "whole cell into the bin holding its center Q"),
        ("B", b_share, "equal share to every bin whose center lies in the corner interval"),
        ("C", c_share, "uniform over the corner footprint, shared by overlap"),
        ("C0", c0_share, "uniform over the λ footprint at the row angle, shared by overlap"),
    ]
    fig, axes = plt.subplots(len(panels), 1, figsize=(8.2, 8.6), sharex=True)
    x = (q_centers - cell.q_center) / cell.q_center * 100
    bar_w = widths / cell.q_center * 100 * 0.86
    span_lam = [(v - cell.q_center) / cell.q_center * 100 for v in cell.lam_footprint]
    span_corner = [(v - cell.q_center) / cell.q_center * 100 for v in cell.corner_footprint]
    for ax, (key, share, text) in zip(axes, panels):
        ax.bar(x, share, width=bar_w, color=COLORS[key], alpha=0.85 if key != "C0" else 0.55, zorder=2)
        # Reference: the true Q density of a cell with angular width, as fraction per bin.
        ax.step(
            np.r_[x - bar_w / 2 / 0.86, x[-1] + bar_w[-1] / 2 / 0.86],
            np.r_[e_share, e_share[-1]],
            where="post",
            color=INK,
            linewidth=1.0,
            zorder=3,
        )
        ax.axvspan(*span_lam, color="#1baf7a", alpha=0.07, zorder=0)
        for v in span_corner:
            ax.axvline(v, color=MUTED, linestyle=":", linewidth=0.9)
        ax.axvline(0.0, color="#2a78d6", linewidth=0.8, alpha=0.6)
        ax.set_ylim(0, 1.05 if key == "A" else max(0.45, share.max() * 1.15))
        ax.set_ylabel("share")
        ax.text(0.01, 0.93, f"{LABELS[key]}: {text}", transform=ax.transAxes, fontsize=8.5, color=INK, va="top")
        ax.grid(axis="x", visible=False)
    axes[-1].set_xlabel("Q relative to the cell center (%)    —    Q bins of ΔQ/Q = 0.5%")
    axes[0].set_title(
        "One cell (θ = 0.60°, Δθ = 0.011°, λ = 5.0 Å, Δλ = 0.1 Å) shared among Q bins\n"
        "black step: true Q distribution of the cell's counts; shaded: λ footprint; dotted: corners",
        fontsize=9.5,
        color=INK,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------
# Figure 3: flat R(Q) ≡ 1
# ---------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class FlatSetup(s.SyntheticSetup):
    def reflectivity(self, q: np.ndarray) -> np.ndarray:
        return np.ones_like(np.asarray(q, dtype=float))


def flat_curves(
    setup: s.SyntheticSetup, lam_step: float, dq: float
) -> tuple[np.ndarray, dict[str, m.QCurve], np.ndarray]:
    reflected, direct = s.noise_free_events(setup)
    rows = np.arange(reflected.row.min(), reflected.row.max() + 1)
    angles = m.row_angles(rows, setup.theta_central, setup.y_ref, setup.pixel_width_mm, setup.distance_mm, "half")
    lam_edges = np.arange(setup.lam_min, setup.lam_max + 0.5 * lam_step, lam_step)
    lam_edges[-1] = min(lam_edges[-1], setup.lam_max)
    cells = m.histogram_cells(reflected, direct, rows, lam_edges)
    q_lo = q_of(setup.theta_central, setup.lam_max) * 0.9
    q_hi = q_of(setup.theta_central, setup.lam_min) * 1.1
    q_edges = m.log_q_edges(q_lo, q_hi, dq)
    curves = {
        "E": m.map_events(reflected, cells, angles, q_edges),
        "A": m.map_point(cells, angles, q_edges, trim=0.0),
        "B": m.map_prototype(cells, angles, q_edges, lam_window=(setup.lam_min, setup.lam_max)),
        "C": m.map_overlap(cells, angles, q_edges, include_angular_width=True),
        "C0": m.map_overlap(cells, angles, q_edges, include_angular_width=False),
    }
    q_centers = np.sqrt(q_edges[1:] * q_edges[:-1])
    covered = m.qline_coverage(q_centers, angles, setup.lam_min, setup.lam_max) >= 1.0
    return q_centers, curves, covered


def plot_tof_film(path: Path, lam_step: float = 0.02, dq: float = 0.01) -> None:
    """Constant TOF next to C0 on the noise-free film, against the truth and the divergence-averaged truth."""
    setup = s.SyntheticSetup()
    reflected, direct = s.noise_free_events(setup)
    rows = np.arange(reflected.row.min(), reflected.row.max() + 1)
    angles = m.row_angles(rows, setup.theta_central, setup.y_ref, setup.pixel_width_mm, setup.distance_mm, "half")
    lam_edges = np.arange(setup.lam_min, setup.lam_max + 0.5 * lam_step, lam_step)
    lam_edges[-1] = min(lam_edges[-1], setup.lam_max)
    cells = m.histogram_cells(reflected, direct, rows, lam_edges)
    q_edges = m.log_q_edges(
        q_of(setup.theta_central, setup.lam_max) * 0.9, q_of(setup.theta_central, setup.lam_min) * 1.1, dq
    )
    q = np.sqrt(q_edges[1:] * q_edges[:-1])
    covered = m.qline_coverage(q, angles, setup.lam_min, setup.lam_max) >= 1.0
    truth = s.bin_averaged_truth(setup, q_edges)
    smeared = s.divergence_averaged_truth(setup, q_edges)
    c0 = m.map_overlap(cells, angles, q_edges, include_angular_width=False)
    tof = m.map_constant_tof(cells, setup.theta_central, q_edges)
    show = covered & np.isfinite(c0.r) & np.isfinite(tof.r)
    tof_color = "#e87ba4"

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8.2, 6.4), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    top.plot(q[show], truth[show], color=INK, linewidth=1.4, label="truth R(Q)")
    top.plot(
        q[show],
        smeared[show],
        color=MUTED,
        linewidth=1.4,
        linestyle="--",
        label="truth averaged over the beam divergence",
    )
    top.plot(q[show], c0.r[show], color=COLORS["C0"], linestyle="none", marker="v", markersize=4, label=LABELS["C0"])
    top.plot(q[show], tof.r[show], color=tof_color, linestyle="none", marker="D", markersize=4, label="T constant TOF")
    top.set_yscale("log")
    top.set_ylabel("R")
    top.legend(fontsize=8, frameon=False, loc="lower left")
    top.set_title(
        f"Synthetic film, θ = {setup.theta_central_deg}°, divergence FWHM ≈ 18% of θ, "
        f"Δλ = {lam_step} Å, ΔQ/Q = {dq:.0%}",
        color=INK,
        fontsize=10,
    )
    bottom.axhline(0.0, color=INK, linewidth=0.8)
    bottom.plot(
        q[show], 100 * (c0.r[show] / truth[show] - 1), color=COLORS["C0"], marker="v", markersize=3, linewidth=1
    )
    bottom.plot(q[show], 100 * (tof.r[show] / truth[show] - 1), color=tof_color, marker="D", markersize=3, linewidth=1)
    bottom.plot(q[show], 100 * (smeared[show] / truth[show] - 1), color=MUTED, linestyle="--", linewidth=1.2)
    bottom.set_ylabel("deviation from truth (%)")
    bottom.set_xlabel("Q (1/Å)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_flat(path: Path) -> None:
    """Small multiples: one column per method, one row per λ grid; the truth is 0 everywhere."""
    setup = FlatSetup()
    steps = (0.02, 0.1)
    keys = ("A", "B", "C", "C0", "E")
    fig, axes = plt.subplots(len(steps), len(keys), figsize=(11.0, 4.8), sharex=True, sharey=True)
    for i, step in enumerate(steps):
        q, curves, covered = flat_curves(setup, step, dq=0.005)
        for j, key in enumerate(keys):
            ax = axes[i, j]
            r = np.where(covered, curves[key].r, np.nan)
            ax.axhline(0.0, color=MUTED, linewidth=0.8)
            ax.plot(q, (r - 1.0) * 100, color=COLORS[key], linewidth=1.0)
            ax.set_ylim(-8, 8)
            if i == 0:
                ax.set_title(LABELS[key].replace(", ", "\n"), fontsize=9, color=INK)
            if j == 0:
                ax.set_ylabel(f"Δλ = {step} Å\nR − 1 (%)")
            if i == len(steps) - 1:
                ax.set_xlabel("Q (1/Å)")
            ax.tick_params(labelsize=8)
    fig.suptitle(
        "Sample with R(Q) ≡ 1: every deviation is the method's own (noise-free, half convention, ΔQ/Q = 0.5%)",
        color=INK,
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------
# Figure 4: A's per-row weights
# ---------------------------------------------------------------------------------------------------


def plot_row_weights(path: Path) -> None:
    setup = s.SyntheticSetup()
    reflected, direct = s.noise_free_events(setup)
    rows = np.arange(reflected.row.min(), reflected.row.max() + 1)
    angles = m.row_angles(rows, setup.theta_central, setup.y_ref, setup.pixel_width_mm, setup.distance_mm, "half")
    lam_edges = np.arange(setup.lam_min, setup.lam_max + 0.0025, 0.005)
    cells = m.histogram_cells(reflected, direct, rows, lam_edges)
    q_edges = m.log_q_edges(0.012, 0.03, 0.01)
    q = q_of(angles.center[:, None], cells.lam_centers[None, :])
    k = np.digitize(q, q_edges) - 1
    target = len(q_edges) // 2
    in_bin = k == target
    n_by_row = in_bin.sum(axis=1)
    weight = n_by_row * len(rows) / n_by_row.sum()
    # Each row's share of the signal at this Q: the mean of its cells in the bin.
    row_mean = np.array([cells.r[i][in_bin[i]].mean() for i in range(len(rows))])
    share = row_mean / row_mean.sum()
    bias = np.sum(weight * share) - 1.0
    centroid = np.sum(share * rows)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.6, 5.4), sharex=True)
    ax1.bar(rows, share * 100, color=MUTED, width=0.8)
    ax1.axvline(centroid, color=INK, linewidth=1.0, linestyle="--")
    ax1.axvline(rows.mean(), color=MUTED, linewidth=1.0, linestyle=":")
    ax1.text(centroid + 0.15, share.max() * 100 * 0.95, f"beam centroid {centroid:.2f}", fontsize=8, color=INK)
    ax1.text(
        rows.mean() - 0.15, share.max() * 100 * 0.95, f"row middle {rows.mean():.1f}", fontsize=8, color=INK, ha="right"
    )
    ax1.set_ylabel("share of the signal\nat this Q (%)")
    ax1.set_title(
        f"Q bin {q_edges[target]:.4f}–{q_edges[target + 1]:.4f} 1/Å, Δλ = 0.005 Å", fontsize=9.5, color=INK, loc="left"
    )
    ax2.bar(rows, weight, color="#2a78d6", width=0.8)
    ax2.axhline(1.0, color=INK, linewidth=1.0)
    ax2.text(
        rows[0] - 0.4,
        1.03,
        f"weight every row should get; A's weighted sum is {100 * bias:+.2f}%",
        fontsize=8,
        color=INK,
    )
    ax2.set_ylim(0.6, 1.35)
    ax2.set_ylabel("weight in A\n(N_rows · n_y / Σ n_y)")
    ax2.set_xlabel("detector row y (higher rows: larger θ, more cells per Q bin)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------
# Figure 5: Poisson study
# ---------------------------------------------------------------------------------------------------


def plot_poisson(path: Path, summary_csv: Path) -> None:
    with summary_csv.open() as f:
        rows = list(csv.DictReader(f))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    for key in ("E", "A", "B", "C", "C0"):
        sel = sorted((r for r in rows if r["method"] == key), key=lambda r: float(r["lam_step"]))
        x = [float(r["lam_step"]) for r in sel]
        style = dict(color=COLORS[key], linestyle=STYLES[key], marker=MARKERS[key] or "x", label=LABELS[key])
        ax1.plot(x, [float(r["median_bias_pct"]) for r in sel], **style)
        ax2.plot(x, [float(r["pull_std"]) for r in sel], **style)
    for ax in (ax1, ax2):
        ax.set_xscale("log")
        ax.set_xlabel("λ bin width (Å)")
    ax1.axhline(0.0, color=MUTED, linewidth=0.8)
    ax1.set_ylabel("median bias from truth (%)")
    ax2.set_yscale("log")
    ax2.axhline(1.0, color=MUTED, linewidth=0.8)
    ax2.set_ylabel("pull width, std of (R − truth)/dR")
    ax2.legend(fontsize=8, frameon=False)
    fig.suptitle("Poisson study: 20 data sets, ~2.7 M reflected events each", color=INK)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------
# Figure 7: the three angle models
# ---------------------------------------------------------------------------------------------------


def plot_angle_models(path: Path) -> None:
    """θ(y) − θ_central under `half`, `full` and the prototype's `meantheta`, at two slit settings."""
    distance, pixel = 1830.0, 0.7
    settings = [
        ("REF_L_198415, 2.34°: S1 1.52 mm, Si 0.98 mm", 1.5205, 0.976),
        ("REF_L_198409, 0.60°: S1 0.39 mm, Si 0.25 mm", 0.394, 0.250),
    ]
    dy = np.linspace(-8.0, 8.0, 321)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9), sharey=True)
    for ax, (title, s1, si) in zip(axes, settings):
        geometry = m.SlitGeometry(s1, si, 1350.0, 135.0, distance, det_sigma=0.8, det_res_fn="rectangular")
        bare = m.SlitGeometry(s1, si, 1350.0, 135.0, distance)
        series = [
            ("full (constantQ)", m.angle_offset(dy, pixel, distance, "full"), "#eb6834", "-"),
            ("half (legacy)", m.angle_offset(dy, pixel, distance, "half"), "#2a78d6", "--"),
            ("meanTheta, no detector response", m.mean_theta_offset(dy * pixel, bare)[0], "#1baf7a", ":"),
            ("meanTheta, DetSigma 0.8 mm", m.mean_theta_offset(dy * pixel, geometry)[0], "#1baf7a", "-"),
        ]
        a = np.arctan(np.array([s1 + si, s1 - si, -(s1 + si), -(s1 - si)]) / (2 * 1350.0))
        y_beam = (np.array([si, -si, -si, si]) / 2 + (135.0 + distance) * np.tan(a)) / pixel
        ax.axvspan(y_beam.min(), y_beam.max(), color=GRID, alpha=0.8, linewidth=0)
        for label, offset, color, style in series:
            ax.plot(dy, np.degrees(offset), color=color, linestyle=style, linewidth=1.6, label=label)
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
        ax.set_xlabel("row offset from the specular pixel, y − y_ref (pixels)")
    axes[0].set_ylabel("row angle − θ_central (degrees)")
    axes[0].legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=Path("figures"))
    parser.add_argument("--poisson-csv", type=Path, help="synthetic-poisson_summary.csv from compare_q_methods.py")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plot_plane(args.output / "q_methods_plane.png")
    plot_plane_tof(args.output / "q_methods_plane_tof.png")
    plot_tof_film(args.output / "q_methods_tof_film.png")
    plot_deposit(args.output / "q_methods_deposit.png")
    plot_flat(args.output / "q_methods_flat.png")
    plot_row_weights(args.output / "q_methods_row_weights.png")
    plot_angle_models(args.output / "q_methods_angle_models.png")
    if args.poisson_csv is not None:
        plot_poisson(args.output / "q_methods_poisson.png", args.poisson_csv)


if __name__ == "__main__":
    main()
