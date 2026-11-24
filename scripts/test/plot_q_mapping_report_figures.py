"""
The figure added for REPORT_Q_mapping_comparison.md (the other figures come from plot_q_method_diagrams.py).

`report_q_sweep.png` — median deviation against λ bin width, `half` row angles, for the noise-free
synthetic film and run REF_L_198415, from the `*_summary.csv` outputs of `compare_q_methods.py`.

    python plot_q_mapping_report_figures.py --output figures \
        --synthetic-csv <dir>/synthetic_summary.csv --real-csv <dir>/REF_L_198415_summary.csv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

INK = "#2b2b2b"
MUTED = "#8a8a85"
GRID = "#e4e4e0"
COLORS = {"E": "#8a8a85", "A": "#2a78d6", "B": "#eb6834", "C": "#1baf7a", "C0": "#1baf7a", "legacy": "#e8a317"}
LABELS = {
    "E": "E events",
    "A": "A point",
    "B": "B prototype",
    "C": "C overlap, corner footprint",
    "C0": "C0 overlap, row angle",
    "legacy": "legacy specular_weighted",
}
MARKERS = {"E": None, "A": "o", "B": "s", "C": "^", "C0": "v", "legacy": "D"}
STYLES = {"E": "--", "A": "-", "B": "-", "C": "-", "C0": ":", "legacy": "-."}

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def _read_summary(path: Path, convention: str) -> dict[str, list[tuple[float, float]]]:
    out: dict[str, list[tuple[float, float]]] = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            if row["convention"] != convention or not row["method"] or not row["median_abs_dev_pct"]:
                continue
            out.setdefault(row["method"], []).append((float(row["lam_step"]), float(row["median_abs_dev_pct"])))
    return {k: sorted(v) for k, v in out.items()}


def plot_sweep(path: Path, synthetic_csv: Path, real_csv: Path) -> None:
    panels = (
        (synthetic_csv, "Noise-free synthetic film, reference: the truth"),
        (real_csv, "REF_L_198415, reference: E"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharey=True)
    for ax, (csv_path, title) in zip(axes, panels):
        data = _read_summary(csv_path, "half")
        for key in ("E", "A", "B", "C", "C0", "legacy"):
            if key not in data:
                continue
            x, y = np.array(data[key]).T
            ax.plot(x, y, color=COLORS[key], ls=STYLES[key], marker=MARKERS[key], label=LABELS[key], lw=1.5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("λ bin width (Å)")
        ax.set_title(title, color=INK, fontsize=10, loc="left")
        ax.grid(True, color=GRID, lw=0.8, which="both")
    axes[0].set_ylabel("median |deviation| (%)")
    handles, labels = [], []
    for ax in axes:
        for h, label in zip(*ax.get_legend_handles_labels()):
            if label not in labels:
                handles.append(h)
                labels.append(label)
    fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=Path("figures"))
    parser.add_argument("--synthetic-csv", type=Path, required=True)
    parser.add_argument("--real-csv", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plot_sweep(args.output / "report_q_sweep.png", args.synthetic_csv, args.real_csv)


if __name__ == "__main__":
    main()
