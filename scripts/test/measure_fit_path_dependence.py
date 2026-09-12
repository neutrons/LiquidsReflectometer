#!/usr/bin/env python
"""Measure how far the scaling-factor fit moves when only the minimizer changes.

`tests/test_scaling_factors_workflow.py` sizes the tolerance on `b` against
this "path-dependence floor": the shift in a converged fit parameter when the
input data is identical and only Mantid's route to the optimum differs.

    pixi run python scripts/test/measure_fit_path_dependence.py [--verbose]

Needs tests/data/liquidsreflectometer-data (the suite's fixture); runs the
scaling-factor workflow once per minimizer, a few minutes in total. Workflow
and Mantid chatter is suppressed unless --verbose is given.
"""

import argparse
import contextlib
import os
import tempfile

import mantid.simpleapi as mtd_api

mtd_api.config["default.facility"] = "SNS"
mtd_api.config["default.instrument"] = "REF_L"

from lr_reduction.scaling_factors import LRScalingFactors  # noqa: E402
from lr_reduction.scaling_factors import workflow as sf_workflow  # noqa: E402
from lr_reduction.utils import amend_config  # noqa: E402

# The floor is measured within the Levenberg-Marquardt family: LM-MD (Mantid's
# default) and LM reach the same optimum by different routes. Simplex is a
# control only — a different algorithm with a looser stopping criterion — and
# is excluded from the floor.
LM_FAMILY = ("Levenberg-MarquardtMD", "Levenberg-Marquardt")
CONTROL = ("Simplex",)
MINIMIZERS = LM_FAMILY + CONTROL

FIELDS = ("a", "b", "error_a", "error_b")


def _rows(path):
    rows = []
    with open(path) as fd:
        for line in fd:
            if line.startswith("#") or not line.strip():
                continue
            rows.append(dict(tok.split("=", 1) for tok in line.split() if "=" in tok))
    return rows


def _run(workspace, minimizer, out_dir):
    """Run the workflow with `Fit` forced onto one minimizer."""
    original = LRScalingFactors.Fit

    def patched(*args, **kwargs):
        kwargs["Minimizer"] = minimizer
        return original(*args, **kwargs)

    LRScalingFactors.Fit = patched
    try:
        sf_workflow.process_scaling_factors(
            workspace, out_dir, use_deadtime=True, deadtime=4.6,
            deadtime_tof_step=200, paralyzable=False, wait=False, postfix="_probe",
        )
    finally:
        LRScalingFactors.Fit = original
    return _rows(os.path.join(out_dir, "sf_197912_Si_probe.cfg"))


def _worst_shift(baseline_rows, other_rows):
    worst = dict.fromkeys(FIELDS, 0.0)
    for base_row, other_row in zip(baseline_rows, other_rows):
        for field in FIELDS:
            ref = float(base_row[field])
            if ref == 0.0:
                continue
            worst[field] = max(worst[field], abs((float(other_row[field]) - ref) / ref))
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--nexus-dir", default="tests/data/liquidsreflectometer-data/nexus")
    parser.add_argument("--verbose", action="store_true", help="show workflow and Mantid output")
    args = parser.parse_args()

    if not args.verbose:
        mtd_api.config.setLogLevel(3)  # errors only
    quiet = contextlib.nullcontext() if args.verbose else contextlib.redirect_stdout(open(os.devnull, "w"))

    with quiet, amend_config(data_dir=os.path.abspath(args.nexus_dir)):
        workspace = mtd_api.Load("REF_L_197912")
        results = {}
        for minimizer in MINIMIZERS:
            with tempfile.TemporaryDirectory() as out_dir:
                results[minimizer] = _run(workspace, minimizer, out_dir)

    baseline = MINIMIZERS[0]
    print(f"worst relative shift vs {baseline}, over {len(results[baseline])} rows")
    print(f"{'minimizer':<26}" + "".join(f"{f:>13}" for f in FIELDS))
    for minimizer in MINIMIZERS[1:]:
        worst = _worst_shift(results[baseline], results[minimizer])
        print(f"{minimizer:<26}" + "".join(f"{worst[f]:>13.3e}" for f in FIELDS))

    floor = dict.fromkeys(FIELDS, 0.0)
    for minimizer in LM_FAMILY[1:]:
        for field, shift in _worst_shift(results[LM_FAMILY[0]], results[minimizer]).items():
            floor[field] = max(floor[field], shift)

    print("\npath-dependence floor (Levenberg-Marquardt family only; Simplex is a control):")
    for field in FIELDS:
        print(f"  {field:<10} {floor[field]:.3e}")


if __name__ == "__main__":
    main()
