# `scripts/test`

Measurements behind claims made in the test suite. Anything a committed file
cites has to live in the repository so it can be re-run.

| Script | Answers | Cited from |
|---|---|---|
| `measure_fit_path_dependence.py` | How far do the scaling-factor fit parameters move when only Mantid's minimizer changes? | `tests/test_scaling_factors_workflow.py`, the `_TOL` bar for `b` |
| `roi_estimate_mutations.py` | Does each guard and rule in `src/lr_reduction/roi_estimate.py` have a test that fails without it? One mutation per row, the tests that went red, and the measured table at the file's end. | `tests/unit/lr_reduction/test_roi_estimate.py` (I5–I7 test the battery itself); the roi-estimate and roi-popout-data commit bodies |

Run from the repository root:

```sh
pixi run python scripts/test/measure_fit_path_dependence.py   # add --verbose for workflow/Mantid output
pixi run python scripts/test/roi_estimate_mutations.py --rows 1-25   # then --rows 26-51
```

Both need the test-data submodule at `tests/data/liquidsreflectometer-data`.
The fit measurement takes a few minutes. The battery takes about three minutes
per 25 rows, so run it in chunks under a 600 s limit. It refuses to start
unless `roi_estimate.py` matches `HEAD`, and it restores the file after every
row.
