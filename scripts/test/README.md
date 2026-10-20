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
pixi run python scripts/test/roi_estimate_mutations.py --rows 1-30   # then --rows 31-62, and --rows 52 --with-slow
```

Both need the test-data submodule at `tests/data/liquidsreflectometer-data`.
The fit measurement takes a few minutes. The battery takes about four and a half
minutes per 30 rows, so run it in chunks under a 600 s limit. Its slow census
test (T1b, ~45 s a run) is left out unless `--with-slow`. It refuses to start
unless `roi_estimate.py` matches `HEAD`, and it restores the file after every
row.
