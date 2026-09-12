# `scripts/test`

Measurements behind claims made in the test suite. Anything a committed file
cites has to live in the repository so it can be re-run.

| Script | Answers | Cited from |
|---|---|---|
| `measure_fit_path_dependence.py` | How far do the scaling-factor fit parameters move when only Mantid's minimizer changes? | `tests/test_scaling_factors_workflow.py`, the `_TOL` bar for `b` |

Run from the repository root:

```sh
pixi run python scripts/test/measure_fit_path_dependence.py   # add --verbose for workflow/Mantid output
```

Needs the test-data submodule at `tests/data/liquidsreflectometer-data` and
takes a few minutes.
