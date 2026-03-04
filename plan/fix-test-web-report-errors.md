# Plan: Fix Pre-existing `test_web_report.py` Errors

## Context

Branch `bvacaliuc/new_workflow_ui_plan` (and prior branches) show 7 errors in
`tests/unit/lr_reduction/test_web_report.py`.  These are **not regressions** —
they fail because the tests require real NEXUS instrument data from a git-lfs
submodule that is not initialized locally.

## Root Cause

All 7 test errors trace to the same fixture chain:

```
workspace_sc (fixture)  →  LoadEventNexus("REF_L_201288")  →  File not found
workspace_db (fixture)  →  LoadEventNexus("REF_L_198399")  →  File not found
template_data           depends on  workspace_sc
meta_data               depends on  workspace_sc + workspace_db
```

The NEXUS files live in a git-lfs submodule at
`tests/data/liquidsreflectometer-data/`.  Current submodule status:

```
$ git submodule status
-872ef741...  tests/data/liquidsreflectometer-data   (prefix `-` = not initialized)
```

The `nexus/` subdirectory inside it does not exist locally, so Mantid's
`LoadEventNexus` raises `ValueError: File "REF_L_201288" not found`.

The tests pass in CI because the GitHub Actions workflow initializes submodules
before running pytest.

## Error Messages

All 7 tests produce the same error type at fixture setup:

```
ERROR at setup of test_generate_report_section_reduction_parameters
E   ValueError: Invalid value for property Filename (string) from string "REF_L_201288":
    When setting value of property "Filename": File "REF_L_201288" not found
```

## Affected Tests

| Test function | Fixture dependency |
|---|---|
| `test_generate_report_section_reduction_parameters` | `workspace_sc`, `template_data`, `meta_data` |
| `test_generate_report_plots_reflected_beam` | `workspace_sc`, `template_data` |
| `test_generate_report_plots_direct_beam` | `workspace_db`, `template_data` |
| `test_generate_report_section_run_meta_data` | `workspace_sc` |
| `test_generate_report_sections` | `workspace_sc`, `template_data`, `meta_data` |
| `test_generate_report_section_direct_beam` | `workspace_db`, `template_data` |
| `test_assemble_report` | `workspace_sc`, `template_data`, `meta_data` |

## Fix Options

### Option A: Initialize the Submodule Locally

**Effort:** Minimal (one command)
**Trade-off:** Downloads ~100+ MB of NEXUS files via git-lfs

```bash
git submodule update --init tests/data/liquidsreflectometer-data
```

This is the intended workflow for developers who want to run the full test suite
locally, matching what CI does.  No code changes required.

**When to choose:** If you want to run the full integration test suite locally
and have git-lfs available.

---

### Option B: Mark Tests with `@pytest.mark.datarepo` and Skip by Default

**Effort:** Low (~10 lines changed)
**Trade-off:** Tests are skipped locally unless explicitly requested

The `datarepo` marker is already defined in `pyproject.toml`:
```toml
markers = [
    "datarepo: mark a test as using LiquidsReflectometer-data repository",
]
```

But the `test_web_report.py` tests don't use it.

**Implementation:**

1. Add the marker to the test file:

```python
# tests/unit/lr_reduction/test_web_report.py  (top of file, after imports)
pytestmark = pytest.mark.datarepo
```

2. Run locally with skip:
```bash
pixi run python -m pytest -m "not datarepo"
```

3. CI continues to run all tests (no `-m` filter in the CI command).

**When to choose:** If you want clean local test output but don't want to
change the tests' behavior when data IS available.  Requires remembering (or
configuring) the `-m` flag locally.

---

### Option C: Auto-skip When Data Submodule is Missing (Recommended)

**Effort:** Low (~5 lines added)
**Trade-off:** None — tests run when data exists, auto-skip when it doesn't

**Implementation:**

Add a module-level skip condition at the top of `test_web_report.py`:

```python
import os
import pytest

_DATA_SUBMODULE = os.path.join(
    os.path.dirname(__file__), "..", "..", "data", "liquidsreflectometer-data", "nexus"
)
pytestmark = pytest.mark.skipif(
    not os.path.isdir(_DATA_SUBMODULE),
    reason="NEXUS data submodule not initialized (run: git submodule update --init)"
)
```

**Result:**
- Locally (no submodule): all 7 tests show `SKIPPED` with a clear reason message
- CI (submodule initialized): all 7 tests run normally
- No changes to CI config or pytest invocation needed

**When to choose:** This is the cleanest option.  It makes local test runs
green without requiring LFS data, while keeping the tests fully active in CI.
The skip reason also tells developers how to enable the tests if they want to.

---

### Option D: Leave As-Is

**Effort:** None
**Trade-off:** 7 `ERROR` lines in local test output

The errors are clear about the cause (`File not found`) and don't mask real
problems.  CI passes because the submodule is initialized there.

**When to choose:** If the errors don't bother you and you don't want to touch
the file at all.  This is the simplest approach — the errors are informational,
not harmful.

---

## Recommendation

**Option C** is recommended because:

1. Zero impact on CI behavior
2. Local test output is clean (SKIPPED instead of ERROR)
3. Skip message tells developers exactly how to enable the tests
4. No need to remember special pytest flags
5. Only 5 lines of code added to one file

## Additional Observations

These tests also have a few secondary characteristics worth noting for future
maintenance:

- **Fixture scope mismatch:** `conftest.py` defines `nexus_dir` as
  `scope="session"` while `test_web_report.py` fixtures use `scope="module"`.
  This works but means workspaces are reloaded per-module rather than shared.
- **Hardcoded assertion lengths:** `test_generate_report_section_reduction_parameters`
  asserts `len(report) == 942` and `test_generate_report_section_run_meta_data`
  asserts `len(report) == 132`.  These are brittle — any change to HTML
  formatting or report content will break them.  Consider switching to
  `assert len(report) > 0` or checking for specific content strings.
- **No `conftest.py` in `tests/unit/lr_reduction/`:** All fixtures come from
  `tests/conftest.py` via pytest's directory traversal.  This works but could be
  made more explicit.

## Files to Modify

| File | Action | Lines Changed |
|------|--------|---------------|
| `tests/unit/lr_reduction/test_web_report.py` | **Modify** | Add ~5 lines (module-level skipif) |

## Verification

After applying Option C:

```bash
# Without submodule — should show 7 SKIPPED
pixi run python -m pytest tests/unit/lr_reduction/test_web_report.py -v

# With submodule — should show 7 PASSED
git submodule update --init tests/data/liquidsreflectometer-data
pixi run python -m pytest tests/unit/lr_reduction/test_web_report.py -v
```
