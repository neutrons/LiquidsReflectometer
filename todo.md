# Integrator review gate — `scaling-factor-path-anchor` v1 REJECTED

**Disposition: blocking, per the plan's own declaration.** `test-reviewer` is
declared blocking for this bug-fix phase and returned three blocking findings;
no Integrator departure was needed. Retry budget: attempt 1 of N = 3.

## The fix itself is right, and one judgment call in it is exactly right

State this first, because it is the larger part of the truth:

- **Both halves work.** Verified independently: the two new tests go **RED**
  against the base implementation — *while running from `tests/`*, which is the
  cwd where the original bug is invisible. The `monkeypatch.chdir` design
  genuinely defeats the gate's structural blindness. That is a good answer to a
  hard problem.
- **`raise` over `return 1, 0, 0, 0` was the correct call.** The tidy
  alternative (matching the sibling exit 80 lines below) would silently produce
  **unscaled reflectivity** — plausible-looking R(Q) that a reader takes at face
  value. The design reviewer traced every caller chain and confirmed nothing
  regresses: the single call site is guarded by `scaling_factor_flag`, the
  sanctioned no-scaling path (`flag = False`) is untouched, and autoreduce
  already exited non-zero on the old `TypeError`.
- **Blast radius verified safe** for absolute paths, empty strings, and
  cwd-relative paths that already resolve (as-given is tried first). The
  round-trip is contained: `workflow.write_template` re-reads the raw XML, so a
  machine-specific absolute path is not persisted into
  `REF_L_*_auto_template.xml`.

## B1 — BLOCKING: "109 passed from the repo root" is real but *accidental*

The acceptance claim does not hold for the reason it appears to. I measured
this from both sides and then found the mechanism.

**My measurements, same clone, same commit:**

| invocation | result |
|---|---|
| full suite from repo root | **109 passed** |
| `tests/test_scaling_factors_workflow.py` alone, from root | **5 failed** |
| the 5 tests by node id, from root | **5 failed** |
| `tests/test_reduction.py` alone, from root | 17 passed |

**Cause — an unrestored `chdir` in an unrelated test.**
`tests/test_reduction.py:230`, inside `test_reduce_functional_bck`:

```python
os.chdir(Path(template_dir).parent)     # -> <repo>/tests ; never restored
```

In a full-suite run from the repo root, every test that executes **after** that
line inherits `cwd = tests/`, so every cwd-relative path in the rest of the
session resolves. The suite's root-invocation green is test pollution, not a
property of this fix.

**This accounts for the original nine exactly.** The failing tests were
`test_dead_time.py::test_full_reduction` (runs before `test_reduction.py`) plus
five functions in `test_reduction.py` defined at lines 53, 102, 126, 180 and 202
— **all before line 229** — with `test_q_summing_as_option` parametrized ×4, so
1 + 8 = 9. The two `test_reduction.py` tests that *passed* in the human's
transcript, `test_reduce_workflow_201282` (`:365`) and
`test_background_subtraction` (`:401`), are the two defined **after** the chdir.

**Consequence beyond this slug:** the suite's cwd-independence is currently
unverifiable by running it. Anything after line 229 is exercised from `tests/`
no matter how pytest is invoked, so no full-suite run can tell you whether those
tests are cwd-clean.

**Correcting the blocking reviewer's hypothesis, and my own first reading.**
The review proposed that a green root run implied "a tree that had a stray
root-level `data/`". I checked: there is **no** `data/` at my repo root, the
working tree is clean, and nothing in the suite calls `makedirs`. The real
mechanism is the chdir above. Separately, my clone does carry **16 gitignored
generated artifacts** under `tests/data/` (`REFL_*_partial.txt`,
`*_combined_data_auto.txt`, `REF_L_*_auto_template.xml`, `options.json`) that a
`git archive` copy lacks — a genuine reason pristine and working clones can
disagree, and worth knowing, but not what produced this result.

My own first report of "109 passed from the repo root" was therefore misleading,
and I am correcting it here rather than letting it stand.

**Fix:** restore the cwd (`monkeypatch.chdir`, or a `try/finally`) at
`test_reduction.py:230`. Then re-measure the root invocation honestly — expect
the residual failures below to surface, and decide their scope deliberately.

## B2 — BLOCKING: the headline test does not pin the feature

Verified by the reviewer: absolutize the path in `tests/data/template.xml` **and**
delete the `_resolve_scaling_factor_file(...)` call from `read_template` → **2
passed**. The whole anchoring feature is deletable provided someone edits the
fixture, because `test_template_scaling_factor_path_is_cwd_independent` asserts
the *outcome* (`os.path.isfile`) but never the *precondition* that the template
still declares a relative path.

Fix: assert the fixture still declares a relative path, and use
`os.path.samefile` against the expected `tests/data/sf_197912_Si_auto.cfg` rather
than bare `isfile` — which also closes the decoy hole that a 3-level upward walk
opens.

## B3 — BLOCKING (cheap): the stated invariant is untested

`candidates = [scaling_factor_file]` → `candidates = []` — deleting as-given-first
precedence entirely — leaves both tests green. That precedence is the documented
protection for facility `/SNS/...` templates and for the `cd tests/` gate. A
six-line test with a decoy at both the cwd-relative and anchored locations pins it.

## Advisory — from the design domain, and worth acting on

1. **The anchoring fallback can silently bind the wrong file, unlogged.**
   Probed against a simulated `/SNS/REF_L/IPTS-1234/shared/autoreduce/` layout:
   a stale same-named `sf_*.cfg` two levels up is found and used, R(Q) is
   multiplied by a scaling factor from a file the template author never named,
   and **the user sees nothing** — no log, and the choice is not in the output
   metadata. This is the same silent-wrong-number hazard the `raise` half was
   added to close; the two halves pull in opposite directions on the same axis.
   **One line fixes it:** when the winning candidate differs from the declared
   string, log both. That turns a guess into an auditable fallback.
2. **Record the resolved file in the output metadata**
   (`meta_data["scaling_factor_file"]`). Pre-existing gap, but load-bearing now
   that the path used may differ from the path declared.
3. **The larger silent hazard on the same axis is untouched:**
   `template.py:207` returns `1, 0, 0, 0` when the file parses but no entry
   matches the measured slits/wavelength — no scaling, no log. That is the
   *routine* REF_L failure, far likelier than an absent file, and it produces
   the identical unscaled-but-plausible R(Q). At minimum log it.
4. `_SCALING_FACTOR_ANCHOR_DEPTH = 3` is off-by-one against its own comment
   (iteration 1 anchors *at* the template dir, not above it); the suite pins only
   "≥ 2", and depth 4 would leave the repository entirely. Depth 2 with a clearer
   name, or a test that exercises level 3.
5. `new_reduction_from_template.py:460` is a verbatim fork of `read_template`
   carrying its own `# TODO: Fix to use the one inside template.py`, and did not
   receive the anchoring — so the same template now resolves differently
   depending on which reader you call. Same defect class this branch fixes.
6. `time_resolved.py:257` and `:374` are bare `except:` blocks that swallow the
   new `FileNotFoundError` into a print and continue. Not a regression, but
   fail-loud is not global — worth one line in the PR body.

## Residual cwd-relative paths (enumerated, for the scope decision)

Reads that fail from the root: `test_time_resolved.py:19-20, 56-57`;
`test_scaling_factors_workflow.py:72, 91, 119, 147, 176`. Writes that fail at
`output.py:185`: `test_reduction.py:370, 406`; `test_time_resolved.py:18, 55`.
All are the same one-line `template_dir`/`tmp_path` substitution already applied
once in this diff at `test_reduction.py:159`.

**Either fix them or renegotiate the acceptance criterion explicitly** — but do
not claim root-pytest green while a leaked chdir is what produces it.

## Gate result

`pixi run test-reduction` (charter gate, cwd `tests/`): `test-launcher` **10
passed**; reduction **108 passed, 1 failed**. The single failure is
`test_scaling_factors_workflow.py::test_compute_sf` — the known cross-clone
shared-`/tmp` race, **6th occurrence**. Confirmed not attributable to this slug:
the `scaling_factors` package never imports `template`, the test imports only
`sf_workflow`, and the slug's only source change is `template.py`. Verified
passing in isolation (46 s). No retry consumed.

## Suggested v2 order

1. **B1's mechanism** — restore cwd at `test_reduction.py:230`, then re-measure
   the root invocation and state the real number.
2. **B2 + B3** — pin the fixture precondition and the as-given-first precedence.
3. **Advisory 1** — log declared → resolved. This is the one that keeps the fix
   from re-opening the hazard it was written to close.
4. Decide the residual-path scope deliberately, and correct the commit/PR body.
