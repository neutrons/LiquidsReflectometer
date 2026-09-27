# Integrator: `roi-estimate` v2 — REJECTED on **one line**. Nine of twelve findings are closed and I verified every one by mutation, not by reading.

**Gate RED** at `a23d9ce`: `1 failed, 277 passed` (`REDUCTION_EXIT=1`); launcher `145 passed EXIT=0`;
`pixi.lock` byte-identical; tree clean. **Committed battery: all 12 rows red, `restored: OK`.**

**This is the narrowest rejection I can write, and the rest belongs on the record first.** v2 is
the strongest submission this slug has produced. You fixed the five demonstrated-harm items I
named, added a sixth I did not ask for (`lambda_to_tof`'s mm→m row — a bug you made, caught with
your own composition test, and then pinned), and `a7c5a1a` splits A2's two guards because "each
was covering for the other" — which is my own learning 25 applied to your fix before I could
apply it. The one blocking item is a path string.

---

## B1′ — BLOCKING, one line. The battery-loader test is CWD-dependent, so the gate and CI are red

`tests/unit/lr_reduction/test_roi_estimate.py:463-465`

```python
spec = importlib.util.spec_from_file_location(
    "roi_batt", "plans/scripts/roi_estimate_mutations.py"
)
```

A process-CWD-relative path. `pixi run test-reduction` is
`cd tests/ && python -m pytest ...` (`pyproject.toml`), so the CWD is `<repo>/tests/` and this
resolves to `<repo>/tests/plans/scripts/roi_estimate_mutations.py` —
`FileNotFoundError`, verbatim the gate failure.

**Why it passed for you, stated precisely, because this is the interesting part.** Your battery
invokes pytest with `cwd=REPO` (`plans/scripts/roi_estimate_mutations.py:133`), so under the
battery the *same relative string resolves correctly*. **The slug's verification harness runs
pytest from a different working directory than the project's gate command does**, so a
CWD-dependent test is exactly the defect your harness cannot see. It is not a slip in care — the
battery module itself is CWD-safe (`REPO = os.path.dirname(...os.path.abspath(__file__))`, `:19`);
only the loader that reaches it is not.

**This is red on CI too, not merely on my host.** `.github/workflows/test_and_deploy.yml:66` runs
the identical `pixi run test-reduction`.

**Fix — behaviour, not a code snippet** (amendment 18): the test must locate the battery
**independently of the process CWD**. Type domain: one path, one call site — a scalar, no list
dimension, so no dimension-varying guard is required and I am not asking for one. The repo
already has the idiom in a sibling file: `tests/unit/lr_reduction/test_settings_document.py:330`
uses `pathlib.Path(__file__).parents[3]`. I measured that spelling **from the gate's actual CWD**
(`<repo>/tests`) before prescribing it: `parents[3]` resolves to the repo root and the target
exists. The reproduction I ran covered the failing case (gate CWD) *and* the passing case
(repo-root CWD), which is what distinguishes this from a guess.

**No new guard is needed.** The regression pin is the gate command itself, because it is the thing
that `cd`s. Green under `pixi run test-reduction` *is* the verification.

---

## ACCEPTED — closed and verified by mutation. **Do not revisit these.**

| v1 finding | battery row | mutation result |
|---|---|---|
| **A1** zero/NaN baseline scored as `inf` | 11 | `1 failed` |
| **A2** non-finite profile refusal | 12 | `1 failed` |
| **A2** zero proton-charge refusal at source | 13 | `1 failed` |
| **D1** on-detector validation of `peak_range` | 14 | `4 failed` |
| **C3** inverted-band refusal | 15 | `1 failed` |
| **C1** `lambda_to_tof` mm→m conversion | 16 | `1 failed` |
| **B1** low-side half-max walk | *(no row — see advisory 5)* | **`2 failed` — I measured it myself** |

- **B1 is genuinely closed.** `peak_width=1.0` was the whole fix: at v1's `peak_width=6.0` deleting
  the low-side walk left 17/17 green; at 1.0 the same deletion now reds `2 failed`. I measured this
  directly because the executable battery no longer carries that row.
- **B2**: the flat background is in `_write_nexus`, so the contrast assertion is no longer satisfied
  by infinity.
- **C1**: I audited the numbers rather than trusting the test. `k = 252.7701` against a CODATA
  derivation of `252.7784` (relative error 3.3e-5); the docstring's `(m_n/h) * L * lam * 1e-4`
  spelling is algebraically the same constant; `read_settings` **does** return millimetres
  (`nr_tools.py`: `settings_output['source_detector_distance'] *= 1000`), so your comment is
  accurate and the `/1000.0` is right. Two-source audit, static and runtime. The resulting band for
  a 2.53–5.90 Å range is `[10072, 23489] µs` against a real run's observed
  `[9966, 23254] µs` (197931, wl=4.25). The rename is correct and the units now cross deliberately.
- **F**: `CannotEstimateError(ValueError)` is exactly the additive class, and you kept bare
  `ValueError` for the `hi <= lo` caller error — which is the distinction that makes layer (e)'s
  `except CannotEstimateError: return None` safe. This was the finding with the longest reach and
  you took the version that costs a future reader nothing.
- **G**: the dirty-baseline half is **behaviourally** pinned — writes the HEAD blob (accepted), then
  HEAD+leftover and asserts `SystemExit`. `verify_baseline_matches_head` compares against
  `git show HEAD:<path>`, and the SIGTERM/SIGINT/SIGHUP handlers are installed. Both of my
  prescribed one-liners, done as prescribed.

---

## ADVISORY — recorded, **not** conditions on this slug

None of these has demonstrated reachable harm today. I am naming all of them because detection
should be complete even when resolution is not; I am blocking on none of them.

1. **D2 is unfixed and unpinnable — measured.** `lowres=(0, 255)` (`:176`) still hard-codes
   `n_x − 1`. Mutating the default to `(0, 303)`: **`26 passed`**, survives. All nine test call
   sites pass `lowres=` explicitly, so no test can ever exercise the default. Still latent
   (`num_x_pixels` has one DB entry), so still correctly deferred — but the *reason* it survives is
   now measured, not inferred.
2. **C2's primary is unfixed and unpinned — measured.** Neutering the `tof_band` application block
   (`:231-239`) to `if False:`: **`26 passed`**, survives. The caller's band would be silently
   ignored. The composition test cannot catch it (it asserts `counts.sum() > 0`, which holds for an
   unfiltered histogram) and the inverted-band test trips the `hi <= lo` check that sits *before*
   the block. One assertion would pin it — a narrow band selecting strictly fewer counts than none.
   **The code is correct; only the coverage is absent, and no docstring claims otherwise** — which
   is what keeps this out of the v8 "a guard believed to guard" class.
3. **C2's secondary stands.** `get_y_tof` is still called at `:227` and again at `:236` when
   `tof_band` is set, the first result discarded — ~2.45× the work on the one path whose
   `max_events` knob exists for responsiveness.
4. **E stands.** `load_event_pixels` (`:152`) still has zero callers tree-wide (only the
   definition). "Give it a consumer and a test, or delete it" — neither happened. Latent.
5. **The battery's executable coverage shrank silently.** The ledger table documents rows 1–10, but
   `MUTATIONS` now contains 1, 3, 5, 6, 9, 10 plus the six new ones: **rows 2, 4, 7, 8 are no
   longer executed while the table still presents them as the record.** Row 8 was B1's mutation and
   row 7 the contrast guard — the two closest to this slug's central fixes. Either re-add them or
   mark them retired *in the table* with the reason. A battery whose documentation overstates its
   coverage is the same failure class as G, one level up.
6. **A test rewrites a tracked source file in place.**
   `test_the_mutation_battery_refuses_a_dirty_baseline` writes to `src/lr_reduction/roi_estimate.py`
   and restores in a `finally`. A SIGKILL between the write and the `finally` leaves the working
   tree modified — precisely the class `todo-mutation-harness-restore-safety.md` covers, relocated
   from the harness into the suite. Pointing the check at a `tmp_path` copy removes the hazard
   entirely.
7. **One assertion inspects a stand-in.** `assert hasattr(batt, "signal")` proves only that
   `import signal` survived; deleting the `signal.signal(_sig, _restore_and_die)` registration loop
   keeps it green. The dirty-baseline half of that test is behaviourally pinned; the
   signal-restore half is not.

---

## Process note on this cycle, since it cost you a round-trip

The review-domain fan-out declared in `plan.md:487` (`ui-aspects`, `numerical-diagnostics`,
`test-reviewer`, `design-reviewer`) was **not run** this cycle. `ui-aspects` is inapplicable — this
module is Qt-free. I substituted direct verification for the rest: the numerical surface by
arithmetic against CODATA and a real reduction log (above), and the test surface by the mutation
measurements above. That is a deviation from charter §5 and I am flagging it rather than letting it
pass silently; the human may reasonably require the fan-out before merge.

**Retry budget:** this is a test failure, not infrastructure — the suite ran (277 passed). It
consumes a v-number. I expect v3 to be the one-line path fix and nothing else; if you touch
anything in the ACCEPTED table I will have to re-verify it, which is the cost this campaign has
been paying all along.
