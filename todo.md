# todo.md — Integrator rejection, `roi-popout-dialog` v4 @ e964955 (attempt 4 of 4: the human's one extension; review gate: ui-aspects, design, test)

**Verdict: REJECT, tests only, one block.** v4 used the human's single extension, so this goes to the Analyst's
escalation and the human. There is no v5 on any seat's authority (plan v4 "Authority").

**Production is right.** On real data (IPTS-36119) the deployment-shaped acceptance has now passed four times
(v1–v4), each with identical numbers. The audit table works: every row this seat re-ran reproduces its recorded count
exactly. One credited row (C1) is credited to a test that cannot see the clause's composition.

## What passed (do not redo)

- **Gate:** `pixi run test-reduction`, analysis clone 2 (analysis-node01), 18:17–18:26 EDT: launcher **761 passed**
  (757 + the four V14 legs), reduction **874 passed**, exit 0, clean.
- **The battery's v3/v4 rows, re-run here** in a `git archive` copy of e964955 (`__file__` printed; ledger
  `scripts/mutations-roi-popout-dialog.py` @ 8e6c878, run with the gate env's python): every count equals the commit
  body's. Restore check 2/2.

  | Rows | Result |
  |---|---|
  | M30, M31 | 2, 2 |
  | M32 | 8 |
  | M33, M33b | 1, 1 |
  | M34, M34b | 1, 1 |
  | M35 | 1 |
  | M36 | 2 |
  | M37, M37b | 2, 1 |
  | M38 | 3 |
  | M39 | 1 |
  | M40 | 1 |
  | M41 | 2 |
  | M42 | 1 |
- **Group 1 (P1–P5), A-i, A-ii, A-iii, A-iv:** done as the plan wrote them.
  - Design PASS: every claim in `select_roi`'s new docstring is asserted (E4: save, panel, document, cwd, run folder
    by size/mtime/content, QSettings; E9: panel, document, cwd). Production is byte-identical to v3 except that docstring.
  - ui-aspects PASS on B3/B12/B5: 7 axes, bbox placement, LogNorm at open and after a draw, aspect after a drag on
    both axes and both branches. Its own probe of the Y-TOF branch: 1 failed.
- **Group 2 C2, C3, C4:** credited correctly. The test reviewer's own drops of `refresh_report` alone and of
  `refresh_scalars` alone each red E8 (1 failed). C4 → M40 red (see A-3).
- **Rule (c):** the 174 / 170 / 4-abort and 632 counts are consistent.
- **§8.7 acceptance re-run at e964955** (ledger `scripts/roi-popout-acceptance.py`; reduce_settings.json, REF_L_231801.html,
  the 231800 template; report row 1, edit row 2): **PASS, 53 ok, exit 0**.
  - Dead-time medians 1.0099 and 1.0086.
  - Edited save sha256 `ee0475b3…`, identical to v1–v3.

## BLOCKING (rules a, b; plan v4 Group 2 C1; reproduced here)

### B-1: C1, "the view filter starts at the chopper band when the run has a chopper log", is unpinned where the slot composes the band

- **The credit does not meet C1's own prescription.** C1 was credited to E10
  (`test_the_run_is_titled_from_its_metadata_and_filtered_at_its_chopper_band[from its metadata]`). E10 stubs both
  data-layer calls with lambdas that **ignore their arguments**:
  - `chopper_lambda_range = lambda _path: (2.5, 9.5)`;
  - `lambda_to_tof = lambda _band, _start: (12000.4, 31000.6)`.

  So E10 sees whether a band arrives, never which band. The plan's C1 test asked for "`tof_spins` read
  `lambda_to_tof(chopper_lambda_range(...))` (the data layer's numbers, not literals)".
- **Mutation** (`launcher/apps/settings_editor.py:1217`, one line; an argument-order slip):

  ```
  band = roi_estimate.lambda_to_tof(roi_estimate.chopper_lambda_range(path), meta["start_time"])
  → band = roi_estimate.lambda_to_tof(meta["start_time"], roi_estimate.chopper_lambda_range(path))
  ```

  With it, `test_settings_editor.py` + `test_roi_dialog.py` give **632 passed, exit 0**. That is the unmutated count:
  the mutant survives. Reproduced in a fresh `git archive` copy of e964955, with `settings_editor.__file__` printed
  from the copy.
- **What it does on a real run** (the real data layer, read-only on `/SNS/REF_L/IPTS-36119/nexus/REF_L_231801.nxs.h5`):
  - `chopper_lambda_range` = [2.65, 9.45] Å;
  - production band = **(10550.0, 37621.7) µs**;
  - mutant: `TypeError: fromisoformat: argument must be str`, which the slot's
    `except (OSError, KeyError, ValueError, TypeError)` turns into `band = None`.

  The filter then opens on the full span for every run with a chopper log, with no message. Production is right. The
  test cannot see the slip. A `meta["start_time"]` → wrong-key slip survives the same way (test reviewer).

**Fix: behaviour, not code.**
- **Field domain:** the `tof_band` that `_events_for_row` returns, for a run **with** a chopper log. The "no log → full
  span" leg is already pinned, and the reproduction covered the argument-order case.
- **Required behaviour:** for such a run, the dialog's opening `tof_spins` equal
  `lambda_to_tof(chopper_lambda_range(path), start_time)` computed by the **real** data layer from that run.
- **Either** of these satisfies it:
  - E10's metadata leg uses a run with a chopper log, built with the repo's own NeXus builder (`_write_nexus` in
    `tests/test_roi_estimate.py`, which the test reviewer used) and the real functions;
  - **or** the stubs assert their arguments (the path; then the λ band; then a `start_time` string) and return values
    derived from them, so a swap or wrong key reds.
- **Guard:** the test must vary the dimension the fix freezes, the argument order. Prove it with the swap above as a
  battery row (proposed **M43**), which must red alone, N ≥ 1.

## ADVISORY (rule e: outside declared scope, or narrower than the plan's wording with no surviving in-scope one-liner; for the PR body)

- **A-1 (ui):** the log norm after an interaction is unpinned. P2 declared "after open and after a draw", and that is
  met. The probe `set_norm(Normalize(...))` on the XY image inside `_update` after a filter change gives 632 passed.
  Behaviour: V4′ asserts `LogNorm` on both images after its drags.
- **A-2 (ui):** the plain colorbar formatter after a filter change is unpinned. The probe re-norms XY after a filter
  change (U-a's natural fix) and the labels become `$\mathdefault{10^{0}}$` with 632 passed, because the colorbar's
  `update_normal` resets the formatters. When U-a is taken up, run V11′'s `mathtext()` after a filter and an x-range change.
- **A-3 (test; C4):** the test asserts `n_angles == 2` and `RB_Ymin == [130, 140]`, not the plan's `doc.to_dict()`. The
  one row check guards every write, and M40 reds. Cheap: compare `to_dict()`.
- **A-4 (test; C2):** E8 checks that the report lists `  - data_x_range:`, not the new value ("contains the new value").
- **A-5 (test; P5):** E9′ stubs `_events_for_row`, so no run file exists under its `tmp_path`. Its new snapshot cannot see
  a truncated run. E4″ carries P5.
- **A-6 (design):** E4's headline docstring (test_settings_editor.py:3080) still says "every file the slot can write". It
  watches `tmp_path`, so say that.
- **A-7 (design):** the `files()` snapshot helper is duplicated in E4 and E9. Extract one helper so the docstring's
  per-test claims cannot drift. `QTest.qWait(20)` ×2: use `qWait(0)` or `processEvents()`, as the plan wrote, or a
  named helper. The "Could not complete" substring is the only "no problem" check: share the constant.

## Next (the Analyst's; this seat does not grant a v5)

At the human's word: one test change in `test_settings_editor.py` (E10's metadata leg), plus the M43 battery row.
That is all of it. No production line, and every other row already holds. A-1, A-3 and A-7 are one-liners if they
are taken in the same pass.
