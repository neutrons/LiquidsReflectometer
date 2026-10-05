# todo.md — Integrator rejection, `roi-popout-dialog` v2 @ e017a97 (attempt 2 of 3; review gate: ui-aspects, design, test)

**Verdict: REJECT — two missing pins, both test-only (plus one docstring scope).** v2 closed every v1 item.
- The gate is green.
- The deployment-shaped acceptance passes again on real data.
- All seven v1 survivors red, with the counts v2 records, and so do M27, M28 and M29.

Two declared behaviours still have no test that can fail; both were already gaps in v1 and were not raised then (the
Integrator's miss). Each was reproduced here in a `git archive` copy of e017a97, each surviving with **628 passed** (the
slug's two test files; the unmutated count is also 628 = 79 + 549, as recorded).

v3 is the **last attempt** (N = 3): fix exactly these two, plus the cheap advisory A-1 (it is in E9, which B-1 touches
anyway). Nothing else is asked.

## What passed (do not redo)

- **Gate:** `pixi run test-reduction` from the subject root, analysis clone 2. Launcher **757 passed**, reduction
  **874 passed**, exit 0, clean. The recorded test_roi_dialog 79, test_settings_editor 549 and launcher 757 reproduce,
  and so do d269b5b's RED counts (3 failed / 76; 2 failed / 547).
- **v1's blockers closed.** Each mutation below is applied to e017a97 and reds, each with v2's recorded count and tests:
  - M17 (minor formatter) → 2;
  - M18 (peak `None` overlay) → 1;
  - M19 (`useBS` `None`) → 2;
  - M20–M22 (X/TOF drags dead) → 1;
  - M23 (Cancel → accept) → 2;
  - M24 (OK unwired) → 2;
  - M25 (log toggle) → 1;
  - M27 (Y-TOF limits) → 1;
  - M28 (a file in the cwd) → 1;
  - M29 (a second QSettings key) → 2.
  The test reviewer also reproduced 16 further battery rows exactly, M13's six names among them.
- **§8.7 acceptance re-run at e017a97** (ledger `scripts/roi-popout-acceptance.py`, launcher main() offscreen, IPTS-36119
  read only): **PASS**, as in v1.
  - Against `REF_L_231801.html` with the template's ROI, the peak, background and x lines sit on the report's lines.
  - The XY and Y-TOF images have the same cells as the report, with report = raw × dead time (1.008–1.012 and
    1.002–1.018).
  - Drag + OK + save + reload hold; sha256 `ee0475b3…` (identical to v1); no traceback.
  - B3′: the Y-TOF image now opens on its data (x limits = the TOF edges; PNG checked). The ui reviewer confirmed the
    limits hold through Estimate, a TOF drag and an X drag, a zoom then a nudge, and Home.
- **Design.** D1 is closed: `git blame -w -M -C` gives 3f74d41 0 lifted lines, and the docstring now says so. §8.4's
  grep prints nothing; there is one F9 name; the two #197 comments are marked inferred; `VIEW_MARGIN`/`TOF_STEP` are
  named; `ruff check` is clean.

## BLOCKING

### B-1 (a, d): a slot that saves the settings document passes every test — #197's F2, the bug E4 exists for

- **Mutation:** after `self.refresh_report()` at the end of `select_roi` (`settings_editor.py`), add
  `self.document.save(Path.cwd() / "roi-settings.json")`.
- **Reproduced here:** the two files give **628 passed**, and the mutant **wrote a 1 484-byte settings JSON**
  (`roi-settings.json`) into the working directory. The test and design reviewers each reproduced it independently;
  E9[Ok] alone writes the file and passes.
- **Cause:** E4′'s stub `no_save` *raises* `AssertionError` inside the slot. The slot is `@guarded` (`except Exception`
  → `report_problem`), so the error goes to the panel: no file is written, and the test passes. The design reviewer
  printed the panel holding `'AssertionError: SettingsDocument.save called'`. The tests that leave `save` unpatched
  (E2, E9, …) do not watch the disk.
- **Rule (d):** v2's new docstring sentence "E4 watches every file and key" (`settings_editor.py:1143-1145`) is
  falsified by this command.
- **Fix (tests; domain = every way the slot can write a settings or data file):**
  - `no_save` **records** calls (`calls.append(...)`); E4′ asserts `calls == []` **and** that the panel reports no
    problem after the slot.
  - E9 snapshots the files under `tmp_path` and the working directory before and after each leg, and asserts nothing
    new or changed (apart from the QSettings store E4′ already allows).
  - A battery row M30, "the slot saves the document", must red alone.
  - Reword the docstring to the domain E4′ actually enumerates (e.g. "E4 watches the working directory, the run's
    folder, the QSettings store and `SettingsDocument.save`") — or keep "every file" only if E9's snapshot makes it
    true.

### B-2 (a): B3's "aspect left to the data" has no test

- **Mutation:** `roi_dialog.py:404`, `"aspect": "auto"` → `"aspect": "equal"`.
- **Reproduced here:** **628 passed**.
- **Effect on run 231801** (200 k events, offscreen, 1000×1000):
  - the Y-TOF axes box goes from **273 × 221 px to 273 × 1 px** (the image is a one-pixel sliver);
  - the XY box goes from 273 × 221 to 186 × 221.
- **The pin is required:** §3 B3 ("aspect left to the data") and L1 ("do **not** force `set_aspect('equal')`") are the
  faithful-axes pin §8.6 hands to V2. V2 asserts origin and extent, not aspect. (ui reviewer, reproduced here.)
- **Fix (test):** V2 asserts `get_aspect() == "auto"` on `xy_axis` and `ytof_axis`. Add a battery row M31,
  "`aspect` equal".

## Cheap advisory to fix in v3 (B-1 touches E9 anyway)

- **A-1 (test):** E9 leaves its 8 s `give_up` timer (`QtCore.QTimer.singleShot(8000, give_up)`,
  `test_settings_editor.py:3289`) running after the test.
  - The test reviewer's probe: a later test that shows a dialog and waits 9 s → the dialog is rejected, giving
    `AssertionError: (False, 0)` after E9[Ok] or E9[Cancel].
  - Harmless in the recorded order, but a shuffled order would flake. The design reviewer once saw "1 failed, 627
    passed" on an unmutated copy (not reproduced in 4 reruns, test name not captured).
  - **Fix:** a `QTimer` object stopped in a `finally` (or parented to the tab).

## Advisories (non-blocking; to the PR body with v1's A2–A8, D2–D8, T2–T3)

### ui

- **U-a.** The colour scale is fixed at open. Widening the x range from [120, 130] to the detector leaves Y-TOF's vmax at
  56 against data at 285 (5 907 of 110 842 non-zero cells saturated). The fix would update `norm.vmax` with `set_data`.
  This is v1's A4, made visible.
- **U-b.** The legend is built once: "peak"/"background" show swatches with no band drawn while the peak is unset, and
  "reduction TOF window" always shows (v1's A5).
- **U-c.** The toolbar's "Customize" changes the scale through `set_yscale`, which bypasses `_plain_log_ticks`: switching
  linear → log there gives `$\mathdefault{10^{-2}}$` ticks. The recursion is font-set dependent and was not reproduced
  (F7). The fix would drop "Customize" from `toolitems`, or re-apply the formatters on a scale change.

### design

- **D-a.** `_guarded`'s docstring cites "the battery's M5 aborts the run". M5 is the editor slot's `@guarded`; the
  dialog's rows are G-guard-values/-estimate/-log (true claim, wrong row).
- **D-b (plan, for the Analyst).** §6′ lists M4 under E9, but E9's docstring says M4 cannot red E9 (Cancel's `reject()`
  restores the opening values; E3 reds it). The Developer's statement looks right; the plan row is stale.

### test

- **T-a.** V18 checks Y limits as `0 <= low <= 133, 149 <= high <= 303`, not "within ±VIEW_MARGIN". A margin of 400 is
  still caught by `test_estimate_brings_its_peak_into_view`.
