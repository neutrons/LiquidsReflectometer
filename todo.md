# todo.md — Integrator rejection, `roi-popout-dialog` v3 @ 21a7367 (attempt 3 of 3 — the retry cap; review gate: ui-aspects, design, test)

**Verdict: REJECT, tests only. The cap is reached, so this goes to the Analyst's `review/roi-popout-dialog-escalate` and the
human; there is no v4 without the human.**

**Production is right.** On real data (IPTS-36119) the deployment-shaped acceptance has passed three times, at v1, v2
and v3, with identical numbers. v3 closed both v2 items, every earlier survivor reds, and the gate is green.

What remains are **five declared clauses with no test that can fail**. All were latent since v1, and none was named in
the v1 or v2 work orders. Each was reproduced here in a `git archive` copy of 21a7367: the slug's two test files give
**628 passed** under each, the same as the unmutated count. The fix is roughly six assertions and one byte of fixture,
with no production change.

## What passed (do not redo)

- **Gate:** `pixi run test-reduction`, analysis clone 2: launcher **757 passed**, reduction **874 passed**, exit 0,
  clean.
- **v2's two blockers are closed:**
  - M30 (the slot saves the document) → 2 failed (E4″[accepted], E9′[Ok]);
  - M31 (aspect `"equal"`) → 1 failed (V2′).

  Both match v3's record.
- **Every earlier survivor still reds:** M17–M29, plus the file-write and second-QSettings-key mutants (12 mutants, this
  seat).
- **26 other battery rows** reproduce exactly (test reviewer).
- **Suite counts** match: 79 / 549 / 757.
- **No flakes:** 5 consecutive runs of the two files gave 628 passed each.
- **E9′'s timers** no longer fire into a later test. The probe gives 4 passed; the positive control at e017a97 failed.
- **Design PASS:** both v3 docstrings hold. Removing each of the three `@_guarded` decorators aborts its test (rc 134).
  E4″ covers save, the panel, files and keys.
- **§8.7 acceptance re-run at 21a7367** (ledger `scripts/roi-popout-acceptance.py`): **PASS**, 53 checks.
  - Every report line and threshold cell matches; dead-time factors 1.008–1.012 (XY) and 1.002–1.018 (Y-TOF).
  - Drag, OK, save and reload all hold; sha256 `ee0475b3…`, identical to v1 and v2. No traceback.

## BLOCKING (rule a; each survives with 628 passed — reproduced here)

### B-1: "Cancel writes nothing" / §5 "authored file, no RBnum — must not: RBnum written", on the real lookup path

- **Mutation:** in `_events_for_row`'s asked-for branch, after
  `self.settings.setValue("roi_nexus_dir", str(path.parent))`, add
  `self.document.set_angle_field(row, "RBnum", 999)` → **628 passed**.
- **Effect** (test reviewer's probe): a chosen file, then the dialog **cancelled**, leaves the document changed:
  `RBnum [None, None, None] → [999, None, None]`. The panel says "No problems found".
- **Why it survives:** E4″ and E6 drive this branch but never compare `doc.to_dict()`. E3 and E9 stub
  `_events_for_row` out.
- **Fix:**
  - E4″ (both legs) and E6 (each leg) assert `doc.to_dict()` is unchanged, except the edited `RB_Ymin` in E4″'s
    accepted leg.
  - Add a battery row: "the lookup writes RBnum".

### B-2: B3's "log norm" has no test

- **Mutation:** `_norm` returns `Normalize(vmin=1, …)` instead of `LogNorm(…)` → **628 passed**. The ui reviewer's
  `Normalize(vmin=0, …)` survives too.
- **Why it survives:** no test reads `image.norm`. This is the same class as v2's aspect.
- **Fix:** V2 asserts `isinstance(image.norm, LogNorm)` for both images.

### B-3: B3's "colorbar" (and B12's "colorbars use a plain-text formatter") have no test; V11′'s colorbar leg is vacuous without them

- **Mutation:** the body of the colorbar loop in `_build_images` replaced with `pass` → **628 passed**. The figure then
  has 5 axes instead of 7.
- **Also survives** (ui reviewer): the loop's two `(image, axis)` pairs swapped, which puts the Y-TOF scale beside the
  XY image.
- **Fix:** V2 asserts:
  - each image has a colorbar whose `ax` is in `figure.axes` (7 axes);
  - the colorbar sits beside its own image: the XY colorbar's `x0` lies between `xy_axis.x1` and `ytof_axis.x0`, and
    the Y-TOF colorbar's lies right of `ytof_axis.x1`.

  Add battery rows for "colorbars removed" and "colorbars swapped".

### B-4: B9's "(labelled 'all angles')" has no test (lowest weight)

- **Mutation:** `"x range, all angles (data_x_range):"` → `"x range (data_x_range):"` → **628 passed**.
- **Why it matters:** it is the one place the dialog tells the scientist an edit there changes every angle.
- **Fix:** E8 (or V7′) asserts the x-range label contains "all angles".

### B-5: E4″ cannot see a truncating write to the chosen run, the exact shape of #197's F2

- **Mutation:** in `_events_for_row`, after `load_event_pixels`, add `open(path, "w").close()` → **628 passed**
  (ui reviewer).
- **Why it survives:** E4″'s run fixture is `run.write_bytes(b"")`, and a 0 → 0-byte truncation leaves the content
  snapshot equal.
- **Fix:** give the fixture at least one byte, or snapshot `(size, mtime_ns)` as plan E9′ wrote.

## Advisories (to the PR body, with v1's A2–A8, D2–D8, T2–T3 and v2's U-a..U-c, D-a, D-b, T-a)

- **Deferred save** (design): the slot posting `QTimer.singleShot(0, lambda: self.document.save(...))` survives. E4″
  and E9′ read the disk before processing posted events. Fix: `QTest.qWait(0)` / `processEvents()` before the asserts.
- **Docstring scope** (design): "E4 and E9 watch that: save, the panel, cwd, the run's folder, QSettings" is true only
  as a union. E9′ stubs no `save`, has no run folder and diffs no keys. "E4 watches all five; E9 the panel and the
  working directory" would be exact.
- **Aspect after an interaction** (ui): `set_aspect("equal")` inside `_update` after a filter change survives, because
  V2′ checks it only at open and after a draw. Fix: repeat the check after one TOF drag in V4′.
- **Types-table legs** (ui/test): V14 has no legs for `BkgROI` `[]`, 3-number and 3-/4-zero entries. They share the raise
  path of the tested legs, so this is coverage of the stated cell, not a survivor.
- **Battery provenance** (test): 21a7367's body cites the ledger script at 87206c4. The reviewers' local ledger clones
  lacked it until fetched; it is on the bus.

## For the escalation (the Integrator's view; the decision is the Analyst's and the human's)

- **Why the loop did not converge.** No round's work order listed every declared clause together with its pin. Each
  gate sampled clauses, so each found new latent gaps: v1 had 7, v2 had 2, v3 has 5, and all of v3's existed at v1.
  This seat should have run a clause-by-clause pin audit at v1. Filed as ledger
  `todo-gate-declared-clause-pin-audit.md`.
- **What the PR would carry if merged as-is.** A dialog that works on real data: acceptance matched the stored
  monitor.sns.gov report three times, and the gate has been green in every version. Five pins would be missing, and a
  later regression in the colour scale, the colorbars, the lookup's writes or the label could land unseen.
- **Smallest completion.** The six assertions and one fixture byte above. That is one commit in the two test files; no
  production line moves.
