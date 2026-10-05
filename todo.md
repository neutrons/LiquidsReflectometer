# todo.md — Integrator rejection, `roi-popout-data` v1 @ 64766a8 (attempt 1 of 3; review gate: numerical + test)

**Verdict: REJECT — the numbers are right; one declared equality is false and five declared behaviours have no test that
can fail.** No module arithmetic needs to change. Not stacked (base `exp-review` @ 5a8742d + PR #31's `feature/roi-estimate`
@ 655a67d). Not infrastructure.

## What passed (do not redo)

- **Gate** `pixi run test-reduction` from the subject root, analysis clone 2, 05:19–05:29 EDT: launcher **558 passed**,
  reduction **792 passed**, 699 warnings (= the base), **exit 0**; `git status` clean afterwards.
- Scope (§8.2): `git diff 5a8742d..64766a8` = exactly the four files of §4; no `plans/`, no `todo.md`; `9bc0049` is R100.
- **Numerical reviewer, point-wise:** B6 `background_bands` vs `NR_Reduction._background_roi_sorter` — 754 784 grid cases
  (exhaustive n_y=8; an edge set at n_y=304; 300 000 random), **0 mismatches**, refusal exactly where the sorter returns
  `None`; the real `background_subtract` agrees with the pop-out's bands in 3 000 / 3 000. B7 `default_bkg_roi` — 231 800
  cases: four ascending ints, no 0, inside the detector, bands of `width` rows `gap` from the peak, refused (never clamped)
  exactly when there is no room; the sorter returns it unchanged. B4 — exactly one cell differs from `RefRoi`, the last-edge
  event, on all 8 runs. B5 — `profile_y/∑pc` vs `counts_vs_y` ≤ 3.0e-16; every marginal equal. B1 — packing and TOF
  bit-equal to `get_y_tof`; stride sampling is every k-th on-detector event; detector shape (256, 304) from the database =
  Mantid's IDF. I1 — `lowres` default from `detector_shape`.
- **Test reviewer:** the battery's 49 rows reproduce **verbatim** (names and `-> N failed`) against the file and 5f3d09b's
  body; `restored: OK`, 0 ANCHOR MISS, `git status` empty; every commit-body count reproduced at its own tree; I1–I3, I5–I7,
  B1, B3, B4, B7, B8 pinned; I6 verified over every tracked file's mtime.

## BLOCKING — N-1: "the web report's XY array, cell for cell" is false on most real runs (rule d)

`roi_estimate.py` (`xy_image` docstring, ~:245) and plan B2/F4. **Reproduced here** (`web_report.py:579-583`'s exact
`Integration` of `LoadEventNexus`, vs `xy_image(load_event_pixels(f))`): `201288` equal; **`179932` 2 cells differ** —
(72, 14), (191, 217), sums 300 088 vs 300 086; **`201282` 1 cell** — (216, 202), 6 050 vs 6 049. Numerical reviewer: 45 of
the 63 fixture runs and the IPTS spot check `REF_L_235251` differ, by 1–2 counts in 1–2 cells; the missing events are always
the run's minimum- and/or maximum-TOF event, which Mantid's default `Integration` range (min/max TOF, apparently rounded)
excludes (mechanism inferred from the printed ranges). T1 tests only `201288`, one of the 18 agreeing runs. The module's
arrays are correct (its "sum equals the events held" holds); the claim is not.
**Fix (wording + test; domain = B2's reference claim):** restate B2/F4 (docstring, plan, PR body) as the F4 Y-TOF note
already is — equal except the event(s) at the run's extreme TOFs, which the report's `Integration` drops when its default
range excludes them; add a T1 leg on `179932` asserting the difference is exactly those events' pixels.

## BLOCKING — T-1…T-5: declared behaviours with no test that can fail (rule a; K1 and K4 also rule d)

Test reviewer; K1 and K4 **reproduced here** in an archive copy (baseline 3 environment failures; each mutant adds 0).
- **K1 — B6 "values and types":** `ordered.tolist()` → `ordered.astype(float).tolist()` survives (lists compare `136 ==
  136.0`); falsifies d98ac59's body, the test docstring and the module's "values and types". **Fix:** also compare
  `[type(v) for v in got]`.
- **K2 — B9 called-wrong is `ValueError`, not `CannotEstimateError`:** `CannotEstimateError` subclasses `ValueError`, every
  called-wrong test uses `pytest.raises(ValueError)`; raising `CannotEstimateError` in `_inclusive` reds only I5's META
  anchor check. **Fix:** assert `type(exc) is ValueError` (or `not isinstance(exc, CannotEstimateError)`) for the called-wrong
  and `background_bands` refusals.
- **K3 — §5 / types table "no `bank1_events` → h5py `KeyError`, uncaught":** returning empty arrays instead survives (72
  passed). **Fix:** a builder file without the group, `pytest.raises(KeyError)`.
- **K4 — `bkg_roi is None` → "no background is set for this angle":** deleting the guard (`roi_estimate.py:489-490`)
  survives; the `ndim` refusal fires with the wrong reason and the test matches only "background"; falsifies 5f3d09b's "a
  row for each guard GREEN added" and the README's guard claim. **Fix:** match "no background is set"; add a battery row.
- **K5 — ranges accept `low == high`:** `high < low` → `high <= low` reds only I5's META check (a one-pixel range,
  `profile_y(ev, (120, 120))` = 1174 counts today, would become a `ValueError`). **Fix:** one equal-bounds case per range kind.

## Advisories (non-blocking; carried to the PR body)

Numerical: `default_bkg_roi` silently truncates a fractional `peak_range` and returns non-int bounds for a fractional
`gap`/`width` (docstring promises ints); accepts `gap=True`; `background_bands` refuses entries longer than four that the
reducer would read (conservative; say so); `y_min`/`y_max` non-finite not validated (as the sorter); `counts_vs_y` strides
before dropping off-detector ids, `load_event_pixels` after (0 off-detector ids in 63 fixtures — unreachable here);
`tof_edges` adds a whole last bin where Mantid's `Rebin` narrows it (335 vs 334 edges on 201288 — the pop-out's last TOF
column will not match the report's; nothing claims it does).
Test: **A1** I8's subprocess imports whichever `lr_reduction` the interpreter finds (the editable install), so from another
checkout it tests the wrong tree — pass `PYTHONPATH` from `re_mod.__file__` or assert the child's `__file__`; **A2** §7 says M5
reds T7, row 21 reds only T5 (note it); **A3** the 3-zeros `bkg_roi` case is never named, refusal tests match "background"
never "sentinel"; **A4** a one-row peak `(150, 150)` in `default_bkg_roi` unpinned; **A5** the empty-run cell checks shapes not
all-zero values, `profile_tof` never run on empty events; **A6** M3 uses a percentile window as the crop stand-in (fine); **A7**
"well under a second" and "a second file read per redraw" are unmeasured — mark inferred; **A8** I6's in-test mtime check
covers the module only.

— Integrator, Claude Opus 5.5
