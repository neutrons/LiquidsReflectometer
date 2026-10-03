# todo.md — Integrator rejection, `editor-angle-count` v1 @ 09ac86e (review gate: design domain; ui-aspects harm)

**Verdict: REJECT — one blocking defect (B-1, the edit path), two unguarded declared behaviours (B-3), one falsified docstring claim (B-2).** Gate cycle: v1 →
`review/editor-angle-count` (retry 1 of N=3). Not infrastructure. The load / count / notes / Add / file-boundary
behaviours pass on real files; the defect is in the **edit** path (`set_angle_field`), which the plan's tests never
exercised on a file that has surplus rows.

## What passed (do not redo)

- **Gate** `pixi run test-reduction` from the subject root, analysis clone 2, 23:27–23:38 EDT: launcher **161 passed**,
  reduction **379 passed**, 699 warnings (= the base), **exit 0** — matches the GREEN commit.
- **Scope**: the plan's five files; no ledger-shaped path; base = `exp-review` @ c34c8c5 (PR #33 merged).
- **§8.5, deployment-shaped (analysis node; ledger `scripts/editor-real-file-roundtrip.py` sha256 3fddcb8b11144334 and
  the new `scripts/editor-angle-count-acceptance.py`):**
  - The 13 real files of the seed todo (`IPTS-36119/shared/autoreduce/reduce_settings.json` + the twelve
    `shared/reduced/Aug2026/REFL_*_settings.json`): **13/13, no count line at all** (exit 0). Notes as measured, e.g.
    `reduce_settings.json`: rows 4, reduction 3, "Subtract background (useBS) has 1 extra entry beyond the 3 angles the
    reduction uses; …"; `Aug2026/REFL_231105`: rows 7, reduction 3, notes for `method_per_run` (4 extra) and `useBS` (3).
  - **No-edit harm check:** load → save at c34c8c5 vs at 09ac86e (worktree envs) → **identical files for all 13**.
  - **Tab Add (offscreen):** Load the real file → "Add angle" → type DBname/RB_Ymin/RB_Ymax/RBnum in the new row → Save
    → the typed values share one index (4) in every list (G6 as planned; see Q-1).
  - **First end-to-end proof that an editor-authored file reduces:** a three-angle file built from scratch through
    `add_angle` (required columns, λ band, `tof_max`; compact lists left unset) saves `[]` for `ThetaShift`,
    `ScaleFactor`, `tof_min`, `useBS`, `method_per_run`; reducing run 221472 with it through the harness
    (`reduction_scripts`, `go.sh` → `/SNS/REF_L/shared/autoreduce/reduce_REF_L.py`) → no traceback, `check_headers.py`
    PASS, data columns **byte-identical** to the same run reduced with the real `reduce_settings.json`. At c34c8c5 the
    same flow writes `[null, null, null]` lists and the reduction raises `AttributeError: 'NoneType' object has no
    attribute 'lower'` (F5) — RED proven, GREEN proven.
- **Test reviewer** (it also blocks; see B-3): T9 stubs nothing (real `NR_Reduction.__init__`/`_validate_config`); T10 tells
  `[]` from `[None, …]`; the `DBname` all-None pin works; re-executed C1 6, C8 2, C11 4, C11n 2, C12 3 → exact counts.

## BLOCKING — B-1: an ordinary cell edit on a file with surplus rows pads lists into the surplus rows

**Where.** `settings_document.py:279-283` — `set_angle_field` pads the edited list to `n_angles` (the **table's** rows,
surplus included), not to the edited index or the reduction's count; combined with G7's "some unset" check measured
over the whole list (`:335`, `:346-347`), and the view not re-marking rows after an edit
(`settings_editor.py:341-359`, `_on_cell_changed` → `refresh_report()` only).

**Reproductions** (Integrator, at 09ac86e, independently of both reviewers):
- `Aug2026/REFL_231105_settings.json` (rows 7, reduction 3, `validate() == []`, 2 notes) → `set_angle_field(0, "DBname",
  "A2_fixed.txt")` — an edit to a **real** angle → `DBname` becomes `[…3 values…, None, None, None, None]`,
  **reduction_angles 3 → 7**, notes vanish, and four false problems: `RB_Ymin`/`RB_Ymax`/`BkgROI` "has 3 entries for 7
  angles", `useBS` "has 6 entries for 7 angles". Saved and reloaded, it stays that way (ui-aspects reviewer: rows 1–7, no
  surplus mark, same four problems). The reducer still counts by `RBnum` (3), so the file reduces — the panel is wrong.
- `autoreduce/reduce_settings.json` (rows 4, reduction 3, `validate() == []`) → `set_angle_field(0, "ThetaShift", 0.01)`
  → `ThetaShift [0.01, 0, 0, None]` → `validate()` = "Theta shift (deg) (ThetaShift) is set for some angles but not
  angles [3]: …" — index 3 is the surplus row the reducer never reads; saved `[0.01, 0, 0, null]`, which
  `json_to_config` + the reducer read as `[0.01, 0, 0]`. Same for `ScaleFactor` and `method_per_run`; editing `useBS`
  on row 0 pads it 6 → 7 on the Aug2026 file → "set for some angles but not angles [6]".
- The view: after the edit the rows still read "N (surplus)" while the panel counts them as angles (stale mark).

**Why it blocks.** G4: "a reducible file still reads 'No problems found.'" — after one ordinary edit it does not; G1/G2/
G3: the count and the surplus rows are the slug's subject, and the edit path silently converts surplus rows into
angles and persists that into the file. Reachable by the first gesture a scientist makes on the files this slug exists
for (the seed's Aug2026 files). Not reduction harm (the reducer's `RBnum` count is unaffected) — a panel-truth defect
inside the declared scope.

## Fix — behaviour, not code (amendment 18)

**Type domain.** Per-angle lists of four kinds (angle-defining, default-if-empty, broadcast, optional; §3 table) × an
edit at an index **below** the reduction's count `m`, **at/after** `m` (a surplus row), or beyond every list. The
reproductions covered: angle-defining (`DBname`) and default-if-empty (`ThetaShift`, `useBS`) edits below `m` on files
with surplus rows. Not covered: broadcast and optional lists; an edit **in** a surplus row; a file whose lists are ragged
short (not surplus).

**Required behaviour (G8, new):** *an edit changes exactly the entry edited.* Editing index `i` never adds entries at
any other index of that list beyond what is needed to reach `i` (no growth into surplus rows); `reduction_angles` and
the set of surplus rows change only if the edit itself defines a new angle (an angle-defining value at an index ≥ `m`);
the "some unset" rule (G7) and every count look at entries below `m` only; the view re-derives the surplus marks after
every edit (as after load/Add/Remove).

**Guard (vary what the fix freezes).** On the T1 surplus document **and** on an Aug2026-shaped one (`useBS` and
`method_per_run` longer than the defining lists): for each list kind × edit at index 0, at index `m-1`, and inside a
surplus row → assert `validate()`, `notes()`, `reduction_angles`, the surplus marks (view) and the saved file are as
G8 says; then save → reload → same. A guard that edits only `DBname[0]` freezes the kind dimension.

## BLOCKING — B-3 (test reviewer): two declared G6/G7 behaviours have no test that can fail

- **G7, broadcast "some unset → problem".** `validate()` `settings_document.py:346` — `fills_itself = field.default_if_empty
  or field.broadcast_ok`; T11 drives only `ThetaShift`. Mutant `fills_itself = field.default_if_empty` → the model suite
  **233 passed** (Integrator, scratch copy, import verified; reviewer: all 281 incl. view). Reachable from the tab: typing a
  method into one row → `set_angle_field(2, "method_per_run", "constantQ")` → held `[None, None, 'constantQ']`; under the
  mutant `validate() == []` ("No problems found.") while `NR_Reduction(json_to_config(normalize()))` raises
  `AttributeError: 'NoneType' object has no attribute 'lower'` — F5 silenced. Fix: a T11 leg for `method_per_run` (or
  parametrize T11 over the default-if-empty and broadcast kinds).
- **G6, "compact and a value supplied → expanded with unset entries" for an EMPTY default/broadcast list when n > 0.**
  `add_angle`'s `[None] * n` head; mutant `else []` survives (all 281). At GREEN `ThetaShift: []` on three angles +
  `add_angle(DBname="d", ThetaShift=0.1)` → `[None, None, None, 0.1]` (correct); under the mutant `[0.1]` — the new angle's
  value on angle 0, the misalignment G6 exists to prevent. API-only today (the Add button passes no values). Fix: a T8 leg.

## BLOCKING — B-2: a prescriptive docstring claim is false at 09ac86e (verify-prose rule)

`settings_document.py:15-18`: "*Angles grow together.* `add_angle` mutates **every** per-angle field in one operation."
Falsified (Integrator): on `{DBname×3, …, ThetaShift: [], method_per_run: ["meanTheta"], LambdaMin: None}`,
`add_angle(DBname="d")` leaves `ThetaShift []`, `method_per_run ['meanTheta']`, `LambdaMin None` — by design (G6).
Restate it (e.g. "one index per added angle; compact lists stay compact"). The cited reducer lines `:76-78`/`:41-42`
are now `:77-79`/`:42-43` (pre-existing drift).

## Question for the Analyst (not a rejection ground)

**Q-1.** G6 as planned puts a new angle **after** the surplus rows; on the real file the surplus row then becomes an
angle with unset required entries (saved `RBnum [231800, 231801, 231802, null, 999999]`, `DBname […, null,
"DB_new.dat"]`). The panel reports the `method_per_run`/`useBS` partly-unset lines but not the missing `DBname` (A3).
Both reviewers advise a signal here ("rows 4–7 have no DBname"); whether that is this slug's or A3's follow-on is the
plan's call.

## Advisories (non-blocking; carried to the PR body)

Design reviewer:
- **A2's harm, plan-sanctioned:** a `useBS: [null, null, null]` file reduces with subtraction OFF as loaded (falsy entries)
  and is written `[]` (reducer default ON) by a no-edit save; `changed_vs_seed()` is `{}`, only the note says so. 0 of 111
  real settings files (IPTS-35/36/37) hold an all-null or partly-null default list — not demonstrated on real data.
  Suggest the note say "will be written as [] (on)" when the loaded list was all-null.
- The A2 note's "on (1)" comes from `element_type == "bool"`, not a declared reducer default; the note shows
  `(nr_reduction_calc.py:102-103)` to scientists (`:418`) — a line number that will go stale.
- "default_if_empty or broadcast_ok" appears in `validate()`, `_encode_for_file` and `add_angle`; one helper
  (`Field.fills_itself`) would hold it. The two "set for some angles but not angles" strings are near-duplicates.
- `_length_is_allowed`'s docstring (`:424`) still says "differs from n_angles".
- `set_angle_field` expands a single broadcast entry to `n_angles` (`:268`), surplus included (part of B-1's fix).
- An empty broadcast list given a value expands with `None`, not the reducer's `meanTheta` (plan wording; consider the default).
- Measured on real files (the plan's §5 Aug2026 row differs): `reduction_angles` is 3, 4 or 6 across the twelve; the notes
  are mostly `method_per_run` and `useBS` extras — quote measured notes in the PR body.

UI-aspects reviewer (advise):
- The mark is the row-header text "N (surplus)" + a tooltip; no cell colouring. At the default splitter size about 2 rows
  show, so the surplus rows start out of view and the panel note is the first signal; greying the empty defining cells
  would help (A4: the Developer's call).
- Note wording: "removing the surplus angle" (singular) with 4 surplus rows; "1 extra entry … ignores them"; the note
  does not name the rows ("rows 4–7" would tie panel to mark).
- `settings_editor.py:374` (pre-existing): Remove uses `currentRow()`; after `clearSelection()` it still removes row 4
  with nothing visibly selected — use `selectionModel().selectedRows()`.
- Verified OK: no `.destroy()`; header items replaced via `setVerticalHeaderItem` (no leak); marks re-derived on every
  `refresh_angles` (no hidden per-row state); no stale mark across loads; removing surplus rows one by one ends at
  `['1','2','3']`, "No problems found.".

Test reviewer:
- G1 fallback "longest" not distinguished (T3's fallback lists are equal-length; "shortest" survives).
- T13/V4 `== [True, True, True]` would accept ints (trim pins; harmless); T13 checks only `DBname` for "nothing else
  trimmed".
- The new non-list guard in `add_angle` (clicking Add on a file holding `tof_min: 5` raised `TypeError` at c34c8c5) has
  no test; `_row_header`'s tooltip has none.
- T9 does not red under "rule removed from normalize only" (the plan predicted it would); the commit body explains why
  (G6 keeps lists compact, G7 collapses them; either alone keeps the file reducible). Not a count error.

— Integrator, Claude Opus 5.5
