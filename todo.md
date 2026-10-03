# todo.md — Integrator rejection, `editor-angle-count` v2 @ 77b5a35 (review gate: design domain + Integrator acceptance; harm clause)

**Verdict: REJECT — one defect family (an edit of a COMPACT list: C-1 silent, C-2 the acceptance criterion, C-3), plus two test-only guard gaps (B-T).** Gate cycle: v2 → `review/editor-angle-count`
(retry **2 of N=3** — see "Decompose?" below). Not infrastructure. v1's findings are fixed (B-1, the stale marks, Q-1,
B-2, B-3) — confirmed through the real tab gestures by the ui-aspects reviewer and on real files here. What fails is a
family of cells the v2 matrix did not specify: **editing one angle of a list the reducer fills or broadcasts itself**
(a compact list: an empty default-if-empty list, an empty or single-entry broadcast list, an unset `None` optional list).

## What passed (do not redo)

- **Gate** `pixi run test-reduction` from the subject root, analysis clone 2, 00:10–00:20 EDT: launcher **163 passed**,
  reduction **421 passed**, 699 warnings (= the base), **exit 0** — matches the GREEN commit.
- Scope: the plan's five files; no ledger-shaped path; base c34c8c5; `todo.md` removed in its own commit.
- v1's B-1 fixed (Integrator, `ledger/scripts/editor-angle-count-acceptance.py edit`, 13 real files): `DBname[0]` edit →
  13/13 PASS (no count line, notes and count unchanged, the reducer's reading of every other entry unchanged, the saved
  file accepted by `NR_Reduction`); `useBS` surplus-row edits PASS where the list already reaches the surplus row.
- ui-aspects reviewer (real gestures on `Aug2026/REFL_231105`): edits in real and surplus rows, marks re-derived, Add at
  index m with surplus shifted and still surplus, saved `RBnum` with no null before the new run, Add→Remove restores the
  file exactly, table order = file order.
- design reviewer: G6 revised holds on three real files (reducer's reading of existing angles unchanged after Add);
  B-2's docstring citations correct.
- **Test reviewer:** v1's B-3 closed (N7 → R7, N8 → R8, 1 each); re-executed N1 24, N2 1, N3 4, N4 4, N5 1, N6 17,
  C1 23, C8 2, C11 6, C12 3, V1–V4 1/2/1/2 → exact counts; R1 generates 30 cases spanning every kind × position × document.

## BLOCKING — the compact-list edit family (three reproductions, one rule missing)

**C-1 (design BL-1; reproduced here) — silent.** An editor-authored file holds `method_per_run: []` (G7 writes it so).
`set_angle_field(0, "method_per_run", "constantQ")` → held `['constantQ']` → the reducer broadcasts a length-1 list
(`nr_reduction_calc.py:77-79`): **every angle** is reduced with `constantq` instead of `meantheta`; `validate() == []`,
"No problems found."; the table shows rows 1 and 2 empty. Integrator's run at 77b5a35: authored → reducer
`['meantheta']*3`; after the one-angle edit → reducer `['constantq']*3`, `validate() == []`. v1 reported this case
("set for some angles but not angles [1, 2]"); v2 removed the report (`settings_document.py:290-292` pads to `i + 1`).

**C-2 (Integrator; the v2 acceptance criterion §8.5 itself).** Plan v2: "re-run on the thirteen real files with one edit
each (a `DBname` at index 0; **a `ThetaShift` at index 0**; a surplus-row edit) — no count line, notes unchanged, saved
file reduces where it reduced before." The twelve Aug2026 files hold no `ThetaShift` (unset → the reducer's default 0).
`set_angle_field(0, "ThetaShift", 0.01)` → `[0.01]` → "Theta shift (deg) (ThetaShift) has 1 entries for N angles" on
**12 of 13** files (`editor-angle-count-acceptance.py edit`: 20 failing edits; 12 are this, the rest a probe artefact
noted below). Harm, reduced end to end (harness, run 221472 + 221473 of IPTS-36119): the real `reduce_settings.json`
with `ThetaShift: []` reduces both runs; after that one edit and Save it holds `[0.01]` and run 221473 raises
`IndexError: list index out of range` (`ThCen = theta_motor + self.config.ThetaShift[i]`, `:413` — the reducer never
length-checks `ThetaShift`). The panel did warn; nothing stopped the save.

**C-3 (design BL-2).** `LambdaMin`/`LambdaMax` unset (`None`, derived from the choppers) in all twelve Aug2026 files. A
value typed into a **surplus** row (`set_angle_field(4, "LambdaMin", 3.0)` on `REFL_231105`) → held `[None]*4 + [3.0]` →
false "Lambda min (LambdaMin) is set for some angles but not angles [0, 1, 2]" while `notes()` calls it surplus; saved, the
reducer takes `LambdaMin[0] = None` and the `:452` mask raises `TypeError: '>=' not supported between 'float' and
'NoneType'`. The v2 matrix declared this cell ("for i ≥ m in any other kind … nothing else changes"). (Pre-existing at
v1 and c34c8c5; it blocks because v2 declared the cell.)

## BLOCKING — B-T (test reviewer; one mutant re-run here): two declared v2 behaviours are guarded only where v1 and v2 agree

- **B-T1 — G6 revised, "a value supplied for a compact list lands at index m".** `add_angle`'s two `m` lines
  (`settings_document.py:236`, `:244`) revert to v1's `n` (captured before the loop) with **every** test green
  (Integrator: the model suite 275 passed under that mutant, import verified; reviewer: 325 incl. view). Every with-value
  test runs where `m == n` (no surplus rows). Reachable (API): on the surplus document,
  `add_angle(DBname="d", method_per_run="constantQ")` would put `constantQ` on surplus row 4, not the new angle at 3;
  same for `ThetaShift=`, `LambdaMin=`. Fix: run T7's with-value leg, R8 and the optional-materialize test on
  `_surplus_document()` and assert the value's index is `m`.
- **B-T2 — the optional-list half of "the some-unset rule looks below m"** (`validate()`, `angles = list(value[:count])`).
  Mutant: the optional rule reads the whole list again → 325 passed. R1's optional cells start from `None`, never a
  surplus-length optional list. Reachable under the mutant: `{…useBS×4, LambdaMin: [2.5]*4}`, clear surplus cell 3 →
  false "LambdaMin is set for some angles but not angles [3]". Fix: an R1 leg starting from a surplus-length optional
  list. (C-3 above is the other face of the same rule.)

## Fix — behaviour, not code (amendment 18)

**Type domain (the one v2 left open).** List kind × **compact** state × edit position:

| kind | compact state(s) | what the reducer reads for every angle while compact |
|---|---|---|
| default-if-empty (`ThetaShift`, `ScaleFactor`, `tof_min`, `tof_max`, `useBS`) | `[]` | its own default (`:99-110`: 0, 1, 0, 100000, 1) |
| broadcast (`method_per_run`) | `[]`, `[x]` | `meanTheta` (`:42-43`) / `x` (`:77-79`) |
| optional (`LambdaMin`, `LambdaMax`) | `None` | derived per angle |

× edit at `i < m` (m = 1, m > 1), at `i ≥ m` (a surplus row). The reproductions covered: broadcast `[]` at i=0 (C-1),
default `[]` at i=0 (C-2), optional `None` at i ≥ m (C-3). Not covered: `[x]` broadcast edits (v2's "expand to m" —
keep), default `[]` at `i ≥ m`, optional `None` at `i < m`, m = 1.

**Required behaviour (G9, new — the plan words it):** *editing angle i of a compact list never changes what the reducer
reads for any other existing angle, and never leaves a file the reducer refuses or crashes on without saying so.*
Concretely the plan must choose, per kind, one of: (a) expand to `m` filling the other angles with **the value the
reducer would have used** (default-if-empty: its default; broadcast: `meanTheta` or `x`) — the file reduces, no line; or
(b) expand to `m` with unset entries and report them by name (G7's wording) — not silent, but the file is not reducible
until they are filled. (a) is what §8.5's criterion requires for `ThetaShift`; (b) is v1's behaviour for broadcast. For
optional `None` (no explicit equivalent of "derived"): (b), and the some-unset rule and the file encoder must look below
`m` only, so a surplus-row value is a surplus value (C-3). Where (a) needs the reducer's defaults, declare them once on
the field (the design reviewer's `Field.reducer_default` advisory), not per call site.

**Guard (vary what the fix freezes).** Extend R1 with the compact states: every kind × every compact state × edit at 0,
m-1, and a surplus row, on m = 1 and m > 1 documents and on an Aug2026-shaped file. Assert the invariant directly:
`json_to_config(saved)` read the reducer's way (after `NR_Reduction(config)`, per angle, truncated to `len(RBnum)`) is
unchanged for every angle ≠ i — or a problem line names the angle. Mutations: C-1's `i + 1` padding for an empty broadcast
→ red; C-3's whole-list some-unset check → red; (a)'s fill value replaced by `None` → red.

**Decompose?** This is the second plan-level gap on `set_angle_field` (v1: edits never enumerated; v2: compact cells
never enumerated). Retry 3 of 3 is the last before escalation. The Integrator's recommendation: v3 states the outcome
**per cell** of the table above (kind × state × position) before code, and the acceptance probe's `edit` mode is
extended with those cells (it already runs the reducer's validation and compares readings). Decomposing "edits of compact
lists" into a child leaf is not recommended: C-1 is a regression from v1 and cannot ship.

## Probe artefact (not a finding)

`editor-angle-count-acceptance.py edit`'s third edit (`useBS` at index m) "notes changed" on eight Aug2026 files: there
`useBS` is exactly m long, so writing a surplus value correctly *creates* a surplus entry and its note (G8 as written).
The probe will pick a list that already reaches the surplus row.

## Advisories (non-blocking; carried to the PR body)

Design reviewer:
- A value typed into a surplus row of an unset default list is held but written `[]` (G7 + G8; declared in the commit). The
  note shown before saving does not say the value will not be written. Same mechanism: a loaded `useBS [null, null,
  null, 1]` is written `[]` (on) by a no-edit save (A2 extended); no real file of that shape among the 13.
- G1's count differs from the reducer's on 9 of the 12 Aug2026 files (`DBname` 4/6/8 entries, `RBnum` 3): Add puts the
  new angle at index 6 and pads `RBnum` with nulls; Add→Remove leaves three trailing nulls. Pre-existing (v1), inside G1
  as declared; the reducer's reading of existing angles is unchanged.
- Plan inconsistency: §3's "Add (G6), by the list's state…" paragraph still carries v1's "the new angle is the row after
  them", contradicting revised G6.
- Wording: an edit at index 0 of an empty list reads "has 1 entries for 3 angles" while the same edit at m-1 names the
  unset angles; "set for some angles but not angles [0, 1, 2]" when no real angle is set.
- Weight unchanged: `Field.fills_itself` (`:356`, `:490`, `add_angle`); `_length_is_allowed` docstring (`:434`); A2
  wording; line numbers shown to users; `settings_document.py` is 572 lines (pre-existing).

UI-aspects reviewer (advise):
- **Selection does not follow Add** (`settings_editor.py:363-366`): with real row 3 selected, Add then Remove — the
  natural undo — deletes real angle 3 (`A6_div10_Cd.txt`, run 229199, its peak and BkgROI) and keeps the empty new
  angle. Not silent (panel problems; "Changed from the seed" lists it). Pre-existing (Remove uses `currentRow()`, a v1
  advisory), now interacting with Add. Highest-priority advisory: select and scroll to the new row.
- Inserting at m shifts surplus rows under a fixed selection (Add then Remove removes the former row 4's entries).
- The panel numbers angles from 0 ("not angles [3]", "at angle {i}") while row headers number from 1 — every Add on a
  surplus file points at the wrong row.
- When every `useBS` entry below the count is unset, Save writes `[]` and drops non-null surplus entries; the note
  "useBS has 3 extra entries" disappears after reload (reducer unaffected).
- `_defining_length` counts trailing `None`: type a `DBname` into surplus row 4, then clear it → row 4 stays an angle with
  three "3 entries for 4 angles" problems; only Remove undoes it (and drops that row's surplus entries).
- An angle created by Add or by a defining edit in a surplus row gets no line for its missing `DBname`/`RBnum`/peak (A3).
- Still open from v1: the mark is header text only and surplus rows start out of view; note wording ("the surplus
  angle", rows not named).

Test reviewer:
- R1 checks the saved file only by round-trip, not against the expected content; checks `notes()` only below m.
- R7 weakened v1's T11 save check (`== [0.1, None, None]` → `saved[1:] == [None, None]`); still caught elsewhere (X3: 4).
- R5's `edit-an-angle` leg cannot detect N5; R5 checks only row 3's label.
- Ragged-short and non-list states have no R1 cell (guarded elsewhere).

— Integrator, Claude Opus 5.5
