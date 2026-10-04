# todo.md — Integrator rejection, `editor-combos` v2 @ 6e97703 (attempt 2 of 3; review gate: ui-aspects harm + test)

**Verdict: REJECT — one harm (U-1, the direct-beam cell), two guards that cannot catch what they name (T-1, T-2).**
Gate cycle: v2 → `review/editor-combos`; **v3 is attempt 3 of 3 — a rejection then escalates.** Not infrastructure. The
human's three PR #36 findings are FIXED on the Q method and background cells and on the scalars (below). What fails is
the **direct-beam** column, which v2's tests never drove with a held name outside the listed folder — the normal state of
a reducer-written file.

## What passed (do not redo)

- **Gate** `pixi run test-reduction` from the subject root, analysis clone 2, 22:02–22:13 EDT: launcher **223 passed**,
  reduction **666 passed**, 699 warnings (= the base), **exit 0** (no teardown abort).
- Scope: two files (`launcher/apps/settings_editor.py`, its test); no ledger-shaped path; base e313f38.
- **The human's findings, through real gestures on a shown tab** (Integrator: ledger `scripts/editor-combos-gestures.py`
  on `IPTS-36119/shared/autoreduce/reduce_settings.json` and `…/shared/reduced/Aug2026/REFL_231105_settings.json`;
  ui-aspects reviewer: the full first-use checklist on both files):
  - one `QTest.mouseClick` opens the Q method / background / direct-beam list in every row incl. surplus; the arrow is
    drawn at rest; Return / Enter / Space / Alt+Down / F2 / F4 open a focused cell; scalars open on one click;
  - re-choosing the shown Q method / background value: document byte-identical, "Changed from the seed" unchanged;
  - C9: the Aug2026 file's list reads `meantheta / constantq / constanttof` and a choice is written `constantq`; the
    declared-spelling file offers `meanTheta / constantQ / constantTOF`; unchosen entries are byte-identical after Save;
  - C10: after a choice focus is on the table (cells) / the scroll area (scalars); Down, Up, wheel, Tab change nothing.
- v1's acceptance still holds at v2 (`scripts/editor-combos-acceptance.py`): wheel over 3 scalar combos × 2 focus states +
  12 table drop-downs → nothing changed; the direct-beam list offers the 39 `*.txt`/`*.dat` of `shared/transmission`.
- C6 active-row rule, C1 on closed cells, C7 no-write-on-display, Add/Remove, 500-row cost (22–24 ms, = v1): pass.
- Test reviewer: the three faithful v1 reverts you were told to kill are killed — v1 double-click (3 failed), v1 C9
  (4), v1 commit-on-close with focus kept (11); C1/C6 rows still red; RED 21/605 and 629 tests confirmed.

## BLOCKING — U-1 (ui-aspects; reproduced by the Integrator): click + Return on a direct-beam cell writes the folder's first file

**Reproduction** (Integrator, 6e97703, `Aug2026/REFL_231105_settings.json`, shown tab): held `DBname[0] =
A2_div10_Cd.txt`, **not** among the 39 names the editor lists (the file carries no `_DBpath_override`; the names it holds
live in `shared/transmission/Aug2026/`). One click on the cell opens the list with **current row 0** (nothing visibly
selected; the edit text still reads `A2_div10_Cd.txt`); `Return` → `DBname[0] = 176.txt`, `changed_vs_seed() =
{'DBname': (['A2_div10_Cd.txt', …], ['176.txt', …])}`. Nothing was chosen; the reduction would use another direct beam
for a real angle. The same on an empty cell and on a new row after Add (ui-aspects, both files). At 7452201 the same
gesture changed nothing — a v2 regression. Mechanism: `settings_editor.py:327` opens the list (`_later(editor,
editor.showPopup)`); `:335` `setCurrentIndex(findText(text))` = -1 for a name not in the list, and the list's current
row is still 0; Return in the list emits `activated(0)` → `_choose` → `:339-340` writes `currentText()`.

**Why it blocks:** demonstrated reachable harm (a silent change of a reduction input by the most natural keyboard
gesture), and it fails the human's finding-2 rule as C9/C10 state it ("choosing the value a cell already holds is the
identity"; "a value changes only by a choice made with the list open").

## BLOCKING — T-1 (test reviewer): the C1 wheel test on table cells can no longer fail

`test_the_wheel_never_changes_a_table_drop_down` (V3/V15) opens the cell, and at v2 the **list is already open**
(auto-shown), so Qt's `QComboBox` ignores the wheel whatever the guard does — the test's own docstring ("the open editor
is a closed drop-down until its list is shown") is now false. Mutant P11 (the two delegate editors' `wheelEvent` back to
`QComboBox.wheelEvent`, class guard kept) → **111 passed**. Reachable under P11: open, Escape (the closed drop-down C10
leaves), one wheel notch → writes in all three columns. Battery F2/F3 now red only via the focus-policy test.
**Fix:** in V3, Escape in `editor.view()` first; assert the list is hidden and the editor still open; wheel; assert
`currentText` and `changed_vs_seed() == {}`. Re-describe F2/F3.

## BLOCKING — T-2 (test reviewer; re-run by the Integrator): v1's commit-on-close with focus kept survives at the direct-beam cell

Mutant R3d — remove only `_CandidatesDelegate`'s `editor.activated.connect(… self._choose(editor))` (v1's shape) →
**111 passed** (Integrator, scratch copy, import verified). V18 drives only the `method_per_run` cell; no test picks a
direct-beam name from the list. Under R3d a pick writes nothing and focus stays on the combo — finding 3 on that column.
**Fix:** a V18 leg for the direct-beam cell: choose a name in the open list → the write, and `focusWidget() is
tab.angle_table`.

## Fix — behaviour, not code (amendment 18)

**Domain.** Three drop-down columns (Q method, background, direct beam) × held state {a listed item; not a listed item
(direct beam: a name outside the folder — the reducer-written norm); empty / `None`; implied (C7)} × gesture {click the
shown item; Return / Enter with no item deliberately moved to; Escape; click elsewhere; Tab; arrow to another item +
Return; click another item}. U-1 covered: direct beam × not-listed and empty × Return. Not covered by v2's tests: any
direct-beam held state outside the fixture folder.

**Required behaviour (C11, the plan words it):** *only a deliberate choice writes.* Opening a list never makes an item
current that the user did not move to; Return / Enter / Escape / Tab / click-away with no deliberate move leave the cell
exactly as held (any column, any held state) — `changed_vs_seed()` unchanged. A held direct-beam name not in the folder
stays shown and kept until the user picks or types another.

**Guard (vary what the fix freezes):** the gesture matrix above on all three columns, with the direct-beam column's held
name **outside** the candidate folder (the Aug2026 shape: no override, names in a subfolder) **and** inside it, plus an
empty cell and a new row after Add; assert `changed_vs_seed()` and the saved text. Mutations: the list's current row
left at 0 → red; `setModelData` writing without a choice → red.

## Decisions for the Analyst / the human (not rejection grounds; v3 must settle them, they are the human's gate)

- **D-a — plain Down.** C8/V16 say Down opens a focused cell; v2 makes Down move to the next row (the APG **grid**
  convention: arrows navigate cells; Enter / F2 / Alt+Down open), disclosed at RED and GREEN and pinned by
  `test_down_in_the_table_moves_to_the_next_row`. Strictly a declared behaviour that fails; argued and disclosed. The
  human cited APG; the grid pattern is the APG pattern for a table. Amend V16 or the code — the plan's call, stated in
  the PR body.
- **D-b — re-choosing an implied value.** In a compact column a cell shows the implied value (C7); re-choosing it
  materialises the list (`['meantheta']` → `['meantheta'] * 3`; `useBS []` → `[True] * 3`) with an entry in "Changed from
  the seed". C7 specifies that; the human's "choose the shown value, nothing changes" reads it as a change. Neither real
  file has a compact column.

## Advisories (non-blocking; carried to the PR body)

UI-aspects: with the direct-beam list open only printable keys reach the name and replace it (Backspace/End/Left go to
the list) — correcting a typo needs Escape first; selecting a row by clicking a drop-down cell opens its list, so the first
click on "Remove angle" only closes it (emulated; X11 replay unverified); a mixed-case column lists both `meanTheta` and
the held `meantheta` (by spec).
Test reviewer: V18 asserts focus "not a combo" rather than `is tab.angle_table` / `is tab.scalar_panel` (P7/P8 survive);
V17's fixture is hand-written, not reducer-written (assertions sound); a mouse click on a pop-up item is never exercised
(same `activated` path); unrowed hunks: case-insensitive completion (P4), `selectAll` (P5), `NoWheelComboBox.keyPressEvent`
editable branch (P2), F4 in `_opens_a_cell` (P3, beyond C8), `_in_file_spelling`'s filter (P6); `_list_shown`'s
completer branch is dead (InlineCompletion).

— Integrator, Claude Opus 5.5
