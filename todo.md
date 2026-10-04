# todo.md — Integrator rejection, `editor-sections` v1 @ b9ad10c (attempt 1 of 3; review gate: test domain only)

**Verdict: REJECT — tests only.** The behaviour passes every reviewer and the deployment-shaped acceptance; four declared
items have no test that fails when they break (two §7 rows survive in faithful forms, one S3 clause and one States case
are unasserted). No production change is asked for. Stacked slug: v2 continues on `feature/editor-sections` (base
`feature/editor-paths-header` @ 0cc96e1, merged forward per the posture before `qa/`). Not infrastructure.

## What passed (do not redo)

- **Gate** `pixi run test-reduction` from the subject root, analysis clone 2, 07:48–07:59 EDT: launcher **628 passed**,
  reduction **725 passed**, 699 warnings (= the base), **exit 0**.
- Scope: the plan's four files; no ledger-shaped path; the branch contains 0cc96e1.
- **ui-aspects reviewer PASS:** S1 order on screen at 1200×900 and 700×900, expanded and collapsed; click, Space,
  Return and keypad Enter toggle, focus stays on the heading; collapsed fields skipped by Tab, the scroll follows focus;
  an edit typed then its heading clicked commits; wheel never toggles or changes a value; three collapsed → a second
  tab on the store shows the same three (raw INI `sections\Dead%20time=false`, keyed by name); a problem in a
  collapsed section is listed and the section stays collapsed; Load while collapsed refreshes the hidden editors; a
  read-only store (dir mode 500): toggle works, nothing raised, `_last_error` None; no leaked top-levels.
- **Test reviewer:** all 18 battery rows reproduce exactly as recorded (bb4dd02); RED 72/1003, GREEN 1075, tip 1077;
  54 V9 cells named on a shown tab, panel cells read `tab.report.toPlainText()`, Tab cells send real Tab; V4 writes the
  raw INI string; the QSettings return-type claims hold (same process: bool; new process: `'false'`/`'true'`).
- **Integrator acceptance (§8.4, §8.6)** — ledger `scripts/editor-sections-acceptance.py`, the launcher's own `main()`
  sequence (`ensure_identity`, `migrate_legacy_settings`, `ReductionInterface`), offscreen, `XDG_CONFIG_HOME` pointed at
  scratch: three headings clicked → stored `false`; **a new process started from a fresh login shell**
  (`env -i … bash -lc`, the profile hook ran) → the same three collapsed, the other seven expanded. The real store is
  `/SNS/users/6ov/.config/ORNL/lr_reduction_new_launcher.conf` (`QSettings().fileName()`, constructed, not written) —
  under `$HOME/.config`, which the login hook's cache wipe does not touch. Drive on `IPTS-36574`: ten sections in S1's
  order; Tab from the DB Browse reaches "Runs and angles" via the Angles table, Add/Remove and the scroll area; Space
  then Return toggle it; `mmpix = "not-a-number"` loaded while Instrument geometry is collapsed → the panel names it,
  the section stays collapsed; all ten collapsed → `changed_vs_seed()` and the panel unchanged.

## BLOCKING — B-1: deleting the import-time call survives (§7 "import-time check removed", rule b; S5 / §5 "import fails naming the group", rule a)

U2 calls `fs._check_section_order` directly, so the module-level call can go. **Reproduction** (Integrator, archive copy of
b9ad10c with a resolution test): delete `field_spec.py:795` `_check_section_order(FIELD_SPEC, SECTION_ORDER)` →
**1078 passed** (1077 + the resolution test). The wiring itself works (test reviewer: a copy with an added field in a new
group raises `ValueError: … missing ['Brand new group'] …` at `import lr_reduction.field_spec`).
**Fix (tests; domain = the one import-time call site):** a test that imports a modified copy of the module source (exec
the edited source, or a subprocess) with a field in a group absent from `SECTION_ORDER`, and expects the `ValueError`
naming that group — so removing the call reds it.

## BLOCKING — B-2: "collapse implemented by clearing values" survives (§7, rule b)

The battery mutated only the disabling half (M9). **Mutant** (test reviewer): in `_Section._show_body`, after
`self.body.setVisible(expanded)`, `if not expanded:` clear every `QLineEdit` under the body → **1077 passed**. Editors write
on `editingFinished`, so a programmatic clear leaves the document intact and Save, `changed_vs_seed()` and the panel stay
green; V7 and the load-collapsed cells reload before expanding. No test does collapse → expand → read the editor.
**Fix (tests; domain = every editor kind under a section — line edits, combos, check boxes):** after a C → E toggle,
assert each of the section's editors shows the document's value (text / current entry / check state), on a section that
holds a line edit and one that holds a combo.

## BLOCKING — B-3: S3 "collapsed means its fields take no space" has no test (rule a)

**Reproduction** (Integrator, archive copy): in `_show_body`, set `setRetainSizeWhenHidden(True)` on the body's size policy
→ **1078 passed**; the test reviewer's probe confirms the mutant keeps the space (the next heading does not move up).
**Fix (tests):** in V2 or V9's C cells, assert the space is gone — the next section's heading `y()` decreases by at least
the body's height on collapse, and returns on expand.

## BLOCKING — B-4: G's "the bad value left in the store" is asserted for `"maybe"` only (rule a; rule d on V5's docstring)

V9's G column uses only `"maybe"`; V5 covers `"maybe"`, `[1]` and `""` for "expanded, nothing raised", and its docstring says
"the bad values stay in the store as they were", but it asserts only the orphaned `"Paths"` key. **Mutants** (test
reviewer): at build, remove a garbage value only when it is a non-string → **1077 passed**; only when it is `""` →
**1077 passed** (removing any unparseable value is caught, 2 failures).
**Fix (tests; domain = the three garbage forms the plan names):** V5 asserts the store still holds each garbage value as
written (`[1]` and `""` included), so V5's docstring is true.

## Advisories (non-blocking; carried to the PR body)

Test: **A-1** keypad Enter is untested (dropping `Key_Enter` survives; it works — probe); **A-2** the §7 heading-text-keyed
rename variant was not written (equivalent today — the heading text is the declared name; a rename probe passes); **A-3**
V10's docstring says the header "has no toggle" but V10 asserts neither placement nor checkability — making the header
`QGroupBox` checkable survives (wrapping it in a `_Section` is caught); `tab.paths_header not in tab.sections.values()` is
vacuous; **A-4** a case-sensitive parse survives (Qt writes lower case — contrived); **A-5** the duplicate-group message
prints the whole order, not the repeated group, though the GREEN body says every case names the group; **A-6** the build
cells' "reachable by Tab" legs are proven through the tab-through cells only.
ui-aspects: **A1** an expanded heading is a checked auto-raise `QToolButton` and draws as a pressed button under Fusion —
can read as "this section is on"; the ~7 px arrow is the only collapse cue; **A2** only the heading text (145 px of a
1144 px row, horizontal policy Fixed) is clickable — Expanding + left-aligned would make the row the target; **A3** no cue
on a collapsed heading that holds a problem (not required); **A4** pre-existing Tab stops between the header and the
first heading (the table walks its cells; the scroll area is a stop with no focus indicator); **A5** a read-only store
forgets the state silently (`status()` AccessError — declared).

— Integrator, Claude Opus 5.5
