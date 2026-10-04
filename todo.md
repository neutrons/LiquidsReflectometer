# todo.md — Integrator rejection, `editor-paths-header` v1 @ 1594685 (attempt 1 of 3; review gate: test domain + a reachable-harm UI finding)

**Verdict: REJECT — the report panel is never asserted after a header edit (three faithful mutants survive, and a test
docstring plus the RED commit claim the opposite), and Browse → Choose on the folder the dialog opens at freezes the
derived path into the file, which reads another experiment's direct beams when the file is reused.** Stacked slug: v2
continues on `feature/editor-paths-header` (base `feature/editor-defaults-and-theta` @ 3c4ec39, merged forward per the
posture before `qa/`). Not infrastructure.

## What passed (do not redo)

- **Gate** `pixi run test-reduction` from the subject root, analysis clone 2, 05:38–05:49 EDT: launcher **519 passed**,
  reduction **720 passed**, 699 warnings (= the base), **exit 0**.
- Scope: the plan's five files; no ledger-shaped path; the branch contains 3c4ec39.
- **ui-aspects reviewer PASS on every gesture driven** (QTest on a shown tab): one click reaches each field; Tab order
  IPTS → NeXus → Browse → DB → Browse → table; Browse by click and Space; focus through derived / loaded / malformed
  fields, Return in untouched fields, retyping the held IPTS, typing then backspacing → no write; placeholders follow on
  Tab, Return and click-away; placeholders greyed, `text() == ""`, Ctrl+A Ctrl+C copies `''`; empty IPTS says "set an
  IPTS or type a path", a reported IPTS says it is not a folder name; Ctrl+A Delete Tab → `None`; clear IPTS on real
  IPTS-36574 → `''`, no error, `candidates() == ([], 0)`, Save writes `""`; Load A (override) then real B → nothing of A;
  a typed DB folder → one click on a DBname cell lists it; 700 px layout; lifecycle (no new top-levels, no `destroy()`).
- **Test reviewer:** the 21-row battery — every §7 row and frame row red except F4 (equivalent, agreed); the 40 cells × 2
  controls are named by V12 with value **and type**; V6 runs through the real focus path on an active window; V10 checks
  `== ""` and `type is str`; RED 93/867, GREEN 960, tip 963 and the learnings' "M2 reddened 73" reproduced.
- **Integrator acceptance (§8.4), the real tab, QTest gestures** (ledger `scripts/editor-paths-header-acceptance.py`) on
  IPTS-36119 (frozen overrides), 36517 and 38016 (no `experiment_id`), 36574 (`experiment_id` IPTS-38511): Load shows the
  file's state, seed diff `{}`, Load → Save leaves `experiment_id` and both overrides as the source holds them; a
  no-typing focus pass writes nothing; typing `38016` + Tab stores `IPTS-38016`, unset overrides stay `None` with
  following placeholders, set overrides keep their text; one click on a direct-beam cell lists IPTS-38016's 7 files; a
  typed folder is the override, listed, and survives an IPTS change; clearing it returns to derived; clearing the IPTS
  stores `""` with no slot problem and Save writes `""` / `null`.

## BLOCKING — B-1: no test reads the report panel after a header edit (rules a, b, d)

**Declared:** §3's D × "type an IPTS" cell — "one 'Changed' line (`experiment_id`)" — is the panel's "Changed from the
seed" section (`refresh_report` → `tab.report`). V12's docstring says it checks "whether the report names the field";
the RED commit body (141c1ef) says each cell asserts "… "Changed from the seed", and the report". V12 calls only
`document.changed_vs_seed()` and `document.validate()`; no new test reads `tab.report.toPlainText()`.

**Faithful mutants that survive** (test reviewer; N2+N3 re-run by the Integrator in an archive copy with a resolution
test: **964 passed**): delete `self.refresh_report()` from `_on_ipts_edited` (N2), from `_on_path_edited` (N3), from
`_browse_path` (N4).

**Reproduction of the harm** (Integrator, N2+N3 vs 1594685): load `{"experiment_id": "IPTS-1", "_DBpath_override": 5}`,
type `/data/typed` + Return in the direct-beam path, `36119` + Return in the IPTS. `validate()` is `[]` in both. The
mutant's panel still reads `Problems: - Direct-beam path (_DBpath_override): expected text, got int 5` with no Changed
lines; 1594685's reads `No problems found.` and both Changed lines.

**Fix (tests; domain = the header's three write slots — `_on_ipts_edited`, `_on_path_edited`, `_browse_path` — each a
scalar write):** every header write path asserts the panel text itself: the "Changed from the seed" line naming the edited
field after type-IPTS, type-path, clear-path and Browse, and that an X state's problem line is gone from the panel after
type-path, clear-path and Browse. The guard must vary the slot (each of the three, so N2, N3 and N4 each red) and correct
V12's docstring to what it asserts.

## BLOCKING — U-1: Browse → Choose on the folder the dialog opens at writes the derived path (harm clause)

**Q5** (`[human, 2026-10-02]`): "display only; write an override only when the user edits the field, **because written
overrides freeze absolute paths**." `_browse_path` opens the dialog at `editor.text() or self.document.derived_path(name)`;
pressing Choose without navigating — the dialog's own default, and what a scientist does who clicks Browse to *look* —
writes that derived folder as the override and adds a Changed line. This is the no-op re-choose gesture (choose the value
already shown), the class the human rejected `editor-combos` v1 for (PR #36, finding 2).

**Reproduction** (Integrator, ledger `scripts/editor-paths-browse-noop-probe.py`, real
`/SNS/REF_L/IPTS-36574/shared/autoreduce/reduce_settings.json`, read-only): one click on the direct-beam Browse, dialog
returns its start folder → override `None` → `/SNS/REF_L/IPTS-38511/shared/transmission`; Save writes it. Reusing the
saved file for the next experiment (`reduce_from_file` replaces `experiment_id` with the run's,
`new_reduction_from_file.py:62`; this very file was copied in from IPTS-38511): the **source** resolves `DBpath` to
`/SNS/REF_L/IPTS-38016/shared/transmission`, the **saved** file to `/SNS/REF_L/IPTS-38511/shared/transmission` — another
experiment's direct beams, silently.

**Fix (behaviour; domain = the two path overrides, `_NEXUSpathRB_override` and `_DBpath_override`, each `str | None`, plus
a malformed non-string in X; the reproduction covered `_DBpath_override` in state D):** a Browse whose chosen folder is
the folder the reduction derives now (`derived_path(name)`, compared as paths — a trailing separator is not a different
folder) writes no override in state D: held `None`, placeholder intact, `changed_vs_seed()` unchanged, panel unchanged.
What the same choice does from S / X (the override returns to `None`, i.e. derived, is the natural reading of "a derived
path is never written") is the plan's to state; whichever it states is a cell with a test. The guard must vary the field
(both paths) and the state (D and S at least), and drive the Browse button on the shown tab.

## Advisories (non-blocking; carried to the PR body)

ui-aspects: **A2** whitespace is stored as a path — `"   "` or `"/some/dir "` is held as typed, `validate()` silent, and a
whitespace-only value hides the placeholder (the field looks empty but holds an override; P4's `Path("")` reason applies
to `Path("   ")`): `widget.text().strip() or None` would close it; **A3** a long path shows only its tail after leaving the
field and the tooltip is fixed help text — put the value / derived folder in the tooltip or `setCursorPosition(0)`;
**A4** the header takes ~140 px; at 1200×900 the Angles table shows 2 rows (for `editor-sections`).
Test: **A1** M2 does not red V3/V6 and M4 does not red V3 as §7 names (each row reds elsewhere); no test clears then
Saves to pin the `null` of the `""`-after-clear row; **A2** "Load a second file … with an override" is untested (only
"none"); **A3** P6's IPTS text after Load is asserted only indirectly; **A4** V12's type-path / Browse cells skip the
placeholder; **A5** V10's "no TypeError" rests on the value check (`@guarded` swallows; no `_last_error` assertion);
**A6** V6's own legs pass at the base — the no-typing gate is caught by V12's X × focus-through (F3) and F2.
Integrator: the header shows the file's `experiment_id`; `reduce_from_file` replaces it with the run's IPTS, so for a
file copied between experiments (IPTS-36574's holds IPTS-38511) the header's derived paths are the copied IPTS's — true
to the file, worth a line in the help text. Real autoreduce files without `experiment_id` (IPTS-36517, 38016) show "set
an IPTS or type a path".

— Integrator, Claude Opus 5.5
