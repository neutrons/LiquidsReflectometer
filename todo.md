# todo.md — Integrator rejection, `editor-ipts-inference` v1 @ 4374f89 (attempt 1 of 3; review gate: test, two findings)

**Verdict: REJECT — tests only.** The behaviour passes every reviewer and the deployment-shaped acceptance on the real tree;
two declared behaviours have no test that fails when they break, and one of them guards the human's own flow. Stacked on
`feature/editor-notes-and-report-spelling` (#43 @ 779f787, which contains #40). Not infrastructure.

## What passed (do not redo)

- **Gate** `pixi run test-reduction` from the subject root, analysis clone 2, 07:05–07:18 EDT: launcher **650 passed**,
  reduction **789 passed**, 699 warnings, **exit 0**.
- **Integrator acceptance (§8.3)** — ledger `scripts/editor-ipts-acceptance.py`, the real tab offscreen, the real
  `/SNS/REF_L` tree read (nothing written): the from-scratch file → `IPTS-36119`, NeXus placeholder
  `/SNS/REF_L/IPTS-36119/nexus`, `experiment_id: "" -> "IPTS-36119"` under Changed, Load 0.09 s; the real
  `IPTS-36574/shared/autoreduce/reduce_settings.json` → its `IPTS-38511` kept with a Note that runs 228313–228316 resolve under
  **IPTS-36970** (**correction to plan §8.3**, which expected IPTS-36574 — measured with `ipts_of_run`, confirmed by the ui
  reviewer); a scratch copy with runs from two IPTSs → the first run's (`IPTS-36970`) and a Note naming `229197` under
  `IPTS-36119`; ~10 ms per run warm; the Load dialog with `IPTS-36119` → starts at `/SNS/REF_L/IPTS-36119/shared`, sidebar
  `shared`, `shared/reduced`, `shared/autoreduce`.
- **ui-aspects PASS:** I1 2a–2d, I3, I4, I5, I6 driven on the shown tab; the user's saved sidebar (QtProject.conf
  `shortcuts`) byte-identical after six rejected and one accepted dialog (4374f89's claim holds); keyboard path entry works;
  latency 0.04–0.33 s (20-run cap), the cap's Note shows; an unmounted root → the "not available here" Note, no exception.
- **design PASS:** one resolution site; the file's value never overridden; a reported value held with no lookup; failed and
  cancelled Loads resolve nothing; inference shown under Changed and clearable. **security (advise):** run numbers reach no
  path unless `int`; no glob engine (scandir + `isfile`); the inferred value is always `check()`-clean.
- **test:** every §7 row killed (16, two split); I1 swaps file↔runs, field↔path, path↔"" killed; I2, I4, I5, I6 frames killed;
  recorded counts RED 39/1115, GREEN 1154, tip 1163 reproduced.

## BLOCKING — B-1: I1's order "runs (2a) before the field (2b)" has no test; a faithful swap survives (rule a/b; reachable)

**Reproduced here** (archive copy of 4374f89): moving the `held = _clean_ipts(field_value)` block above the `if lookup is not
None: for run in distinct:` block in `resolve_experiment_id` → the two modules **1164 passed** (1163 + a resolution test).
On a fabricated tree with runs 229197–229199 under `IPTS-36119`: `SettingsDocument.from_dict({"experiment_id": "", "RBnum":
[229197, 229198, 229199]}).resolve_ipts("IPTS-1", root=tree)` → **`IPTS-1`** under the mutant (`IPTS-36119` at 4374f89). No
test has runs that resolve **and** a header that held another IPTS (U2 row 5's runs have no hit; V1/U2 row 3 have an empty
field). **This is the human's flow:** a previous file leaves `IPTS-1` in the header, the from-scratch file is loaded — the
mutant keeps `IPTS-1` and the reducer cannot find the NeXus files.
**Fix (tests; domain = I1's five steps):** a U2 case "runs resolve, field held `IPTS-1` → the runs' IPTS" and a V2 leg through
Load; a battery row for the runs↔field swap.

## BLOCKING — B-2: "the remembered folder is still written after a Load and a Save" (I5) has no test (rule a)

Deleting `self.settings.setValue("settings_editor_dir", str(Path(path).parent))` survives in `load_settings` (1163 passed)
and in `save_settings` (1163 passed) (test reviewer); no test reads the value back after a Load or a Save — V5 only presets
it. I5's "a remembered folder under the IPTS wins" (cell 16) never arises in use if the write is lost.
**Fix (tests):** assert `settings_editor_dir` after `_load` and after `save_settings`.

## Advisories (non-blocking; carried to the PR body)

Design: **A1** a no-edit Load → Save now writes the inferred IPTS into a file that held `""`/`null`/no key — the declared
intent (I1, V1, §5, I7), shown and clearable; the editor's load→save identity no longer holds for such files (B8 as written is
about booleans); `reduce_from_file` overwrites `experiment_id` anyway — say this in the PR body; **A2** 2b carries a previous
*inference* into an unrelated file (the human's rule (b), literally); **A3** resolution sits in `set_document`, so any
re-adoption re-infers (would break I6 after a clear; no caller today) — resolve in `load_settings` or a `resolve=` flag; a tab
built with an injected document that names runs does a facility lookup at construction; **A4** 571 `IPTS-*` folders now (F6:
286); first contact ~1.35 s cold, then ~13 ms/run; `MAX_RUN_LOOKUPS`'s comment ("~0.1 s a run … ~2 s at most") overstates
the warm cost ~8× and omits the cold touch; `cap=MAX_RUN_LOOKUPS` bound at definition; **A5** Qt's widget dialog writes
`history`/`lastVisited`/`sidebarWidth` (normal); a sidebar entry added during the dialog is dropped; the filter is a global
event callback while modal; the header's Browse is still native (inconsistent within the tab); **A6** `settings_document.py`
1112 lines, `settings_editor.py` 1307 — the ~180-line IPTS block is a natural `ipts_lookup.py`; two "is under the IPTS" rules
(relpath vs commonpath); `ipts_of_run` has no production caller; `_facility_layout` builds a config per call; the
`REF_L_{run}.nxs.h5` spelling now in four places; **A7** an `experiment_id` of only spaces gets no I4 problem and a malformed
Note (two definitions of "empty").
ui-aspects: **A1** nothing says *why* the IPTS appeared (2a/2b/2c) and the Changed line is below the fold at 1000×700 — a
one-line Note, or the IPTS/Changed lines first; **A2** the file-wins Note sits below "No problems found." though a launcher
reduction would not find the runs — add "where they are not"; **A3** the cleared-field Problem could name the IPTS the runs
are under; **A4** the IPTS sidebar replaces the user's bookmarks while open — consider appending; **A5** the lookup is on the
GUI thread with no busy cursor or timeout (a hung mount would freeze Load); **A6** the two-IPTS Note does not say the first
was chosen.
Security (advise): **F1** `_distinct_runs` dedups with `in` on a list over all of `RBnum`, and `notes()` repeats it on every
refresh — quadratic: 30k runs → 3.5–4 s per Load and 1.6 s per refresh (a deliberately oversized file; the editor treats >500
angles as a mistake) — `dict.fromkeys`, test `skipped` against a set; **F2** `IPTS-(\d+)` matches non-ASCII digits (use
`[0-9]`/`re.ASCII`, pass hits through `_clean_ipts`; not reachable on the root-owned real tree); **F3** no time budget on the
lookup (as ui A5); **F4** `load_start_folder` raises `TypeError` on a non-str remembered folder (`QSettings` can return
`None`) — `@guarded` catches it but Load/Save then fail.
Test: rows 1, 2, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14 each lack one column's assertion (header / Changed / Problems) — only
contrived row-specific mutants survive; **row 14's cell looks wrong** (after loading `""` and clearing, nothing shows as
Changed — the value equals the seed) — the Analyst's; U2 row 10's problem check is `any("experiment_id" in …)`; I4 "first"
unpinned (matters past `MAX_REPORTED_PROBLEMS`); `experiment_id: null` with unresolvable runs → `validate() == []` though
`NEXUSpathRB` raises (outside the plan's states; as at the base) — the Analyst's; the cap constant unpinned; the Notes print
IPTS names and runs raw though §2 says `file_spelling` (design/ui).

— Integrator, Claude Opus 5.5
