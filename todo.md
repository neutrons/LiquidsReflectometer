# todo.md — Integrator rejection, `launcher-test-teardown` v1 @ 8585197 (attempt 1 of 3; review gate: test domain)

**Verdict: REJECT — two declared clauses with no test that can fail. Test-only fixes; the fix itself is sound.**
Gate cycle: v1 → `review/launcher-test-teardown`. Not infrastructure.

## What passed (do not redo)

- **Gate** `pixi run test-reduction` from the subject root, analysis clone 2, 01:02–01:14 EDT: launcher **321 passed**
  (116 s — the subprocess self-tests), reduction **666 passed**, 699 warnings (= the base), **exit 0**.
- Scope: `launcher/tests/conftest.py`, `launcher/tests/test_harness.py` only; nothing under `src/` or `launcher/apps/`.
- **The defect is gone, independently, on this host (the deploy-target class):** a test that opens a `QCompleter` pop-up
  under the real `isolated_qapp`, widgets kept alive past the test body, in 12 fresh `pytest` processes — base b86237b's
  conftest **7 of 12 failed** (abort / "pure virtual method called"); 8585197's conftest **0 of 12, and 0 of 12 again**.
  The Advisor's standalone probe here: `fixture` 7 of 10 failed, `skip-popups` 0 of 10.
- **Test reviewer:** the faithful revert (b86237b's teardown, byte for byte) is caught every time — T1 red in 6 of 6
  invocations (72 subprocess runs, every failure exit −6); P1–P5, F1–F3, F5, F6 red as recorded; counts exact.
- **Design reviewer (advisory):** the window-type allowlist is exact for Qt 5.15 (`windowType()` masks the hint bits;
  Sheet 0x5 / Drawer 0x7 matched exactly; Tool 0xb = Popup|Dialog); the gdb site reproduced
  (`sipQCompleter::~sipQCompleter` ← `deleteChildren` ← `~QLineEdit` … under DeferredDelete).

## BLOCKING — B1: "hidden or visible (both handled)" (§3 types and states) — the hidden half has no test

Every widget that reaches `_drain_test_windows` in the tests is `show()`n (T3, the freed-during-drain module,
`test_the_drain_frees_the_test_windows_now_and_leaves_the_owned_ones`); T4's hidden widgets test `_is_test_window` only,
not the loop. Mutant: `conftest.py:58` `if _is_test_window(widget):` → `if widget.isVisible() and _is_test_window(widget):`
→ **25 passed** (Integrator, scratch copy, the copy's conftest verified in force; the reviewer: 27 incl. the subprocess
tests). **Fix:** in the drain test leave one or more of the `freed` windows un-shown and keep the `isdeleted` assertion.

## BLOCKING — B2: H2 "identity restored first" has no assertion

Deleting the three `setOrganizationName/Domain/ApplicationName(prev_*)` lines in `isolated_qapp` → 25 passed on the
harness file; 272 passed with `test_settings_persistence.py` and `test_settings_editor.py` (test reviewer). H2 is declared
as delivered by GREEN; the gap pre-dates this slug, but the plan states the clause. **Fix (pins "first" too):** in
`_LEFT_OPEN_WINDOWS_MODULE`, let `_Window.closeEvent` record `QCoreApplication.organizationName()` and assert it does not
start with `test-org-` (restored before the drain closed the window); optionally a fixture-free test that the identity
equals the module-import value afterwards.

## Advisories (non-blocking; carried to the PR body)

Test reviewer:
- Power: `_TEARDOWN_RUNS = 12`, fails on any nonzero exit and requires `1 passed` per run. At the measured ~0.48 abort
  rate the false-pass probability per invocation is 0.52^12 ≈ 4e-4; at the worst observed 4/12 it is 0.67^12 ≈ 0.8 % —
  the `(1/2)**12` comment / plan A1 overstate the margin slightly; N = 16–20 would bound it under 1e-4.
- Widget lifetime correct (module-level `_LEFT_OPEN`); crashed runs abort inside `_drain_test_windows` — red for the
  right reason. No leak across tests (`topLevelWidgets() == []` after the harness and editor files).
- F4 (processEvents removed) recorded equivalent — true for H2's declared outcomes, not strictly (events a `closeEvent`
  posts to a survivor would go undelivered).
- T2 (QMenu) is a pin, and its menu is parented. `--basetemp` in `_run_inner_pytest` has no row (concurrent inner runs
  would otherwise share pytest's rotating basetemp — flaky, not falsely green).

Design reviewer:
- The rule now leaves a parentless `QMenu` a test keeps alive visible and active into the next test **when the
  `QApplication` is shared** (repro in the reviewer's scratch; with a fresh application per test both base and tip free
  it). No test depends on it today.
- `conftest.py:42-43` comment too pessimistic ("left to die with the process") — a Python-owned parentless popup dies
  with its last reference; `:37-38` names `_drain_test_windows` as the crash site, which did not exist when measured.
- Parented windows (`QDialog(owner)`) are no longer closed — deleted with their parent, no `closeEvent`; leak if the
  parent is a left-alone popup. Not named in the comment or H1.
- `_repeat_inner_pytest` runs in parallel and is slower than serial here (12 parallel 29–35 s vs serial 13.5 s);
  serial gave a lower crash rate (weaker bound) — keep parallel and say why, or serial with larger N. A run killed at
  120 s on a slow host reads as a failure unrelated to the crash; T2 (a pin) could use fewer runs.
- `test_harness.py` is 442 lines; the teardown block would split into `test_harness_teardown.py`.
- `editor-combos`' inline completer (`settings_editor.py:358-384`) no longer has a test-side reason; whether to move to
  `PopupCompletion` is a UX decision for a successor slug (with its own subprocess repeat test); its docstring's past
  tense should be updated if it stays.

— Integrator, Claude Opus 5.5
