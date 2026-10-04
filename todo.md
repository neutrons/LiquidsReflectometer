# todo.md — Integrator rejection, `editor-defaults-and-theta` v1 @ 060834d (attempt 1 of 3; review gate: test domain + a D3 defect)

**Verdict: REJECT — one small code defect (D3 after a second load) and the plan's declared Operation × state cells
left untested.** Stacked slug: v2 continues on `feature/editor-defaults-and-theta` (base `feature/launcher-test-teardown`
@ ed663f7, merged forward per the posture before `qa/`). Not infrastructure.

## What passed (do not redo)

- **Gate** `pixi run test-reduction` from the subject root, analysis clone 2, 03:10–03:23 EDT: launcher **353 passed**,
  reduction **712 passed**, 699 warnings (= the base), **exit 0**.
- Scope: the plan's five files; no ledger-shaped path; the branch contains ed663f7.
- **Design reviewer PASS.** D7 table over every loaded state (absent / falsy forms / `True` / `"detector_angle"` any case /
  `"sample_angle"` any case / rejected strings / `1` / lists): the reducer reads the saved file exactly as the source,
  incl. a real `useCalcTheta: true` file (`IPTS-36517/…/reduce_settings_air.json`). D1/D2/D3/D4 as declared; the
  import-time label check survives `python -O` and raises on a missing label.
- **Numerical-diagnostics reviewer PASS.** Two-source audit: `nr_reduction_config.py:98-99` rectangular/0.8; editor
  start gaussian/1.0 on `Field.editor_start` (`field_spec.py:684-686`); four real files with sha256 — none read
  differently by the reducer after load → save; a keys-absent file saves rectangular/0.8. `DetSigma` is the kernel's
  standard deviation in mm in every branch (`nr_tools.py:387-402`, `nr_reduction_calc.py:681` `_calc_detector_convolution`,
  constantq at `:660-661`). Reducer files hash-identical; a `calc_beam_on_detector` fingerprint unchanged.
- **Test reviewer:** every §7 row red (21/21), the faithful reverts (a)–(e) red; Types-and-states rows assert type and
  value; U7 asserts the reducer's own reading; RED 46/731 and 789 confirmed.
- **Integrator acceptance (§8.4), the real tab** (ledger `scripts/editor-defaults-acceptance.py`): fresh tab holds and
  shows gaussian/1.0, seed diff `{}`; `IPTS-36574` (keys absent) saves rectangular/0.8; `IPTS-36517` (gaussian/0.65)
  byte-identical on the three keys; `IPTS-38016` (`useCalcTheta: true`) → `"detector_angle"`, reducer reads
  `detector_angle` from both; theta combo: items `False / True / trust sample angle`, one click opens, re-choose shown →
  no write, a choice → canonical `sample_angle`, focus released, Down/Up/wheel afterwards → no change.

## BLOCKING — D-1: after a second Load the theta control keeps a stale entry (D3 "exactly three entries")

**Reproduction** (Integrator, 060834d): load a file holding `useCalcTheta: "true"` (a rejected spelling — shown as an
extra entry, reported: D6, correct), then Load `IPTS-36574/shared/autoreduce/reduce_settings.json` → held `False`,
shown `False`, but items `['False', 'True', 'trust sample angle', 'true']` — four entries; choosing the stale one would
write a value no loaded file holds. **Fix (behaviour):** the extra entry exists only while the held value needs it;
every refresh (load, set_document, a choice that replaces the raw value) leaves exactly the three entries plus, at most,
the raw entry of the value **now** held.

## BLOCKING — T-1: the plan's Operation × state cells (§3: "every cell is a required outcome") — 17 have no test

Test reviewer's audit, confirmed against the test module (probe: all 54 gesture × state combinations hold in production —
this is missing guards, not a defect, apart from D-1):

| Operation | False | detector_angle | sample_angle | rejected |
|---|---|---|---|---|
| Tab out with no deliberate move | — | — | — | — |
| Click-away with no deliberate move | — | — | — | — |
| Wheel on the closed combo | fresh tab only | — | — | — |
| Up / Down on the closed combo | — | — | — | — |
| Load a second file omitting the key | trivial | — | covered (V5) | — (D-1) |

**Verify-prose:** the docstring of `test_each_gesture…` says the arrow keys on the closed combo are covered for this
field; only V18 does it, for `DetResFn`. A claim that fails its own check is blocking. **Fix:** parametrize the gesture
test over the four held states × {Tab out, click-away, wheel, Up, Down} on a shown tab, and add load-second-file legs from
each state asserting the held value, the shown entry **and the item list** (D-1's guard); correct the docstring.

## Advisories (non-blocking; carried to the PR body)

Test reviewer: **A1** V6 cannot tell which entry is selected when the raw value spells a label (`"True"`, `"False"`,
`"trust sample angle"` as strings): a mutant selecting by `findText` survives — assert `currentIndex() == 3` or
`type(currentData()) is type(loaded)` and add a held-`"True"`-string column; **A2** removing the import-time loop
(`field_spec.py`, after `EDITOR_START_NAMES`) survives (U6 still catches the §5 scenario), and a check ignoring the
`False` label survives; **A4** theta combo: Enter / Alt+Down open and reach-by-Tab untested (they work), `DetSigma` has no
first-use QTest, `test_the_theta_control_is_a_choice_not_a_checkbox` uses `setCurrentText`; **A6** no view test loads a
list (shown as itself, reported, no exception); `Field.value_for` has no production caller.
Design reviewer: the falsy-forms docstring overstates ("no reader tells those apart" — `web_report.py:432` and the `.dat`
`# Config:` header print them); `_check_choice_labels` compares by set equality (`0 == False`) while `label_for` matches by
type; stale docstrings ("asserted at import", "(a copy)"); `True → detector_angle` rides `CALC_THETA_CHOICES[0]` (name it);
the theta combo writes on `currentTextChanged` — `activated`/`currentIndexChanged` is the robust hook; file sizes.
Numerical reviewer (for `library-defaults-detres` and the scientists): no real file uses gaussian/1.0 (gaussian files hold
0.8 or 0.65); since σ is the standard deviation in both shapes, the requested 1.0 changes the shape **and** widens the
kernel by 25 % (the convolved profile ~19 % in a test geometry) — explicitly requested, worth confirming at sign-off;
`DetSigma`'s help text says "Width" — "standard deviation (mm)" would prevent full-width entries; the §8.3 statement could
name `_calc_detector_convolution` and the gaussian half-width clip.

— Integrator, Claude Opus 5.5
