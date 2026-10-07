# todo.md — Integrator rejection, `roi-popout-dialog` v5 @ d35e631 (attempt 5 of 5: the human's second extension; review gate: ui-aspects, design, test)

**Verdict: REJECT, tests only, one block: C1's start-time dimension.**

**The human's v6 conditions (A-94, M-46), read against this verdict. All three hold:**
- **(a) No production line.** The fix is in `launcher/tests/test_settings_editor.py` and the ledger battery only.
  `git diff dc28b9a d35e631 -- launcher/apps src` is empty.
- **(b) §8.7 acceptance passing with identical numbers.** PASS, 53 ok, exit 0. Dead-time medians 1.0099 / 1.0086 and edited-save
  sha256 `ee0475b3…`, the same as v1–v4.
- **(c) Only test-side pins of clauses already in §4's table.** The one block is **C1** (§4: "B8 the filter starts at the
  chopper band when the run has one, else the full span"; escalate §0: "C1's clause fixes … **which** band (the data
  layer's composition of λ-range and start time)"). The advisories below are PR-body items, not part of the block.

## What passed (do not redo)

- **Gate:** `pixi run test-reduction`, analysis clone 2 (analysis-node01), 00:02–00:10 EDT, 2026-10-07. Launcher **761 passed**,
  reduction **874 passed**, exit 0, clean.
- **Battery, re-run here.** M30–M44 with M33b, M34b, M37b and M43b: 19 rows, each red with the commit body's count. Run in a `git archive` copy
  of d35e631, ledger table @ 8dacf78, the gate env's python, `__file__` printed; restore 2/2.
  - **M43** is I-54's argument swap, byte for byte: 1 failed. **I-54's block is closed.**
  - M43b: 1 failed. M44: 1 failed.
  - M41 still reds with `qWait(0)`: 2 failed.
- **ui-aspects PASS.**
  - Q2: a `set_norm(Normalize)` on either image in either branch of `_update` reds (two probes + M44).
  - B3/B12 pins intact (+2 lines, 0 removed).
- **design:** production byte-identical to v4. Q4 done as written. Every docstring in the diff holds but one sentence (A-4 below,
  ruled advisory).
- **test:**
  - Q3/C4: `doc.to_dict() == gone` can fail (probe: 1 failed).
  - The `_nexus_builder` load has no vacuity or order risk.
  - 3 baseline runs: 632 each, no flake.

## BLOCKING (rules a, b, d; C1; reproduced here)

### B-1: E10 cannot tell the run's own start time from another date in a same-geometry epoch

- **Mechanism.** `lambda_to_tof` uses `start_time` only to choose the instrument geometry epoch (`source-det-distance`).
  - The epochs are: 15.75 m from 2014-10-10, 15.282 m from 2024-08-26, 15.75 m from 2025-01-01.
  - `_write_nexus`'s default start, 2025-03-04, has the same distance as before 2024-08-26 and as today.
  - So the run's own start time and any other date in those epochs give the same band. This contradicts E10's
    docstring ("the band the data layer composes from that run's file and **start time**") and escalate §0 ("every
    dimension a clause fixes needs a mutation that varies it").
- **Mutation** (`launcher/apps/settings_editor.py:1217`, one line; the start time becomes a literal, standing for a refactor
  that takes "now" or a document date):
  ```
  lambda_to_tof(roi_estimate.chopper_lambda_range(path), meta["start_time"])
  → lambda_to_tof(roi_estimate.chopper_lambda_range(path), "2026-10-07")
  ```
  With it, `test_settings_editor.py` + `test_roi_dialog.py` give **632 passed, exit 0**: the mutant survives.
  Reproduced in a fresh `git archive` copy of d35e631, with `__file__` printed from the copy.
- **What it would do** (the real builder and data layer):

  | Run's start (builder default `chopper_lam=4.25`) | Band, own start time (µs) | Band, literal date (µs) |
  |---|---|---|
  | builder default 2025-03-04 | (9554.7, 23090.5) | (9554.7, 23090.5), identical, so invisible |
  | 2024-10-01 (the 15.282 m epoch) | (9270.8, 22404.4) | (9554.7, 23090.5), ~3 % off, silently |

  Production is right, so no harm is reachable. The test cannot see the slip.

**Fix: behaviour, not code. Complete the dimensions so that v6 is the last round.**

C1's band is composed from **(i)** the run's file → λ-range, **(ii)** the run's start time → geometry epoch, and **(iii)** the
conversion of the band to the opening spins (floor low, ceil high). Each needs a mutation that varies it:

| Dimension | Status at d35e631 | Required in v6 |
|---|---|---|
| (i) λ-range from the run's file | **pinned**: M43 (swap) reds; a literal λ range would change the band from `chopper_lam` | nothing |
| (ii) start time from the run's metadata | **unpinned**: the literal date survives (above); M43b reds only through KeyError | E10's metadata run gets a start time in the **15.282 m epoch** (e.g. `start_time="2024-10-01T12:00:00-04:00"`); new battery row **M43c** (the literal date above) must red alone, N ≥ 1 |
| (iii) floor/ceil to the spins | pinned **only by the no-metadata leg**: E10's run (`chopper_lam=4.6`, default date) has band (10948.105, 24483.944), whose fractions round the same either way (probe P2: `round` survives the metadata leg) | with (ii)'s 2024-10-01 start the band is (10622.79, 23756.42): floor 10622 ≠ round 10623, ceil 23757 ≠ round 23756, so the same test change pins (iii); assert it, and add a `round` battery row (**M43d**) that must red the metadata leg alone |
| the file's identity (which run) | **pinned** by E6 (8 legs) | nothing |

The domain is the `tof_band` that `_events_for_row` returns for a run **with** a chopper log; the no-log leg is already
pinned. No production line.

## ADVISORY (rule e; for the PR body; **not** part of the block or of condition (c))

- **A-1 (test; C4):** `gone` / `before_doc` are shallow (`to_dict()` is `dict(self._config.__dict__)`), so an in-place list
  write would also change the snapshot. The probe `…get("RB_Ymin").__setitem__(0, 999)` in the guard gives 632 passed.
  Production writes rebind (`set`, `set_angle_field`), so nothing is reachable. One-liner: `copy.deepcopy(doc.to_dict())`.
- **A-2 (test; A-i):** `_run_posted_events` (`qWait(0)`) runs zero-delay timers, chains included. A `singleShot(1)` would
  not run (Qt 5.15.15 probe). M41's declared shape is covered. The docstring's "what the slot posted" overstates it.
- **A-3 (test):** E10's expected band uses the same `roi_estimate` functions as the slot, so a mutation inside the data layer
  cannot red E10. That is data-layer scope (#44's tests), outside this slug.
- **A-4 (design; rule d, ruled ADVISORY here):** `_nexus_builder`'s docstring says it loads by path because "launcher/tests is
  not a package beside it". `launcher/tests/__init__.py` is tracked, and `import tests.unit.lr_reduction.test_roi_estimate`
  works from the repo root (reproduced here). Its operative reason (one builder definition) holds, and the sentence
  states no declared behaviour, pin or guarantee. That makes it a prose finding outside declared scope, consistent with
  I-41. Reword or import by name.
- **A-5 (design):** the helper re-executes a test module by path to reach the private `_write_nexus`. Longer term, move it to a
  shared test-support module. One "Could not complete" literal remains at `test_settings_editor.py:767` (outside the slug).
- **A-6 (ui; I-54 A-2, carried):** the plain colorbar formatter after a re-norm is unpinned. Bind it to U-a.
- **Carried from v4's PR-body list:** I-54 A-2 (above), A-4 (E8 checks the report line, not the value), A-5 (E9′ has no run
  file, so it cannot see P5).

## Next (the Analyst's, under the human's pre-authorisation)

Conditions (a), (b) and (c) hold, so v6 is dispatchable without a further line from the human.
- **Scope:** dimension (ii) (a 15.282 m-epoch start time + M43c) and (iii) in the same test (+ M43d). That is one test change and
  one battery row.
- **Optional:** A-1 (`deepcopy`) and A-4 (one docstring sentence) are one-liners. Including them is the Analyst's call; they are
  not required by this verdict.
- **A v7 returns to the human.**
