# todo.md — Integrator rejection, `roi-popout-dialog` v1 @ 791bdbe (attempt 1 of 3; review gate: ui-aspects, design, test)

**Verdict: REJECT, for tests and one docstring.** The pop-out works on real data: the gate is green, and the
deployment-shaped acceptance matched the run's stored web report and found no traceback (below). But six declared
behaviours or §7 rows have no test that fails when they break (rules a/b), and one docstring claim is false (rule d).
Every blocking item below was **reproduced here** in a `git archive` copy of 791bdbe. No production behaviour has to
change except the docstring wording. Stacked on `feature/roi-popout-data` @ f424aec, with the editor stack merged in
(K2 @ 8ca43ce, which contains #43 @ 779f787). Not infrastructure.

## What passed (do not redo)

- **Gate:** `pixi run test-reduction` from the subject root, analysis clone 2, 14:34–14:50 EDT. Launcher **747 passed**,
  reduction **874 passed**, exit 0; the tree was clean afterwards.
- **Scope:** the slug's own diff, taken against the merged predecessors' tree (`git merge-tree --write-tree
  agentic/feature/roi-popout-data agentic/feature/editor-ipts-inference` = 672e986), touches exactly §4's four files.
  The §8.4 grep prints nothing. No `.destroy(` appears in either module.
- **Deployment-shaped acceptance (§8.7):** ledger `scripts/roi-popout-acceptance.py`. It builds the launcher with its
  own main() sequence, offscreen, with QSettings under scratch, then runs these phases:
  - **Real file:** loaded `/SNS/REF_L/IPTS-36119/shared/autoreduce/reduce_settings.json`, read-only (RBnum 231800–231802;
    4 rows, because `useBS` has 4 entries).
    - Each row opened with one click on its cell and one on Select ROI, on the reducer's own file name for its run.
    - The title names the run, and "1 event in 2/3" for the two runs over 2 M events.
    - The peak and background overlays sit on the row's values. Cancel leaves `to_dict()` equal.
    - Row 4 (no run) asked "The NeXus file of angle 4"; cancelled, nothing opened.
  - **Against the stored web report** (`REF_L_231801.html`, the monitor.sns.gov content). The row carried the ROI from
    the reduction-time template (`REF_L_231800_auto_template.xml`: peak 145–155, background [142, 158, 0, 0], x 65–220).
    - The peak band, the background bands' outer edges and the x range sit exactly on the report's lines.
    - XY (full TOF span typed into the view filter): the same 11 088 cells pass the report's 1 % threshold. Every report
      cell is the dialog's raw count × 1.0078–1.0120 (the report's counts are dead-time corrected; template
      `dead_time_correction True`, 4.2 µs). No cell is below the raw count.
    - Y-TOF (x range typed to the whole detector): the report's bins are the dialog's 50 µs bins (residual 0). The same
      69 286 cells pass the threshold, with a factor of 1.0016–1.0180.
    - The log toggle went off and on, with a full draw each time: no error, and no traceback on stderr.
  - **Drag and save:** a drag on the Y profile of row 2 (peak 130–167 → 127–169), then OK. Only RB_Ymin/RB_Ymax were
    reported and written; the table shows 127/169. Saved to scratch, the file differs from an unedited save only in
    RB_Ymin[2] and RB_Ymax[2], and a reload holds them. sha256 `ee0475b339d5c076d39bc38abb4e6778ddff5afa46987c88db227400fdb37254`.
- **Reviewers confirmed:**
  - The lift is byte-identical to 65c83d9 :396–690.
  - A5 credit is in the docstring and the commit body, with no personal trailer.
  - One writer per field, through `set_angle_field` / `set`. Every slot is guarded.
  - V2/V3/V13/E2/E7 pin what §8.6 names.
  - Counts reproduced: test_roi_dialog **72**, launcher **747**.
  - These recorded battery rows reproduced exactly: M1, M2a, M3, M4, M6, M7, M10, M11a, M14a, M15, M16, Fspan-v, Flookup.

## BLOCKING

Each item: what fails, the command, what was seen, and the fix. "Survives" means `launcher/tests/test_roi_dialog.py`
plus `launcher/tests/test_settings_editor.py` gave **618 passed**, the same as the unmutated copy. The full
`launcher/tests` run gives only the 2 `test_overplot_axes` failures, which are an archive-copy artefact (no data) and
appear on the unmutated copy too.

### B-1 (d): `select_roi`'s docstring says "no file is opened for writing (E4)"; the launcher's settings file is written

- **Where:** `launcher/apps/settings_editor.py:1143`.
- **Reproduced:** on a tab with one row, the file dialog returns `REF_L_231801.nxs.h5` and `exec_` returns Rejected.
  After one click on Select ROI and `tab.settings.sync()`, `$XDG_CONFIG_HOME/ORNL/lr_reduction_new_launcher.conf`
  **changed** and now holds `roi_nexus_dir=/SNS/REF_L/IPTS-36119/nexus`.
- **Not the defect:** the write itself is B2's "the last one used", and it is benign.
- **The defect:** the claim is false, and E4 (which pins `SettingsDocument.save` and the working directory) does not
  test what the sentence says. §8.5 names exactly this class of claim ("cannot write a file").
- **Source:** the wording comes from plan §3 ("F2: the dialog and the slot open no file for writing; pinned by E4"),
  which contradicts B2. The Analyst is told so too.
- **Fix (docstring + E4; domain = every file the slot can write):** say what is true, e.g. "no settings or data file is
  written; the folder a chosen file came from is remembered in the launcher's QSettings (B2)". Then make E4 assert it:
  the only change after a Select ROI is that one QSettings key, and no file appears in `tmp_path` or the folder of the
  chosen run.

### B-2 (b): the minor-tick plain formatter can be removed and V11 stays green; mathtext then reaches the screen

- **Where:** `launcher/apps/roi_dialog.py:87`, `axis.set_minor_formatter(LogFormatter(labelOnlyBase=True))`.
- **Reproduced:** deleting the line → **618 passed**.
- **Effect on a real run:** run 231801 (200 k events), with each profile's y axis zoomed to 20–80 (what the toolbar's zoom
  does), then a full draw. **48 minor labels** read `'$\\mathdefault{2\\times10^{0}}$'` and similar; the unmutated copy
  has 0. This is F7/L6's tick-label route, and M11's ("the plain formatter removed from the log axes"). V11 collects only
  the major `get_xticklabels()`/`get_yticklabels()`.
- **Fix (test; domain = every tick and offset label on the three log profiles and both colorbars, major and minor):**
  - V11 draws once at the opening limits and once with a profile zoomed inside one decade.
  - It then asserts that no `get_*ticklabels(minor=True)`, major label or offset text contains `$`.
  - The colorbar minor leg should red on its own. Today it reds only through a pyparsing warning in the empty-run test.

### B-3 (a): the §3 types cell "RB_Ymin/RB_Ymax `None` → no peak overlay" has no test

- **Mutation:** `roi_dialog.py:538`, `self._move("peak", peak)` → `self._move("peak", self._spin_values(self.peak_spins))`.
- **Reproduced:** **618 passed**.
- **Effect:** on run 231801 with `RB_Ymin=None`, the peak overlay is visible on `y_axis`, `xy_axis` and `ytof_axis` (a
  band from −1 to 155 across the lower detector). At 791bdbe it is visible on none.
- **Fix (test):** in V15's (or V14's) `RB_Ymin=None` leg, assert every `dialog.overlays["peak"][…]` artist is hidden,
  and stays hidden until both edges are set.

### B-4 (a): the §3 types cell "`useBS` `[]`/`None` → treated as on" has no test

- **Mutation:** `roi_dialog.py:216`, `values.get("useBS") != 0` → `bool(values.get("useBS"))`.
- **Reproduced:** **618 passed**.
- **Effect:** with `useBS=None`, the status line reads "Background drawn but not subtracted: useBS is off for this
  angle". At 791bdbe it is empty.
- **Reachable:** `doc.set("useBS", [])` makes `angle_row(0)["useBS"]` `None` (ui reviewer, measured). The reducer's
  default for `[]` is on.
- **Fix (test):** `useBS=None` and `useBS=[]`-row legs that assert "not subtracted" is absent, beside the existing 0
  and False legs.

### B-5 (a): B6's drags on the X and TOF profiles have no test

- **Mutation:** both `_x_range_selected` and `_tof_range_selected` (`roi_dialog.py:317-321`) made `pass`.
- **Reproduced:** **618 passed**. The test reviewer also ran each alone, and the swap of the two: each survives.
- **Cause:** V4 drags only on `y_axis`, and the battery has no row for these two callbacks.
- **Production works:** the test reviewer's probe, with the existing `shown`/`drag` helpers:
  - `drag(dialog, dialog.x_axis, 70, 180)` → x spins [70, 180] and `changes() == {"data_x_range": [70, 180]}`;
  - a TOF drag moves the view filter, with `changes() == {}`.
- **Fix (test):** V4 legs for the X and TOF profiles that assert those spins and those `changes()`. Add battery rows for
  each callback made dead, and for the two swapped.

### B-6 (a): "Cancel writes nothing" / "OK writes what changed" are never pressed as gestures; Cancel wired to accept survives

- **Mutation:** `roi_dialog.py:286`, `buttons.rejected.connect(self.reject)` → `…connect(self.accept)`.
- **Reproduced:** **618 passed**. In the launcher, Cancel would then write the edits.
- **Also survives (test reviewer):** leaving OK unconnected (`:285`).
- **Cause:** V6/V7 call `dialog.reject()`/`accept()` directly, and every E-test replaces `exec_`. Plan §6: gestures go
  through QTest, and a direct call is not a gesture. (The acceptance above pressed the real buttons, and production is
  right.)
- **Fix (test):** V6 and V7 press the `QDialogButtonBox` buttons with `QTest.mouseClick`. One E-test runs the real
  modal `exec_()`, driven from a `QTimer` that edits and presses Cancel (and, in a second leg, OK), as
  `scripts/roi-popout-acceptance.py` does.

### B-7 (a): B6's "log toggle" can stop switching to linear and no test fails

- **Mutation:** `roi_dialog.py:562`, `axis.set_yscale("log" if log_scale else "linear")` → `axis.set_yscale("log")`.
- **Reproduced:** **618 passed**.
- **Cause:** the only toggle test is the `_set_log_scale` injection leg, which checks the guard, not the scale.
- **Fix (test):** click `log_check` (a gesture), then assert `get_yscale() == "linear"` on all three profiles. Click it
  again and assert `"log"` and the plain formatters (B-2).

## To show (not counted as a block): M13's recorded count

The commit records "M13 … -> 5 failed", naming `each_plot_is_the_data_layers_for_the_ranges_shown`. The test reviewer
tried three forms of M13 and got **4 failed** each time, with the same set: `drag_past_the_detector_edge`,
`a_nudge…`, and `a_reversed_range` ×2. The three forms were:
- rebuild in `_values_changed`;
- `cla()` + rebuild in `_update`;
- the plan's literal `figure.clear()` + recreate.

V13 reds under every form, so §7 holds. The battery that defines the row's exact text (ledger
`scripts/mutations-roi-popout-dialog.py`, 42c48dd) could not be fetched here, because code.ornl.gov was returning 502.
**In v2, quote the M13 row's code and its observed count.** If it gives 5, the record stands.

## Advisories (non-blocking; carried to the PR body; fix in v2 where cheap)

### ui-aspects

- **A1. The Y-TOF image fills only part of its panel.**
  - Cause: the overlays are made as `axvspan(0, 1)`, which pulls x = 0 into the limits, and `_reset_limits` never sets
    `ytof_axis`.
  - Seen here on run 231801: the axis starts near −3 000 µs while the events start near 8 000 µs. On run 220050 the
    image fills 32 % of the panel.
  - Every real run shows this, and a scientist comparing with monitor.sns.gov will see it.
  - **Recommended in v2:** set the Y-TOF x limits to the TOF edges, and pin them with `get_xlim()`.
- **A2. The background spins keep the bands of the opening peak** for an adjacent `[a, b, 0, 0]` entry. Moving the peak
  redraws the bands, but the spins keep the old four bounds. A one-pixel nudge then makes three stale values the
  written background. What is drawn still equals what is written after the nudge. (The design reviewer flagged this
  too.)
- **A3. A log axis zoomed below one decade shows no numbers**, because `labelOnlyBase=True` leaves the minor labels
  empty. `LogFormatter(labelOnlyBase=False, minor_thresholds=…)` stays plain text.
- **A4. The `LogNorm` is fixed at open.** A narrowed view filter dims the XY image (max 11 against vmax 85); the
  colorbar stays accurate.
- **A5.** The TOF legend lists "reduction TOF window" when the row has none.
- **A6.** A one-row peak (140–140) leaves OK enabled but draws a zero-width overlay. Spans run from pixel centre to pixel
  centre, so each edge pixel is half-covered (open question A8).
- **A7.** A constructor failure (e.g. a reversed `tof_band`) is reported, but leaves a hidden, half-built dialog child
  (§5 "half-built dialog"; the trigger is contrived).
- **A8.** A chopper band wholly outside the event span clamps both TOF spins to one edge, giving empty plots with no
  note.

### design

- **D1 (attribution, a factual claim).** The module docstring and b041aa2's body name **welbournR (3f74d41)** as an
  author of the lifted lines.
  - `git log -L 396,690:launcher/apps/json_settings_builder.py agentic/feature/harden-review-branch` lists only 65c83d9,
    ab22307, 1e692c7, 8191e49 and f5513c7.
  - `git blame -w -M -C` attributes 0 lifted lines to 3f74d41, whose hunks are the constants, `read_nexus_metadata` and
    the tab.
  - The name was carried over from plan F2 (the Analyst is told).
  - Correct the docstring, or say "authors of PR #197's dialog path". 65c83d9's 8 lines (the `parse_math` title) are
    credited as the source SHA but not as an author.
- **D2. `_adjacent` (`:219-223`) re-implements the data layer's background rule.** It treats any zero as "adjacent",
  while `background_bands` needs exactly two. With `[120, 125, 160, 0]` and no peak, the status line gives the wrong
  reason until a peak is set. `_pixel`/`_number` mirror `roi_estimate._whole`.
  - Follow-up for roi-popout-data: a Qt-free validator for the written form.
- **D3.** `f"REF_L_{run}.nxs.h5"` appears twice in `_events_for_row`, and a third copy is in `settings_document.py:159`.
  Use one name per slot, plus a data-layer `run_file(config, run)` follow-up.
  - When RBnum is set but the file is missing, the file dialog gives no "looked in <folder>" hint.
- **D4.** `_guarded` (dialog) nearly duplicates `guarded` (editor).
- **D5.** Size: `roi_dialog.py` is 640 lines; `settings_editor.py` grew from 1308 to 1412.
  - A Qt-free `roi_dialog_state.py` (opening values, background state, `_valid`, `changes()`) would be testable without
    a QApplication.
- **D6. Magic numbers and drift:**
  - magic numbers: the ±30-row zoom margin (`:572`) and the TOF spin step 100 (`:268`);
  - unreachable: `_mode`'s fallback (`:307`);
  - naming drift: `BACKGROUND_LEFT/RIGHT` against `bkg_low/high`.
- **D7 (§8.5).** Two prescriptive comments carried over from #197 cite no test and are not marked inferred: "A layout
  engine would run on every redraw" (`:140`) and "useblit keeps a drag from redrawing" (`:170`).
- **D8.** The QSettings key `roi_nexus_dir` is shared across IPTSs. It is written before the run is read, so a failed
  read still remembers the folder (benign; see B-1).

### test

- **T1.** `_reset_limits`'s TOF and X x-limits have no test (deleting `:573-574` survives).
- **T2.** V15 sets "not set" with `setValue(UNSET)` rather than typing it.
- **T3.** E6/E10 replace the whole read chain (fine under A4). The editor-side chain meets a real file only in §8.7
  (above).

### Analyst (plan corrections, correct-and-flag)

- **P1.** §3 "F2: … open no file for writing; pinned by E4" contradicts B2's "the last one used" (B-1).
- **P2.** F2's author list includes 3f74d41, which wrote none of the lifted lines (D1).
- **P3.** The V9 sparse-run row is the estimator's to decide (the Developer's e0bba12, recorded as an advisory for
  roi-popout-data).
