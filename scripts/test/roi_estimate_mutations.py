#!/usr/bin/env python3
"""Mutation battery for roi_estimate (roi-estimate, then roi-popout-data). Run from the repo root:

    pixi run python scripts/test/roi_estimate_mutations.py [--rows 1-25,32]

Each row applies one mutation to ``src/lr_reduction/roi_estimate.py``, runs the module's tests, and records
which tests failed. ``--rows`` runs a subset, so a long battery can be run in chunks under the 600 s harness
ceiling; every chunk ends with the same baseline run.

What a row's result means:

- A kill is a test that fails with the mutation in the file and passes on the restored module. The suite
  runs once more, unmutated, after the rows (the baseline). A test failing there is subtracted from every
  row, and the run says so.
- The battery's own tests (``META``: I5-I7) are deselected. I5 reads the module for its anchors, so with a
  mutation in the file it fails for every row: a kill that says nothing about the module. A ``META`` name in
  any row's failures is reported, never counted.
- ``SURVIVED`` (no kill) and ``ANCHOR MISS`` (the text to mutate is not in the module exactly once) make
  the battery exit non-zero.

Restore safety per the Developer contract: originals held in memory AND a mode-600 backup, restore in a
``finally`` and in a signal handler, sha256 compared after every mutation with an abort on mismatch, and a
per-invocation timeout far below the 600 s harness ceiling.
"""

import argparse
import hashlib
import os
import re
import signal
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MOD = os.path.join(REPO, "src/lr_reduction/roi_estimate.py")
T = "tests/unit/lr_reduction/test_roi_estimate.py"

# The battery's own tests (I5, I6, I7): deselected from every run, never counted as a kill.
META = (
    "test_the_battery_table_lists_the_rows_it_runs",
    "test_the_mutation_battery_refuses_a_dirty_baseline",
    "test_a_sigterm_during_the_battery_leaves_the_module_restored",
)

MUTATIONS = [
    # ---- roi-estimate (PR #31): v1 rows 1-10, v2 rows 11-16 ----
    (1, 'motors read [0] (the angle BEFORE the move) instead of [-1]',
     'meta[name] = float(f[dpath][-1])',
     'meta[name] = float(f[dpath][0])'),
    (2, 'a missing chopper log defaults instead of refusing',
     '                raise KeyError(\n'
     '                    f"{path} has no chopper log at {dpath!r} — the wavelength "\n'
     '                    f"band cannot be derived for this run, and must not be guessed"\n'
     '                )',
     '                return [2.4, 5.8]'),
    (3, 'inline the band maths (a fifth copy) instead of calling the library',
     '        return nr_tools.get_lam_range(chopper_lam, chopper_speed)',
     '        return [chopper_lam - 1.75 * 60.0 / chopper_speed - 0.15,\n                chopper_lam + 1.75 * 60.0 / chopper_speed - 0.15]'),
    (4, 'fork the id unpacking instead of calling get_y_tof',
     '    _, y_tof, _ = binary_processing.get_y_tof(\n'
     '        tof_array, event_id, e_offset, list(lowres), pcharge, n_y=n_y, n_x=n_x\n'
     '    )',
     '    y_tof = np.zeros((n_y, len(tof_array)))\n'
     '    np.add.at(y_tof, (event_id % n_y, np.zeros(len(event_id), dtype=int)), 1)'),
    (5, 'hard-code the detector shape instead of the instrument DB',
     '    settings = nr_tools.read_settings(start_time)\n    return int(settings["num_x_pixels"]), int(settings["num_y_pixels"])',
     '    return 256, 304'),
    (6, 'drop the no-counts guard',
     '    if counts.size == 0 or not np.any(np.isfinite(counts) & (counts > 0)):',
     '    if False:'),
    (7, 'drop the contrast guard',
     '    if contrast < min_contrast:\n        raise CannotEstimateError(',
     '    if False:\n        raise CannotEstimateError('),
    (9, 'drop the no-room refusal (re-anchored on the two-sided default_bkg_roi, B7)',
     '    if b0 < 1 or b3 > n_y - 1:',
     '    if False:'),
    (10, 'normalise by the DASlogs series instead of the run total',
     '        pcharge = np.asarray(f["entry/proton_charge"][:])',
     '        pcharge = np.asarray(f["entry/DASlogs/proton_charge/value"][:])'),
    (11, 'A1: let a zero/NaN baseline score as infinite contrast again',
     '    if not np.isfinite(baseline) or baseline <= 0:',
     '    if False:'),
    (12, 'A2: drop the non-finite profile refusal',
     '    if not np.all(np.isfinite(counts)):',
     '    if False:'),
    (13, 'A2: drop the zero proton-charge refusal at source',
     '    if not np.isfinite(total_charge) or total_charge <= 0:',
     '    if False:'),
    (14, 'D1: drop the on-detector validation of peak_range',
     '    if not 0 <= peak_low <= peak_high <= n_y - 1:',
     '    if False:'),
    (15, 'C3: drop the inverted-band refusal',
     '    if hi <= lo:',
     '    if False:'),
    (16, 'C1: lambda_to_tof forgets the mm->m conversion',
     '    flight_m = float(settings["source_detector_distance"]) / 1000.0',
     '    flight_m = float(settings["source_detector_distance"])'),
    # ---- roi-popout-data: the plan's section 7, M1-M15 ----
    (17, 'M1: xy_image counts in id order (x * n_y + y), reshaped (n_y, n_x) without the transpose',
     '    flat = np.bincount(events.y[keep] * events.n_x + events.x[keep], minlength=events.n_x * events.n_y)',
     '    flat = np.bincount(events.x[keep] * events.n_y + events.y[keep], minlength=events.n_x * events.n_y)'),
    (18, 'M2: the off-detector filter removed',
     '    on_detector = (event_id >= 0) & (event_id < n_x * n_y)',
     '    on_detector = np.ones(len(event_id), dtype=bool)'),
    (19, 'M3: tof_edges cropped to a window (5th-95th percentile, standing in for the chopper band), not the span',
     '    low, high = float(events.tof.min()), float(events.tof.max())',
     '    low, high = (float(v) for v in np.percentile(events.tof, [5, 95]))'),
    (20, 'M4: y_tof_image ignores x_range',
     '    selected = (events.x >= x_low) & (events.x <= x_high)',
     '    selected = np.ones(len(events.x), dtype=bool)'),
    (21, 'M5: the last TOF bin made half-open (an event on the last edge dropped)',
     '    selected = (events.x >= x_low) & (events.x <= x_high)',
     '    selected = (events.x >= x_low) & (events.x <= x_high) & (events.tof < edges[-1])'),
    (22, 'M6: the packing swapped (x = id % n_y, y = id // n_y)',
     '    return RunEvents(x=event_id // n_y, y=event_id % n_y,',
     '    return RunEvents(x=event_id % n_y, y=event_id // n_y,'),
    (23, 'M7: background_bands sorts the entry without putting the peak in place of the two zeros',
     '        ordered[0], ordered[1] = y_min, y_max\n        ordered = np.sort(ordered)',
     '        pass'),
    (24, "M8: background_bands' zero-count refusal -> if False",
     '    elif zeros != 0:',
     '    elif False:'),
    (25, 'M9: default_bkg_roi clamps to the detector instead of refusing',
     '    if b0 < 1 or b3 > n_y - 1:',
     '    b0, b3 = max(b0, 0), min(b3, n_y - 1)\n    if False:'),
    (26, 'M10: default_bkg_roi returns only the low-side band',
     '    return b0, b1, b2, b3',
     '    return b0, b1'),
    (27, 'M11: stride sampling -> the head slice [:max_events]',
     '        event_id, tof = event_id[::stride], tof[::stride]',
     '        event_id, tof = event_id[:max_events], tof[:max_events]'),
    (28, "M12: counts_vs_y's default lowres -> the literal (0, 255)",
     '        lowres = (0, n_x - 1)',
     '        lowres = (0, 255)'),
    (29, "M13: the tof_band selection in counts_vs_y -> if False",
     '    if tof_band is not None:\n        # get_y_tof clips',
     '    if False:\n        # get_y_tof clips'),
    (30, 'M14: import qtpy at the module top',
     'from dataclasses import dataclass\n',
     'from dataclasses import dataclass\n\nimport qtpy  # noqa: F401\n'),
    (31, "M15: tof_edges' no-events refusal removed",
     '    if len(events.tof) == 0:\n        raise CannotEstimateError(',
     '    if False:\n        raise CannotEstimateError('),
    # ---- roi-popout-data frame: each guard GREEN added beyond section 7 ----
    (32, "background_bands accepts a NaN or infinite bound",
     '    if not numeric or not np.all(np.isfinite(bounds)):',
     '    if not numeric:'),
    (33, "background_bands accepts bounds that are not numbers",
     '    numeric = np.issubdtype(bounds.dtype, np.integer) or np.issubdtype(bounds.dtype, np.floating)',
     '    numeric = True'),
    (34, "background_bands accepts a nested list (the per-angle BkgROI)",
     '    if bounds.ndim != 1:',
     '    if False:'),
    (36, "background_bands lets numpy's error for a ragged entry through",
     '    try:\n        bounds = np.asarray(bkg_roi)\n    except ValueError as exc:  # ragged nesting\n'
     '        raise ValueError(f"a background is one angle\'s four pixel bounds, not {bkg_roi!r}") from exc',
     '    bounds = np.asarray(bkg_roi)'),
    (37, "default_bkg_roi accepts gap < 0 or width < 1",
     '    if gap < 0 or width < 1:',
     '    if False:'),
    (38, "load_event_pixels accepts max_events <= 0",
     '    if max_events is not None and max_events <= 0:',
     '    if False:'),
    (39, "the shared edges check -> if False (one edge, a repeated edge)",
     '    if edges.ndim != 1 or len(edges) < 2 or np.any(np.diff(edges) <= 0):',
     '    if False:'),
    (40, "a reversed range is an empty selection, not a ValueError",
     '    if high < low:\n        raise ValueError(f"{name}',
     '    if False:\n        raise ValueError(f"{name}'),
    (41, "tof_edges trusts ceil(): the rounding guard -> if False",
     '    if edges[-1] < high:',
     '    if False:'),
    (42, "tof_edges without max(1, ...): a zero span gives a single edge",
     '    edges = low + bin_width * np.arange(max(1, int(np.ceil((high - low) / bin_width))) + 1)',
     '    edges = low + bin_width * np.arange(int(np.ceil((high - low) / bin_width)) + 1)'),
    (43, "tof_edges accepts bin_width <= 0",
     '    if not bin_width > 0:',
     '    if False:'),
    (44, "n_off_detector not counted",
     '    n_off_detector = int(np.count_nonzero(~on_detector))',
     '    n_off_detector = 0'),
    (45, "the stride not recorded",
     'stride=stride, n_off_detector=n_off_detector)',
     'stride=1, n_off_detector=n_off_detector)'),
    (46, "the TOF band selects every event",
     '    return (tof >= low) & (tof <= high)',
     '    return np.ones(len(tof), dtype=bool)'),
    (47, "xy_image drops its band",
     '    keep = _in_band(events.tof, tof_band)\n    flat',
     '    keep = _in_band(events.tof, None)\n    flat'),
    (48, "profile_y ignores x_range",
     '    keep = (events.x >= x_low) & (events.x <= x_high) & _in_band(events.tof, tof_band)',
     '    keep = _in_band(events.tof, tof_band)'),
    (49, "profile_x ignores y_range",
     '        keep &= (events.y >= y_low) & (events.y <= y_high)\n    return np.bincount(events.x[keep]',
     '        pass\n    return np.bincount(events.x[keep]'),
    (50, "profile_tof ignores x_range",
     '        keep &= (events.x >= x_low) & (events.x <= x_high)',
     '        pass'),
    (51, "profile_tof ignores y_range",
     '        keep &= (events.y >= y_low) & (events.y <= y_high)\n    counts, _',
     '        pass\n    counts, _'),
]

# Table rows that no longer run, with the reason (I5: the table lists MUTATIONS + RETIRED).
RETIRED = {
    8: "its target, default_bkg_roi's low-side fit check, is gone: B7 (roi-popout-data) returns a band on "
       "each side, and the no-room refusal that replaced both fit checks is row 9",
    35: "SURVIVED its first run: background_bands' text refusal was unreachable, because a string is a 0-d "
        "array, which the ndim refusal (row 34) catches; removed from the module in dfe1683",
}

# "FAILED path::test[id] - msg" and "ERROR path::test - msg"; a collection error has no "::test".
_FAILED = re.compile(r"^(?:FAILED|ERROR) (\S+)")


def sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def verify_baseline_matches_head():
    """Refuse to start unless the target matches its committed blob.

    The first version took `sha(MOD)` of whatever was on disk and called that
    "clean". With zero git references, a leftover from a killed run BECAME the
    baseline: the battery then restored to the mutated text and printed
    `restored: OK` with the mutation still in the file. That is precisely the
    failure `todo-mutation-harness-restore-safety` was written about, and since
    this is the campaign's first committed battery it is the reference
    implementation — so the bug propagates by being copied.

    Detection complete: compare against `git show HEAD:<path>`, not against the
    working tree. A battery that cannot tell dirty from clean cannot honestly
    report anything.
    """
    rel = os.path.relpath(MOD, REPO)
    try:
        blob = subprocess.run(
            ["git", "show", f"HEAD:{rel}"], cwd=REPO,
            capture_output=True, check=True,
        ).stdout
    except (subprocess.CalledProcessError, OSError) as exc:
        raise SystemExit(f"ABORT: cannot read the HEAD blob for {rel}: {exc}") from exc
    if hashlib.sha256(blob).hexdigest() != sha(MOD):
        raise SystemExit(
            f"ABORT: {rel} differs from HEAD — the working tree is dirty, so it "
            f"cannot serve as the mutation baseline. Commit or restore it first. "
            f"(A leftover from a killed run looks exactly like this.)"
        )


def parse_rows(spec):
    """``"1-16,32"`` -> {1, ..., 16, 32}; None -> every row."""
    if spec is None:
        return None
    rows = set()
    for part in spec.split(","):
        low, _, high = part.partition("-")
        rows.update(range(int(low), int(high or low) + 1))
    return rows


def run_tests():
    """The module's tests, META deselected: (failed test ids, summary line)."""
    source = open(os.path.join(REPO, T), encoding="utf-8").read()
    deselect = [arg for name in META if f"def {name}(" in source for arg in ("--deselect", f"{T}::{name}")]
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", "-rfE", "--no-header", "-p", "no:cacheprovider",
             "--timeout=90", "--timeout-method=thread", *deselect, T],
            cwd=REPO, capture_output=True, text=True, timeout=240,
        )
    except subprocess.TimeoutExpired:
        return None, "HUNG (harness timeout)"
    failed = {m.group(1).partition("::")[2] or "<collection error>"
              for m in map(_FAILED.match, proc.stdout.splitlines()) if m}
    lines = [x for x in proc.stdout.strip().splitlines() if x.strip()]
    return failed, (lines[-1] if lines else f"exit {proc.returncode}")


def describe(failed):
    """Failed test ids, parametrized cases folded: "test_a (x3), test_b"."""
    counts = {}
    for test in sorted(failed):
        base = test.split("[", 1)[0]
        counts[base] = counts.get(base, 0) + 1
    return ", ".join(name if n == 1 else f"{name} (x{n})" for name, n in counts.items())


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--rows", help='rows to run, e.g. "1-25" or "9,14,32-36" (default: every row)')
    wanted = parse_rows(parser.parse_args().rows)

    # Baseline first: never adopt the working tree sight-unseen.
    verify_baseline_matches_head()

    orig = open(MOD, encoding="utf-8").read()
    clean = sha(MOD)
    fd, bak = tempfile.mkstemp(prefix="roi-", suffix=".bak")
    os.close(fd)
    os.chmod(bak, 0o600)
    open(bak, "w", encoding="utf-8").write(orig)

    # A `finally` does not run when the process is SIGTERM'd, which is exactly
    # how the 600 s harness ceiling kills a long battery — demonstrated: killed
    # at the limit, the target was left mutated. Restore from the handler too.
    def _restore_and_die(signum, _frame):
        with open(MOD, "w", encoding="utf-8") as fh:
            fh.write(orig)
        print(f"\nsignal {signum}: restored {MOD} before exiting", file=sys.stderr)
        raise SystemExit(128 + signum)

    for _sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(_sig, _restore_and_die)

    results = []
    try:
        for row, desc, old, new in MUTATIONS:
            if wanted is not None and row not in wanted:
                continue
            if orig.count(old) != 1:
                results.append((row, desc, None, f"ANCHOR MISS (x{orig.count(old)})"))
                print(f"[{row}] ANCHOR MISS — {desc}", flush=True)
                continue
            try:
                open(MOD, "w", encoding="utf-8").write(orig.replace(old, new))
                failed, summary = run_tests()
            finally:
                open(MOD, "w", encoding="utf-8").write(orig)
                if sha(MOD) != clean:
                    raise SystemExit(f"ABORT: {MOD} did not restore cleanly")
            results.append((row, desc, failed, summary))
            print(f"[{row}] {desc}\n     {summary}", flush=True)
    finally:
        open(MOD, "w", encoding="utf-8").write(orig)
        print("restored:", "OK" if sha(MOD) == clean else "*** DIRTY ***")
        os.unlink(bak)

    baseline, summary = run_tests()
    print(f"baseline (unmutated): {summary}")
    if baseline is None:
        raise SystemExit("ABORT: the baseline run hung, so no row's failures can be read as kills")
    if baseline:
        print(f"*** the baseline fails {describe(baseline)}: subtracted from every row below ***")

    bad = 0
    print("\n=== ledger rows ===")
    for row, desc, failed, summary in results:
        if failed is None:
            bad += 1
            print(f"| {row} | {desc} | {summary} |")
            continue
        meta = {test for test in failed if test.split("[", 1)[0] in META}
        kills = failed - baseline - meta
        if not kills:
            bad += 1
        observed = f"{describe(kills)} -> {len(kills)} failed" if kills else f"SURVIVED ({summary})"
        if meta:
            observed += f"; META (not counted): {describe(meta)}"
        print(f"| {row} | {desc} | {observed} |")
    if bad:
        raise SystemExit(f"{bad} row(s) SURVIVED, hung or missed their anchor")


if __name__ == "__main__":
    main()


# Measured 2026-10-05 on feature/roi-popout-data at dfe1683, in two chunks: --rows 1-25 (186 s) and
# --rows 26-51 (194 s). 49 rows ran and all 49 are red; rows 8 and 35 are retired (RETIRED). Each chunk
# printed "restored: OK" and the baseline "69 passed, 3 deselected" (META); `git status --porcelain` showed
# only this file, uncommitted at the time.
#
# | # | mutation | observed (kills: tests red with the mutation, green without) |
# |---|---|---|
# | 1 | motors read [0] (the angle BEFORE the move) instead of [-1] | test_read_nexus_metadata_reads_the_run_identity_and_angles, test_read_nexus_metadata_takes_the_LAST_motor_sample_not_the_first -> 2 failed |
# | 2 | a missing chopper log defaults instead of refusing | test_chopper_lambda_range_does_not_cache_across_files, test_chopper_lambda_range_refuses_when_the_chopper_log_is_ABSENT -> 2 failed |
# | 3 | inline the band maths (a fifth copy) instead of calling the library | test_chopper_lambda_range_uses_the_library_not_a_private_copy -> 1 failed |
# | 4 | fork the id unpacking instead of calling get_y_tof | test_counts_vs_y_calls_the_histogrammer_once_with_or_without_a_band, test_counts_vs_y_default_lowres_follows_the_detector_database, test_counts_vs_y_excludes_events_outside_the_x_range, test_counts_vs_y_goes_through_the_library_histogrammer, test_profile_y_agrees_with_the_library_histogrammer -> 5 failed |
# | 5 | hard-code the detector shape instead of the instrument DB | test_counts_vs_y_default_lowres_follows_the_detector_database, test_the_event_geometry_follows_the_instrument_database, test_the_geometry_comes_from_the_instrument_database_not_a_literal -> 3 failed |
# | 6 | drop the no-counts guard | test_estimate_peak_range_refuses_an_empty_detector -> 1 failed |
# | 7 | drop the contrast guard | test_estimate_peak_range_refuses_a_featureless_detector -> 1 failed |
# | 8 | (retired) drop the low-side FIT CHECK | not run: its target is gone with B7 (RETIRED) |
# | 9 | drop the no-room refusal (re-anchored on the two-sided default_bkg_roi, B7) | test_default_bkg_roi_never_returns_a_band_off_the_detector, test_default_bkg_roi_refuses_a_peak_that_leaves_no_room, test_default_bkg_roi_survives_the_reducer -> 3 failed |
# | 10 | normalise by the DASlogs series instead of the run total | test_a_zero_proton_charge_is_refused_rather_than_bracketing_pixel_zero, test_counts_vs_y_calls_the_histogrammer_once_with_or_without_a_band, test_counts_vs_y_default_lowres_follows_the_detector_database, test_counts_vs_y_excludes_events_outside_the_x_range, test_counts_vs_y_finds_the_injected_peak, test_counts_vs_y_goes_through_the_library_histogrammer, test_counts_vs_y_with_a_band_counts_exactly_the_events_inside_it, test_estimate_peak_range_brackets_the_injected_peak, test_estimate_peak_range_reports_its_contrast_when_asked, test_profile_y_agrees_with_the_library_histogrammer, test_the_lambda_range_and_the_tof_band_compose -> 11 failed |
# | 11 | A1: let a zero/NaN baseline score as infinite contrast again | test_a_sparse_profile_is_refused_not_scored_as_infinite_contrast, test_a_sparse_real_run_gives_images_and_a_refused_estimate -> 2 failed |
# | 12 | A2: drop the non-finite profile refusal | test_a_non_finite_profile_is_refused_even_if_it_reaches_the_estimator -> 1 failed |
# | 13 | A2: drop the zero proton-charge refusal at source | test_a_zero_proton_charge_is_refused_rather_than_bracketing_pixel_zero -> 1 failed |
# | 14 | D1: drop the on-detector validation of peak_range | test_default_bkg_roi_refuses_a_peak_that_is_not_on_the_detector (x5) -> 5 failed |
# | 15 | C3: drop the inverted-band refusal | test_counts_vs_y_refuses_an_inverted_band -> 1 failed |
# | 16 | C1: lambda_to_tof forgets the mm->m conversion | test_the_lambda_range_and_the_tof_band_compose -> 1 failed |
# | 17 | M1: xy_image counts in id order (x * n_y + y), reshaped (n_y, n_x) without the transpose | test_profiles_are_marginals_of_the_images, test_xy_image_is_what_the_web_report_plots, test_xy_image_puts_the_injected_peak_at_its_row_and_columns -> 3 failed |
# | 18 | M2: the off-detector filter removed | test_off_detector_ids_are_dropped_and_counted -> 1 failed |
# | 19 | M3: tof_edges cropped to a window (5th-95th percentile, standing in for the chopper band), not the span | test_profiles_are_marginals_of_the_images, test_tof_edges_hold_the_latest_event_when_rounding_falls_short, test_tof_edges_span_every_event_not_the_chopper_band, test_y_tof_image_counts_only_the_x_range -> 4 failed |
# | 20 | M4: y_tof_image ignores x_range | test_profiles_are_marginals_of_the_images, test_y_tof_image_counts_only_the_x_range -> 2 failed |
# | 21 | M5: the last TOF bin made half-open (an event on the last edge dropped) | test_y_tof_image_keeps_the_event_at_the_last_edge -> 1 failed |
# | 22 | M6: the packing swapped (x = id % n_y, y = id // n_y) | test_load_event_pixels_holds_the_events_and_the_detector_shape, test_profile_y_agrees_with_the_library_histogrammer, test_profiles_are_marginals_of_the_images, test_stride_sampling_is_not_a_time_slice, test_tof_edges_span_every_event_not_the_chopper_band, test_xy_image_is_what_the_web_report_plots, test_xy_image_puts_the_injected_peak_at_its_row_and_columns, test_y_tof_image_counts_only_the_x_range, test_y_tof_image_keeps_the_event_at_the_last_edge -> 9 failed |
# | 23 | M7: background_bands sorts the entry without putting the peak in place of the two zeros | test_background_bands_are_the_rows_the_reducer_averages (x3) -> 3 failed |
# | 24 | M8: background_bands' zero-count refusal -> if False | test_background_bands_refuses_what_the_reducer_cannot_use (x3) -> 3 failed |
# | 25 | M9: default_bkg_roi clamps to the detector instead of refusing | test_default_bkg_roi_never_returns_a_band_off_the_detector, test_default_bkg_roi_refuses_a_peak_that_leaves_no_room, test_default_bkg_roi_survives_the_reducer -> 3 failed |
# | 26 | M10: default_bkg_roi returns only the low-side band | test_default_bkg_roi_defaults_are_the_reviewed_three_and_five, test_default_bkg_roi_sits_outside_the_peak_with_a_gap, test_default_bkg_roi_survives_the_reducer -> 3 failed |
# | 27 | M11: stride sampling -> the head slice [:max_events] | test_stride_sampling_is_not_a_time_slice -> 1 failed |
# | 28 | M12: counts_vs_y's default lowres -> the literal (0, 255) | test_counts_vs_y_default_lowres_follows_the_detector_database -> 1 failed |
# | 29 | M13: the tof_band selection in counts_vs_y -> if False | test_counts_vs_y_with_a_band_counts_exactly_the_events_inside_it -> 1 failed |
# | 30 | M14: import qtpy at the module top | test_the_module_imports_without_any_gui_package -> 1 failed |
# | 31 | M15: tof_edges' no-events refusal removed | test_an_empty_run_gives_zero_images_and_refuses_edges -> 1 failed |
# | 32 | background_bands accepts a NaN or infinite bound | test_background_bands_refuses_what_the_reducer_cannot_use (x2) -> 2 failed |
# | 33 | background_bands accepts bounds that are not numbers | test_background_bands_refuses_what_the_reducer_cannot_use -> 1 failed |
# | 34 | background_bands accepts a nested list (the per-angle BkgROI) | test_background_bands_refuses_what_the_reducer_cannot_use (x2) -> 2 failed |
# | 35 | (retired) background_bands' text refusal -> if False | SURVIVED its first run (unreachable: the ndim refusal catches text); removed in dfe1683 (RETIRED) |
# | 36 | background_bands lets numpy's error for a ragged entry through | test_background_bands_refuses_what_the_reducer_cannot_use -> 1 failed |
# | 37 | default_bkg_roi accepts gap < 0 or width < 1 | test_called_wrong_names_what_is_wrong -> 1 failed |
# | 38 | load_event_pixels accepts max_events <= 0 | test_called_wrong_is_a_value_error_not_an_empty_answer -> 1 failed |
# | 39 | the shared edges check -> if False (one edge, a repeated edge) | test_called_wrong_names_what_is_wrong -> 1 failed |
# | 40 | a reversed range is an empty selection, not a ValueError | test_called_wrong_is_a_value_error_not_an_empty_answer, test_called_wrong_names_what_is_wrong -> 2 failed |
# | 41 | tof_edges trusts ceil(): the rounding guard -> if False | test_tof_edges_hold_the_latest_event_when_rounding_falls_short -> 1 failed |
# | 42 | tof_edges without max(1, ...): a zero span gives a single edge | test_tof_edges_give_one_bin_when_every_event_shares_a_tof -> 1 failed |
# | 43 | tof_edges accepts bin_width <= 0 | test_called_wrong_is_a_value_error_not_an_empty_answer -> 1 failed |
# | 44 | n_off_detector not counted | test_off_detector_ids_are_dropped_and_counted -> 1 failed |
# | 45 | the stride not recorded | test_stride_sampling_is_not_a_time_slice -> 1 failed |
# | 46 | the TOF band selects every event | test_profiles_are_marginals_of_the_images -> 1 failed |
# | 47 | xy_image drops its band | test_called_wrong_is_a_value_error_not_an_empty_answer, test_profiles_are_marginals_of_the_images -> 2 failed |
# | 48 | profile_y ignores x_range | test_profile_y_agrees_with_the_library_histogrammer -> 1 failed |
# | 49 | profile_x ignores y_range | test_profiles_are_marginals_of_the_images -> 1 failed |
# | 50 | profile_tof ignores x_range | test_profiles_are_marginals_of_the_images -> 1 failed |
# | 51 | profile_tof ignores y_range | test_profiles_are_marginals_of_the_images -> 1 failed |
#
# This battery's earlier runs are in `git log --follow` of this file. Two of them taught something that
# still holds:
#
# - roi-estimate v1 (2026-09-20, 10 rows): rows 5 and 8 survived the first pass.
#   - Row 5: the database holds one entry for each pixel count, so the literal 256/304 equals the
#     database's value, and no value assertion can tell them apart. The guard asserts PROVENANCE instead:
#     move the database under the function, and the answer must move.
#   - Row 8's clamps were dead code, each inside a branch whose own condition forbade the out-of-range
#     case.
# - roi-popout-data's first pass (2026-10-05; module at 33ce53b, this file uncommitted): 46 of 51 red.
#   - Row 14 (D1): B7's no-room refusal also says "detector".
#   - Rows 33, 34, 36: no test reached the numeric, ndim or ragged-entry guard.
#   - Row 35: the text refusal was unreachable.
#   525f07c pinned rows 14, 33, 34 and 36, and dfe1683 removed row 35's guard. A guard that no row can
#   red reads as protection and gives none.
