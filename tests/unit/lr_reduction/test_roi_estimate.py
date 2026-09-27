"""Tests for the Qt-free ROI/metadata estimation module (T1 slug 1).

The fixture is BUILT, not committed. The plan describes a synthetic
`tests/data/*.nxs.h5`; generating it from a committed builder gives the same
file without putting a binary in git, and the builder is the part a reader needs
in order to know what the fixture actually asserts. It also lets each test vary
one thing — notably the PRESENCE of a log, which amendment 21 requires.
"""

import os

import h5py
import numpy as np
import pytest

from lr_reduction import roi_estimate as re_mod

N_Y = 304
N_X = 256


def _write_nexus(
    path,
    *,
    peak_y=150,
    peak_width=1.0,
    n_events=60000,
    title="synthetic round-trip run",
    run_number=213628,
    seq_num=1,
    seq_id=7,
    ths=0.6,
    thi=0.6,
    tthd=1.2,
    chopper_lam=4.25,
    chopper_speed=60.0,
    with_chopper=True,
    proton_charge=1.0e12,
    start_time="2025-03-04T11:22:33-05:00",
):
    """Write a minimal REF_L-shaped NeXus file.

    Only the groups this module reads. Event ids follow the instrument's own
    packing — `x = id // n_y`, `y = id % n_y` (`binary_processing.get_y_tof`) —
    so a test that asserts a peak position is asserting about the same
    convention production uses, not about a convention invented here.
    """
    rng = np.random.default_rng(1234)
    # A flat background under the peak. Without it the profile's median is 0,
    # which is A1's defective regime: the contrast score returned inf and the
    # assertions below were satisfied BY the bug. Real detectors are not
    # backgroundless, and `peak_width=1.0` matches the real bracket width
    # (median 3 px across 63 REF_L files; the old 6.0 gave 13 and left the
    # low-side half-max walk deletable with every test still green).
    n_bkg = max(1, n_events // 5)
    y_peak = np.clip(rng.normal(peak_y, peak_width, n_events), 0, N_Y - 1)
    y_bkg = rng.integers(0, N_Y, n_bkg)
    y = np.concatenate([y_peak, y_bkg]).astype(np.int64)
    x = rng.integers(100, 160, len(y))
    event_id = x * N_Y + y
    tof = rng.uniform(10000.0, 40000.0, len(y))

    with h5py.File(path, "w") as f:
        entry = f.create_group("entry")
        entry.create_dataset("title", data=[title.encode()])
        entry.create_dataset("start_time", data=[start_time.encode()])
        entry.create_dataset("run_number", data=[str(run_number).encode()])
        # The TOTAL accumulated charge, which is what `load_and_extract` passes
        # to `get_y_tof` as `pcharge` — not the DASlogs time series beside it.
        entry.create_dataset("proton_charge", data=np.array([proton_charge]))

        events = entry.create_group("bank1_events")
        events.create_dataset("event_id", data=event_id)
        events.create_dataset("event_time_offset", data=tof)

        logs = entry.create_group("DASlogs")

        def log(name, values):
            g = logs.create_group(name)
            g.create_dataset("value", data=np.asarray(values))
            g.create_dataset("average_value", data=np.asarray([np.mean(values)]))

        log("BL4B:Mot:ths.RBV", [ths - 0.01, ths])
        log("BL4B:Mot:thi.RBV", [thi - 0.01, thi])
        log("BL4B:Mot:tthd.RBV", [tthd - 0.01, tthd])
        log("BL4B:CS:Autoreduce:Sequence:Num", [seq_num])
        log("BL4B:CS:Autoreduce:Sequence:Id", [seq_id])
        log("proton_charge", [1.0e12, 1.0e12])
        if with_chopper:
            log("LambdaRequest", [chopper_lam])
            log("SpeedRequest1", [chopper_speed])
    return path


@pytest.fixture
def nexus(tmp_path):
    return _write_nexus(tmp_path / "REF_L_213628.nxs.h5")


# -- metadata ---------------------------------------------------------------


def test_read_nexus_metadata_reads_the_run_identity_and_angles(nexus):
    meta = re_mod.read_nexus_metadata(nexus)

    assert meta["run_number"] == 213628
    assert meta["seq_num"] == 1
    assert meta["seq_id"] == 7
    assert meta["title"] == "synthetic round-trip run"
    assert meta["ths"] == pytest.approx(0.6)
    assert meta["thi"] == pytest.approx(0.6)
    assert meta["tthd"] == pytest.approx(1.2)
    assert meta["start_time"].startswith("2025-03-04")


def test_read_nexus_metadata_takes_the_LAST_motor_sample_not_the_first(tmp_path):
    """`binary_processing.get_log_values` uses [-1] for the motors, and so must this.

    A motor log's first sample is where the axis was before it moved. Taking it
    would report the previous run's angle for every run, silently and plausibly.
    """
    path = _write_nexus(tmp_path / "m.nxs.h5", ths=0.9)
    meta = re_mod.read_nexus_metadata(path)
    assert meta["ths"] == pytest.approx(0.9)


# -- chopper window (amendment 21: vary PRESENCE, not only value) ------------


def test_chopper_lambda_range_uses_the_library_not_a_private_copy(nexus):
    """#197 grew a third copy of this maths at 3.5 A; the library default is 3.4.

    Asserting equality with `nr_tools.get_lam_range` rather than with a literal
    is the point: a literal here would BE a fourth copy, and would keep passing
    if the library were corrected.
    """
    from lr_reduction.nr_tools import get_lam_range

    window = re_mod.chopper_lambda_range(nexus)

    assert window == pytest.approx(get_lam_range(4.25, 60.0))


def test_chopper_lambda_range_refuses_when_the_chopper_log_is_ABSENT(tmp_path):
    """Amendment 21: absent is a state, not a value.

    The window must be derived or refused — never replayed from a previous
    file and never quietly defaulted, which would hand the caller a band that
    belongs to a different measurement.
    """
    path = _write_nexus(tmp_path / "nochop.nxs.h5", with_chopper=False)
    with pytest.raises(KeyError, match="chopper"):
        re_mod.chopper_lambda_range(path)


def test_chopper_lambda_range_does_not_cache_across_files(tmp_path):
    """The stale-window failure, stated as a test rather than trusted."""
    a = _write_nexus(tmp_path / "a.nxs.h5", chopper_lam=4.25, chopper_speed=60.0)
    b = _write_nexus(tmp_path / "b.nxs.h5", chopper_lam=9.0, chopper_speed=30.0)

    first = re_mod.chopper_lambda_range(a)
    second = re_mod.chopper_lambda_range(b)

    assert first != second
    absent = _write_nexus(tmp_path / "c.nxs.h5", with_chopper=False)
    with pytest.raises(KeyError):
        re_mod.chopper_lambda_range(absent)


# -- counts vs y ------------------------------------------------------------


def test_counts_vs_y_finds_the_injected_peak(nexus):
    """One array per detector row, peaking where the events were put."""
    counts = re_mod.counts_vs_y(nexus, lowres=(100, 160))

    assert counts.shape == (N_Y,)
    assert int(np.argmax(counts)) == pytest.approx(150, abs=3)
    assert counts.sum() > 0


def test_counts_vs_y_goes_through_the_library_histogrammer(nexus, monkeypatch):
    """`binary_processing.get_y_tof` owns the id unpacking and the x filter.

    Pinned by observation rather than by trust: forking that maths is exactly
    how the bandwidth constant reached four values, and a private copy here
    would drift the same way.
    """
    from lr_reduction import binary_processing

    called = {}
    real = binary_processing.get_y_tof

    def spy(*args, **kwargs):
        called["yes"] = True
        return real(*args, **kwargs)

    monkeypatch.setattr(binary_processing, "get_y_tof", spy)
    re_mod.counts_vs_y(nexus, lowres=(100, 160))
    assert called.get("yes"), "counts_vs_y did not use the library histogrammer"


def test_counts_vs_y_excludes_events_outside_the_x_range(nexus):
    """The x filter is the library's; this asserts it is actually reaching it."""
    inside = re_mod.counts_vs_y(nexus, lowres=(100, 160)).sum()
    outside = re_mod.counts_vs_y(nexus, lowres=(0, 10)).sum()
    assert outside < inside
    assert outside == 0


# -- peak estimate (amendment 21: an empty detector is a STATE) -------------


def test_estimate_peak_range_brackets_the_injected_peak(nexus):
    counts = re_mod.counts_vs_y(nexus, lowres=(100, 160))
    low, high = re_mod.estimate_peak_range(counts)

    assert low < 150 < high
    # <= 5, not < 40: real REF_L brackets are median 3 px. At the old width the
    # low-side half-max walk was deletable with all tests green.
    assert high - low <= 5, "the bracket is far wider than a real specular peak"


def test_estimate_peak_range_refuses_an_empty_detector():
    """No counts is a state, not a peak at pixel 0.

    `argmax` of an all-zero array is 0, so an unguarded estimator returns a
    confident ROI at the edge of the detector for a run that recorded nothing.
    """
    with pytest.raises(ValueError, match="no counts"):
        re_mod.estimate_peak_range(np.zeros(N_Y))


def test_estimate_peak_range_refuses_a_featureless_detector():
    """Flat illumination has an argmax too, and it means nothing.

    The contrast score is what separates "a peak" from "the largest sample of
    noise", and without it the caller cannot tell the two apart.
    """
    rng = np.random.default_rng(0)
    flat = rng.normal(100.0, 1.0, N_Y)
    with pytest.raises(ValueError, match="contrast"):
        re_mod.estimate_peak_range(flat)


def test_estimate_peak_range_reports_its_contrast_when_asked(nexus):
    counts = re_mod.counts_vs_y(nexus, lowres=(100, 160))
    low, high, contrast = re_mod.estimate_peak_range(counts, with_contrast=True)
    assert low < 150 < high
    assert np.isfinite(contrast), "an infinite contrast means the baseline was zero"
    assert contrast > 1.0


# -- background ROI ---------------------------------------------------------


def test_default_bkg_roi_sits_outside_the_peak_with_a_gap():
    low, high = re_mod.default_bkg_roi((140, 160), n_y=N_Y, gap=5, width=10)
    assert high < 140 - 5 or low > 160 + 5


def test_default_bkg_roi_never_returns_a_band_off_the_detector():
    """Swept, not sampled: every peak position must give an on-detector band.

    A single position tested only the branch that position happens to take. The
    sweep is what shows the low-side and high-side fit checks are both doing
    work — and it is how the redundant `max`/`min` clamps were found to be dead
    code sitting inside branches that already forbade the out-of-range case.
    """
    for peak_low in range(0, N_Y - 1, 7):
        peak = (peak_low, min(peak_low + 12, N_Y - 1))
        try:
            low, high = re_mod.default_bkg_roi(peak, n_y=N_Y, gap=5, width=10)
        except ValueError:
            continue  # refused for want of room, which is the other contract
        assert 0 <= low <= high <= N_Y - 1, f"off-detector band {(low, high)} for peak {peak}"
        assert high < peak[0] or low > peak[1], f"band {(low, high)} overlaps peak {peak}"


def test_default_bkg_roi_refuses_a_peak_that_leaves_no_room():
    with pytest.raises(ValueError, match="no room"):
        re_mod.default_bkg_roi((0, N_Y - 1), n_y=N_Y, gap=5, width=10)


# -- Qt-free (VR-2 / VR-4) --------------------------------------------------


def test_the_module_imports_without_any_gui_package():
    """The whole point of splitting this out of the tab.

    Run in a subprocess with the Qt bindings blocked at import, because this
    process has already imported them — asserting on `sys.modules` here would
    pass for a module that imports Qt eagerly, which is precisely the
    regression worth catching.
    """
    import subprocess
    import sys
    import textwrap

    program = textwrap.dedent(
        """
        import sys

        class _Blocker:
            def find_module(self, name, path=None):
                if name.split(".")[0] in ("qtpy", "PyQt5", "PyQt6", "PySide2", "PySide6"):
                    raise ImportError(f"{name} must not be imported by roi_estimate")
                return None

        sys.meta_path.insert(0, _Blocker())
        import lr_reduction.roi_estimate  # noqa: F401
        print("OK")
        """
    )
    proc = subprocess.run([sys.executable, "-c", program], capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr
    assert "OK" in proc.stdout


def test_the_geometry_comes_from_the_instrument_database_not_a_literal(monkeypatch):
    """256x304 at 15.75 m was true for part of the instrument's life only.

    Asserting `== (256, 304)` cannot show this: `settings.json` currently holds
    exactly one entry for each pixel count, so the database value and the
    literal agree, and a hard-coded `return 256, 304` passed this test. Measured
    — it was mutation row 5, and it survived.

    So the test asserts PROVENANCE instead of value: move the database and the
    answer must move with it. That is the property the slug actually needs,
    because the failure being prevented is a future geometry change the literal
    would not follow.
    """
    from lr_reduction import nr_tools

    assert re_mod.detector_shape("2025-03-04T11:22:33-05:00") == (N_X, N_Y)

    real = nr_tools.read_settings

    def moved(time):
        settings = dict(real(time))
        settings["num_x_pixels"] = 512
        settings["num_y_pixels"] = 608
        return settings

    monkeypatch.setattr(re_mod.nr_tools, "read_settings", moved)
    assert re_mod.detector_shape("2025-03-04T11:22:33-05:00") == (512, 608)


# -- v2: the five demonstrated-harm fixes -----------------------------------


def test_a_sparse_profile_is_refused_not_scored_as_infinite_contrast():
    """A1: `baseline == 0` made the contrast guard a no-op.

    `inf < min_contrast` is False for EVERY threshold, so a zero baseline passed
    unconditionally — and a zero baseline is what a sparse, low-flux run has.
    Reachable on 5 of 63 real REF_L files (70-77% empty rows); those are exactly
    the runs the guard exists for. A profile of zeros with one 3-count row
    returned a confident ROI with infinite contrast.
    """
    sparse = np.zeros(N_Y)
    sparse[200] = 3.0
    with pytest.raises(ValueError, match="baseline|contrast"):
        re_mod.estimate_peak_range(sparse)


def test_a_zero_proton_charge_is_refused_rather_than_bracketing_pixel_zero(tmp_path):
    """A2: `y_tof /= 0` walks every guard and returns (0, 0).

    Empty rows become NaN and occupied rows inf, so `any(counts > 0)` is True on
    the infs; argmax finds the first NaN at index 0; peak is NaN so both walks
    stop; baseline is NaN and `nan > 0` is False, so contrast is inf and that
    guard passes too. The result is (0, 0) — verbatim the failure the no-counts
    guard says it prevents. `get_deadtime_correction` already treats zero charge
    as a real state.
    """
    path = _write_nexus(tmp_path / "zero_pc.nxs.h5", proton_charge=0.0)
    with pytest.raises(ValueError, match="proton charge|no counts"):
        counts = re_mod.counts_vs_y(path, lowres=(100, 160))
        re_mod.estimate_peak_range(counts)


@pytest.mark.parametrize(
    "peak", [(400, 410), (-50, -40), (300, 320), (-5, 5)], ids=["past-end", "negative", "straddles-end", "straddles-zero"]
)
def test_default_bkg_roi_refuses_a_peak_that_is_not_on_the_detector(peak):
    """D1: each branch checks only ONE end, so the other end went unchecked.

    Measured before the fix: `(400, 410)` on a 304-row detector returned
    `(385, 394)` and `(-50, -40)` returned `(-34, -25)` — both the failure the
    docstring says it prevents. I had removed the clamps as "dead code"; they
    were dead only under an unstated precondition (that `peak_range` is on the
    detector) which the function never checked and my sweep never violated.

    Refusing rather than clamping, deliberately: clamping `(400, 410)` yields a
    band from the wrong end of the detector, silently. `RB_Ymin`/`RB_Ymax` reach
    the resolver from layer (c) as unvalidated file input, so this is reachable.
    """
    with pytest.raises(ValueError, match="detector"):
        re_mod.default_bkg_roi(peak, n_y=N_Y, gap=5, width=10)


def test_the_lambda_range_and_the_tof_band_compose(nexus):
    """C1: the two functions disagreed on units and the composition was silent.

    `chopper_tof_window` was named for TOF and returned Angstrom, while
    `counts_vs_y(tof_band=...)` is microseconds. The natural composition
    filtered events to 2.4-5.8 us, produced an all-zero profile, and made
    `estimate_peak_range` blame the run for a unit error.
    """
    lam = re_mod.chopper_lambda_range(nexus)
    meta = re_mod.read_nexus_metadata(nexus)
    band = re_mod.lambda_to_tof(lam, meta["start_time"])

    assert band[1] > band[0]
    assert band[0] > 1000.0, "a TOF band in microseconds, not Angstrom"

    counts = re_mod.counts_vs_y(nexus, lowres=(100, 160), tof_band=band)
    assert counts.sum() > 0, "the composed band selected no events"


def test_counts_vs_y_refuses_an_inverted_band(nexus):
    """C3: nothing pinned the `hi <= lo` refusal.

    Remove it and `np.linspace(hi, lo, n)` descends, `d_tof` goes negative,
    `np.digitize` clips into valid indices, and the function returns a plausible
    profile computed from inverted bins with no error at all.
    """
    with pytest.raises(ValueError, match="band"):
        re_mod.counts_vs_y(nexus, lowres=(100, 160), tof_band=(40000.0, 10000.0))


def test_the_mutation_battery_refuses_a_dirty_baseline():
    """G: the committed battery adopted whatever was on disk as "clean".

    It took `sha(MOD)` of the working tree with zero git references, so a
    leftover from a killed run became the baseline and the battery printed
    `restored: OK` with the mutation still in the file. This is the campaign's
    reference battery and it had the exact failure
    `todo-mutation-harness-restore-safety` was written about.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "roi_batt", "plans/scripts/roi_estimate_mutations.py"
    )
    batt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(batt)

    import subprocess

    assert hasattr(batt, "verify_baseline_matches_head"), "no git-backed baseline check"
    assert hasattr(batt, "signal"), "no signal module — a SIGTERM'd run cannot restore"

    original = open(batt.MOD, encoding="utf-8").read()
    rel = os.path.relpath(batt.MOD, batt.REPO)
    head_blob = subprocess.run(
        ["git", "show", f"HEAD:{rel}"], cwd=batt.REPO, capture_output=True, check=True
    ).stdout.decode()

    try:
        # Matches HEAD -> accepted. Written explicitly rather than relying on the
        # working tree being clean, so the test is deterministic while the slug
        # is mid-edit.
        with open(batt.MOD, "w", encoding="utf-8") as fh:
            fh.write(head_blob)
        batt.verify_baseline_matches_head()

        # Differs from HEAD -> refused. This is the case that matters: a
        # leftover from a killed run must not become the baseline.
        with open(batt.MOD, "w", encoding="utf-8") as fh:
            fh.write(head_blob + "\n# leftover from a killed run\n")
        with pytest.raises(SystemExit, match="dirty|HEAD|baseline"):
            batt.verify_baseline_matches_head()
    finally:
        with open(batt.MOD, "w", encoding="utf-8") as fh:
            fh.write(original)
