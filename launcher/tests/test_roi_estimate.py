"""Tests for the Qt-free ROI/metadata estimation module (T1 slug 1).

The fixture is BUILT, not committed. The plan describes a synthetic
`tests/data/*.nxs.h5`; generating it from a committed builder gives the same
file without putting a binary in git, and the builder is the part a reader needs
in order to know what the fixture actually asserts. It also lets each test vary
one thing — notably the PRESENCE of a log, which amendment 21 requires.
"""

import glob
import os
from pathlib import Path

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
    events=None,
    with_events=True,
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
    if events is not None:
        # roi-popout-data: the caller's own (event_id, tof) arrays, for the cases a random draw cannot pin:
        # off-detector ids, an event order, no events at all.
        event_id = np.asarray(events[0], dtype=np.int64)
        tof = np.asarray(events[1], dtype=float)

    with h5py.File(path, "w") as f:
        entry = f.create_group("entry")
        entry.create_dataset("title", data=[title.encode()])
        entry.create_dataset("start_time", data=[start_time.encode()])
        entry.create_dataset("run_number", data=[str(run_number).encode()])
        # The TOTAL accumulated charge, which is what `load_and_extract` passes
        # to `get_y_tof` as `pcharge` — not the DASlogs time series beside it.
        entry.create_dataset("proton_charge", data=np.array([proton_charge]))

        if with_events:  # roi-popout-data v2, K3: a file without the group
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
    """roi-popout-data B7: a band on each side, in the reducer's four-bound form (F5)."""
    b0, b1, b2, b3 = re_mod.default_bkg_roi((140, 160), n_y=N_Y, gap=5, width=10)
    assert b1 < 140 - 5 and b2 > 160 + 5


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
            b0, b1, b2, b3 = re_mod.default_bkg_roi(peak, n_y=N_Y, gap=5, width=10)
        except ValueError:
            continue  # refused for want of room, which is the other contract
        assert 0 < b0 <= b1 < peak[0] and peak[1] < b2 <= b3 <= N_Y - 1, f"bands {(b0, b1, b2, b3)} for peak {peak}"


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

    # roi-popout-data I8 (F12): asserted on sys.modules in a fresh interpreter, as test_settings_document.py
    # does, not through a find_module finder. Python 3.12 dropped that hook, so there the guard would pass
    # vacuously (inferred from the language's removal notice; not run: this environment is Python 3.11).
    program = textwrap.dedent(
        """
        import sys
        import lr_reduction.roi_estimate
        bindings = ("qtpy", "PyQt5", "PyQt6", "PySide2", "PySide6")
        print(lr_reduction.roi_estimate.__file__)
        print(sorted(m for m in sys.modules if m.split(".")[0] in bindings))
        """
    )
    # roi-popout-data v2, K6 (A1): the child imports this checkout's module, not whichever lr_reduction the
    # environment's editable install points at (from another checkout that would test another tree).
    source = os.path.dirname(os.path.dirname(os.path.abspath(re_mod.__file__)))
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(filter(None, [source, os.environ.get("PYTHONPATH")]))}
    proc = subprocess.run([sys.executable, "-c", program], capture_output=True, text=True, timeout=120, env=env)
    assert proc.returncode == 0, proc.stderr
    imported, bindings = proc.stdout.strip().splitlines()[-2:]
    assert os.path.realpath(imported) == os.path.realpath(re_mod.__file__), imported
    assert bindings == "[]", proc.stdout


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
    with pytest.raises(ValueError, match="proton charge"):
        re_mod.counts_vs_y(path, lowres=(100, 160))


def test_a_non_finite_profile_is_refused_even_if_it_reaches_the_estimator():
    """The SECOND half of A2, pinned independently of the first.

    There are two guards on this path — `counts_vs_y` refuses a zero charge at
    source, and `estimate_peak_range` refuses a non-finite profile downstream.
    A single test that composed them left BOTH mutations green, because each
    guard covered for the other: drop either one and the survivor still raised.
    Measured as battery rows 12 and 13, both surviving.

    So this one hands the estimator the profile directly. A zero-charge divide
    yields NaN on empty rows and inf on occupied ones; unguarded, `any(counts >
    0)` passes on the infs, argmax finds the first NaN at index 0, both walks
    stop, and the function returns (0, 0).
    """
    profile = np.full(N_Y, np.nan)
    profile[150] = np.inf
    with pytest.raises(ValueError, match="NaN|inf|non-finite"):
        re_mod.estimate_peak_range(profile)


@pytest.mark.parametrize(
    "peak",
    [(400, 410), (-50, -40), (300, 320), (-5, 5), (200, 100)],
    ids=["past-end", "negative", "straddles-end", "straddles-zero", "reversed"],
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

    roi-popout-data: B7's two-sided band refuses every off-detector peak above through its no-room check as
    well, whose message also says "detector", so deleting D1 left this test green (battery row 14 SURVIVED).
    The test now matches D1's own words, and adds the case only D1 refuses: a reversed peak, for which the
    no-room check returns (192, 196, 104, 108).
    """
    with pytest.raises(ValueError, match="not an ascending range"):
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


# ===========================================================================
# roi-popout-data (plan v1, triage 5fa056b). PR #31's advisories I1-I8, and the pop-out's data contract
# B1-B9 on this module. Counts are raw event counts (the web report's quantity). Pixel ranges are inclusive
# [low, high], as the reducer's lowres and peak masks are.
# ===========================================================================

_BINDINGS_ROOT = Path(__file__).resolve().parents[3]
_BATTERY = _BINDINGS_ROOT / "scripts" / "test" / "roi_estimate_mutations.py"


def _load_battery():
    """The committed battery, by the path anchored to this file (the gate runs from tests/)."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("roi_batt", str(_BATTERY))
    batt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(batt)
    return batt


def _events_file(tmp_path, ids, tofs, name="ev.nxs.h5", **kwargs):
    return _write_nexus(tmp_path / name, events=(ids, tofs), **kwargs)


def _ids(x, y):
    return np.asarray(x, dtype=np.int64) * N_Y + np.asarray(y, dtype=np.int64)


# -- I1-I3: counts_vs_y ------------------------------------------------------


def _spy_get_y_tof(monkeypatch):
    from lr_reduction import binary_processing

    calls = []
    real = binary_processing.get_y_tof

    def spy(tof_array, event_id, e_offset, lowres, *args, **kwargs):
        calls.append(list(lowres))
        return real(tof_array, event_id, e_offset, lowres, *args, **kwargs)

    monkeypatch.setattr(binary_processing, "get_y_tof", spy)
    return calls


def test_counts_vs_y_default_lowres_follows_the_detector_database(nexus, monkeypatch):
    """I1 (L3): the default means every X pixel of this run's detector. A literal (0, 255) agrees with the
    database today, so the pin moves the database: the default must follow it."""
    calls = _spy_get_y_tof(monkeypatch)
    re_mod.counts_vs_y(nexus)
    assert calls[-1] == [0, N_X - 1]

    real = re_mod.nr_tools.read_settings

    def moved(time):
        settings = dict(real(time))
        settings["num_x_pixels"] = 512
        settings["num_y_pixels"] = 608
        return settings

    monkeypatch.setattr(re_mod.nr_tools, "read_settings", moved)
    re_mod.counts_vs_y(nexus)
    assert calls[-1] == [0, 511]


def test_counts_vs_y_with_a_band_counts_exactly_the_events_inside_it(tmp_path):
    """I2: a narrow band selects strictly fewer counts than no band, and exactly the in-band events (charge 1)."""
    ids = _ids([120] * 4000, [150] * 4000)
    tofs = np.r_[np.full(1000, 15000.0), np.full(3000, 30000.0)]
    path = _events_file(tmp_path, ids, tofs, proton_charge=1.0)
    total = re_mod.counts_vs_y(path, lowres=(100, 160)).sum()
    in_band = re_mod.counts_vs_y(path, lowres=(100, 160), tof_band=(14000.0, 16000.0)).sum()
    assert in_band < total
    assert in_band == pytest.approx(1000)


def test_counts_vs_y_calls_the_histogrammer_once_with_or_without_a_band(nexus, monkeypatch):
    """I3: one get_y_tof call per counts_vs_y call (the band is applied to the events first)."""
    calls = _spy_get_y_tof(monkeypatch)
    re_mod.counts_vs_y(nexus, lowres=(100, 160))
    assert len(calls) == 1
    re_mod.counts_vs_y(nexus, lowres=(100, 160), tof_band=(12000.0, 30000.0))
    assert len(calls) == 2


# -- I5-I7: the battery ------------------------------------------------------


def test_the_battery_table_lists_the_rows_it_runs():
    """I5: the battery's documented table and MUTATIONS name the same rows. A row in the table that no longer
    runs is in RETIRED with its reason, and every anchor occurs exactly once in the module (no ANCHOR MISS)."""
    import re

    batt = _load_battery()
    source = _BATTERY.read_text(encoding="utf-8")
    table = {int(m) for m in re.findall(r"^# \|\s*(\d+)\s*\|", source, flags=re.MULTILINE)}
    running = {row for row, *_ in batt.MUTATIONS}
    retired = getattr(batt, "RETIRED", {})
    assert table == running | set(retired), (sorted(table), sorted(running), sorted(retired))
    assert all(isinstance(reason, str) and reason.strip() for reason in retired.values())
    module = open(batt.MOD, encoding="utf-8").read()
    assert [row for row, _desc, old, _new in batt.MUTATIONS if module.count(old) != 1] == []


def _battery_copy(tmp_path, slow_test=None):
    """A git repository in tmp_path holding copies of the module and the battery (I6, I7)."""
    import shutil
    import subprocess

    repo = tmp_path / "repo"
    (repo / "src" / "lr_reduction").mkdir(parents=True)
    (repo / "scripts" / "test").mkdir(parents=True)
    (repo / "tests" / "unit" / "lr_reduction").mkdir(parents=True)
    shutil.copy(re_mod.__file__, repo / "src" / "lr_reduction" / "roi_estimate.py")
    shutil.copy(_BATTERY, repo / "scripts" / "test" / "roi_estimate_mutations.py")
    (repo / "tests" / "unit" / "lr_reduction" / "test_roi_estimate.py").write_text(
        slow_test or "def test_nothing():\n    pass\n")
    git = ["git", "-c", "user.name=test", "-c", "user.email=test@example.invalid"]
    subprocess.run([*git, "init", "-q"], cwd=repo, check=True)
    subprocess.run([*git, "add", "-A"], cwd=repo, check=True)
    subprocess.run([*git, "commit", "-q", "-m", "copy"], cwd=repo, check=True)
    return repo


def test_the_mutation_battery_refuses_a_dirty_baseline(tmp_path, monkeypatch):
    """G (v2), and I6: on a tmp_path copy, never on the tracked module. The battery adopted whatever was on
    disk as "clean", so a leftover from a killed run became the baseline. It must compare against HEAD. The
    tracked module's mtime is unchanged across this test."""
    tracked = Path(re_mod.__file__)
    mtime = os.stat(tracked).st_mtime_ns
    repo = _battery_copy(tmp_path)
    batt = _load_battery()
    module = repo / "src" / "lr_reduction" / "roi_estimate.py"
    monkeypatch.setattr(batt, "REPO", str(repo))
    monkeypatch.setattr(batt, "MOD", str(module))
    batt.verify_baseline_matches_head()  # matches HEAD: accepted
    module.write_text(module.read_text(encoding="utf-8") + "\n# leftover from a killed run\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="dirty|HEAD|baseline"):
        batt.verify_baseline_matches_head()
    assert os.stat(tracked).st_mtime_ns == mtime


def test_a_sigterm_during_the_battery_leaves_the_module_restored(tmp_path):
    """I7: the behaviour the hasattr(batt, "signal") stand-in only implied. On a tmp_path copy, the battery is
    stopped by SIGTERM while a mutation is in the file (its pytest is sleeping), and the module must be as
    committed afterwards."""
    import signal
    import subprocess
    import sys
    import time

    repo = _battery_copy(tmp_path, slow_test="import time\n\n\ndef test_slow():\n    time.sleep(60)\n")
    module = repo / "src" / "lr_reduction" / "roi_estimate.py"
    original = module.read_text(encoding="utf-8")
    proc = subprocess.Popen([sys.executable, str(repo / "scripts" / "test" / "roi_estimate_mutations.py")],
                            cwd=repo, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
    try:
        deadline = time.monotonic() + 60
        while module.read_text(encoding="utf-8") == original and time.monotonic() < deadline:
            time.sleep(0.05)
        assert module.read_text(encoding="utf-8") != original, "the battery never wrote a mutation"
        os.kill(proc.pid, signal.SIGTERM)
        assert proc.wait(timeout=30) == 128 + signal.SIGTERM
        assert module.read_text(encoding="utf-8") == original
    finally:
        try:
            os.killpg(proc.pid, signal.SIGKILL)  # the battery's sleeping pytest, left behind by its exit
        except ProcessLookupError:
            pass


# -- B1, B8: one read of the file, into RunEvents ----------------------------


def test_load_event_pixels_holds_the_events_and_the_detector_shape(tmp_path):
    """B1: x = id // n_y, y = id % n_y (get_y_tof's packing); n_x, n_y from the database; stride 1."""
    ids = _ids([100, 120, 159], [150, 3, 300])
    path = _events_file(tmp_path, ids, [11000.0, 12000.0, 13000.0])
    events = re_mod.load_event_pixels(path)
    assert (events.n_x, events.n_y) == (N_X, N_Y)
    assert list(events.x) == [100, 120, 159] and list(events.y) == [150, 3, 300]
    assert list(events.tof) == [11000.0, 12000.0, 13000.0]
    assert events.stride == 1 and events.n_off_detector == 0


def test_load_event_pixels_reads_the_file_once(nexus, monkeypatch):
    """B1 and acceptance 4: one read of the file per load (an h5py.File spy)."""
    opened = []
    real = re_mod.h5py.File

    def spy(*args, **kwargs):
        opened.append(args[0])
        return real(*args, **kwargs)

    monkeypatch.setattr(re_mod.h5py, "File", spy)
    re_mod.load_event_pixels(nexus)
    assert len(opened) == 1


def test_the_event_geometry_follows_the_instrument_database(nexus, monkeypatch):
    """B8: no geometry literal. Move the database, and the shapes follow."""
    real = re_mod.nr_tools.read_settings

    def moved(time):
        settings = dict(real(time))
        settings["num_x_pixels"] = 512
        settings["num_y_pixels"] = 608
        return settings

    monkeypatch.setattr(re_mod.nr_tools, "read_settings", moved)
    events = re_mod.load_event_pixels(nexus)
    assert (events.n_x, events.n_y) == (512, 608)
    assert re_mod.xy_image(events).shape == (608, 512)


# -- B2: the XY image --------------------------------------------------------


def _report_xy(path, n_x, n_y):
    """The web report's XY array for a run (web_report.py:579-583): Integration of LoadEventNexus."""
    from mantid.simpleapi import Integration, LoadEventNexus

    workspace = LoadEventNexus(Filename=path, OutputWorkspace="xy_cross_check")
    signal = Integration(InputWorkspace=workspace, OutputWorkspace="xy_cross_check_sum").extractY()
    return np.reshape(signal, (n_x, n_y)).T


def _extreme_tof_pixels(events):
    """``{(y, x): count}`` of the events at the run's minimum and maximum TOF."""
    pixels = {}
    at = (events.tof == events.tof.min()) | (events.tof == events.tof.max())
    for y, x in zip(events.y[at], events.x[at]):
        pixels[(int(y), int(x))] = pixels.get((int(y), int(x)), 0) + 1
    return pixels


def _differing(image, report):
    """``{(y, x): image - report}`` where the two differ."""
    difference = image - report
    return {(int(y), int(x)): int(difference[y, x]) for y, x in zip(*np.nonzero(difference))}


def test_xy_image_is_the_web_reports_array_but_for_the_extreme_tof_events(nexus_dir):
    """T1 (F4, corrected at v2: N-1). The report's array is the event image less the event(s) at the run's minimum
    and/or maximum TOF, which Mantid's Integration drops when its default range excludes them. 201288 agrees
    everywhere; 179932 differs at exactly its two extreme events' pixels, (y, x) = (72, 14) and (191, 217), by one
    count each (the Integrator's measurement, reproduced)."""
    agree = os.path.join(nexus_dir, "REF_L_201288.nxs.h5")
    events = re_mod.load_event_pixels(agree)
    assert _differing(re_mod.xy_image(events), _report_xy(agree, events.n_x, events.n_y)) == {}

    differ = os.path.join(nexus_dir, "REF_L_179932.nxs.h5")
    events = re_mod.load_event_pixels(differ)
    image, report = re_mod.xy_image(events), _report_xy(differ, events.n_x, events.n_y)
    extreme = _extreme_tof_pixels(events)
    assert _differing(image, report) == extreme == {(72, 14): 1, (191, 217): 1}
    assert image.sum() - report.sum() == sum(extreme.values())


def test_xy_image_against_the_web_report_over_every_fixture_run(nexus_dir):
    """T1b (v2, the census v1 lacked; F4's population claim, not one sample). Over every REF_L run in the fixture
    data, xy_image minus the report's array is non-negative, non-zero only at the extreme-TOF events' pixels, and
    at most their count there. Measured at v2: 63 runs, 18 equal, 45 differing, no exception, ~45 s."""
    pytest.importorskip("mantid.simpleapi", reason="the report's array is Mantid's Integration; no Mantid here")
    paths = sorted(glob.glob(os.path.join(nexus_dir, "REF_L_*.nxs.h5")))
    assert len(paths) >= 60, "precondition: the fixture data is present"
    differing_runs = 0
    for path in paths:
        events = re_mod.load_event_pixels(path)
        image = re_mod.xy_image(events)
        differing = _differing(image, _report_xy(path, events.n_x, events.n_y))
        extreme = _extreme_tof_pixels(events)
        assert all(0 < count <= extreme.get(pixel, 0) for pixel, count in differing.items()), (path, differing)
        differing_runs += bool(differing)
    assert differing_runs > 0, "precondition: some run differs, or the census proves nothing about the drop"


def test_xy_image_puts_the_injected_peak_at_its_row_and_columns(nexus):
    """T2: image[y, x], shape (n_y, n_x); the builder's peak is centred on y = 150 at x 100-159.

    The builder truncates N(150, 1) to ints, so the peak's centre sits on the 149/150 boundary and the two rows
    share it: measured, 20647 and 20330 of 72000 counts. Row 150 alone never "dominates", as the plan had it.
    Asserted instead: rows 148-152 hold most of the counts, their centroid is 149.5 (a one-row shift moves it
    by 1), and nothing lands outside x 100-159.
    """
    image = re_mod.xy_image(re_mod.load_event_pixels(nexus))
    assert image.shape == (N_Y, N_X)
    assert image[148:153, 100:160].sum() > 0.75 * image.sum()
    rows = np.arange(145, 156)
    profile = image[rows].sum(axis=1)
    assert abs((profile * rows).sum() / profile.sum() - 149.5) < 0.1
    assert image[:, :100].sum() == 0 and image[:, 160:].sum() == 0


def test_off_detector_ids_are_dropped_and_counted(tmp_path):
    """T3: ids outside [0, n_x*n_y) are dropped and counted, never reshaped into an edge pixel."""
    ids = np.r_[_ids([120] * 10, [150] * 10), [N_X * N_Y, N_X * N_Y + 5, -1]]
    path = _events_file(tmp_path, ids, np.full(13, 20000.0))
    events = re_mod.load_event_pixels(path)
    assert events.n_off_detector == 3
    assert re_mod.xy_image(events).sum() == 10


# -- B3, B4: TOF edges and the Y-TOF image -----------------------------------


def test_tof_edges_span_every_event_not_the_chopper_band(nexus):
    """T4 (R14): the edges cover the events' full TOF span (50 us bins, the web report's); the chopper band is
    an overlay, not a crop. An event outside the band is inside the edges and counted."""
    events = re_mod.load_event_pixels(nexus)
    band = re_mod.lambda_to_tof(re_mod.chopper_lambda_range(nexus), re_mod.read_nexus_metadata(nexus)["start_time"])
    outside = (events.tof < band[0]) | (events.tof > band[1])
    assert outside.any(), "precondition: the builder has events outside the chopper band"
    edges = re_mod.tof_edges(events)
    assert edges[0] <= events.tof.min() and edges[-1] >= events.tof.max()
    assert np.allclose(np.diff(edges), 50.0)
    image = re_mod.y_tof_image(events, (0, N_X - 1), edges)
    assert image.shape == (N_Y, len(edges) - 1)
    assert image.sum() == len(events.tof)


def test_y_tof_image_counts_only_the_x_range(nexus):
    """T5: events outside the inclusive x range contribute nothing; inside, every event is counted."""
    events = re_mod.load_event_pixels(nexus)
    edges = re_mod.tof_edges(events)
    assert re_mod.y_tof_image(events, (0, 99), edges).sum() == 0
    assert re_mod.y_tof_image(events, (100, 159), edges).sum() == len(events.tof)
    inside = (events.x >= 100) & (events.x <= 129)
    assert re_mod.y_tof_image(events, (100, 129), edges).sum() == inside.sum()


def test_y_tof_image_keeps_the_event_at_the_last_edge(tmp_path):
    """T5 (F4): an event at exactly the last edge is counted. Mantid's half-open last bin drops it."""
    path = _events_file(tmp_path, _ids([120] * 3, [150] * 3), [1000.0, 1050.0, 1100.0])
    events = re_mod.load_event_pixels(path)
    image = re_mod.y_tof_image(events, (0, N_X - 1), np.array([1000.0, 1050.0, 1100.0]))
    assert image.sum() == 3 and image[150, 1] == 2


# -- B5: the three profiles --------------------------------------------------


def test_profile_y_agrees_with_the_library_histogrammer(nexus_dir):
    """T6 (F8): profile_y over the run's proton charge equals counts_vs_y, through get_y_tof."""
    path = os.path.join(nexus_dir, "REF_L_201288.nxs.h5")
    with h5py.File(path, "r") as f:
        charge = float(np.sum(f["entry/proton_charge"][:]))
    events = re_mod.load_event_pixels(path)
    expected = re_mod.counts_vs_y(path, lowres=(50, 200))
    assert np.allclose(re_mod.profile_y(events, (50, 200)) / charge, expected, rtol=1e-12, atol=0)


def test_profiles_are_marginals_of_the_images(nexus):
    """T7: profile_y is the Y-TOF image summed over TOF; profile_x is the XY image summed over Y; profile_tof
    counts every event. Each restriction is pinned the same way: profile_x's rows and band, xy_image's band,
    profile_tof's columns and rows."""
    events = re_mod.load_event_pixels(nexus)
    edges = re_mod.tof_edges(events)
    assert np.array_equal(re_mod.profile_y(events, (100, 159)), re_mod.y_tof_image(events, (100, 159), edges).sum(axis=1))
    assert np.array_equal(re_mod.profile_x(events), re_mod.xy_image(events).sum(axis=0))
    assert re_mod.profile_tof(events, edges).sum() == len(events.tof)
    band = (15000.0, 25000.0)
    in_band = (events.tof >= band[0]) & (events.tof <= band[1])
    assert re_mod.profile_x(events, tof_band=band).sum() == in_band.sum()
    assert re_mod.profile_y(events, (100, 159), tof_band=band).sum() == in_band.sum()
    assert re_mod.xy_image(events, tof_band=band).sum() == in_band.sum()
    assert np.array_equal(re_mod.profile_x(events, tof_band=band), re_mod.xy_image(events, tof_band=band).sum(axis=0))
    assert np.array_equal(re_mod.profile_x(events, y_range=(148, 152)), re_mod.xy_image(events)[148:153].sum(axis=0))
    assert np.array_equal(re_mod.profile_tof(events, edges, x_range=(100, 129), y_range=(148, 152)),
                          re_mod.y_tof_image(events, (100, 129), edges)[148:153].sum(axis=0))


def test_tof_edges_give_one_bin_when_every_event_shares_a_tof(tmp_path):
    """B3: a span of zero is one bin, [t, t + bin_width], never a single edge that bins nothing."""
    path = _events_file(tmp_path, _ids([120] * 3, [150] * 3), np.full(3, 20000.0))
    events = re_mod.load_event_pixels(path)
    edges = re_mod.tof_edges(events)
    assert np.array_equal(edges, [20000.0, 20050.0])
    assert re_mod.y_tof_image(events, (0, N_X - 1), edges).sum() == 3


def test_tof_edges_hold_the_latest_event_when_rounding_falls_short(tmp_path):
    """B3: this span is 1532 bins of 50 us plus one rounding ulp. ceil((high - low) / 50) is 1532, and low +
    50 * 1532 lands one ulp short of the latest event (a measured float64 pair). The edges are extended, so
    the event is binned, not dropped past the last edge."""
    tofs = np.array([34429.573096232525, 111029.57309623253])
    assert tofs[0] + 50.0 * 1532 < tofs[1], "precondition: the pair still falls short"
    path = _events_file(tmp_path, _ids([120, 120], [150, 150]), tofs)
    events = re_mod.load_event_pixels(path)
    edges = re_mod.tof_edges(events)
    assert edges[-1] >= tofs[1]
    assert re_mod.y_tof_image(events, (0, N_X - 1), edges).sum() == 2


# -- B6: the reducer's background bands --------------------------------------

_ACCEPTED = [[133, 149, 0, 0], [120, 130, 150, 160], [160, 150, 130, 120], [0, 0, 10, 300], [133.5, 149.5, 0, 0]]
_REFUSED = [  # (entry, the words of its own refusal) -- v2, K4: each reason matched, not only "background"
    ([0, 10, 150, 160], "sentinel"), ([0, 0, 0, 0], "sentinel"), ([0, 140, 150, 160], "sentinel"),
    ([0, 0, 0, 150], "sentinel"),
    ([121, 130], "needs four bounds"), ([133, 149, 0], "needs four bounds"), ([], "needs four bounds"),
    (None, "no background is set"), ("120, 130", "flat list"),
    ([float("nan"), 149, 0, 0], "finite numbers"), ([133, float("inf"), 0, 0], "finite numbers"),
    ([[133, 149, 0, 0], [120, 130, 150, 160]], "flat list"), (["120", "130", "150", "160"], "finite numbers"),
    ([[120, 130, 150, 160]] * 4, "flat list"), ([[133, 149], [0, 0, 0]], "one angle's four pixel bounds, not"),
]


@pytest.mark.parametrize("bkg", _ACCEPTED)
def test_background_bands_are_the_rows_the_reducer_averages(bkg):
    """T8 (F6): point-wise equal to the reducer's _background_roi_sorter where it returns four bounds, values
    and types: a fractional bound stays as the sorter keeps it (the reducer's mask then starts at the next row),
    never truncated to an int."""
    from lr_reduction.nr_reduction_calc import NR_Reduction

    expected = NR_Reduction._background_roi_sorter(None, bkg, 136, 146).tolist()
    (b0, b1), (b2, b3) = re_mod.background_bands(bkg, 136, 146)
    assert [b0, b1, b2, b3] == expected
    # v2, K1: "values and types" -- 136 == 136.0, so == alone passes a float for an int.
    assert [type(v) for v in (b0, b1, b2, b3)] == [type(v) for v in expected]


def _against_the_sorter(bkg, y_min, y_max):
    """None when background_bands agrees with the sorter on ``bkg`` (values and types; a refusal exactly where the
    sorter returns no four bounds), else a description of the disagreement."""
    from lr_reduction.nr_reduction_calc import NR_Reduction

    sorted_ = NR_Reduction._background_roi_sorter(None, list(bkg), y_min, y_max)
    try:
        got = [v for band in re_mod.background_bands(list(bkg), y_min, y_max) for v in band]
    except ValueError as exc:
        if type(exc) is not ValueError:
            return f"{bkg}: {type(exc).__name__}, not ValueError"
        return None if sorted_ is None or len(sorted_) != 4 else f"{bkg}: refused, the sorter gives {sorted_}"
    if sorted_ is None or len(sorted_) != 4:
        return f"{bkg}: {got}, the sorter gives {sorted_}"
    want = sorted_.tolist()
    return None if got == want and [type(v) for v in got] == [type(v) for v in want] else f"{bkg}: {got} != {want}"


def test_background_bands_is_the_sorter_on_an_exhaustive_grid_and_random_entries():
    """v2, K1 (after the numerical reviewer's grid): every four-bound entry from 0..7 (4096, every zero count, peak
    3-4), and 2000 seeded entries on the detector with zeros mixed in. Python ints come back for int entries, and
    a refusal where the sorter returns None."""
    import itertools

    misses = [m for bkg in itertools.product(range(8), repeat=4) if (m := _against_the_sorter(bkg, 3, 4))]
    rng = np.random.default_rng(2026)
    for _ in range(2000):
        bkg = [int(v) for v in rng.integers(0, N_Y, 4)]
        for i in rng.choice(4, size=int(rng.integers(0, 3)), replace=False):
            bkg[i] = 0
        low = int(rng.integers(1, N_Y - 20))
        if (m := _against_the_sorter(bkg, low, low + int(rng.integers(0, 15)))):
            misses.append(m)
    assert misses == [], misses[:5]


@pytest.mark.parametrize("bkg, reason", _REFUSED)
def test_background_bands_refuses_what_the_reducer_cannot_use(bkg, reason):
    """T8 (F6): one, three or four zeros (pixel 0 is the reducer's sentinel), a length other than four, no
    entry, not numbers: ValueError naming the background, never None or two bounds. So is a bound that is not
    finite, which the sorter passes through: a NaN band selects no rows and its centre is NaN, and an infinite
    one's centre is infinite. So is the per-angle list of entries in place of one angle's entry.

    Three cases each reach one guard that nothing else catches (battery rows 33, 34, 36 SURVIVED without them):
    - four numbers as text: no numeric check means np.isfinite raises TypeError;
    - a four-angle list without zeros: shape (4, 4) passes the length and zero checks, and without the ndim
      check four lists come back as "bounds";
    - a ragged entry: numpy's own ValueError does not name the background.

    v2: each refusal is matched by its own reason (K4: a deleted None guard let the ndim refusal answer with the
    wrong one), and is exactly ValueError (K2: CannotEstimateError subclasses it)."""
    with pytest.raises(ValueError, match="background") as refused:
        re_mod.background_bands(bkg, 136, 146)
    assert refused.type is ValueError and reason in str(refused.value), refused.value


# -- B7: a default background the reducer accepts ----------------------------


def test_default_bkg_roi_survives_the_reducer():
    """T9 (F5): swept over every peak position, the result is four ascending ints, inside the detector and
    never 0 (the sentinel). The reducer's sorter returns it unchanged, and a side with no room is refused."""
    from lr_reduction.nr_reduction_calc import NR_Reduction

    gap, width = 3, 5
    for peak_low in range(0, N_Y - 1):
        peak = (peak_low, min(peak_low + 8, N_Y - 1))
        room = peak[0] - gap - width >= 1 and peak[1] + gap + width <= N_Y - 1
        if not room:
            with pytest.raises(ValueError, match="room"):
                re_mod.default_bkg_roi(peak, n_y=N_Y, gap=gap, width=width)
            continue
        bounds = re_mod.default_bkg_roi(peak, n_y=N_Y, gap=gap, width=width)
        assert len(bounds) == 4 and all(isinstance(b, int) for b in bounds)
        assert 0 < bounds[0] <= bounds[1] < peak[0] and peak[1] < bounds[2] <= bounds[3] <= N_Y - 1
        assert [int(v) for v in NR_Reduction._background_roi_sorter(None, list(bounds), *peak)] == list(bounds)


def test_default_bkg_roi_defaults_are_the_reviewed_three_and_five():
    """A3: gap 3 and width 5 either side, #197's values that the scientists reviewed."""
    assert re_mod.default_bkg_roi((140, 160), n_y=N_Y) == (132, 136, 164, 168)


def test_default_bkg_roi_takes_a_one_row_peak_and_refuses_what_is_not_a_whole_row():
    """v2, K7 (A4, the numerical advisory). A one-row peak is legitimate: four ints around it. A fractional peak
    edge, gap or width, or a bool, is refused rather than truncated: the result promises whole rows. An integral
    float (a JSON 150.0) is a whole row."""
    assert re_mod.default_bkg_roi((150, 150), n_y=N_Y) == (142, 146, 154, 158)
    bounds = re_mod.default_bkg_roi((150.0, 160.0), n_y=N_Y)
    assert bounds == (142, 146, 164, 168) and all(type(b) is int for b in bounds)
    for kwargs in ({"peak_range": (150.5, 160)}, {"peak_range": (150, 160.5)}, {"gap": 2.5}, {"width": 4.5},
                   {"gap": True}, {"width": True}, {"peak_range": (True, 160)}):
        call = {"peak_range": (150, 160), "n_y": N_Y, **kwargs}
        with pytest.raises(ValueError, match="whole") as refused:
            re_mod.default_bkg_roi(**call)
        assert refused.type is ValueError, kwargs


# -- B9, T10-T12: refusals, sparse and empty runs, sampling -------------------


def test_a_sparse_real_run_gives_images_and_a_refused_estimate(nexus_dir):
    """T10 (F8): 201284 has ~6 000 events. The images and profiles are returned, and the estimator refuses."""
    events = re_mod.load_event_pixels(os.path.join(nexus_dir, "REF_L_201284.nxs.h5"))
    assert re_mod.xy_image(events).sum() > 0
    with pytest.raises(re_mod.CannotEstimateError):
        re_mod.estimate_peak_range(re_mod.profile_y(events, (50, 200)))


def test_stride_sampling_is_not_a_time_slice(tmp_path):
    """T11: a peak only in the second half of the event list is present at max_events = n // 4 (stride 4). A
    head slice would hold none of it."""
    rng = np.random.default_rng(7)
    half = 4000
    ids = np.r_[_ids(rng.integers(100, 160, half), rng.integers(0, N_Y, half)), _ids([120] * half, [200] * half)]
    path = _events_file(tmp_path, ids, rng.uniform(10000.0, 40000.0, 2 * half))
    events = re_mod.load_event_pixels(path, max_events=len(ids) // 4)
    assert events.stride == 4
    profile = re_mod.profile_y(events, (100, 159))
    assert int(np.argmax(profile)) == 200


def test_an_empty_run_gives_zero_images_and_refuses_edges(tmp_path):
    """T12 (§3 states): no events gives all-zero arrays of the right shape; tof_edges refuses, with
    CannotEstimateError, never a bare min() error."""
    path = _events_file(tmp_path, np.array([], dtype=np.int64), np.array([], dtype=float))
    events = re_mod.load_event_pixels(path)
    assert len(events.tof) == 0
    # v2, A5: all-zero values, not only shapes, and profile_tof as well.
    zeros = {
        "xy_image": (re_mod.xy_image(events), (N_Y, N_X)),
        "profile_y": (re_mod.profile_y(events, (0, N_X - 1)), (N_Y,)),
        "profile_x": (re_mod.profile_x(events), (N_X,)),
        "y_tof_image": (re_mod.y_tof_image(events, (0, N_X - 1), np.array([0.0, 50.0])), (N_Y, 1)),
        "profile_tof": (re_mod.profile_tof(events, np.array([0.0, 50.0, 100.0])), (2,)),
    }
    for name, (array, shape) in zeros.items():
        assert array.shape == shape and not array.any(), name
    with pytest.raises(re_mod.CannotEstimateError):
        re_mod.tof_edges(events)


def test_called_wrong_is_a_value_error_not_an_empty_answer(nexus):
    """B9: a reversed range or band, a bin width <= 0, max_events <= 0: ValueError, never an empty selection
    reported as "no counts"."""
    events = re_mod.load_event_pixels(nexus)
    for call in (lambda: re_mod.profile_y(events, (160, 100)),
                 lambda: re_mod.xy_image(events, tof_band=(30000.0, 10000.0)),
                 lambda: re_mod.profile_x(events, y_range=(200, 100)),
                 lambda: re_mod.y_tof_image(events, (160, 100), re_mod.tof_edges(events)),
                 lambda: re_mod.tof_edges(events, bin_width=0),
                 lambda: re_mod.load_event_pixels(nexus, max_events=0)):
        with pytest.raises(ValueError) as refused:
            call()
        # v2, K2: exactly ValueError. CannotEstimateError ("looked, nothing to offer") subclasses it, so
        # pytest.raises(ValueError) alone passes a called-wrong site that raises the other.
        assert refused.type is ValueError, refused.value


def test_called_wrong_names_what_is_wrong(nexus):
    """B9: the refusals numpy would not make. One TOF edge bins nothing (numpy returns an empty histogram), and
    an edge repeated is a zero-width bin that counts an event on it; a reversed profile_tof range selects
    nothing; a negative gap or a zero width puts the default band inside the peak or reverses it."""
    events = re_mod.load_event_pixels(nexus)
    edges = re_mod.tof_edges(events)
    for call, match in ((lambda: re_mod.y_tof_image(events, (100, 159), [20000.0]), "edges"),
                        (lambda: re_mod.y_tof_image(events, (100, 159), [20000.0, 20000.0, 20050.0]), "edges"),
                        (lambda: re_mod.profile_tof(events, [20000.0]), "edges"),
                        (lambda: re_mod.profile_tof(events, [20000.0, 20000.0, 20050.0]), "edges"),
                        (lambda: re_mod.profile_tof(events, edges, x_range=(160, 100)), "reversed"),
                        (lambda: re_mod.profile_tof(events, edges, y_range=(200, 100)), "reversed"),
                        (lambda: re_mod.default_bkg_roi((140, 160), n_y=N_Y, gap=-1), "gap"),
                        (lambda: re_mod.default_bkg_roi((140, 160), n_y=N_Y, width=0), "width")):
        with pytest.raises(ValueError, match=match) as refused:
            call()
        assert refused.type is ValueError, refused.value  # v2, K2


# -- roi-popout-data v2: K3, K5 ------------------------------------------------


def test_a_file_without_bank1_events_is_a_key_error_not_an_empty_run(tmp_path):
    """v2, K3 (the plan's types table): a file without the event group is h5py's KeyError, uncaught. It is not
    read as a run with no events, which would draw empty images for a file that is not a run at all."""
    path = _write_nexus(tmp_path / "no_events.nxs.h5", with_events=False)
    with pytest.raises(KeyError):
        re_mod.load_event_pixels(path)


def test_a_range_whose_bounds_are_equal_is_one_pixel_or_one_tof(nexus):
    """v2, K5: low == high is a legal inclusive range, one pixel (or one TOF value), for every range kind; only
    high < low is refused. Each selection counts exactly the events at that pixel or TOF."""
    events = re_mod.load_event_pixels(nexus)
    edges = re_mod.tof_edges(events)
    t = float(events.tof[0])
    at_x, at_y, at_t = events.x == 120, events.y == 150, events.tof == t
    assert at_x.any() and at_y.any() and at_t.any(), "precondition: events at the chosen pixel and TOF"
    assert re_mod.profile_y(events, (120, 120)).sum() == at_x.sum()
    assert re_mod.profile_x(events, y_range=(150, 150)).sum() == at_y.sum()
    assert re_mod.xy_image(events, tof_band=(t, t)).sum() == at_t.sum()
    assert re_mod.y_tof_image(events, (120, 120), edges).sum() == at_x.sum()
    assert re_mod.profile_tof(events, edges, x_range=(120, 120), y_range=(150, 150)).sum() == (at_x & at_y).sum()
