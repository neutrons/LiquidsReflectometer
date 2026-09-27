"""Qt-free ROI and metadata estimation for a REF_L run.

Pure functions over a NeXus file: no Qt, no widget state, no module-level
caches. They exist so three callers can share one implementation — the settings
editor's ROI guesser, the resolver's layer-(e) ``dataset_probe``, and the #197
comparison — instead of the three growing their own, which is how the chopper
bandwidth ended up with four different values in one tree (see
:func:`chopper_tof_window`).

Nothing here re-derives instrument maths. Where the library already has a
function, this module calls it; where the geometry is time-dependent, it comes
from the instrument database (``nr_tools.read_settings``) and never from a
literal, because the detector was 256x304 at 15.75 m only for part of the
instrument's life.
"""

import h5py
import numpy as np

from lr_reduction import nr_tools


class CannotEstimateError(ValueError):
    """No estimate is available for this run — a genuine refusal, not a caller error.

    Layer (e)'s `dataset_probe` needs to tell "I looked and there is nothing to
    guess" apart from "you called me wrong". Both were plain `ValueError`, so a
    probe wrapper written as `except ValueError: return None` would swallow
    programming errors as "no guess" and fall through to layer (f) wearing a
    badge that reads "default" — the silent-wrong-value class. Introduced now so
    layer (e) inherits it; introducing it later changes what the probe catches.
    """


#: DASlogs paths, matching `binary_processing.get_log_values` exactly. Motors
#: are read at [-1] (where the axis ended up) and the autoreduce counters at [0]
#: (they are set once per run).
_MOTOR_LOGS = {
    "ths": "entry/DASlogs/BL4B:Mot:ths.RBV/value",
    "thi": "entry/DASlogs/BL4B:Mot:thi.RBV/value",
    "tthd": "entry/DASlogs/BL4B:Mot:tthd.RBV/value",
}
_SEQUENCE_LOGS = {
    "seq_num": "entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Num/value",
    "seq_id": "entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Id/value",
}
_CHOPPER_LAMBDA_LOG = "entry/DASlogs/LambdaRequest/value"
_CHOPPER_SPEED_LOG = "entry/DASlogs/SpeedRequest1/value"


def _text(dataset):
    """Decode a NeXus string dataset, which is a length-1 array of bytes."""
    value = dataset[0]
    return value.decode() if isinstance(value, bytes) else str(value)


def read_nexus_metadata(path):
    """Run identity and the three angles, as plain Python values.

    Returns a dict rather than a dataclass because the callers merge it into
    settings dictionaries; a dataclass here would be unpacked at every call
    site.
    """
    meta = {}
    with h5py.File(str(path), "r") as f:
        meta["title"] = _text(f["entry/title"])
        meta["start_time"] = _text(f["entry/start_time"])
        meta["run_number"] = int(_text(f["entry/run_number"]))
        for name, dpath in _MOTOR_LOGS.items():
            # [-1], not [0]: a motor log's first sample is where the axis was
            # BEFORE it moved, so [0] reports the previous run's angle for every
            # run — plausibly, and silently.
            meta[name] = float(f[dpath][-1])
        for name, dpath in _SEQUENCE_LOGS.items():
            meta[name] = int(f[dpath][0])
    return meta


def chopper_lambda_range(path, scaled_width=None):
    """The wavelength band this run actually measured, in **Angstrom**.

    Named for what it returns. It was `chopper_tof_window`, which promised TOF
    and returned wavelength, so composing it with `counts_vs_y(tof_band=...)` —
    documented in microseconds — filtered events to 2.4-5.8 us, produced an
    all-zero profile, and made `estimate_peak_range` blame the run for a unit
    error. Use :func:`lambda_to_tof` to cross the units deliberately.

    Delegates to :func:`lr_reduction.nr_tools.get_lam_range`. **Do not inline
    this maths.** It already exists at four values in this tree — 3.4 in
    ``get_lam_range``'s signature, "3.5" in that function's own docstring,
    ``1.3 * 60 / speed`` (equivalent to 2.6, and without the -0.15 shift) at
    ``event_reduction.py:54-55``, and 3.5 again in the #197 settings builder.
    A fifth copy here would make the eventual reconciliation harder, so the
    default is deliberately ``None`` -> whatever the library's default is, and
    a caller who needs a different band passes it explicitly.

    Raises ``KeyError`` when the run has no chopper log. It does **not** fall
    back to a default band: the band belongs to one measurement, so a guessed
    one would silently describe a different run's configuration. Absent is a
    state of its own, not a missing value to paper over.
    """
    with h5py.File(str(path), "r") as f:
        for dpath in (_CHOPPER_LAMBDA_LOG, _CHOPPER_SPEED_LOG):
            if dpath not in f:
                raise KeyError(
                    f"{path} has no chopper log at {dpath!r} — the wavelength "
                    f"band cannot be derived for this run, and must not be guessed"
                )
        chopper_lam = float(f[_CHOPPER_LAMBDA_LOG][0])
        chopper_speed = float(f[_CHOPPER_SPEED_LOG][0])

    if chopper_speed == 0:
        raise ValueError(f"{path} reports chopper speed 0 — no band can be derived")

    if scaled_width is None:
        return nr_tools.get_lam_range(chopper_lam, chopper_speed)
    return nr_tools.get_lam_range(chopper_lam, chopper_speed, scaled_width=scaled_width)


def lambda_to_tof(lam_range, start_time):
    """Convert an Angstrom band to a TOF band in microseconds for this run.

    de Broglie, with the moderator-to-detector flight path from the
    time-indexed instrument database rather than a literal — the distance has
    three entries in `settings.json` and has genuinely changed.

    tof[us] = (m_n / h) * L[m] * lam[A] * 1e-4
    """
    settings = nr_tools.read_settings(start_time)
    # read_settings reports this in MILLIMETRES (15750.0 for the 15.75 m
    # path) — the reducer works in mm. Writing the conversion without checking
    # gave a band of 9.5e6 us, ~1000x the run's whole TOF span; the composition
    # test is what caught it, which is the case for having written one.
    flight_m = float(settings["source_detector_distance"]) / 1000.0
    # m_n/h in us/(m*A): 252.7701 (CODATA), the standard neutron TOF constant.
    k = 252.7701
    lo, hi = float(lam_range[0]), float(lam_range[1])
    return k * flight_m * lo, k * flight_m * hi


def detector_shape(start_time):
    """(n_x, n_y) for the given run time, from the instrument database.

    Never the literals 256/304: those were right for part of the instrument's
    life only, and `settings.json` is time-indexed precisely because they
    changed.
    """
    settings = nr_tools.read_settings(start_time)
    return int(settings["num_x_pixels"]), int(settings["num_y_pixels"])


def load_event_pixels(path, max_events=None):
    """Event ``(x, y, tof)`` arrays for bank 1, optionally sub-sampled.

    The id packing is the instrument's, taken from
    ``binary_processing.get_y_tof``: ``x = id // n_y``, ``y = id % n_y``. ``n_y``
    comes from the instrument database for this run's start time, so a run from
    a period with a different detector is unpacked with that period's geometry.
    """
    with h5py.File(str(path), "r") as f:
        event_id = np.asarray(f["entry/bank1_events/event_id"][:])
        tof = np.asarray(f["entry/bank1_events/event_time_offset"][:])
        start_time = _text(f["entry/start_time"])

    if max_events is not None and len(event_id) > max_events:
        # Stride rather than head: the first N events are the first N
        # microseconds of the run, which is a time slice, not a sample of it.
        step = int(np.ceil(len(event_id) / max_events))
        event_id = event_id[::step]
        tof = tof[::step]

    _, n_y = detector_shape(start_time)
    return event_id // n_y, event_id % n_y, tof


def counts_vs_y(path, lowres=(0, 255), max_events=None, n_tof_bins=200, tof_band=None):
    """Counts per detector row, through the library histogrammer.

    Calls :func:`lr_reduction.binary_processing.get_y_tof` rather than
    unpacking event ids here. That function already owns the id packing, the
    x-range filter and the proton-charge normalisation; a private copy would be
    a second place for the packing convention to drift, which is the defect
    class this module exists to avoid.

    ``tof_band`` optionally restricts the sum to ``(tof_min, tof_max)`` in
    microseconds — the caller's way of looking only at the band the chopper
    actually delivered.
    """
    from lr_reduction import binary_processing

    with h5py.File(str(path), "r") as f:
        event_id = np.asarray(f["entry/bank1_events/event_id"][:])
        e_offset = np.asarray(f["entry/bank1_events/event_time_offset"][:])
        # `entry/proton_charge` (the run total), matching what
        # `binary_processing.load_and_extract` hands `get_y_tof`. The DASlogs
        # entry beside it is the time SERIES, and passing that makes the
        # normalisation a broadcast against the histogram rather than a divide.
        pcharge = np.asarray(f["entry/proton_charge"][:])
        start_time = _text(f["entry/start_time"])

    if max_events is not None and len(event_id) > max_events:
        step = int(np.ceil(len(event_id) / max_events))
        event_id = event_id[::step]
        e_offset = e_offset[::step]

    total_charge = float(np.sum(pcharge))
    if not np.isfinite(total_charge) or total_charge <= 0:
        raise CannotEstimateError(
            f"proton charge is {total_charge!r} — dividing by it yields a NaN/inf "
            f"profile that walks every downstream guard and brackets pixel 0"
        )

    n_x, n_y = detector_shape(start_time)

    if len(e_offset) == 0:
        return np.zeros(n_y)

    lo = float(np.min(e_offset)) if tof_band is None else float(tof_band[0])
    hi = float(np.max(e_offset)) if tof_band is None else float(tof_band[1])
    if hi <= lo:
        # Not cosmetic: without this `np.linspace(hi, lo, n)` descends, `d_tof`
        # goes negative, `np.digitize` clips into valid indices, and the
        # function returns a plausible profile from inverted bins with no error.
        raise ValueError(f"empty or inverted TOF band {(lo, hi)!r}")
    tof_array = np.linspace(lo, hi, n_tof_bins)

    _, y_tof, _ = binary_processing.get_y_tof(
        tof_array, event_id, e_offset, list(lowres), pcharge, n_y=n_y, n_x=n_x
    )

    if tof_band is not None:
        # get_y_tof clips out-of-range events into the edge bins rather than
        # dropping them, so a band has to be applied to the events, not to the
        # histogram: trimming columns here would keep the clipped strays.
        in_band = (e_offset >= lo) & (e_offset <= hi)
        _, y_tof, _ = binary_processing.get_y_tof(
            tof_array, event_id[in_band], e_offset[in_band], list(lowres),
            pcharge, n_y=n_y, n_x=n_x,
        )

    return y_tof.sum(axis=1)


def estimate_peak_range(counts, min_contrast=1.5, smooth=3, with_contrast=False):
    """Bracket the specular peak: smooth, take the max, walk out to half-max.

    Two states are refused rather than answered, because both have an
    ``argmax`` and neither has a peak:

    * **no counts at all** — ``argmax`` of zeros is 0, so an unguarded
      estimator returns a confident ROI at the detector edge for a run that
      recorded nothing;
    * **no contrast** — flat illumination has a largest sample, and reporting
      it as a peak is reporting noise with a straight face.

    ``min_contrast`` is peak height over the median of everything outside the
    bracket. Returns ``(low, high)``, or ``(low, high, contrast)`` when
    ``with_contrast``.
    """
    counts = np.asarray(counts, dtype=float)
    if not np.all(np.isfinite(counts)):
        # Checked BEFORE the no-counts guard: a zero proton-charge divide makes
        # empty rows NaN and occupied rows inf, and `isfinite & > 0` is then
        # False everywhere — so the no-counts guard would fire and report the
        # wrong cause. Unguarded entirely, `any(counts > 0)` passes on the infs,
        # argmax finds the first NaN at index 0, both walks stop, and the
        # function returns (0, 0).
        raise CannotEstimateError(
            "the row profile contains NaN or inf — normalisation produced a "
            "non-finite profile (a zero proton charge does this), so no peak "
            "can be estimated"
        )
    if counts.size == 0 or not np.any(np.isfinite(counts) & (counts > 0)):
        raise CannotEstimateError("no counts on the detector — no peak can be estimated")

    if smooth and smooth > 1:
        kernel = np.ones(int(smooth)) / float(smooth)
        smoothed = np.convolve(counts, kernel, mode="same")
    else:
        smoothed = counts

    centre = int(np.argmax(smoothed))
    peak = float(smoothed[centre])
    half = peak / 2.0

    low = centre
    while low > 0 and smoothed[low - 1] >= half:
        low -= 1
    high = centre
    while high < len(smoothed) - 1 and smoothed[high + 1] >= half:
        high += 1

    # Baseline is the median of the WHOLE row profile, not of what falls
    # outside the bracket. On a featureless detector the half-max walk runs
    # almost edge to edge — every bin is above half of a flat maximum — so
    # "outside" is empty, and an empty outside was reading as infinite
    # contrast: the one input the contrast score exists to reject scored best.
    baseline = float(np.median(smoothed))

    # Refuse, do NOT score as infinite. `inf < min_contrast` is False for every
    # threshold, so the previous `else float("inf")` made this guard a no-op on
    # exactly the runs it exists for: a sparse, low-flux profile has median 0.
    # Reachable on 5 of 63 real REF_L files (70-77% empty rows) — a profile of
    # zeros with one 3-count row scored (199, 201, inf).
    if not np.isfinite(baseline) or baseline <= 0:
        raise CannotEstimateError(
            f"the row profile has a non-positive baseline ({baseline!r}) — it is "
            f"too sparse to separate a peak from nothing, so no contrast can be "
            f"computed and no ROI is offered"
        )
    contrast = peak / baseline

    if contrast < min_contrast:
        raise CannotEstimateError(
            f"contrast {contrast:.2f} is below {min_contrast} — the detector is "
            f"featureless here, so the largest bin is noise, not a peak"
        )

    if with_contrast:
        return low, high, contrast
    return low, high


def default_bkg_roi(peak_range, n_y, gap=5, width=10):
    """A background band beside the peak, on whichever side has room.

    Prefers the low side, matching the tab's habit, and falls back to the high
    side. The **fit checks are the protection**: a band is only returned when it
    lies wholly on the detector, because a negative row or one past the last
    pixel indexes silently in numpy and yields a background taken from the wrong
    end of the detector.

    An earlier draft also clamped the result with ``max(0, ...)`` / ``min(n_y-1,
    ...)``. Those were **dead code** — each sat inside a branch whose own
    condition already forbids the out-of-range case, so neither could ever
    change a value, and a mutation deleting them passed every test. Removed
    rather than kept as reassurance: an unreachable guard reads as protection
    and provides none.
    """
    peak_low, peak_high = int(peak_range[0]), int(peak_range[1])

    # Each branch below checks only ONE end, so the other end went unchecked and
    # the function returned bands off the detector: (400, 410) on a 304-row
    # detector gave (385, 394), and (-50, -40) gave (-34, -25). I had removed
    # clamps here as "dead code"; they were dead only under an unstated
    # precondition — that `peak_range` is on the detector — which nothing
    # checked. Validating and refusing rather than clamping, because clamping
    # (400, 410) returns a band from the wrong END of the detector, silently,
    # which is the failure this docstring names. RB_Ymin/RB_Ymax reach the
    # resolver from layer (c) as unvalidated file input.
    if not 0 <= peak_low <= peak_high <= n_y - 1:
        raise ValueError(
            f"peak {peak_range} is not on a {n_y}-pixel detector (rows 0-{n_y - 1}); "
            f"refusing rather than clamping, which would return a band from the "
            f"wrong end"
        )

    high_edge = peak_low - gap - 1
    if high_edge - width + 1 >= 0:
        return high_edge - width + 1, high_edge

    low_edge = peak_high + gap + 1
    if low_edge + width - 1 <= n_y - 1:
        return low_edge, low_edge + width - 1

    raise ValueError(
        f"no room for a {width}-pixel background with a {gap}-pixel gap beside "
        f"peak {peak_range} on a {n_y}-pixel detector"
    )
