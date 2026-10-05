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

The ROI pop-out's data (roi-popout-data) is here too, Qt-free: one read of a
run's events (:func:`load_event_pixels`), the two detector images the web report
draws (:func:`xy_image`, :func:`y_tof_image` over :func:`tof_edges`), the three
profiles (:func:`profile_y`, :func:`profile_x`, :func:`profile_tof`), the
background bands as the reducer will use them (:func:`background_bands`), and a
default background in the reducer's shape (:func:`default_bkg_roi`). Counts are
raw event counts, the web report's quantity; pixel ranges are inclusive
``[low, high]``, as the reducer's ``lowres`` and peak masks are.
"""

from dataclasses import dataclass

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


@dataclass(frozen=True)
class RunEvents:
    """One run's bank-1 events, read once, for the pop-out's images and profiles.

    ``x`` and ``y`` are pixel indices, ``x = id // n_y`` and ``y = id % n_y``
    (the packing ``binary_processing.get_y_tof`` uses); ``tof`` is the event time
    offset in microseconds. ``n_x`` and ``n_y`` are the detector shape for the
    run's start time, from the instrument database. ``stride`` is 1 when every
    event is held and k when every k-th is. ``n_off_detector`` counts the ids
    outside ``[0, n_x * n_y)``, which are dropped.
    """

    x: np.ndarray
    y: np.ndarray
    tof: np.ndarray
    n_x: int
    n_y: int
    stride: int = 1
    n_off_detector: int = 0


def load_event_pixels(path, max_events=None):
    """Bank 1's events as :class:`RunEvents`: one read of the file.

    The id packing is the instrument's, taken from
    ``binary_processing.get_y_tof``: ``x = id // n_y``, ``y = id % n_y``. The
    shape comes from the instrument database for this run's start time, so a
    run from a period with a different detector is unpacked with that period's
    geometry. Ids outside the detector are dropped and counted
    (``n_off_detector``), never reshaped into an edge pixel. With
    ``max_events``, every k-th event is held (``stride`` = k): a sample of the
    whole run, never its first N events, which would be a time slice of it.
    """
    if max_events is not None and max_events <= 0:
        raise ValueError(f"max_events must be positive, or None for every event, not {max_events!r}")
    with h5py.File(str(path), "r") as f:
        event_id = np.asarray(f["entry/bank1_events/event_id"][:], dtype=np.int64)
        tof = np.asarray(f["entry/bank1_events/event_time_offset"][:], dtype=float)
        start_time = _text(f["entry/start_time"])

    n_x, n_y = detector_shape(start_time)
    on_detector = (event_id >= 0) & (event_id < n_x * n_y)
    n_off_detector = int(np.count_nonzero(~on_detector))
    event_id, tof = event_id[on_detector], tof[on_detector]

    stride = 1
    if max_events is not None and len(event_id) > max_events:
        stride = int(np.ceil(len(event_id) / max_events))
        event_id, tof = event_id[::stride], tof[::stride]

    return RunEvents(x=event_id // n_y, y=event_id % n_y, tof=tof, n_x=n_x, n_y=n_y,
                     stride=stride, n_off_detector=n_off_detector)


def _inclusive(name, bounds):
    """``(low, high)`` of an inclusive range; reversed is a caller's error, never an empty selection."""
    low, high = bounds
    if high < low:
        raise ValueError(f"{name} {tuple(bounds)!r} is reversed: ranges are inclusive [low, high]")
    return low, high


def _in_band(tof, tof_band):
    """The events inside an inclusive TOF band (all of them when ``tof_band`` is None)."""
    if tof_band is None:
        return np.ones(len(tof), dtype=bool)
    low, high = _inclusive("tof_band", tof_band)
    return (tof >= low) & (tof <= high)


def _edges(edges):
    """``edges`` as floats. Fewer than two, or not strictly increasing, is a caller's error: numpy bins
    nothing for one edge, and counts an event on a repeated edge in a zero-width bin."""
    edges = np.asarray(edges, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or np.any(np.diff(edges) <= 0):
        raise ValueError(f"edges must be at least two strictly increasing TOF values ({edges.size} given)")
    return edges


def xy_image(events, tof_band=None):
    """Counts per detector pixel, ``image[y, x]``, shape ``(n_y, n_x)``.

    The web report's XY array, cell for cell (``web_report.py:579-583``), with no
    band. With ``tof_band`` (inclusive, microseconds), only the events inside it.
    """
    keep = _in_band(events.tof, tof_band)
    flat = np.bincount(events.y[keep] * events.n_x + events.x[keep], minlength=events.n_x * events.n_y)
    return flat.reshape(events.n_y, events.n_x)


def tof_edges(events, bin_width=50.0):
    """TOF bin edges over the events' full span, ``bin_width`` microseconds apart.

    The full span, never the chopper window: a window is an overlay on the
    image, not a crop of it. 50 us is the web report's bin (``web_report.py:602``).
    The last edge is at or past the latest event, which :func:`y_tof_image`
    then counts (its last bin is closed). No events, nothing to bin:
    :class:`CannotEstimateError`.
    """
    if not bin_width > 0:
        raise ValueError(f"bin_width must be positive, not {bin_width!r}")
    if len(events.tof) == 0:
        raise CannotEstimateError("the run has no events, so there is no TOF span to bin")
    low, high = float(events.tof.min()), float(events.tof.max())
    edges = low + bin_width * np.arange(max(1, int(np.ceil((high - low) / bin_width))) + 1)
    if edges[-1] < high:  # floating-point: never leave the latest event outside
        edges = np.append(edges, edges[-1] + bin_width)
    return edges


def y_tof_image(events, x_range, edges):
    """Counts per (Y pixel, TOF bin) for the events inside the inclusive ``x_range``.

    Shape ``(n_y, len(edges) - 1)``. An event inside the edges is counted once,
    including one at exactly the last edge (the last bin is closed). Mantid's
    histogram drops that event, so this array and ``RefRoi``'s can differ in that
    one cell (F4).
    """
    x_low, x_high = _inclusive("x_range", x_range)
    edges = _edges(edges)
    selected = (events.x >= x_low) & (events.x <= x_high)
    image, _, _ = np.histogram2d(events.y[selected], events.tof[selected],
                                 bins=[np.arange(events.n_y + 1), edges])
    return image.astype(np.int64)


def profile_y(events, x_range, tof_band=None):
    """Counts per detector row for the events inside the inclusive ``x_range`` (and ``tof_band``).

    With no band, the Y-TOF image summed over TOF, and, divided by the run's
    proton charge, :func:`counts_vs_y` (F8).
    """
    x_low, x_high = _inclusive("x_range", x_range)
    keep = (events.x >= x_low) & (events.x <= x_high) & _in_band(events.tof, tof_band)
    return np.bincount(events.y[keep], minlength=events.n_y)


def profile_x(events, y_range=None, tof_band=None):
    """Counts per detector column, for the events inside ``y_range`` and ``tof_band`` (each optional)."""
    keep = _in_band(events.tof, tof_band)
    if y_range is not None:
        y_low, y_high = _inclusive("y_range", y_range)
        keep &= (events.y >= y_low) & (events.y <= y_high)
    return np.bincount(events.x[keep], minlength=events.n_x)


def profile_tof(events, edges, x_range=None, y_range=None):
    """Counts per TOF bin over ``edges``, for the events inside ``x_range`` and ``y_range`` (each optional)."""
    keep = np.ones(len(events.tof), dtype=bool)
    if x_range is not None:
        x_low, x_high = _inclusive("x_range", x_range)
        keep &= (events.x >= x_low) & (events.x <= x_high)
    if y_range is not None:
        y_low, y_high = _inclusive("y_range", y_range)
        keep &= (events.y >= y_low) & (events.y <= y_high)
    counts, _ = np.histogram(events.tof[keep], bins=_edges(edges))
    return counts


def counts_vs_y(path, lowres=None, max_events=None, n_tof_bins=200, tof_band=None):
    """Counts per detector row, through the library histogrammer.

    ``lowres`` is the inclusive X pixel range. ``None`` means every X pixel of this
    run's detector, ``(0, n_x - 1)``, from the instrument database. A literal
    such as ``(0, 255)`` agrees with the database today, and would not follow it.

    Calls :func:`lr_reduction.binary_processing.get_y_tof` rather than
    unpacking event ids here. That function already owns the id packing, the
    x-range filter and the proton-charge normalisation; a private copy would be
    a second place for the packing convention to drift, which is the defect
    class this module exists to avoid.

    ``tof_band`` optionally restricts the sum to ``(tof_min, tof_max)`` in
    microseconds — the caller's way of looking only at the band the chopper
    actually delivered. The band is applied to the events, before the one
    ``get_y_tof`` call.
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
    if lowres is None:
        lowres = (0, n_x - 1)

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

    if tof_band is not None:
        # get_y_tof clips out-of-range events into the edge bins rather than
        # dropping them, so a band has to be applied to the events, not to the
        # histogram: trimming columns would keep the clipped strays.
        in_band = (e_offset >= lo) & (e_offset <= hi)
        event_id, e_offset = event_id[in_band], e_offset[in_band]

    _, y_tof, _ = binary_processing.get_y_tof(
        tof_array, event_id, e_offset, list(lowres), pcharge, n_y=n_y, n_x=n_x
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


def background_bands(bkg_roi, y_min, y_max):
    """The two background bands the reducer averages for one angle's ``BkgROI`` entry.

    ``((b0, b1), (b2, b3))``, inclusive: the four bounds ``NR_Reduction._background_roi_sorter``
    (``nr_reduction_calc.py:827-840``) returns and the reducer's masks use (``:865-866``), computed the
    same way, values and types. A fractional bound stays as the sorter keeps it, and the reducer's mask
    then starts at the next row. The bounds are sorted. Exactly two zeros mean "adjacent to the peak": the
    first two sorted bounds become ``y_min`` and ``y_max`` (the peak's own edges), and the four are sorted
    again. Mirrored, not shared: extracting the sorter would put this slug on the reduction path, which
    belongs with the reducer's own fix for pixel 0; T8 pins the two point-wise.

    Every entry on which the reducer would fail is a ``ValueError`` naming the reason. Its sorter returns
    ``None`` for one, three or four zeros, because pixel 0 is its sentinel. It returns a short entry
    unchanged, which the reducer then indexes at ``[3]``. A NaN or infinite bound passes through, and the
    band's centre is then NaN or infinite.
    """
    if bkg_roi is None:
        raise ValueError("no background is set for this angle (its BkgROI entry is None)")
    if isinstance(bkg_roi, (str, bytes)):
        raise ValueError(f"a background is four pixel bounds, not the text {bkg_roi!r}")
    try:
        bounds = np.asarray(bkg_roi)
    except ValueError as exc:  # ragged nesting
        raise ValueError(f"a background is one angle's four pixel bounds, not {bkg_roi!r}") from exc
    if bounds.ndim != 1:
        raise ValueError(
            f"a background is one angle's four bounds, not a {bounds.shape} array (the per-angle "
            f"BkgROI list?): {bkg_roi!r}"
        )
    if len(bounds) != 4:
        raise ValueError(
            f"a background needs four bounds, two per band; {bkg_roi!r} has {len(bounds)}, "
            f"and the reducer reads the fourth"
        )
    numeric = np.issubdtype(bounds.dtype, np.integer) or np.issubdtype(bounds.dtype, np.floating)
    if not numeric or not np.all(np.isfinite(bounds)):
        raise ValueError(f"a background's four bounds must be finite numbers: {bkg_roi!r}")
    ordered = np.sort(bounds)
    zeros = int(np.sum(ordered == 0))
    if zeros == 2:
        ordered[0], ordered[1] = y_min, y_max
        ordered = np.sort(ordered)
    elif zeros != 0:
        raise ValueError(
            f"a background with {zeros} zero bound(s): pixel 0 is the reducer's sentinel for "
            f"'adjacent to the peak', which takes exactly two, so it returns None for {bkg_roi!r}"
        )
    b0, b1, b2, b3 = ordered.tolist()
    return (b0, b1), (b2, b3)


def default_bkg_roi(peak_range, n_y, gap=3, width=5):
    """A default background in the reducer's form: a band on each side of the peak.

    Four ascending ints ``(b0, b1, b2, b3)``: ``[b0, b1]`` ends ``gap`` rows below
    the peak and ``[b2, b3]`` starts ``gap`` rows above it, each ``width`` rows
    wide. That is the entry ``BkgROI`` holds, and ``_background_roi_sorter``
    returns it unchanged. The defaults, 3 and 5, are #197's, which the
    scientists reviewed.

    Refused, never clamped. A peak off the detector is refused, because clamping
    returns a band from the wrong end (``RB_Ymin``/``RB_Ymax`` arrive as unvalidated
    file input). So is a side with no room: a band may not reach row 0, which is
    the reducer's sentinel and turns its sorter's answer into ``None``, nor run
    past the last row.
    """
    peak_low, peak_high = int(peak_range[0]), int(peak_range[1])
    if gap < 0 or width < 1:
        raise ValueError(f"gap must be >= 0 and width >= 1, not gap={gap!r}, width={width!r}")
    if not 0 <= peak_low <= peak_high <= n_y - 1:
        raise ValueError(
            f"peak {peak_range} is not on a {n_y}-pixel detector (rows 0-{n_y - 1}); "
            f"refusing rather than clamping, which would return a band from the "
            f"wrong end"
        )
    b1 = peak_low - gap - 1
    b0 = b1 - width + 1
    b2 = peak_high + gap + 1
    b3 = b2 + width - 1
    if b0 < 1 or b3 > n_y - 1:
        raise ValueError(
            f"no room for a {width}-pixel background with a {gap}-pixel gap on each side of "
            f"peak {peak_range} on a {n_y}-pixel detector (a band starts at row 1: row 0 is "
            f"the reducer's sentinel)"
        )
    return b0, b1, b2, b3
