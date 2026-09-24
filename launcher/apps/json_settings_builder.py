#!/usr/bin/env python
"""
Small Qt interface to build the JSON settings file read by the new reduction
workflow (``lr_reduction.new_reduction_from_file``) and by autoreduction
(``shared/autoreduce/reduce_settings.json``).

A settings file holds two kinds of entries:

* per-run arrays, with one entry per angle setting of a sequence, ordered by
  sequence number: ``method_per_run``, ``DBname``, ``RB_Ymin``, ``RB_Ymax``,
  ``BkgROI``, ``useBS``, ``ThetaShift`` and ``ScaleFactor``.
* scalar options that apply to the whole reduction.

Each row of the table corresponds to one entry of those arrays and can be tied
to a NeXus file. Loading the NeXus file fills in the run number, the sequence
number and the title, and optionally estimates the specular peak and the
background ROIs from the events.
"""

import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np
from qtpy import QtCore
from qtpy.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

try:
    from matplotlib.figure import Figure
    from matplotlib.widgets import SpanSelector

    try:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
    except ImportError:  # matplotlib < 3.5
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
except ImportError:  # the interface works without the ROI plots
    Figure = None

from launcher.apps.file_batch import parse_run_list
from lr_reduction import nr_tools as tools

# Detector geometry: event_id = x_pixel * N_Y + y_pixel
N_Y = 304
N_X = 256

# Maximum number of events read from a run. The profiles only serve to choose
# pixel ranges, so a sub-sample is plenty, and it keeps the ROI dialog quick on
# long runs.
MAX_EVENTS = 2000000

# Chopper band used to select the useful TOF range (see the ROI selector)
#CHOPPER_BANDWIDTH = 3.5
MODERATOR_DETECTOR_DISTANCE = 15.75  # meters

METHODS = ["meantheta", "constantq", "constanttof"]

# Per-run arrays written to the settings file, in the order they are written
PER_RUN_KEYS = [
    "method_per_run",
    "DBname",
    "RB_Ymin",
    "RB_Ymax",
    "BkgROI",
    "useBS",
    "ThetaShift",
    "ScaleFactor",
]

# Scalar options: (key, label, kind, default, extra)
#   int_pair   extra = (minimum, maximum)
#   int/float  extra = (minimum, maximum[, decimals, step])
#   choice     extra = [(label, value), ...]
#   float_opt  same as float plus an initial value, with a checkbox to leave
#              the key out of the file and take the value from the run logs
GLOBAL_FIELDS = [
    ("data_x_range", "Data x range [pixels]", "int_pair", [50, 220], (0, N_X - 1)),
    ("Normalize", "Normalize to critical edge", "bool", False, None),
    ("AutoScale", "Auto-scale between angles", "bool", False, None),
    (
        "useCalcTheta",
        "Theta from fitted peak",
        "choice",
        False,
        [("off", False), ("on", True)]
        # can go back to full version in future once properly implemented
        #[("off", False), ("detector_angle", "detector_angle"), ("sample_angle", "sample_angle")],
    ),
    ("useGravity", "Gravity correction", "bool", True, None),
    ("save8col", "Save 8-column output", "bool", False, None),
    ("qmin", "q min [1/A]", "float", 0.001, (0.0, 10.0, 5, 0.001)),
    ("qmax", "q max [1/A]", "float", 0.5, (0.0, 10.0, 5, 0.01)),
    ("dqbin", "dq bin [1/A]", "float", 0.005, (0.0, 1.0, 5, 0.001)),
    ("Qnorm", "q for normalization [1/A]", "float", 0.015, (0.0, 1.0, 5, 0.001)),
    ("Qline_threshold", "q-line threshold", "float", 1.0, (0.0, 1.0, 3, 0.05)),
    ("tof_bin", "TOF bin [us]", "int", 50, (1, 10000)),
    ("IncidentTheta", "Incident theta [deg]", "float_opt", None, (-10.0, 10.0, 3, 0.01, 4.0)),
    ("dead_time", "Dead time [us]", "float", 4.2, (0.0, 100.0, 3, 0.1)),
    ("dead_time_tof_step", "Dead time TOF step [us]", "int", 50, (1, 10000)),
    ("DetResFn", "Detector resolution", "choice", "gaussian",
     [("rectangular", "rectangular"), ("gaussian", "gaussian"), ("none", "none")]),
    ("DetSigma", "Detector sigma [pixel]", "float", 1.0, (0.0, 20.0, 3, 0.1)),
    ("peak_pad", "Peak fit padding [pixels]", "int", 1, (0, 50)),
    ("peak_type", "Peak shape", "choice", "supergauss", [("supergauss", "supergauss"), ("gauss", "gauss")]),
]

GLOBAL_KEYS = [item[0] for item in GLOBAL_FIELDS]


@dataclass
class RunRow:
    """Settings for a single angle setting (one entry of each per-run array)."""

    nexus_path: str = ""
    run: str = ""
    seq: int = 0
    seq_id: str = ""
    title: str = ""
    method: str = "meantheta"
    db_name: str = ""
    y_min: int = 0
    y_max: int = 0
    bkg: list = field(default_factory=lambda: [0, 0, 0, 0])
    use_bs: int = 1
    theta_shift: float = 0.0
    scale_factor: float = 1.0
    tof_range: object = None  # TOF range the ROIs were chosen on, None for the chopper band
    profile: object = None  # cached counts vs y pixel, for plotting


def _first_value(h5_file, path, default=None):
    """Return the first element of a NeXus dataset, decoded if it holds bytes."""
    try:
        value = h5_file[path][0]
    except (KeyError, IndexError, OSError, TypeError, ValueError):
        return default
    if isinstance(value, (bytes, bytearray)):
        return value.decode("utf-8", errors="ignore")
    return value


def read_nexus_metadata(file_path):
    """
    Read the information needed to place a run in the sequence.

    Returns a dictionary with the run number, the autoreduction sequence number
    and id, the title and the sample/detector angles.
    """
    with h5py.File(file_path, "r") as h5_file:
        meta = {
            "run": _first_value(h5_file, "entry/run_number", ""),
            "title": _first_value(h5_file, "entry/title", ""),
            "seq_num": _first_value(h5_file, "entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Num/value", 0),
            "seq_id": _first_value(h5_file, "entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Id/value", 0),
            "ths": _first_value(h5_file, "entry/DASlogs/ths/value"),
            "thi": _first_value(h5_file, "entry/DASlogs/thi/value"),
            "tthd": _first_value(h5_file, "entry/DASlogs/tthd/value"),
        }

    if not meta["run"]:
        match = re.search(r"REF_L_(\d+)", os.path.basename(str(file_path)))
        meta["run"] = match.group(1) if match else ""
    meta["run"] = str(meta["run"]).strip()
    try:
        meta["seq_num"] = int(meta["seq_num"])
    except (TypeError, ValueError):
        meta["seq_num"] = 0
    return meta

# change to use the same logic as the reduction, i.e. to use nr_tools.get_lam_range
def _chopper_tof_range(h5_file):
    """TOF range [us] covered by the chopper wavelength band, or None."""
    chopper_lambda = _first_value(h5_file, "entry/DASlogs/BL4B:Det:TH:BL:Lambda/value")
    chopper_speed = _first_value(h5_file, "entry/DASlogs/BL4B:Det:TH:BL:Frequency/value")
    if chopper_lambda is None or chopper_speed is None or float(chopper_speed) == 0:
        return None

    wl_min, wl_max = tools.get_lam_range(chopper_lambda, chopper_speed, scaled_width=3.4)
    #half_band = CHOPPER_BANDWIDTH / 2.0 * 60.0 / float(chopper_speed)
    #wl_min = float(chopper_lambda) - half_band
    #wl_max = float(chopper_lambda) + half_band
    tof_min = 252.78 * wl_min * MODERATOR_DETECTOR_DISTANCE
    tof_max = 252.78 * wl_max * MODERATOR_DETECTOR_DISTANCE
    if tof_max <= tof_min:
        return None
    return tof_min, tof_max


def load_events(file_path):
    """
    The events of a run, sub-sampled for very large files.

    Returns
    -------
    tuple
        The x (low resolution) pixel, the y (reflectivity) pixel and the time
        of flight of each event, and the TOF range of the chopper wavelength
        band, which falls back to the range covered by the events.
    """
    with h5py.File(file_path, "r") as h5_file:
        events = h5_file["entry/bank1_events/event_id"]
        stride = max(1, int(np.ceil(events.shape[0] / MAX_EVENTS)))
        event_id = events[::stride]
        time_offset = h5_file["entry/bank1_events/event_time_offset"][::stride]
        chopper_range = _chopper_tof_range(h5_file)

    keep = (event_id >= 0) & (event_id < N_X * N_Y)
    event_id = event_id[keep]
    time_offset = time_offset[keep].astype(float)

    if time_offset.size == 0:
        raise ValueError("The file holds no event in the detector")
    full_range = (float(time_offset.min()), float(time_offset.max()))
    if chopper_range is None:
        chopper_range = full_range
    else:  # the band can reach outside the frame
        chopper_range = (max(chopper_range[0], full_range[0]), min(chopper_range[1], full_range[1]))
        if chopper_range[1] <= chopper_range[0]:
            chopper_range = full_range

    return event_id // N_Y, event_id % N_Y, time_offset, chopper_range


def counts_in_range(pixels, selection, n_pixels):
    """Counts per pixel of ``pixels``, keeping only the events of ``selection``."""
    return np.bincount(pixels[selection], minlength=n_pixels)[:n_pixels].astype(float)


def counts_vs_y(file_path, x_range=(50, 200), tof_range=None):
    """
    Counts per vertical pixel, summed over the useful x pixels and over a TOF
    range, which defaults to the chopper wavelength band.
    """
    x_pixel, y_pixel, tof, chopper_range = load_events(file_path)
    if tof_range is None:
        tof_range = chopper_range
    inside = (x_pixel >= min(x_range)) & (x_pixel <= max(x_range))
    inside &= (tof >= min(tof_range)) & (tof <= max(tof_range))
    return counts_in_range(y_pixel, inside, N_Y)


def read_db_header(file_path, max_lines=20):
    """
    The ``key = value`` header of a processed direct beam file.

    Only the header is read, so that the file list stays cheap to build; the
    reduction reads the whole file with ``nr_tools.load_db_file``.
    """
    meta = {}
    try:
        with open(file_path, "r", errors="ignore") as handle:
            for _ in range(max_lines):
                line = handle.readline()
                if not line or not line.startswith("#"):
                    break
                if "=" in line:
                    key, value = line.lstrip("#").split("=", 1)
                    meta[key.strip()] = value.strip()
    except OSError:
        return {}
    return meta


def estimate_peak_range(profile, max_half_width=25):
    """
    Estimate the specular peak range from a counts-vs-pixel profile.

    The peak is taken at the maximum of the smoothed profile and its edges at
    half of the background-subtracted maximum.

    Returns
    -------
    tuple
        ``(y_min, y_max, contrast)`` where ``contrast`` is the peak height over
        the background level. A low contrast means the estimate is unreliable
        and should be checked against the profile plot.
    """
    smoothed = np.convolve(profile, np.ones(3) / 3.0, mode="same")
    peak = int(np.argmax(smoothed))

    # Background taken next to the peak: the detector edges are much quieter
    # than the region around the specular peak and would bias the estimate low
    distance = np.abs(np.arange(N_Y) - peak)
    wings = (distance > max_half_width) & (distance <= max_half_width + 35)
    background = float(np.median(smoothed[wings])) if np.any(wings) else 0.0
    half_maximum = background + (smoothed[peak] - background) / 2.0

    y_min = peak
    while y_min > max(0, peak - max_half_width) and smoothed[y_min - 1] > half_maximum:
        y_min -= 1
    y_max = peak
    while y_max < min(N_Y - 1, peak + max_half_width) and smoothed[y_max + 1] > half_maximum:
        y_max += 1

    contrast = smoothed[peak] / background if background > 0 else float("inf")
    return y_min, y_max, contrast


def default_bkg_roi(y_min, y_max, gap=3, width=5):
    """Background ROIs on either side of the peak, as used in the examples."""
    return [
        int(max(0, y_min - gap - width)),
        int(max(0, y_min - gap)),
        int(min(N_Y - 1, y_max + gap)),
        int(min(N_Y - 1, y_max + gap + width)),
    ]


def _move_span(patch, low, high):
    """
    Move a shaded region, whose vertical extent is in axes coordinates.

    ``axvspan`` returns a rectangle from matplotlib 3.9 on and a polygon
    before, and the two do not take the same coordinates.
    """
    if hasattr(patch, "set_bounds"):  # a rectangle
        patch.set_bounds(low, 0, high - low, 1)
    else:  # a polygon
        patch.set_xy([[low, 0], [low, 1], [high, 1], [high, 0], [low, 0]])


class ROISelectionDialog(QDialog):
    """
    Pick the peak, the background, the x pixel range and the TOF range by
    dragging on the counts profiles of a run.

    The peak and the background belong to the angle setting being edited. The
    x pixel range is the ``data_x_range`` option, which is shared by every
    angle setting. The TOF range only selects the events the profiles are
    made of: at the highest angles nearly every count outside the reflected
    signal is background, which flattens the profile, and narrowing the TOF
    range around the signal brings the peak out.
    """

    PEAK, BACKGROUND_LEFT, BACKGROUND_RIGHT = "peak", "left", "right"
    TOF_BIN = 100.0  # microseconds

    def __init__(self, row, x_range, events=None, parent=None):
        QDialog.__init__(self, parent)
        self.setWindowTitle(f"Run {row.run}: peak, background and ranges" if row.run else "Select the ranges")
        self.resize(950, 850)

        self.x_pixel, self.y_pixel, self.tof, self.chopper_range = events or load_events(row.nexus_path)
        self.tof_edges = np.arange(self.tof.min(), self.tof.max() + self.TOF_BIN, self.TOF_BIN)
        self.y_profile = np.zeros(N_Y)
        self._updating = False
        self._tof_label = ""

        layout = QVBoxLayout()
        self.setLayout(layout)

        # A layout engine would run on every redraw, which is too slow while dragging
        self.figure = Figure(figsize=(8, 7))
        self.y_axis, self.tof_axis, self.x_axis = self.figure.subplots(
            3, 1, gridspec_kw={"height_ratios": [2, 1, 1]}
        )
        self.figure.subplots_adjust(left=0.1, right=0.98, top=0.94, bottom=0.08, hspace=0.45)
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(NavigationToolbar(self.canvas, self))
        layout.addWidget(self.canvas, stretch=1)

        layout.addWidget(self._build_controls(row, x_range))

        title = row.title or ""
        self.y_axis.set_title(f"{title} (sequence {row.seq})" if row.seq else title)
        self.y_axis.set_xlabel("y pixel (reflectivity direction)")
        self.x_axis.set_xlabel("x pixel (low resolution direction)")
        for axis in (self.y_axis, self.tof_axis, self.x_axis):
            axis.set_ylabel("counts")
        (self.y_line,) = self.y_axis.plot([], [], drawstyle="steps-mid", color="tab:blue")
        (self.tof_line,) = self.tof_axis.plot([], [], drawstyle="steps-mid", color="tab:blue")
        (self.x_line,) = self.x_axis.plot([], [], drawstyle="steps-mid", color="tab:blue")

        # The shaded regions are made once and moved, so that the legends can
        # be made once as well: rebuilding them on every change is slow
        self.peak_span = self.y_axis.axvspan(0, 1, color="tab:green", alpha=0.25, label="peak")
        self.bkg_spans = [
            self.y_axis.axvspan(0, 1, color="tab:red", alpha=0.2, label="background"),
            self.y_axis.axvspan(0, 1, color="tab:red", alpha=0.2),
        ]
        self.tof_span = self.tof_axis.axvspan(0, 1, color="tab:green", alpha=0.25, label="TOF range")
        self.x_span = self.x_axis.axvspan(0, 1, color="tab:green", alpha=0.25, label="x range")
        for axis in (self.y_axis, self.tof_axis, self.x_axis):
            axis.legend(loc="upper right", fontsize="small")

        self._set_log_scale(True)
        self._update_profiles()
        self._reset_limits()

        # Dragging on a plot sets the range it shows; on the upper plot, the
        # range selected by the radio buttons
        # useblit keeps the drag from redrawing the whole figure at every mouse
        # move, which freezes the dialog on runs with many events
        self.selectors = [
            SpanSelector(
                axis, callback, "horizontal", useblit=True,
                props={"facecolor": "tab:blue", "alpha": 0.2}, drag_from_anywhere=True,
            )
            for axis, callback in (
                (self.y_axis, self._y_range_selected),
                (self.tof_axis, self._tof_range_selected),
                (self.x_axis, self._x_range_selected),
            )
        ]

    def _build_controls(self, row, x_range):
        """Radio buttons to choose what a drag sets, and the values themselves."""
        box = QWidget()
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        box.setLayout(grid)

        grid.addWidget(QLabel("Drag on the upper plot to set:"), 0, 0)
        modes = QHBoxLayout()
        self.mode_buttons = []
        for label, mode in [("the peak", self.PEAK),
                            ("the left background", self.BACKGROUND_LEFT),
                            ("the right background", self.BACKGROUND_RIGHT)]:
            button = QRadioButton(label)
            self.mode_buttons.append((button, mode))
            modes.addWidget(button)
        self.mode_buttons[0][0].setChecked(True)
        modes.addStretch(1)
        grid.addLayout(modes, 0, 1, 1, 4)

        grid.addWidget(QLabel("Peak (RB_Ymin, RB_Ymax):"), 1, 0)
        self.peak_spins = [self._make_spin(row.y_min, N_Y - 1), self._make_spin(row.y_max, N_Y - 1)]
        for column, spin in enumerate(self.peak_spins):
            grid.addWidget(spin, 1, 1 + column)

        grid.addWidget(QLabel("Background (BkgROI):"), 2, 0)
        self.bkg_spins = [self._make_spin(value, N_Y - 1) for value in row.bkg]
        for column, spin in enumerate(self.bkg_spins):
            grid.addWidget(spin, 2, 1 + column)

        grid.addWidget(QLabel("x range, all settings:"), 3, 0)
        self.x_spins = [self._make_spin(x_range[0], N_X - 1), self._make_spin(x_range[1], N_X - 1)]
        for column, spin in enumerate(self.x_spins):
            grid.addWidget(spin, 3, 1 + column)

        tof_range = row.tof_range or self.chopper_range
        grid.addWidget(QLabel("TOF range [us], profiles only:"), 4, 0)
        self.tof_spins = [self._make_spin(value, int(self.tof.max()) + 1) for value in tof_range]
        for column, spin in enumerate(self.tof_spins):
            spin.setSingleStep(100)
            spin.setToolTip("Only the events of this TOF range are counted in the profiles above")
            grid.addWidget(spin, 4, 1 + column)
        band_btn = QPushButton("Whole band")
        band_btn.setToolTip("Go back to the TOF range of the chopper wavelength band")
        band_btn.clicked.connect(self._reset_tof_range)
        grid.addWidget(band_btn, 4, 3)

        estimate_btn = QPushButton("Estimate the peak")
        estimate_btn.setToolTip("Set the peak and background from the profile shown, over the TOF range above")
        estimate_btn.clicked.connect(self._estimate)
        grid.addWidget(estimate_btn, 5, 0)

        self.log_check = QCheckBox("Log scale")
        self.log_check.setChecked(True)
        self.log_check.toggled.connect(self._set_log_scale)
        grid.addWidget(self.log_check, 5, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        grid.addWidget(buttons, 5, 3, 1, 2)

        return box

    def _make_spin(self, value, maximum):
        spin = QSpinBox()
        spin.setRange(0, maximum)
        spin.setValue(int(value))
        spin.valueChanged.connect(self._values_changed)
        return spin

    # --------------------------------------------------------- selections

    def _mode(self):
        for button, mode in self.mode_buttons:
            if button.isChecked():
                return mode
        return self.PEAK

    def _y_range_selected(self, low, high):
        mode = self._mode()
        if mode == self.PEAK:
            spins = self.peak_spins
        elif mode == self.BACKGROUND_LEFT:
            spins = self.bkg_spins[:2]
        else:
            spins = self.bkg_spins[2:]
        self._apply_range(spins, low, high)

    def _x_range_selected(self, low, high):
        self._apply_range(self.x_spins, low, high)

    def _tof_range_selected(self, low, high):
        self._apply_range(self.tof_spins, low, high)

    def _apply_range(self, spins, low, high):
        low, high = sorted((int(round(low)), int(round(high))))
        low = max(spins[0].minimum(), min(spins[0].maximum(), low))
        high = max(low + 1, min(spins[1].maximum(), high))
        self._updating = True
        try:
            spins[0].setValue(low)
            spins[1].setValue(high)
        finally:
            self._updating = False
        self._values_changed()

    def _values_changed(self):
        if not self._updating:
            self._update_profiles()

    def _reset_tof_range(self):
        self._apply_range(self.tof_spins, *self.chopper_range)

    def _estimate(self):
        y_min, y_max, _contrast = estimate_peak_range(self.y_profile)
        self._updating = True
        try:
            self.peak_spins[0].setValue(y_min)
            self.peak_spins[1].setValue(y_max)
            for spin, value in zip(self.bkg_spins, default_bkg_roi(y_min, y_max)):
                spin.setValue(value)
        finally:
            self._updating = False
        self._update_profiles()
        self._reset_limits()

    # ------------------------------------------------------------ plotting

    def _update_profiles(self):
        """Recompute the three profiles for the current ranges and redraw them."""
        y_min, y_max, bkg, x_range, tof_range = self.values()

        inside_x = (self.x_pixel >= x_range[0]) & (self.x_pixel <= x_range[1])
        inside_tof = (self.tof >= tof_range[0]) & (self.tof <= tof_range[1])
        if y_max > y_min:
            inside_peak = (self.y_pixel >= y_min) & (self.y_pixel <= y_max)
            label = "time of flight [us], counts inside the peak"
        else:  # no peak chosen yet, so show every pixel rather than nothing
            inside_peak = np.ones(self.y_pixel.shape, dtype=bool)
            label = "time of flight [us], counts over the whole detector"
        if label != self._tof_label:
            self.tof_axis.set_xlabel(label)
            self._tof_label = label

        self.y_profile = counts_in_range(self.y_pixel, inside_x & inside_tof, N_Y)
        x_profile = counts_in_range(self.x_pixel, inside_peak & inside_tof, N_X)
        # The TOF profile keeps the whole range, so that the signal can be
        # found again after a narrow range has been selected
        tof_counts, _edges = np.histogram(self.tof[inside_x & inside_peak], bins=self.tof_edges)

        self.y_line.set_data(np.arange(N_Y), self.y_profile)
        self.x_line.set_data(np.arange(N_X), x_profile)
        self.tof_line.set_data((self.tof_edges[:-1] + self.tof_edges[1:]) / 2.0, tof_counts)

        _move_span(self.peak_span, y_min, y_max)
        _move_span(self.bkg_spans[0], bkg[0], bkg[1])
        _move_span(self.bkg_spans[1], bkg[2], bkg[3])
        _move_span(self.tof_span, tof_range[0], tof_range[1])
        _move_span(self.x_span, x_range[0], x_range[1])

        for axis in (self.y_axis, self.tof_axis, self.x_axis):
            axis.relim()
            axis.autoscale_view(scalex=False)
        self.canvas.draw_idle()

    def _set_log_scale(self, log_scale):
        for axis in (self.y_axis, self.tof_axis, self.x_axis):
            axis.set_yscale("log" if log_scale else "linear")
        self.canvas.draw_idle()

    def _reset_limits(self):
        """Show the region around the ROIs, or around the peak if there is none yet."""
        y_min, y_max, bkg, _x_range, _tof_range = self.values()
        edges = [value for value in list(bkg) + [y_min, y_max] if value]
        if not edges:
            edges = [int(np.argmax(self.y_profile))]
        self.y_axis.set_xlim(max(0, min(edges) - 30), min(N_Y - 1, max(edges) + 30))
        self.tof_axis.set_xlim(self.tof_edges[0], self.tof_edges[-1])
        self.x_axis.set_xlim(0, N_X - 1)
        self.canvas.draw_idle()

    def values(self):
        """The peak range, the background ROI, the x pixel range and the TOF range."""
        return (
            self.peak_spins[0].value(),
            self.peak_spins[1].value(),
            [spin.value() for spin in self.bkg_spins],
            [spin.value() for spin in self.x_spins],
            [spin.value() for spin in self.tof_spins],
        )


class JSONSettingsBuilderTab(QWidget):
    """Build and edit the JSON settings file used by the new reduction."""

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.setWindowTitle("Reduction settings builder")

        self.rows = []
        self.db_files = []  # direct beam files offered for DBname
        self.db_headers = {}
        self.db_directory = ""
        self.experiment = ""  # the experiment the directories were set from
        self.extra_keys = {}  # settings keys we read but do not edit here
        self._updating = False
        self.original_settings = None  # store original state for reset

        self.settings = QtCore.QSettings()

        layout = QVBoxLayout()
        self.setLayout(layout)
        layout.addWidget(self._build_paths_box())

        splitter = QSplitter(QtCore.Qt.Vertical)
        splitter.addWidget(self._build_runs_box())
        splitter.addWidget(self._build_options_box())
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        layout.addWidget(splitter, stretch=1)

        self.status_label = QLabel("")
        layout.addWidget(self.status_label)

        self._read_user_settings()
        self.scan_db_files()
        self._refresh_table()

    # ------------------------------------------------------------------ UI

    def _build_paths_box(self):
        box = QGroupBox("Files")
        grid = QGridLayout()
        grid.setColumnStretch(1, 1)
        box.setLayout(grid)

        grid.addWidget(QLabel("Experiment (IPTS-...):"), 0, 0)
        self.experiment_edit = QLineEdit()
        self.experiment_edit.setPlaceholderText("IPTS-12345")
        self.experiment_edit.editingFinished.connect(self._experiment_changed)
        grid.addWidget(self.experiment_edit, 0, 1)

        grid.addWidget(QLabel("NeXus directory:"), 1, 0)
        self.nexus_dir_edit = QLineEdit()
        grid.addWidget(self.nexus_dir_edit, 1, 1)
        nexus_btn = QPushButton("Browse")
        nexus_btn.clicked.connect(lambda: self._browse_directory(self.nexus_dir_edit))
        grid.addWidget(nexus_btn, 1, 2)

        grid.addWidget(QLabel("Direct beam directory:"), 2, 0)
        self.db_dir_edit = QLineEdit()
        self.db_dir_edit.setToolTip("Directory holding the processed direct beam files, offered for DBname")
        self.db_dir_edit.editingFinished.connect(self.scan_db_files)
        grid.addWidget(self.db_dir_edit, 2, 1)
        db_btn = QPushButton("Browse")
        db_btn.clicked.connect(lambda: self._browse_directory(self.db_dir_edit, self.scan_db_files))
        grid.addWidget(db_btn, 2, 2)

        grid.addWidget(QLabel("Settings file:"), 3, 0)
        self.settings_file_edit = QLineEdit()
        self.settings_file_edit.setPlaceholderText("<IPTS>/shared/autoreduce/reduce_settings.json")
        grid.addWidget(self.settings_file_edit, 3, 1)

        buttons = QHBoxLayout()
        load_btn = QPushButton("Load...")
        load_btn.clicked.connect(self.load_settings)
        buttons.addWidget(load_btn)
        save_btn = QPushButton("Save")
        save_btn.clicked.connect(self.save_settings)
        buttons.addWidget(save_btn)
        save_as_btn = QPushButton("Save as...")
        save_as_btn.clicked.connect(self.save_settings_as)
        buttons.addWidget(save_as_btn)
        reset_btn = QPushButton("Reset")
        reset_btn.setToolTip("Clear all edited inputs back to the original state")
        reset_btn.clicked.connect(self._reset_settings)
        buttons.addWidget(reset_btn)
        grid.addLayout(buttons, 3, 2)

        return box

    def _build_runs_box(self):
        box = QGroupBox("Angle settings (one row per data file, in sequence order)")
        vbox = QVBoxLayout()
        box.setLayout(vbox)

        self.headers = [
            "NeXus file", "Run", "Seq", "Title", "method", "DBname", "RB_Ymin", "RB_Ymax",
            "Bkg 1", "Bkg 2", "Bkg 3", "Bkg 4", "useBS", "ThetaShift", "ScaleFactor",
        ]
        self.table = QTableWidget(0, len(self.headers))
        self.table.setHorizontalHeaderLabels(self.headers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.verticalHeader().setDefaultSectionSize(26)
        self.table.horizontalHeader().setSectionResizeMode(self.COL_TITLE, QHeaderView.Stretch)
        self.table.itemChanged.connect(self._item_changed)
        vbox.addWidget(self.table)

        runs = QHBoxLayout()
        runs.addWidget(QLabel("Run numbers:"))
        self.runs_edit = QLineEdit()
        self.runs_edit.setPlaceholderText("233005-233008, 233012")
        self.runs_edit.setToolTip("Run numbers of the sequence, as a list and/or ranges")
        self.runs_edit.returnPressed.connect(self.add_rows_from_runs)
        runs.addWidget(self.runs_edit)
        add_runs_btn = QPushButton("Add runs")
        add_runs_btn.setToolTip("Add one row per run, reading each NeXus file from the NeXus directory")
        add_runs_btn.clicked.connect(self.add_rows_from_runs)
        runs.addWidget(add_runs_btn)
        nexus_btn = QPushButton("Add rows from NeXus...")
        nexus_btn.setToolTip("Add one row per selected NeXus file, for data outside the NeXus directory")
        nexus_btn.clicked.connect(self.add_rows_from_nexus)
        runs.addWidget(nexus_btn)
        vbox.addLayout(runs)

        buttons = QHBoxLayout()
        for label, slot, tip in [
            ("Add row", self.add_row, "Add an empty row, to be filled in by hand"),
            ("Duplicate", self.duplicate_row, None),
            ("Remove", self.remove_rows, None),
            ("Up", lambda: self.move_row(-1), None),
            ("Down", lambda: self.move_row(1), None),
            ("Sort by seq", self.sort_by_sequence, "Sort rows by autoreduction sequence number"),
            ("Assign DB files", self.assign_db_files, "Give the rows the available direct beam files, in order"),
            ("Estimate ROI", self.estimate_selected_rois, "Estimate the peak and background ROIs from the NeXus file"),
            ("Select ROI...", self.edit_roi, "Pick the peak, the background and the x range on the counts profiles"),
        ]:
            button = QPushButton(label)
            button.clicked.connect(slot)
            if tip:
                button.setToolTip(tip)
            buttons.addWidget(button)
        buttons.addStretch(1)
        vbox.addLayout(buttons)

        options = QHBoxLayout()
        self.auto_roi_check = QCheckBox("Estimate ROIs when a NeXus file is loaded")
        self.auto_roi_check.setToolTip("Only for rows without a peak range; use 'Estimate ROI' to work them out again")
        self.auto_roi_check.setChecked(True)
        options.addWidget(self.auto_roi_check)
        self.save_runs_check = QCheckBox("Write run numbers (RBnum) to the file")
        self.save_runs_check.setToolTip("RBnum is set from the runs being reduced, so it is only informative here")
        options.addWidget(self.save_runs_check)
        options.addStretch(1)
        vbox.addLayout(options)

        return box

    def _build_options_box(self):
        box = QGroupBox("Reduction options")
        grid = QGridLayout()
        box.setLayout(grid)

        self.global_widgets = {}
        per_column = int(np.ceil(len(GLOBAL_FIELDS) / 3.0))
        for index, (key, label, kind, default, extra) in enumerate(GLOBAL_FIELDS):
            row = index % per_column
            column = 2 * (index // per_column)
            grid.addWidget(QLabel(label + ":"), row, column)
            widget = self._make_global_widget(key, kind, default, extra)
            grid.addWidget(widget, row, column + 1)
        for column in (1, 3, 5):
            grid.setColumnStretch(column, 1)

        return box

    def _make_global_widget(self, key, kind, default, extra):
        """Create the editor for one scalar option and register it."""
        if kind == "bool":
            widget = QCheckBox()
            widget.setChecked(bool(default))
        elif kind == "choice":
            widget = QComboBox()
            for label, value in extra:
                widget.addItem(label, value)
            widget.setCurrentIndex(max(0, [value for _, value in extra].index(default)))
        elif kind == "int":
            widget = QSpinBox()
            widget.setRange(extra[0], extra[1])
            widget.setValue(int(default))
        elif kind == "int_pair":
            widget = QWidget()
            hbox = QHBoxLayout()
            hbox.setContentsMargins(0, 0, 0, 0)
            widget.setLayout(hbox)
            widget.spins = []
            for value in default:
                spin = QSpinBox()
                spin.setRange(extra[0], extra[1])
                spin.setValue(int(value))
                hbox.addWidget(spin)
                widget.spins.append(spin)
        else:  # float and float_opt
            spin = QDoubleSpinBox()
            spin.setRange(extra[0], extra[1])
            spin.setDecimals(extra[2])
            spin.setSingleStep(extra[3])
            spin.setValue(float(default if default is not None else extra[4]))
            if kind == "float":
                widget = spin
            else:
                widget = QWidget()
                hbox = QHBoxLayout()
                hbox.setContentsMargins(0, 0, 0, 0)
                widget.setLayout(hbox)
                auto = QCheckBox("from logs")
                auto.setToolTip("Leave the key out of the settings file and take the value from the run logs")
                auto.toggled.connect(spin.setDisabled)
                auto.setChecked(True)
                hbox.addWidget(spin)
                hbox.addWidget(auto)
                widget.spin = spin
                widget.auto = auto

        self.global_widgets[key] = (kind, widget)
        return widget

    # ------------------------------------------------------- global options

    def _global_value(self, key):
        """Current value of a scalar option, or None if it should not be written."""
        kind, widget = self.global_widgets[key]
        if kind == "bool":
            return widget.isChecked()
        if kind == "choice":
            return widget.currentData()
        if kind == "int":
            return widget.value()
        if kind == "int_pair":
            return [spin.value() for spin in widget.spins]
        if kind == "float":
            return widget.value()
        return None if widget.auto.isChecked() else widget.spin.value()

    def _set_global_value(self, key, value):
        kind, widget = self.global_widgets[key]
        if kind == "bool":
            widget.setChecked(bool(value))
        elif kind == "choice":
            # can put this back in once fully implemented.
            #if key == "useCalcTheta" and value is True:  # legacy value for the detector angle
            #    value = "detector_angle"
            index = widget.findData(value)
            if index < 0 and isinstance(value, str):
                index = widget.findText(value.lower())
            if index >= 0:
                widget.setCurrentIndex(index)
        elif kind == "int":
            widget.setValue(int(value))
        elif kind == "int_pair":
            for spin, item in zip(widget.spins, value):
                spin.setValue(int(item))
        elif kind == "float":
            widget.setValue(float(value))
        else:
            widget.auto.setChecked(value is None)
            if value is not None:
                widget.spin.setValue(float(value))

    # ------------------------------------------------------------ the table

    COL_NEXUS, COL_RUN, COL_SEQ, COL_TITLE, COL_METHOD, COL_DB = range(6)
    COL_YMIN, COL_YMAX = 6, 7
    COL_BKG = 8  # four columns
    COL_BS, COL_TSHIFT, COL_SCALE = 12, 13, 14

    def _refresh_table(self):
        self._updating = True
        try:
            self.table.setRowCount(len(self.rows))
            for index, row in enumerate(self.rows):
                button = QPushButton(os.path.basename(row.nexus_path) if row.nexus_path else "Load...")
                button.setToolTip(row.nexus_path or "Load the NeXus file for this angle setting")
                button.clicked.connect(lambda _checked=False, i=index: self.load_nexus_for_row(i))
                self.table.setCellWidget(index, self.COL_NEXUS, button)

                method = QComboBox()
                method.addItems(METHODS)
                if row.method in METHODS:
                    method.setCurrentIndex(METHODS.index(row.method))
                method.currentTextChanged.connect(lambda text, i=index: self._method_changed(i, text))
                self.table.setCellWidget(index, self.COL_METHOD, method)

                self._set_item(index, self.COL_RUN, row.run)
                self._set_item(index, self.COL_SEQ, str(row.seq) if row.seq else "", editable=False)
                if row.seq_id:
                    self.table.item(index, self.COL_SEQ).setToolTip(f"Sequence id {row.seq_id}")
                self._set_item(index, self.COL_TITLE, row.title, editable=False)
                self.table.item(index, self.COL_TITLE).setToolTip(row.title)
                db_combo = QComboBox()
                db_combo.setEditable(True)
                db_combo.addItems(self.db_files)
                if row.db_name and row.db_name not in self.db_files:
                    db_combo.addItem(row.db_name)
                db_combo.setCurrentText(row.db_name)
                db_combo.setToolTip(self._db_tooltip(row.db_name))
                db_combo.currentTextChanged.connect(lambda text, i=index: self._db_changed(i, text))
                self.table.setCellWidget(index, self.COL_DB, db_combo)
                self._set_item(index, self.COL_YMIN, str(row.y_min))
                self._set_item(index, self.COL_YMAX, str(row.y_max))
                for offset, value in enumerate(row.bkg):
                    self._set_item(index, self.COL_BKG + offset, str(value))
                self._set_item(index, self.COL_TSHIFT, str(row.theta_shift))
                self._set_item(index, self.COL_SCALE, str(row.scale_factor))

                item = QTableWidgetItem("")
                item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
                item.setCheckState(QtCore.Qt.Checked if row.use_bs else QtCore.Qt.Unchecked)
                self.table.setItem(index, self.COL_BS, item)

            self.table.resizeColumnsToContents()
            self.table.horizontalHeader().setSectionResizeMode(self.COL_TITLE, QHeaderView.Stretch)
            self._widen_db_column()
        finally:
            self._updating = False

    def _widen_db_column(self):
        """Fit the direct beam names, which the editable combo boxes hide otherwise."""
        names = self.db_files + [row.db_name for row in self.rows]
        if not names:
            return
        metrics = self.table.fontMetrics()
        width = max(metrics.horizontalAdvance(name) for name in names) + 45
        self.table.setColumnWidth(self.COL_DB, max(self.table.columnWidth(self.COL_DB), width))

    def _set_item(self, row_index, column, text, editable=True):
        item = QTableWidgetItem(str(text))
        if not editable:
            item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        self.table.setItem(row_index, column, item)

    def _method_changed(self, row_index, text):
        if not self._updating:
            self.rows[row_index].method = text

    def _item_changed(self, item):
        """Write an edited cell back to the row, reverting invalid numbers."""
        if self._updating:
            return
        row = self.rows[item.row()]
        column = item.column()
        text = item.text().strip()
        try:
            if column == self.COL_RUN:
                if text != row.run:
                    row.run = text
                    self.attach_run(row, text)
                    self._refresh_table()
            elif column == self.COL_BS:
                row.use_bs = 1 if item.checkState() == QtCore.Qt.Checked else 0
            elif column == self.COL_YMIN:
                row.y_min = int(text)
            elif column == self.COL_YMAX:
                row.y_max = int(text)
            elif self.COL_BKG <= column < self.COL_BKG + 4:
                row.bkg[column - self.COL_BKG] = int(text)
            elif column == self.COL_TSHIFT:
                row.theta_shift = float(text)
            elif column == self.COL_SCALE:
                row.scale_factor = float(text)
        except ValueError:
            self.status_label.setText(f"'{text}' is not a valid number")
            self._refresh_table()

    def scan_db_files(self):
        """List the direct beam files available in the direct beam directory."""
        directory = self.db_dir_edit.text().strip()
        names, headers = [], {}
        if os.path.isdir(directory):
            for name in sorted(os.listdir(directory)):
                path = os.path.join(directory, name)
                if name.startswith(".") or not os.path.isfile(path):
                    continue
                names.append(name)
                headers[name] = read_db_header(path)
        if names != self.db_files or directory != self.db_directory:
            self.db_files = names
            self.db_headers = headers
            self.db_directory = directory
            self._refresh_table()
            self.status_label.setText(
                f"{len(names)} direct beam file(s) in {directory}" if names
                else f"No direct beam file found in {directory or 'the direct beam directory'}"
            )

    def _db_tooltip(self, name):
        header = self.db_headers.get(name)
        if not header:
            return name if name in self.db_files else f"{name} is not in the direct beam directory"
        details = [f"{key} = {header[key]}" for key in ("db_runs", "db_pixel", "tthd") if key in header]
        return "\n".join([name] + details)

    def _db_changed(self, row_index, name):
        if not self._updating:
            self.rows[row_index].db_name = name
            widget = self.table.cellWidget(row_index, self.COL_DB)
            if widget is not None:
                widget.setToolTip(self._db_tooltip(name))

    def assign_db_files(self):
        """Give the rows the available direct beam files, in order."""
        if not self.db_files:
            QMessageBox.warning(self, "Direct beam files", "No direct beam file found in the direct beam directory")
            return
        if len(self.db_files) != len(self.rows):
            answer = QMessageBox.question(
                self,
                "Direct beam files",
                f"There are {len(self.db_files)} direct beam file(s) for {len(self.rows)} angle setting(s).\n"
                "Assign them in order as far as they go?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if answer != QMessageBox.Yes:
                return
        for row, name in zip(self.rows, self.db_files):
            row.db_name = name
        self._refresh_table()
        self.status_label.setText(f"Assigned {min(len(self.db_files), len(self.rows))} direct beam file(s) in order")

    def _selected_indices(self):
        return sorted({index.row() for index in self.table.selectedIndexes()})

    # --------------------------------------------------------- row actions

    def add_row(self):
        self.rows.append(RunRow())
        self._refresh_table()

    def duplicate_row(self):
        indices = self._selected_indices()
        if not indices:
            return
        for index in reversed(indices):
            source = self.rows[index]
            copy = RunRow(
                method=source.method,
                db_name=source.db_name,
                y_min=source.y_min,
                y_max=source.y_max,
                bkg=list(source.bkg),
                use_bs=source.use_bs,
                theta_shift=source.theta_shift,
                scale_factor=source.scale_factor,
            )
            self.rows.insert(index + 1, copy)
        self._refresh_table()

    def remove_rows(self):
        for index in reversed(self._selected_indices()):
            del self.rows[index]
        self._refresh_table()

    def move_row(self, step):
        indices = self._selected_indices()
        if len(indices) != 1:
            return
        index = indices[0]
        target = index + step
        if not 0 <= target < len(self.rows):
            return
        self.rows[index], self.rows[target] = self.rows[target], self.rows[index]
        self._refresh_table()
        self.table.selectRow(target)

    def sort_by_sequence(self):
        # Rows without a NeXus file yet, and so without a sequence, go last
        self.rows.sort(key=lambda row: (not row.seq_id, row.seq_id, row.seq or 10000, row.run))
        self._refresh_table()
        self._check_sequence()

    def _check_sequence(self):
        sequence_ids = sorted({row.seq_id for row in self.rows if row.seq_id})
        if len(sequence_ids) > 1:
            self.status_label.setText(
                f"The runs belong to different sequences ({', '.join(sequence_ids)}): "
                "a settings file describes the angle settings of one sequence"
            )
            return
        sequences = [row.seq for row in self.rows if row.seq]
        if sequences and sequences != list(range(1, len(self.rows) + 1)):
            self.status_label.setText(
                f"Sequence numbers are {sequences}: the arrays are indexed by sequence number, "
                "so check that the rows are in the right order"
            )

    # -------------------------------------------------------- NeXus loading

    def _nexus_directory(self):
        directory = self.nexus_dir_edit.text().strip()
        return directory if os.path.isdir(directory) else ""

    def load_nexus_for_row(self, row_index):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select the NeXus file for this angle setting", self._nexus_directory(), "NeXus (*.nxs.h5);;All files (*)"
        )
        if not file_path:
            return
        self._apply_nexus(self.rows[row_index], file_path)
        self._refresh_table()

    def add_rows_from_nexus(self):
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "Select the NeXus files of the sequence", self._nexus_directory(), "NeXus (*.nxs.h5);;All files (*)"
        )
        if not file_paths:
            return
        for file_path in file_paths:
            row = RunRow()
            self._apply_nexus(row, file_path)
            self.rows.append(row)
        self.sort_by_sequence()

    def add_rows_from_runs(self):
        """Add a row per run number typed in the run number field."""
        text = self.runs_edit.text().strip()
        if not text:
            return
        try:
            run_numbers = parse_run_list(text)
        except ValueError as error:
            QMessageBox.warning(self, "Run numbers", str(error))
            return
        if not self._nexus_directory():
            QMessageBox.warning(self, "NeXus directory", "Set the experiment or a valid NeXus directory first")
            return

        known = {row.run for row in self.rows}
        missing = []
        for run_number in run_numbers:
            run = str(run_number)
            if run in known:
                continue
            row = RunRow(run=run)
            if not self.attach_run(row, run):
                missing.append(run)
            self.rows.append(row)
        self.runs_edit.clear()
        self.sort_by_sequence()
        if missing:
            self.status_label.setText(f"No NeXus file found for run(s) {', '.join(missing)}")

    def attach_run(self, row, run):
        """Attach the NeXus file of a run number, from the NeXus directory."""
        directory = self._nexus_directory()
        if not directory:
            self.status_label.setText("Set the experiment or a valid NeXus directory to load a run")
            return False
        candidate = os.path.join(directory, f"REF_L_{run}.nxs.h5")
        if not os.path.isfile(candidate):
            self.status_label.setText(f"REF_L_{run}.nxs.h5 was not found in {directory}")
            return False
        self._apply_nexus(row, candidate)
        return True

    def _apply_nexus(self, row, file_path):
        """Fill a row from a NeXus file, and estimate its ROIs if requested."""
        try:
            meta = read_nexus_metadata(file_path)
        except (OSError, KeyError, ValueError) as error:
            QMessageBox.critical(self, "Read error", f"Could not read {file_path}:\n{error}")
            return
        row.nexus_path = file_path
        row.run = meta["run"]
        row.seq = meta["seq_num"]
        row.seq_id = str(meta["seq_id"] or "")
        row.title = str(meta["title"])
        row.profile = None
        # A row that already holds a peak range, typically read from a settings
        # file, keeps it: use "Estimate ROI" to work it out again
        if self.auto_roi_check.isChecked() and row.y_min >= row.y_max:
            self._estimate_roi(row)

    def _estimate_roi(self, row):
        if not row.nexus_path or not os.path.isfile(row.nexus_path):
            return False
        x_range = self._global_value("data_x_range")
        try:
            row.profile = counts_vs_y(row.nexus_path, x_range=x_range, tof_range=row.tof_range)
        except (OSError, KeyError, ValueError) as error:
            self.status_label.setText(f"Could not read events from {os.path.basename(row.nexus_path)}: {error}")
            return False
        y_min, y_max, contrast = estimate_peak_range(row.profile)
        row.y_min, row.y_max = y_min, y_max
        row.bkg = default_bkg_roi(y_min, y_max)
        if contrast < 3:
            self.status_label.setText(
                f"Run {row.run}: weak peak (contrast {contrast:.1f}); open 'Select ROI...' and narrow "
                "the TOF range around the signal to bring it out"
            )
        return True

    def estimate_selected_rois(self):
        indices = self._selected_indices() or range(len(self.rows))
        estimated = 0
        for index in indices:
            if self._estimate_roi(self.rows[index]):
                estimated += 1
        self._refresh_table()
        if estimated:
            self.status_label.setText(f"Estimated the ROIs of {estimated} row(s)")

    def edit_roi(self):
        """Pick the ROIs of the selected row on the counts profiles."""
        if Figure is None:
            QMessageBox.warning(self, "Select ROI", "matplotlib is not available")
            return
        indices = self._selected_indices()
        if not indices:
            QMessageBox.information(self, "Select ROI", "Select a row first")
            return
        row = self.rows[indices[0]]
        if not row.nexus_path or not os.path.isfile(row.nexus_path):
            QMessageBox.warning(self, "Select ROI", "Load the NeXus file of this row first")
            return

        x_range = self._global_value("data_x_range")
        # The file is read here rather than in the dialog, so that a failure to
        # read is the only thing reported as one
        self.status_label.setText(f"Reading the events of run {row.run}...")
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        QApplication.processEvents()
        try:
            events = load_events(row.nexus_path)
        except (OSError, KeyError, ValueError) as error:
            QMessageBox.critical(self, "Read error", f"Could not read {row.nexus_path}:\n{error}")
            return
        finally:
            QApplication.restoreOverrideCursor()
            self.status_label.setText("")

        dialog = ROISelectionDialog(row, x_range, events=events, parent=self)
        if dialog.exec_() != QDialog.Accepted:
            return

        row.y_min, row.y_max, row.bkg, new_x_range, tof_range = dialog.values()
        # Keep the TOF range, so that estimating the ROIs again uses the range
        # the peak was chosen on
        row.tof_range = tof_range if list(tof_range) != list(dialog.chopper_range) else None
        if new_x_range != x_range:
            self._set_global_value("data_x_range", new_x_range)
            for other in self.rows:
                other.profile = None  # computed with the previous x range
            self.status_label.setText(f"The x range of every angle setting is now {new_x_range}")
        row.profile = dialog.y_profile
        self._refresh_table()

    # ---------------------------------------------------- load / save JSON

    def to_settings(self):
        """Assemble the settings dictionary from the interface."""
        settings = {
            "method_per_run": [row.method for row in self.rows],
            "DBname": [row.db_name for row in self.rows],
            "RB_Ymin": [row.y_min for row in self.rows],
            "RB_Ymax": [row.y_max for row in self.rows],
            "BkgROI": [list(row.bkg) for row in self.rows],
            "useBS": [int(row.use_bs) for row in self.rows],
        }
        # Only written when they differ from the defaults applied by the reduction
        if any(row.theta_shift for row in self.rows):
            settings["ThetaShift"] = [row.theta_shift for row in self.rows]
        if any(row.scale_factor != 1 for row in self.rows):
            settings["ScaleFactor"] = [row.scale_factor for row in self.rows]
        if self.save_runs_check.isChecked():
            settings["RBnum"] = [int(row.run) if str(row.run).isdigit() else row.run for row in self.rows]

        for key in GLOBAL_KEYS:
            value = self._global_value(key)
            if value is not None:
                settings[key] = value
        settings.update(self.extra_keys)
        return settings

    def from_settings(self, settings):
        """Fill the interface from a settings dictionary."""
        # Store the original settings for the reset button
        self.original_settings = json.loads(json.dumps(settings))
        
        arrays = {key: settings[key] for key in PER_RUN_KEYS if isinstance(settings.get(key), list)}
        n_rows = max((len(value) for value in arrays.values()), default=0)
        run_numbers = settings.get("RBnum") or []

        self.rows = []
        for index in range(n_rows):
            def entry(key, default):
                values = arrays.get(key, [])
                return values[index] if index < len(values) else default

            bkg = list(entry("BkgROI", [0, 0, 0, 0]))
            bkg = [int(value) for value in (bkg + [0, 0, 0, 0])[:4]]
            self.rows.append(
                RunRow(
                    run=str(run_numbers[index]) if index < len(run_numbers) else "",
                    method=str(entry("method_per_run", "meantheta")).lower(),
                    db_name=str(entry("DBname", "")),
                    y_min=int(entry("RB_Ymin", 0)),
                    y_max=int(entry("RB_Ymax", 0)),
                    bkg=bkg,
                    use_bs=int(entry("useBS", 1)),
                    theta_shift=float(entry("ThetaShift", 0.0)),
                    scale_factor=float(entry("ScaleFactor", 1.0)),
                )
            )
        self.save_runs_check.setChecked(bool(run_numbers))

        for key, _label, _kind, default, _extra in GLOBAL_FIELDS:
            self._set_global_value(key, settings.get(key, default))

        handled = set(PER_RUN_KEYS) | set(GLOBAL_KEYS) | {"RBnum"}
        self.extra_keys = {key: value for key, value in settings.items() if key not in handled}
        self._refresh_table()
        if self.extra_keys:
            self.status_label.setText(
                "Kept unchanged: " + ", ".join(sorted(self.extra_keys)) + " (not editable here)"
            )

    def load_settings(self):
        start = self.settings_file_edit.text().strip() or self._default_settings_dir()
        file_path, _ = QFileDialog.getOpenFileName(self, "Load reduction settings", start, "JSON (*.json);;All files (*)")
        if not file_path:
            return
        try:
            with open(file_path, "r") as handle:
                settings = json.load(handle)
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Load error", f"Could not read {file_path}:\n{error}")
            return
        if not isinstance(settings, dict):
            QMessageBox.critical(self, "Load error", "The settings file must hold a dictionary")
            return
        match = re.search(r"IPTS-\d+", file_path)
        if match and not self.experiment_edit.text().strip():
            self.experiment_edit.setText(match.group(0))
            self._experiment_changed()  # sets the directories, and a default settings path
        self.settings_file_edit.setText(file_path)  # which the file just loaded overrides
        self.from_settings(settings)
        attached = sum(self.attach_run(row, row.run) for row in self.rows if row.run and not row.nexus_path)
        self._refresh_table()
        self.status_label.setText(f"Loaded {file_path}" + (f", with {attached} NeXus file(s)" if attached else ""))

    def save_settings_as(self):
        start = self.settings_file_edit.text().strip() or os.path.join(self._default_settings_dir(), "reduce_settings.json")
        file_path, _ = QFileDialog.getSaveFileName(self, "Save reduction settings", start, "JSON (*.json);;All files (*)")
        if file_path:
            self.settings_file_edit.setText(file_path)
            self.save_settings()

    def save_settings(self):
        file_path = self.settings_file_edit.text().strip()
        if not file_path:
            self.save_settings_as()
            return

        errors, warnings = self.validate()
        if errors:
            QMessageBox.critical(self, "Invalid settings", "\n".join(errors))
            return
        if warnings:
            answer = QMessageBox.question(
                self,
                "Check the settings",
                "\n".join(warnings) + "\n\nSave anyway?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if answer != QMessageBox.Yes:
                return

        try:
            Path(file_path).parent.mkdir(parents=True, exist_ok=True)
            with open(file_path, "w") as handle:
                json.dump(self.to_settings(), handle, indent=2)
        except OSError as error:
            QMessageBox.critical(self, "Save error", f"Could not write {file_path}:\n{error}")
            return
        self._write_user_settings()
        self.status_label.setText(f"Saved {file_path}")

    def _reset_settings(self):
        """Clear all inputs and restore to empty/default state."""
        answer = QMessageBox.question(
            self,
            "Clear all?",
            "Clear all loaded/edited inputs and all angle settings?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No,
        )
        if answer != QMessageBox.Yes:
            return
        
        # Clear all rows from the table
        self.rows = []
        self.original_settings = None
        
        # Reset global options to their defaults
        for key, _label, _kind, default, _extra in GLOBAL_FIELDS:
            self._set_global_value(key, default)
        
        # Reset extra keys
        self.extra_keys = {}
        
        # Reset checkboxes
        self.auto_roi_check.setChecked(True)
        self.save_runs_check.setChecked(False)
        
        # Refresh the table to show nothing
        self._refresh_table()
        self.status_label.setText("Cleared all inputs and angle settings")

    def validate(self):
        """Return the errors that prevent saving and the warnings worth a look."""
        errors, warnings = [], []
        if not self.rows:
            errors.append("Add at least one angle setting")
        for index, row in enumerate(self.rows, start=1):
            label = f"Row {index}" + (f" (run {row.run})" if row.run else "")
            if not row.db_name:
                errors.append(f"{label}: DBname is empty")
            if row.y_min >= row.y_max:
                errors.append(f"{label}: RB_Ymin must be smaller than RB_Ymax")
            if row.method not in METHODS:
                errors.append(f"{label}: unknown method '{row.method}'")
            if sorted(row.bkg) != row.bkg:
                warnings.append(f"{label}: the background ROI values are not in increasing order")
            if row.bkg[1] > row.y_min or row.bkg[2] < row.y_max:
                warnings.append(f"{label}: the background ROI overlaps the peak")
            db_dir = self.db_dir_edit.text().strip()
            if row.db_name and db_dir and not os.path.isfile(os.path.join(db_dir, row.db_name)):
                warnings.append(f"{label}: {row.db_name} was not found in {db_dir}")

        sequences = [row.seq for row in self.rows if row.seq]
        if sequences and sequences != list(range(1, len(self.rows) + 1)):
            warnings.append(
                f"The rows are not in sequence order (sequence numbers {sequences}); "
                "each entry of the arrays is used for the matching sequence number"
            )
        sequence_ids = sorted({row.seq_id for row in self.rows if row.seq_id})
        if len(sequence_ids) > 1:
            warnings.append(
                f"The runs belong to different sequences ({', '.join(sequence_ids)}); "
                "a settings file describes the angle settings of one sequence"
            )
        x_range = self._global_value("data_x_range")
        if x_range[0] >= x_range[1]:
            errors.append("The data x range must be increasing")
        if self._global_value("qmin") >= self._global_value("qmax"):
            errors.append("q min must be smaller than q max")
        return errors, warnings

    # -------------------------------------------------------------- paths

    def _default_settings_dir(self):
        experiment = self.experiment_edit.text().strip()
        if experiment:
            return str(Path("/SNS/REF_L") / experiment / "shared" / "autoreduce")
        return ""

    def _experiment_changed(self):
        """
        Point the directories at the experiment.

        The paths follow the experiment whenever it changes, since that is the
        point of the field. Nothing is touched when the experiment is only
        confirmed again, so a directory edited by hand survives.
        """
        experiment = self.experiment_edit.text().strip()
        if not experiment or experiment == self.experiment:
            return
        self.experiment = experiment
        base = Path("/SNS/REF_L") / experiment
        self.nexus_dir_edit.setText(str(base / "nexus"))
        self.db_dir_edit.setText(str(base / "shared" / "transmission"))
        self.settings_file_edit.setText(str(base / "shared" / "autoreduce" / "reduce_settings.json"))
        self.scan_db_files()
        self.status_label.setText(f"Directories set to {base}")

    def _browse_directory(self, line_edit, on_change=None):
        directory = QFileDialog.getExistingDirectory(self, "Select a directory", line_edit.text().strip())
        if directory:
            line_edit.setText(directory)
            if on_change is not None:
                on_change()

    def _read_user_settings(self):
        self.experiment_edit.setText(self.settings.value("json_builder_experiment_id", ""))
        self.experiment = self.experiment_edit.text().strip()
        self.nexus_dir_edit.setText(self.settings.value("json_builder_nexus_dir", ""))
        self.db_dir_edit.setText(self.settings.value("json_builder_db_dir", ""))
        self.settings_file_edit.setText(self.settings.value("json_builder_settings_file", ""))

    def _write_user_settings(self):
        self.settings.setValue("json_builder_experiment_id", self.experiment_edit.text())
        self.settings.setValue("json_builder_nexus_dir", self.nexus_dir_edit.text())
        self.settings.setValue("json_builder_db_dir", self.db_dir_edit.text())
        self.settings.setValue("json_builder_settings_file", self.settings_file_edit.text())


def main():
    import sys

    app = QApplication(sys.argv)
    window = JSONSettingsBuilderTab()
    window.resize(1200, 800)
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()


__all__ = ["JSONSettingsBuilderTab"]
