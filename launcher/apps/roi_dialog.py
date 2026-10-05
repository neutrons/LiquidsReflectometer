"""The ROI pop-out of the settings editor (roi-popout-dialog), lifted from upstream PR #197.

Source: neutrons/LiquidsReflectometer PR #197, ``launcher/apps/json_settings_builder.py`` as hardened at
``agentic/feature/harden-review-branch`` @ 65c83d9 (upstream ``exp-json-settings-builder`` @ 3ce5e20 differs in the
dialog only by 65c83d9's ``parse_math=False`` title): ``_move_span`` (:396-406) and ``ROISelectionDialog``
(:409-690). Authors of the lifted lines: Mathieu Doucet (f5513c7, 8191e49, 1e692c7, ab22307) and welbournR
(3f74d41).

This commit lifts those lines verbatim and wires nothing. The names they take from their own module (``N_Y``,
``N_X``, ``load_events``, ``counts_in_range``, ``estimate_peak_range``, ``default_bkg_roi``) are deliberately not
lifted: the next commits re-seat the dialog on ``lr_reduction.roi_estimate``, the data layer, so no geometry, event
reading or band arithmetic lives in the launcher.
"""
# ruff: noqa: F821 -- the lifted dialog's own-module names; the next commit replaces them with the data layer

import numpy as np
from qtpy.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QRadioButton,
    QSpinBox,
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
        # parse_math=False: `title` is the NeXus run title, arbitrary text from
        # the file. matplotlib parses `$...$` as math, so a title carrying a
        # literal `$` reaches the mathtext parser — and a symbol no font
        # provides sends `_mathtext._get_glyph` into its fallback chain. Plot
        # text built from data should never be parsed as markup.
        self.y_axis.set_title(
            f"{title} (sequence {row.seq})" if row.seq else title, parse_math=False
        )
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
