"""The ROI pop-out of the settings editor (roi-popout-dialog), lifted from upstream PR #197.

Source: neutrons/LiquidsReflectometer PR #197, ``launcher/apps/json_settings_builder.py`` as hardened at
``agentic/feature/harden-review-branch`` @ 65c83d9 (upstream ``exp-json-settings-builder`` @ 3ce5e20 differs in the
dialog only by 65c83d9's ``parse_math=False`` title): ``_move_span`` (:396-406) and ``ROISelectionDialog``
(:409-690). Authors of the lifted lines: Mathieu Doucet (f5513c7, 8191e49, 1e692c7, ab22307) and welbournR
(3f74d41). The lift is commit b041aa2, verbatim; the changes since are this slug's.

The dialog is a view over one run's events (``lr_reduction.roi_estimate.RunEvents``) and one Angles row's values.
It reads no file and writes no document: the tab's slot resolves the run, loads its events once and applies what
:meth:`ROISelectionDialog.changes` reports. Every number it draws comes from ``lr_reduction.roi_estimate``: the
images, the profiles, the background bands the reducer averages, the estimate and the default background. No
detector geometry, distance or band arithmetic lives here.
"""

import functools
import traceback

import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter
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

from lr_reduction import roi_estimate

try:
    from matplotlib.figure import Figure
    from matplotlib.widgets import SpanSelector

    try:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
    except ImportError:  # matplotlib < 3.5
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
except ImportError:  # the editor works without the ROI plots; its button says why it is unavailable
    Figure = None

#: A spin box's "not set": one below every pixel index, shown as the words.
UNSET = -1


def _move_span(patch, low, high, vertical=True):
    """
    Move a shaded region, whose other extent is in axes coordinates.

    ``axvspan`` returns a rectangle from matplotlib 3.9 on and a polygon
    before, and the two do not take the same coordinates. A vertical span
    (``axvspan``) covers ``low``-``high`` along the data x axis; a horizontal
    one (``axhspan``, the bands drawn on the images, whose Y is vertical) along
    the data y axis.
    """
    if vertical:
        if hasattr(patch, "set_bounds"):  # a rectangle
            patch.set_bounds(low, 0, high - low, 1)
        else:  # a polygon
            patch.set_xy([[low, 0], [low, 1], [high, 1], [high, 0], [low, 0]])
    elif hasattr(patch, "set_bounds"):
        patch.set_bounds(0, low, 1, high - low)
    else:
        patch.set_xy([[0, low], [1, low], [1, high], [0, high], [0, low]])


def _guarded(method):
    """A slot that raises reports into the status line instead of leaving the slot: an exception out of a PyQt slot
    reaches qFatal() and aborts the launcher (L3). Typing a value passes through intermediate ones, and none may be
    fatal."""

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        try:
            return method(self, *args, **kwargs)
        except Exception as exc:  # noqa: BLE001 -- the point is to catch everything
            traceback.print_exc()
            self.status.setText(f"{type(exc).__name__}: {exc}")
            return None

    return wrapper


def _plain_log_ticks(axis):
    """Log tick labels as plain text ("1e+02"). The default formatter writes mathtext, which a font without the
    glyphs sends into matplotlib's fallback recursion (F7, L6); no ``parse_math`` setting reaches a tick label."""
    axis.set_major_formatter(LogFormatter())
    axis.set_minor_formatter(LogFormatter(labelOnlyBase=True))


class ROISelectionDialog(QDialog):
    """
    One Angles row's peak, background and the shared x range, shown on the run's detector images and profiles, and
    adjusted by dragging on the profiles or typing.

    The peak and the background belong to the row. The x pixel range is ``data_x_range``, shared by every angle.
    The TOF range is a view filter only: it selects the events the profiles and the XY image are made of (at the
    highest angles nearly every count outside the reflected signal is background, and narrowing the TOF range
    brings the peak out), and it is never reported. The reduction's own TOF window (``tof_min``/``tof_max``) is
    drawn when the row has one, and edited in the table.
    """

    PEAK, BACKGROUND_LEFT, BACKGROUND_RIGHT = "peak", "left", "right"

    def __init__(self, events, values, title="", tof_band=None, parent=None):
        QDialog.__init__(self, parent)
        self.events = events
        self.n_x, self.n_y = int(events.n_x), int(events.n_y)
        self._updating = False
        self._bkg_edited = False
        self._notes = []

        stride = f" (every {events.stride}th event)" if getattr(events, "stride", 1) > 1 else ""
        self.setWindowTitle(f"{title}: peak, background and ranges{stride}" if title else "Select the ranges")
        self.resize(1000, 1000)

        if len(events.tof):
            self.tof_edges = roi_estimate.tof_edges(events)
        else:  # no events: empty plots over a nominal span; Estimate refuses
            self.tof_edges = np.array([0.0, 1.0])
        self._opening = self._opening_values(values, tof_band)

        layout = QVBoxLayout()
        self.setLayout(layout)
        # A layout engine would run on every redraw, which is too slow while dragging
        self.figure = Figure(figsize=(9, 9))
        grid = self.figure.add_gridspec(4, 2, height_ratios=[2.2, 1.6, 1, 1])
        self.xy_axis = self.figure.add_subplot(grid[0, 0])
        self.ytof_axis = self.figure.add_subplot(grid[0, 1])
        self.y_axis = self.figure.add_subplot(grid[1, :])
        self.tof_axis = self.figure.add_subplot(grid[2, :])
        self.x_axis = self.figure.add_subplot(grid[3, :])
        self.figure.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.06, hspace=0.5, wspace=0.35)
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(NavigationToolbar(self.canvas, self))
        layout.addWidget(self.canvas, stretch=1)
        layout.addWidget(self._build_controls())
        self.status = QLabel("")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)

        # parse_math=False: the title is the NeXus run title, arbitrary text from the file. matplotlib parses
        # `$...$` as math, so a title carrying a literal `$` reaches the mathtext parser (B12).
        self.xy_axis.set_title(f"{title}{stride}", parse_math=False, fontsize="small")
        self._build_images()
        self._build_profiles()
        self._build_overlays()

        self._set_log_scale(True)
        self._update(images=False)
        self._reset_limits()
        self._report_states()

        # Dragging on a profile sets the range it shows; on the Y profile, the range the radio buttons choose.
        # useblit keeps a drag from redrawing the whole figure at every mouse move.
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

    # ------------------------------------------------------- opening values

    def _opening_values(self, values, tof_band):
        """The row's values as the spins hold them, and why any of them cannot be shown (B5, types table)."""
        def pixel(value, top):
            return int(value) if isinstance(value, (int, np.integer)) and 0 <= value <= top else UNSET

        y_min, y_max = pixel(values.get("RB_Ymin"), self.n_y - 1), pixel(values.get("RB_Ymax"), self.n_y - 1)
        self._bkg_entry = values.get("BkgROI")
        self._bkg_reason = None
        bkg = [UNSET] * 4
        if self._bkg_entry is None:
            self._bkg_reason = "no background is set for this angle"
        else:
            try:
                if UNSET in (y_min, y_max) and self._adjacent(self._bkg_entry):
                    raise ValueError("the background is adjacent to the peak (two zero bounds), and the peak is not set")
                bands = roi_estimate.background_bands(self._bkg_entry, y_min, y_max)
                bkg = [int(round(v)) for band in bands for v in band]
            except ValueError as exc:
                self._bkg_reason = str(exc)
        x_range = values.get("data_x_range")
        if (isinstance(x_range, (list, tuple)) and len(x_range) == 2
                and all(isinstance(v, (int, np.integer)) and 0 <= v < self.n_x for v in x_range)):
            x_values = [int(x_range[0]), int(x_range[1])]
        else:
            x_values = [0, self.n_x - 1]
            self._notes.append(f"data_x_range {x_range!r} is not two pixels on this detector; showing every X "
                               f"pixel (written only if you change it)")
        band = tof_band if tof_band is not None else (self.tof_edges[0], self.tof_edges[-1])
        tof = [int(np.floor(band[0])), int(np.ceil(band[1]))]
        window = (values.get("tof_min"), values.get("tof_max"))
        self._tof_window = window if None not in window else None
        self._use_bs = values.get("useBS") is not False and values.get("useBS") != 0
        return {"peak": [y_min, y_max], "bkg": bkg, "x": x_values, "tof": tof}

    @staticmethod
    def _adjacent(entry):
        """True for an entry in the "adjacent to the peak" form, whose two zero bounds the reducer replaces with the
        peak's edges."""
        return isinstance(entry, (list, tuple)) and any(isinstance(v, (int, float)) and v == 0 for v in entry)

    # ------------------------------------------------------------ controls

    def _build_controls(self):
        """Radio buttons to choose what a drag sets, and the values themselves."""
        box = QWidget()
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        box.setLayout(grid)

        grid.addWidget(QLabel("Drag on the Y profile to set:"), 0, 0)
        modes = QHBoxLayout()
        self.mode_buttons = []
        for label, mode in [("the peak", self.PEAK),
                            ("the low background", self.BACKGROUND_LEFT),
                            ("the high background", self.BACKGROUND_RIGHT)]:
            button = QRadioButton(label)
            self.mode_buttons.append((button, mode))
            modes.addWidget(button)
        self.mode_buttons[0][0].setChecked(True)
        modes.addStretch(1)
        grid.addLayout(modes, 0, 1, 1, 4)

        opening = self._opening
        grid.addWidget(QLabel("Peak (RB_Ymin, RB_Ymax):"), 1, 0)
        self.peak_spins = [self._make_spin(value, self.n_y - 1) for value in opening["peak"]]
        for column, spin in enumerate(self.peak_spins):
            grid.addWidget(spin, 1, 1 + column)

        grid.addWidget(QLabel("Background (BkgROI):"), 2, 0)
        self.bkg_spins = [self._make_spin(value, self.n_y - 1, background=True) for value in opening["bkg"]]
        for column, spin in enumerate(self.bkg_spins):
            grid.addWidget(spin, 2, 1 + column)

        grid.addWidget(QLabel("x range, all angles (data_x_range):"), 3, 0)
        self.x_spins = [self._make_spin(value, self.n_x - 1, unset=False) for value in opening["x"]]
        for column, spin in enumerate(self.x_spins):
            grid.addWidget(spin, 3, 1 + column)

        grid.addWidget(QLabel("TOF range [us], view filter only:"), 4, 0)
        low, high = int(np.floor(self.tof_edges[0])), int(np.ceil(self.tof_edges[-1]))
        self.tof_spins = []
        for value in opening["tof"]:
            spin = self._make_spin(min(max(value, low), high), high, unset=False, minimum=low)
            spin.setSingleStep(100)
            spin.setToolTip("Only the events of this TOF range make the Y and X profiles and the XY image; "
                            "never written (the reduction's TOF window is edited in the table)")
            self.tof_spins.append(spin)
        for column, spin in enumerate(self.tof_spins):
            grid.addWidget(spin, 4, 1 + column)

        self.estimate_button = QPushButton("Estimate the peak")
        self.estimate_button.setToolTip("Set the peak and a default background from the Y profile shown")
        self.estimate_button.clicked.connect(self._estimate)
        grid.addWidget(self.estimate_button, 5, 0)

        self.log_check = QCheckBox("Log scale")
        self.log_check.setChecked(True)
        self.log_check.toggled.connect(self._set_log_scale)
        grid.addWidget(self.log_check, 5, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        self.ok_button = buttons.button(QDialogButtonBox.Ok)
        grid.addWidget(buttons, 5, 3, 1, 2)
        return box

    def _make_spin(self, value, maximum, unset=True, background=False, minimum=0):
        spin = QSpinBox()
        spin.setRange(UNSET if unset else minimum, maximum)
        if unset:
            spin.setSpecialValueText("not set")
        spin.setValue(int(value))
        spin.valueChanged.connect(self._background_edited if background else self._values_changed)
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
        if mode != self.PEAK:
            self._bkg_edited = True
        self._apply_range(spins, low, high)

    def _x_range_selected(self, low, high):
        self._apply_range(self.x_spins, low, high)

    def _tof_range_selected(self, low, high):
        self._apply_range(self.tof_spins, low, high)

    def _apply_range(self, spins, low, high):
        low, high = sorted((int(round(low)), int(round(high))))
        bottom = 0 if spins[0].minimum() == UNSET else spins[0].minimum()  # a drag never sets "not set"
        low = max(bottom, min(spins[0].maximum(), low))
        high = max(low, min(spins[1].maximum(), high))
        self._updating = True
        try:
            spins[0].setValue(low)
            spins[1].setValue(high)
        finally:
            self._updating = False
        self._values_changed()

    @_guarded
    def _background_edited(self, *_value):
        if not self._updating:
            self._bkg_edited = True
            self._values_changed()

    @_guarded
    def _values_changed(self, *_value):
        if not self._updating:
            self._update()
            self._report_states()

    @_guarded
    def _estimate(self, *_checked):
        """B7: the data layer's estimate on the Y profile shown, and its default background. A refusal is a
        message and changes no value; a peak with no room for a background sets the peak and says so."""
        try:
            y_min, y_max = roi_estimate.estimate_peak_range(self.y_profile)
        except roi_estimate.CannotEstimateError as exc:
            self.status.setText(f"Estimate refused: {exc}")
            return
        try:
            background = list(roi_estimate.default_bkg_roi((y_min, y_max), self.n_y))
            note = None
        except ValueError as exc:
            background = [UNSET] * 4
            note = f"Estimate set the peak; background not set: {exc}"
        self._updating = True
        try:
            self.peak_spins[0].setValue(int(y_min))
            self.peak_spins[1].setValue(int(y_max))
            for spin, value in zip(self.bkg_spins, background):
                spin.setValue(int(value))
            self._bkg_edited = True
        finally:
            self._updating = False
        self._update()
        self._reset_limits()
        self._report_states()
        if note:
            self.status.setText(note)

    # ------------------------------------------------------------ plotting

    def _build_images(self):
        """B3: the web report's two detector images, as the data layer computes them; log colour, colorbars."""
        x_range, tof_filter = self._spin_values(self.x_spins), self._spin_values(self.tof_spins)
        self._xy_filter = tuple(tof_filter)
        self._ytof_x_range = tuple(x_range)
        xy = roi_estimate.xy_image(self.events, tof_band=self._xy_filter)
        ytof = roi_estimate.y_tof_image(self.events, self._ytof_x_range, self.tof_edges)
        common = {"origin": "lower", "aspect": "auto", "interpolation": "nearest", "cmap": "viridis"}
        self.xy_image = self.xy_axis.imshow(
            xy, extent=(-0.5, self.n_x - 0.5, -0.5, self.n_y - 0.5), norm=self._norm(xy), **common)
        self.ytof_image = self.ytof_axis.imshow(
            ytof, extent=(self.tof_edges[0], self.tof_edges[-1], -0.5, self.n_y - 0.5), norm=self._norm(ytof), **common)
        self.xy_axis.set_xlabel("x pixel")
        self.xy_axis.set_ylabel("y pixel")
        self.ytof_axis.set_xlabel("TOF [us]")
        self.ytof_axis.set_ylabel("y pixel")
        self.ytof_axis.set_title("Y vs TOF, inside the x range", fontsize="small")
        for image, axis in ((self.xy_image, self.xy_axis), (self.ytof_image, self.ytof_axis)):
            colorbar = self.figure.colorbar(image, ax=axis)
            # through the colorbar, which re-applies its own formatters when it draws
            colorbar.formatter = LogFormatter()
            colorbar.minorformatter = LogFormatter(labelOnlyBase=True)

    @staticmethod
    def _norm(array):
        return LogNorm(vmin=1, vmax=max(2.0, float(np.max(array)) if array.size else 2.0))

    def _build_profiles(self):
        self.y_axis.set_xlabel("y pixel (reflectivity direction)")
        self.x_axis.set_xlabel("x pixel (low resolution direction)")
        self.tof_axis.set_xlabel("time of flight [us], every event of the run")
        for axis in (self.y_axis, self.tof_axis, self.x_axis):
            axis.set_ylabel("counts")
        (self.y_line,) = self.y_axis.plot([], [], drawstyle="steps-mid", color="tab:blue")
        (self.tof_line,) = self.tof_axis.plot([], [], drawstyle="steps-mid", color="tab:blue")
        (self.x_line,) = self.x_axis.plot([], [], drawstyle="steps-mid", color="tab:blue")

    def _build_overlays(self):
        """The overlays are made once and moved (B11). Each is drawn on every plot that has its axis (B4):
        ``overlays[name][axis]`` is the artist, a span along the data x axis on a profile and along the data y axis
        on an image's Y."""
        ax = {"y_axis": self.y_axis, "xy_axis": self.xy_axis, "ytof_axis": self.ytof_axis,
              "x_axis": self.x_axis, "tof_axis": self.tof_axis}
        style = {
            "peak": {"color": "tab:green", "alpha": 0.25, "label": "peak"},
            "bkg_low": {"color": "tab:red", "alpha": 0.2, "label": "background"},
            "bkg_high": {"color": "tab:red", "alpha": 0.2},
            "x_range": {"color": "tab:green", "alpha": 0.25, "label": "x range"},
            "tof_filter": {"color": "tab:blue", "alpha": 0.12, "label": "TOF view filter"},
            "tof_window": {"color": "tab:orange", "alpha": 0.2, "label": "reduction TOF window"},
        }
        on = {
            "peak": [("y_axis", True), ("xy_axis", False), ("ytof_axis", False)],
            "bkg_low": [("y_axis", True), ("xy_axis", False), ("ytof_axis", False)],
            "bkg_high": [("y_axis", True), ("xy_axis", False), ("ytof_axis", False)],
            "x_range": [("x_axis", True), ("xy_axis", True)],
            "tof_filter": [("tof_axis", True), ("ytof_axis", True)],
            "tof_window": [("tof_axis", True), ("ytof_axis", True)],
        }
        self._vertical = {}
        self.overlays = {}
        for name, legs in on.items():
            self.overlays[name] = {}
            for axis_name, vertical in legs:
                axis = ax[axis_name]
                image = axis_name in ("xy_axis", "ytof_axis")
                kwargs = dict(style[name])
                if image:
                    kwargs.pop("label", None)
                    kwargs["alpha"] = min(0.35, kwargs["alpha"] + 0.1)
                span = axis.axvspan(0, 1, **kwargs) if vertical else axis.axhspan(0, 1, **kwargs)
                self.overlays[name][axis_name] = span
                self._vertical[(name, axis_name)] = vertical
        if self._tof_window is not None:
            for axis_name, artist in self.overlays["tof_window"].items():
                _move_span(artist, float(self._tof_window[0]), float(self._tof_window[1]),
                           self._vertical[("tof_window", axis_name)])
        else:
            for artist in self.overlays["tof_window"].values():
                artist.set_visible(False)
        for axis in (self.y_axis, self.tof_axis, self.x_axis):
            axis.legend(loc="upper right", fontsize="small")

    def _spin_values(self, spins):
        return [spin.value() for spin in spins]

    def _peak(self):
        y_min, y_max = self._spin_values(self.peak_spins)
        return (y_min, y_max) if UNSET not in (y_min, y_max) and y_min <= y_max else None

    def _bands(self):
        """The background bands the reducer will average, or None: from the spins once the background has been
        edited, else from the row's own entry and the current peak (an untouched ``[a, b, 0, 0]`` follows the peak,
        as the reducer's does)."""
        if self._bkg_edited:
            bkg = self._spin_values(self.bkg_spins)
            if UNSET in bkg:
                return None
            entry = bkg
        else:
            if self._bkg_reason is not None:
                return None
            entry = self._bkg_entry
        peak = self._peak()
        if peak is None:
            if not self._bkg_edited and self._adjacent(entry):
                return None  # the adjacent form needs the peak's edges
            peak = (UNSET, UNSET)  # four explicit bounds do not use the peak
        try:
            return roi_estimate.background_bands(entry, *peak)
        except ValueError:
            return None

    def _reversed(self):
        """The ranges that are reversed now. Typing "190" into an x bound passes through 1 and 19: a step on the way
        to a value, never data. The plots keep their last state and OK waits."""
        return [name for name, spins in (("the x range", self.x_spins), ("the TOF filter", self.tof_spins))
                if spins[0].value() > spins[1].value()]

    def _update(self, images=True):
        """Recompute the profiles for the current ranges and move the overlays. An image is recomputed only when
        its own input changed: the XY image's view filter, the Y-TOF image's x range (B11)."""
        if self._reversed():
            return
        x_range = tuple(self._spin_values(self.x_spins))
        tof_filter = tuple(self._spin_values(self.tof_spins))
        peak = self._peak()
        if images and tof_filter != self._xy_filter:
            self._xy_filter = tof_filter
            self.xy_image.set_data(roi_estimate.xy_image(self.events, tof_band=tof_filter))
        if images and x_range != self._ytof_x_range:
            self._ytof_x_range = x_range
            self.ytof_image.set_data(roi_estimate.y_tof_image(self.events, x_range, self.tof_edges))

        self.y_profile = roi_estimate.profile_y(self.events, x_range, tof_band=tof_filter)
        x_profile = roi_estimate.profile_x(self.events, y_range=peak, tof_band=tof_filter)
        # The TOF profile keeps the whole span, so that the signal can be found again after a narrow filter
        tof_profile = roi_estimate.profile_tof(self.events, self.tof_edges, x_range=x_range, y_range=peak)
        self.y_line.set_data(np.arange(self.n_y), self.y_profile)
        self.x_line.set_data(np.arange(self.n_x), x_profile)
        self.tof_line.set_data((self.tof_edges[:-1] + self.tof_edges[1:]) / 2.0, tof_profile)

        self._move("peak", peak)
        bands = self._bands()
        self._move("bkg_low", bands[0] if bands else None)
        self._move("bkg_high", bands[1] if bands else None)
        self._move("x_range", x_range)
        self._move("tof_filter", tof_filter)

        for axis in (self.y_axis, self.tof_axis, self.x_axis):
            axis.relim()
            axis.autoscale_view(scalex=False)
        self.canvas.draw_idle()

    def _move(self, name, edges):
        for axis_name, artist in self.overlays[name].items():
            if edges is None:
                artist.set_visible(False)
                continue
            artist.set_visible(True)
            _move_span(artist, float(edges[0]), float(edges[1]), self._vertical[(name, axis_name)])

    @_guarded
    def _set_log_scale(self, log_scale):
        for axis in (self.y_axis, self.tof_axis, self.x_axis):
            axis.set_yscale("log" if log_scale else "linear")
            if log_scale:
                _plain_log_ticks(axis.yaxis)
        self.canvas.draw_idle()

    def _reset_limits(self):
        """Show the region around the ROIs, or around the peak of the profile if there is none yet."""
        edges = [value for value in self._spin_values(self.peak_spins + self.bkg_spins) if value != UNSET]
        if not edges:
            edges = [int(np.argmax(self.y_profile))] if np.any(self.y_profile) else [self.n_y // 2]
        self.y_axis.set_xlim(max(0, min(edges) - 30), min(self.n_y - 1, max(edges) + 30))
        self.tof_axis.set_xlim(self.tof_edges[0], self.tof_edges[-1])
        self.x_axis.set_xlim(0, self.n_x - 1)
        self.canvas.draw_idle()

    def _report_states(self):
        """The status line (B5) and whether OK is available (B10)."""
        notes = list(self._notes) + [f"{name} is reversed" for name in self._reversed()]
        if self._bands() is None:
            if self._bkg_edited:
                notes.append("Background (BkgROI): not set — give all four bounds, or none")
            elif self._bkg_reason is not None:
                notes.append(f"Background (BkgROI): not set — {self._bkg_reason}")
        elif not self._use_bs:
            notes.append("Background drawn but not subtracted: useBS is off for this angle")
        self.status.setText("; ".join(notes))
        self.ok_button.setEnabled(self._valid())

    def _valid(self):
        """B10: OK needs a peak, low <= high, and a background either untouched or four ascending non-zero
        integers (never pixel 0, the reducer's sentinel), or not set at all."""
        if self._peak() is None or self._reversed():
            return False
        if not self._bkg_edited:
            return True
        bkg = self._spin_values(self.bkg_spins)
        if all(value == UNSET for value in bkg):
            return True
        return UNSET not in bkg and 0 not in bkg and bkg == sorted(bkg)

    # ------------------------------------------------------------- results

    def changes(self):
        """B9: the fields whose values differ from those the dialog opened with, as the document holds them.
        ``RB_Ymin``/``RB_Ymax`` (ints), ``BkgROI`` (four ascending ints, only once the background was edited) and
        ``data_x_range`` (two ints, for every angle). The TOF view filter is never among them (B8)."""
        out = {}
        opening = self._opening
        for name, value, before in zip(("RB_Ymin", "RB_Ymax"), self._spin_values(self.peak_spins), opening["peak"]):
            if value != before and value != UNSET:
                out[name] = int(value)
        bkg = self._spin_values(self.bkg_spins)
        if self._bkg_edited and UNSET not in bkg and bkg != opening["bkg"]:
            out["BkgROI"] = [int(v) for v in bkg]
        x_range = self._spin_values(self.x_spins)
        if x_range != opening["x"]:
            out["data_x_range"] = [int(v) for v in x_range]
        return out

    def reject(self):
        """Cancel: the opening values come back, so a cancelled dialog reports nothing (V6)."""
        self._updating = True
        try:
            for spins, key in ((self.peak_spins, "peak"), (self.bkg_spins, "bkg"), (self.x_spins, "x"),
                               (self.tof_spins, "tof")):
                for spin, value in zip(spins, self._opening[key]):
                    spin.setValue(int(value))
            self._bkg_edited = False
        finally:
            self._updating = False
        QDialog.reject(self)
