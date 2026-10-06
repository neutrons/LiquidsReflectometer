"""The ROI pop-out (roi-popout-dialog): #197's dialog re-seated on ``lr_reduction.roi_estimate``.

Events are built in memory through the real ``RunEvents``, never a mock, and no test reads the data submodule (A4).
Gestures go through QTest: clicks, typed values, and a press-move-release on the canvas at coordinates taken from
the axes' ``transData`` (measured, not reasoned: L2). ``changes()``, what OK reports, is read directly: it is the
dialog's API.
"""

import inspect

import numpy as np
import pytest
from matplotlib.colors import LogNorm
from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtTest import QTest

from launcher.apps import roi_dialog
from launcher.apps.roi_dialog import ROISelectionDialog
from lr_reduction import roi_estimate
from lr_reduction.roi_estimate import RunEvents

pytestmark = pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")

N_X, N_Y = 256, 304
PEAK = (136, 146)
UNSET = roi_dialog.UNSET


def make_events(n_x=N_X, n_y=N_Y, peak=PEAK, n=40000, seed=0):
    """A run: a peak at `peak` rows over x 100-159, a flat background, TOF uniform over 10-40 ms."""
    rng = np.random.default_rng(seed)
    x_lo, x_hi = min(100, n_x // 2), min(160, n_x)
    y = np.r_[rng.integers(peak[0], peak[1] + 1, n), rng.integers(0, n_y, n // 4)]
    x = np.r_[rng.integers(x_lo, x_hi, n), rng.integers(0, n_x, n // 4)]
    tof = rng.uniform(10000.0, 40000.0, len(y))
    return RunEvents(x=x, y=y, tof=tof, n_x=n_x, n_y=n_y)


def featureless_events():
    """Events spread evenly over the detector: no peak, so the estimator refuses. Measured: "contrast 1.09 is below
    1.5 — the detector is featureless here". The plan's sparse run is not refused: 300 events (seed 1) are bracketed
    at (138, 141) with contrast 4.0, the estimator's call (src/ is out of this slug's scope)."""
    rng = np.random.default_rng(2)
    n = 200000
    return RunEvents(x=rng.integers(0, N_X, n), y=rng.integers(0, N_Y, n), tof=rng.uniform(1e4, 4e4, n),
                     n_x=N_X, n_y=N_Y)


def row_values(**changes):
    values = {"RB_Ymin": PEAK[0], "RB_Ymax": PEAK[1], "BkgROI": [133, 149, 0, 0], "data_x_range": [50, 200],
              "tof_min": None, "tof_max": None, "useBS": True}
    return {**values, **changes}


def make_dialog(events=None, title="Si Ir Air-221472-2.", tof_band=None, **changes):
    return ROISelectionDialog(events if events is not None else make_events(), row_values(**changes), title=title,
                              tof_band=tof_band)


def close(dialog):
    dialog.close()
    dialog.deleteLater()


def type_into(spin, value):
    """Type a value into a spin box as a user does: select its text, type, Return."""
    spin.setFocus()
    spin.lineEdit().selectAll()
    QTest.keyClicks(spin.lineEdit(), str(value))
    QTest.keyClick(spin.lineEdit(), QtCore.Qt.Key_Return)


def span_edges(artist, vertical):
    """The data edges a span artist covers: along x for a vertical span (axvspan), along y otherwise."""
    if vertical:
        return artist.get_x(), artist.get_x() + artist.get_width()
    return artist.get_y(), artist.get_y() + artist.get_height()


def mathtext(dialog):
    """Draw, then every label on every axes (the three profiles, the two images and both colorbars) that carries
    mathtext: major and minor tick labels and the offset text, on both axes of each."""
    dialog.canvas.draw()
    texts = []
    for axes in dialog.figure.axes:
        for axis in (axes.xaxis, axes.yaxis):
            texts += [label.get_text() for label in axis.get_ticklabels() + axis.get_ticklabels(minor=True)]
            texts.append(axis.get_offset_text().get_text())
    return [text for text in texts if "$" in text]


def visible_bands(dialog):
    return [axis for name in ("bkg_low", "bkg_high") for axis, artist in dialog.overlays[name].items()
            if artist.get_visible()]


# Which overlays sit on which axes, and whether they span the data x axis there (vertical) or the y axis.
ON = {
    "peak": [("y_axis", True), ("xy_axis", False), ("ytof_axis", False)],
    "bkg_low": [("y_axis", True), ("xy_axis", False), ("ytof_axis", False)],
    "bkg_high": [("y_axis", True), ("xy_axis", False), ("ytof_axis", False)],
    "x_range": [("x_axis", True), ("xy_axis", True)],
    "tof_filter": [("tof_axis", True), ("ytof_axis", True)],
}


def test_the_dialog_draws_two_images_and_three_profiles():
    """V1 (B3): the XY and Y-TOF images, then counts per Y, per TOF and per X."""
    dialog = make_dialog()
    assert len(dialog.xy_axis.images) == 1 and len(dialog.ytof_axis.images) == 1
    for axis in (dialog.y_axis, dialog.tof_axis, dialog.x_axis):
        assert len(axis.lines) == 1
    close(dialog)


def test_the_images_are_the_data_layers_arrays():
    """V2″ (B3, L1, B12): each image's array is the data layer's, unmodified: get_array() is the data, whatever the
    display does.
    - origin is "lower", and the extent puts each pixel's centre on its index.
    - The aspect is left to the data ("auto", after open and after a draw): forced to "equal", the Y-TOF image is a
      one-pixel sliver (I-48, B-2).
    - The colour scale is logarithmic (a LogNorm on each image).
    - Each image has its own colorbar, beside it: seven axes in all, the XY colorbar between the two images and the
      Y-TOF colorbar right of its own (I-50, B-2 and B-3). Without the colorbars, V11′'s colorbar leg reads nothing."""
    events = make_events()
    band = (15000.0, 30000.0)
    dialog = make_dialog(events, tof_band=band)
    xy = dialog.xy_axis.images[0]
    np.testing.assert_array_equal(np.asarray(xy.get_array()), roi_estimate.xy_image(events, tof_band=band))
    assert xy.origin == "lower" and list(xy.get_extent()) == [-0.5, N_X - 0.5, -0.5, N_Y - 0.5]
    edges = roi_estimate.tof_edges(events)
    ytof = dialog.ytof_axis.images[0]
    np.testing.assert_array_equal(np.asarray(ytof.get_array()), roi_estimate.y_tof_image(events, (50, 200), edges))
    assert ytof.origin == "lower" and list(ytof.get_extent()) == [edges[0], edges[-1], -0.5, N_Y - 0.5]
    for drawn in (False, True):
        if drawn:
            dialog.canvas.draw()
        assert [axes.get_aspect() for axes in (dialog.xy_axis, dialog.ytof_axis)] == ["auto", "auto"], drawn
        assert all(isinstance(image.norm, LogNorm) for image in (xy, ytof)), drawn
    assert len(dialog.figure.axes) == 7
    bars = [image.colorbar for image in (xy, ytof)]
    assert all(bar is not None and bar.ax in dialog.figure.axes for bar in bars)
    xy_box, ytof_box = dialog.xy_axis.get_position(), dialog.ytof_axis.get_position()
    assert xy_box.x1 <= bars[0].ax.get_position().x0 <= ytof_box.x0, "the XY colorbar sits between the two images"
    assert ytof_box.x1 <= bars[1].ax.get_position().x0, "the Y-TOF colorbar sits right of its own image"
    close(dialog)


def test_each_plot_is_the_data_layers_for_the_ranges_shown():
    """B3, B8 (acceptance 5: "what is drawn is what the reducer uses"): after typed changes to the x range, the peak and
    the view filter, each plot is the data layer's for the ranges shown. The Y profile is over the x range, through
    the filter; the X profile over the peak rows, through the filter; the TOF profile over the x range and the peak
    rows, across the whole span; the XY image through the filter; the Y-TOF image over the x range, across the whole
    span."""
    events = make_events()
    dialog = make_dialog(events)
    for spin, value in zip(dialog.x_spins + dialog.peak_spins + dialog.tof_spins, (80, 170, 138, 144, 15000, 32000)):
        type_into(spin, value)
    x_range, peak, band = (80, 170), (138, 144), (15000, 32000)
    edges = roi_estimate.tof_edges(events)
    np.testing.assert_array_equal(dialog.y_line.get_ydata(), roi_estimate.profile_y(events, x_range, tof_band=band))
    np.testing.assert_array_equal(dialog.x_line.get_ydata(), roi_estimate.profile_x(events, y_range=peak, tof_band=band))
    np.testing.assert_array_equal(dialog.tof_line.get_ydata(),
                                  roi_estimate.profile_tof(events, edges, x_range=x_range, y_range=peak))
    np.testing.assert_array_equal(dialog.tof_line.get_xdata(), (edges[:-1] + edges[1:]) / 2.0)
    np.testing.assert_array_equal(np.asarray(dialog.xy_axis.images[0].get_array()),
                                  roi_estimate.xy_image(events, tof_band=band))
    np.testing.assert_array_equal(np.asarray(dialog.ytof_axis.images[0].get_array()),
                                  roi_estimate.y_tof_image(events, x_range, edges))
    close(dialog)


def test_the_background_overlay_is_what_the_reducer_averages():
    """V3 (B4, F6): [133, 149, 0, 0] with the peak at 136-146 is two bands, (133, 136) and (146, 149), as
    background_bands (the reducer's sorter) gives them -- never one band across the peak -- on every plot with a Y
    axis."""
    dialog = make_dialog()
    (b0, b1), (b2, b3) = roi_estimate.background_bands([133, 149, 0, 0], *PEAK)
    for name, expected in (("bkg_low", (b0, b1)), ("bkg_high", (b2, b3))):
        for axis, vertical in ON[name]:
            artist = dialog.overlays[name][axis]
            assert artist.get_visible() and span_edges(artist, vertical) == expected, (name, axis)
    close(dialog)


def canvas_point(dialog, axis, x):
    """The canvas widget's point at data x and mid-height of `axis` (Qt's y runs down, matplotlib's up)."""
    canvas = dialog.canvas
    ratio = canvas.devicePixelRatioF() or 1.0
    x_display = axis.transData.transform((x, 0))[0]
    y_display = axis.transAxes.transform((0, 0.5))[1]
    return QtCore.QPoint(int(round(x_display / ratio)), int(round(canvas.height() - y_display / ratio)))


def drag(dialog, axis, x0, x1):
    """Press at data x0, move to x1, release, at mid-height of `axis`, on the canvas widget.

    The press and release are QTest's. The moves are the QMouseEvents the platform delivers, sent to the canvas:
    QTest.mouseMove delivers no move to an offscreen widget (measured: the canvas saw button_press_event and
    button_release_event, no motion_notify_event, and the selector set nothing). Sent this way, matplotlib sees two
    motion_notify_events and the SpanSelector its drag."""
    canvas = dialog.canvas
    canvas.draw()
    QTest.mousePress(canvas, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier, canvas_point(dialog, axis, x0))
    for x in ((x0 + x1) / 2, x1):
        move = QtGui.QMouseEvent(QtCore.QEvent.MouseMove, QtCore.QPointF(canvas_point(dialog, axis, x)),
                                 QtCore.Qt.NoButton, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier)
        QtWidgets.QApplication.sendEvent(canvas, move)
    QTest.mouseRelease(canvas, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier, canvas_point(dialog, axis, x1))
    QtWidgets.QApplication.processEvents()


def click(dialog, axis, x):
    """A press and a release at data x, mid-height of `axis`, with no move between."""
    canvas = dialog.canvas
    canvas.draw()
    point = canvas_point(dialog, axis, x)
    QTest.mousePress(canvas, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier, point)
    QTest.mouseRelease(canvas, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier, point)
    QtWidgets.QApplication.processEvents()


def shown(dialog, mode="peak"):
    """The dialog on screen (a drag needs a laid-out canvas), with the radio button for `mode` chosen."""
    dialog.resize(950, 1000)
    dialog.show()
    QTest.qWaitForWindowExposed(dialog)
    button = next(button for button, button_mode in dialog.mode_buttons if button_mode == mode)
    QTest.mouseClick(button, QtCore.Qt.LeftButton)
    return dialog


@pytest.mark.parametrize("mode, span, reported", [
    ("peak", (120, 128), {"RB_Ymin": 120, "RB_Ymax": 128}),
    ("left", (120, 128), {"BkgROI": [120, 128, 146, 149]}),
    ("right", (150, 160), {"BkgROI": [133, 136, 150, 160]}),
])
def test_dragging_on_the_y_profile_sets_the_chosen_range(mode, span, reported):
    """V4 (B6): with the radio button chosen, a drag on the Y profile sets that range's spins and its overlay, and OK
    reports it. A click afterwards is not a drag: once a span has been drawn, matplotlib reports a click as a span of
    zero width, and the range stays."""
    dialog = shown(make_dialog(), mode)
    spins = {"peak": dialog.peak_spins, "left": dialog.bkg_spins[:2], "right": dialog.bkg_spins[2:]}[mode]
    overlay = {"peak": "peak", "left": "bkg_low", "right": "bkg_high"}[mode]
    drag(dialog, dialog.y_axis, *span)
    assert [spin.value() for spin in spins] == list(span)
    assert span_edges(dialog.overlays[overlay]["y_axis"], True) == span
    click(dialog, dialog.y_axis, 170)
    assert [spin.value() for spin in spins] == list(span)
    assert dialog.changes() == reported
    close(dialog)


def test_dragging_on_the_x_and_tof_profiles_sets_the_range_and_the_filter():
    """V4′ (B6, B8): a drag on the X profile sets data_x_range's spins and its overlay on every axes with an X axis, and
    OK reports it. A drag on the TOF profile sets the view filter, the Y profile and the XY image follow it, and OK
    reports nothing more. A TOF drag lands within a canvas pixel of its ends (here about 30 us), so the filter is read
    back from the spins."""
    events = make_events()
    dialog = shown(make_dialog(events))
    drag(dialog, dialog.x_axis, 70, 180)
    assert [spin.value() for spin in dialog.x_spins] == [70, 180]
    for axis, vertical in ON["x_range"]:
        assert span_edges(dialog.overlays["x_range"][axis], vertical) == (70, 180), axis
    assert dialog.changes() == {"data_x_range": [70, 180]}
    drag(dialog, dialog.tof_axis, 15000, 30000)
    band = tuple(spin.value() for spin in dialog.tof_spins)
    per_pixel = abs(dialog.tof_axis.transData.inverted().transform((1, 0))[0]
                    - dialog.tof_axis.transData.inverted().transform((0, 0))[0])
    assert abs(band[0] - 15000) <= 2 * per_pixel and abs(band[1] - 30000) <= 2 * per_pixel, (band, per_pixel)
    np.testing.assert_array_equal(np.asarray(dialog.xy_axis.images[0].get_array()),
                                  roi_estimate.xy_image(events, tof_band=band))
    np.testing.assert_array_equal(dialog.y_line.get_ydata(), roi_estimate.profile_y(events, (70, 180), tof_band=band))
    assert dialog.changes() == {"data_x_range": [70, 180]}
    dialog.canvas.draw()
    assert [axes.get_aspect() for axes in (dialog.xy_axis, dialog.ytof_axis)] == ["auto", "auto"]  # A-ii: after a drag too
    close(dialog)


def test_a_drag_past_the_detector_edge_stops_at_the_edge():
    """B6: on a profile zoomed out past the detector, a drag that starts or ends beyond its rows stops at the edge
    row, never at "not set"; a drag wholly off the detector changes nothing and says so."""
    dialog = shown(make_dialog())
    dialog.y_axis.set_xlim(-40, 340)
    drag(dialog, dialog.y_axis, -20, 10)
    assert [spin.value() for spin in dialog.peak_spins] == [0, 10]
    drag(dialog, dialog.y_axis, 290, 330)
    assert [spin.value() for spin in dialog.peak_spins] == [290, N_Y - 1]
    drag(dialog, dialog.y_axis, -35, -10)
    assert [spin.value() for spin in dialog.peak_spins] == [290, N_Y - 1]
    assert "off the detector" in dialog.status.text(), dialog.status.text()
    close(dialog)


@pytest.mark.parametrize("name, axis, vertical", [(n, a, v) for n, legs in ON.items() for a, v in legs])
def test_a_typed_value_moves_its_overlay_on_every_plot(name, axis, vertical):
    """V5 (B6): typing into an entry box moves that overlay, on each plot that has its axis, to the typed edges."""
    dialog = make_dialog(BkgROI=[120, 125, 160, 165])
    spins = {"peak": dialog.peak_spins, "bkg_low": dialog.bkg_spins[:2], "bkg_high": dialog.bkg_spins[2:],
             "x_range": dialog.x_spins, "tof_filter": dialog.tof_spins}[name]
    new = {"peak": (130, 150), "bkg_low": (110, 118), "bkg_high": (170, 180), "x_range": (60, 190),
           "tof_filter": (12000, 31000)}[name]
    for spin, value in zip(spins, new):
        type_into(spin, value)
    assert span_edges(dialog.overlays[name][axis], vertical) == pytest.approx(new)
    close(dialog)


def test_the_reductions_tof_window_is_drawn_when_the_row_has_one():
    """B4 (A2): the row's tof_min/tof_max are drawn on the TOF profile and the Y-TOF image, and never reported. A row
    without both draws none; a pair that is not two ascending numbers draws none and says so."""
    dialog = make_dialog(tof_min=12000.0, tof_max=30000)
    for axis in ("tof_axis", "ytof_axis"):
        artist = dialog.overlays["tof_window"][axis]
        assert artist.get_visible() and span_edges(artist, True) == pytest.approx((12000.0, 30000.0)), axis
    assert dialog.changes() == {}
    close(dialog)
    for window, noted in (((None, None), False), ((12000.0, None), False), (("12000", 30000.0), True),
                          ((30000.0, 12000.0), True), ((float("nan"), 30000.0), True), ((12000.0, float("inf")), True)):
        dialog = make_dialog(tof_min=window[0], tof_max=window[1])
        assert not any(artist.get_visible() for artist in dialog.overlays["tof_window"].values()), window
        assert ("TOF window" in dialog.status.text()) is noted, (window, dialog.status.text())
        close(dialog)


def press(dialog, which):
    """Press the dialog's own Ok or Cancel button, as a user does."""
    QTest.mouseClick(dialog.buttons.button(getattr(QtWidgets.QDialogButtonBox, which)), QtCore.Qt.LeftButton)


def test_cancel_reports_nothing():
    """V6′: edited, then Cancel pressed: the dialog is rejected and reports no change. Also for a background whose
    bound is not a whole row, which the spins show rounded."""
    for edits, background in (((("peak_spins", 0, 140), ("x_spins", 1, 180)), [133, 149, 0, 0]),
                              ((("bkg_spins", 3, 150),), [133.5, 149, 0, 0])):
        dialog = make_dialog(BkgROI=background)
        for spins, index, value in edits:
            type_into(getattr(dialog, spins)[index], value)
        press(dialog, "Cancel")
        assert dialog.result() == QtWidgets.QDialog.Rejected, background
        assert dialog.changes() == {}, background
        close(dialog)


def test_ok_reports_only_what_changed():
    """V7′ (B9): OK, pressed, reports only the fields whose values differ from those the dialog opened with. An
    untouched [a, b, 0, 0] background, or an unset one, is not reported, so the table keeps it byte for byte."""
    for background in ([133, 149, 0, 0], None):
        dialog = make_dialog(BkgROI=background)
        press(dialog, "Ok")
        assert dialog.result() == QtWidgets.QDialog.Accepted and dialog.changes() == {}, background
        close(dialog)
    dialog = make_dialog()
    type_into(dialog.peak_spins[0], 140)
    press(dialog, "Ok")
    assert dialog.result() == QtWidgets.QDialog.Accepted and dialog.changes() == {"RB_Ymin": 140}
    close(dialog)
    dialog = make_dialog()
    for spin, value in zip(dialog.bkg_spins, (120, 125, 160, 165)):
        type_into(spin, value)
    type_into(dialog.x_spins[0], 60)
    press(dialog, "Ok")
    assert dialog.result() == QtWidgets.QDialog.Accepted
    assert dialog.changes() == {"BkgROI": [120, 125, 160, 165], "data_x_range": [60, 200]}
    close(dialog)


def test_an_edited_background_is_reported_when_the_rows_own_would_draw_otherwise():
    """V7 (B4, B9): a peak edge at "not set" is never reported. A background touched and put back is not reported; put
    back and then the peak moved, it is, as drawn: the row's [133, 149, 0, 0] would follow the new peak, and the
    dialog drew the bounds edited here."""
    dialog = make_dialog()
    type_into(dialog.peak_spins[1], UNSET)
    assert "RB_Ymax" not in dialog.changes()
    close(dialog)
    dialog = make_dialog()
    type_into(dialog.bkg_spins[0], 130)
    type_into(dialog.bkg_spins[0], 133)
    assert dialog.changes() == {}
    type_into(dialog.peak_spins[0], 138)
    type_into(dialog.peak_spins[1], 144)
    assert span_edges(dialog.overlays["bkg_low"]["y_axis"], True) == (133, 136)
    assert dialog.changes() == {"RB_Ymin": 138, "RB_Ymax": 144, "BkgROI": [133, 136, 146, 149]}
    close(dialog)


def test_the_view_filter_is_never_reported():
    """V8 (B8, A2): the TOF entry is a view filter: narrowing it changes the profiles, and OK reports nothing about
    TOF (tof_min/tof_max are the table's)."""
    dialog = make_dialog()
    before = dialog.y_axis.lines[0].get_ydata().copy()
    type_into(dialog.tof_spins[0], 20000)
    type_into(dialog.tof_spins[1], 25000)
    assert not np.array_equal(dialog.y_axis.lines[0].get_ydata(), before)
    dialog.accept()
    assert dialog.changes() == {}
    close(dialog)


def test_an_estimate_refusal_is_a_message_not_a_guess():
    """V9 (B7): a featureless run's profile is refused (CannotEstimateError). The status line carries the refusal,
    and no value changes: no bracket around noise."""
    dialog = make_dialog(featureless_events())
    before = [spin.value() for spin in dialog.peak_spins + dialog.bkg_spins]
    QTest.mouseClick(dialog.estimate_button, QtCore.Qt.LeftButton)
    assert "refus" in dialog.status.text().lower() or "cannot" in dialog.status.text().lower(), dialog.status.text()
    assert [spin.value() for spin in dialog.peak_spins + dialog.bkg_spins] == before
    close(dialog)


def test_estimate_sets_a_background_the_reducer_accepts():
    """V10 (B7, F5): Estimate sets the peak and a default background from the data layer, which refuses rather than
    clamping to row 0: the reported BkgROI is four ascending non-zero ints that background_bands accepts."""
    dialog = make_dialog(RB_Ymin=None, RB_Ymax=None, BkgROI=None)
    QTest.mouseClick(dialog.estimate_button, QtCore.Qt.LeftButton)
    dialog.accept()
    changes = dialog.changes()
    bkg = changes["BkgROI"]
    assert len(bkg) == 4 and 0 not in bkg and bkg == sorted(bkg)
    roi_estimate.background_bands(bkg, changes["RB_Ymin"], changes["RB_Ymax"])
    close(dialog)


def test_estimate_brings_its_peak_into_view():
    """B7: the Y profile's view moves to the estimated peak, wherever the row's own peak was."""
    dialog = make_dialog(RB_Ymin=20, RB_Ymax=30, BkgROI=None)
    assert dialog.y_axis.get_xlim()[1] < PEAK[0]
    QTest.mouseClick(dialog.estimate_button, QtCore.Qt.LeftButton)
    low, high = dialog.y_axis.get_xlim()
    assert low <= dialog.peak_spins[0].value() and dialog.peak_spins[1].value() <= high, (low, high)
    close(dialog)


@pytest.mark.parametrize("background", [None, [120, 125, 160, 165]], ids=["unset", "four bounds"])
def test_estimate_near_a_detector_edge_sets_the_peak_and_leaves_the_background(background):
    """B7 (F5; failure matrix, edge): the data layer finds no room for a background beside a peak at rows 2-8.
    Estimate sets the peak it found, leaves the background as it was (unset, or the row's four bounds) and says there
    is no room: never a 0 bound."""
    dialog = make_dialog(make_events(peak=(2, 8)), RB_Ymin=None, RB_Ymax=None, BkgROI=background)
    before = [spin.value() for spin in dialog.bkg_spins]
    QTest.mouseClick(dialog.estimate_button, QtCore.Qt.LeftButton)
    found = roi_estimate.estimate_peak_range(dialog.y_profile)
    assert [spin.value() for spin in dialog.peak_spins] == list(found)
    assert [spin.value() for spin in dialog.bkg_spins] == before
    assert "no room" in dialog.status.text(), dialog.status.text()
    assert dialog.changes() == {"RB_Ymin": found[0], "RB_Ymax": found[1]}
    close(dialog)


def test_file_text_and_log_ticks_never_reach_the_math_parser():
    """V11′ (B12, F7, L6): the run title is file text, so it is not parsed as maths. No label on the log profiles or the
    colorbars carries mathtext ("$"): major and minor tick labels and the offset text. That holds after a draw at the
    opening limits, and again with a profile zoomed inside one decade (what the toolbar's zoom does), where minor
    labels exist. A sparse run's colorbars span less than a decade, so their minor labels exist too."""
    dialog = make_dialog(title="$\\foo$ run")
    assert dialog.xy_axis.title.get_parse_math() is False
    assert not mathtext(dialog), mathtext(dialog)[:5]
    dialog.y_axis.set_ylim(20, 80)
    assert not mathtext(dialog), mathtext(dialog)[:5]
    close(dialog)
    sparse = make_events(n=600)
    dialog = make_dialog(sparse)
    assert dialog.xy_axis.images[0].get_array().max() < 10  # the colorbar spans less than a decade
    assert not mathtext(dialog), mathtext(dialog)[:5]
    close(dialog)


def test_geometry_follows_the_events_not_a_literal():
    """V12 (B3, F3): a 128 x 200 detector gives spin maxima 199 (Y) and 127 (X), and images of that shape."""
    events = make_events(n_x=128, n_y=200, peak=(90, 100))
    dialog = make_dialog(events, RB_Ymin=90, RB_Ymax=100, BkgROI=[85, 87, 103, 105], data_x_range=[10, 100])
    assert dialog.peak_spins[0].maximum() == 199 and dialog.bkg_spins[3].maximum() == 199
    assert dialog.x_spins[1].maximum() == 127
    assert dialog.xy_axis.images[0].get_array().shape == (200, 128)
    assert dialog.ytof_axis.images[0].get_array().shape[0] == 200
    close(dialog)


def test_a_sampled_run_says_so_in_its_titles():
    """Failure matrix, long run: a run read one event in N (load_event_pixels' stride) says so in the window title and
    over the XY image, so that its counts are not taken for the run's totals."""
    events = make_events()
    sampled = RunEvents(x=events.x, y=events.y, tof=events.tof, n_x=N_X, n_y=N_Y, stride=7)
    dialog = make_dialog(sampled)
    assert "1 event in 7" in dialog.windowTitle() and "1 event in 7" in dialog.xy_axis.get_title()
    close(dialog)
    dialog = make_dialog(events)
    assert "event in" not in dialog.windowTitle() and "event in" not in dialog.xy_axis.get_title()
    close(dialog)


def test_a_run_without_events_draws_empty_plots_and_estimate_refuses(recwarn):
    """Types table, events: a run without events draws empty plots over a nominal TOF span, and Estimate refuses. An
    empty profile is a state, not a problem: no warning from autoscaling a log axis over zero counts."""
    empty = RunEvents(x=np.array([], dtype=np.int64), y=np.array([], dtype=np.int64), tof=np.array([], dtype=float),
                      n_x=N_X, n_y=N_Y)
    dialog = make_dialog(empty)
    dialog.canvas.draw()
    assert not [str(w.message) for w in recwarn if issubclass(w.category, UserWarning)]
    assert not np.any(dialog.y_line.get_ydata()) and not np.any(dialog.xy_axis.images[0].get_array())
    QTest.mouseClick(dialog.estimate_button, QtCore.Qt.LeftButton)
    assert "Estimate refused: no counts" in dialog.status.text(), dialog.status.text()
    close(dialog)


def test_a_nudge_moves_artists_and_rebuilds_nothing():
    """V13 (B11, R8): changing a value moves artists; the axes, images, lines and overlays are the same objects
    afterwards, and each image axes still holds one image. An image is recomputed only when its own input changes:
    the XY image for the view filter, the Y-TOF image for the x range, neither for the peak."""
    dialog = make_dialog()

    def objects():
        return (list(dialog.figure.axes), [axis.images[0] for axis in (dialog.xy_axis, dialog.ytof_axis)],
                [axis.lines[0] for axis in (dialog.y_axis, dialog.tof_axis, dialog.x_axis)],
                {name: dict(artists) for name, artists in dialog.overlays.items()})

    def arrays():
        return [axis.images[0].get_array() for axis in (dialog.xy_axis, dialog.ytof_axis)]

    before, (xy, ytof) = objects(), arrays()
    for spin, value, recomputed in ((dialog.peak_spins[0], 138, ()), (dialog.tof_spins[0], 20000, ("xy",)),
                                    (dialog.x_spins[0], 60, ("ytof",))):
        type_into(spin, value)
        after, (new_xy, new_ytof) = objects(), arrays()
        assert len(after[0]) == len(before[0]) and all(a is b for a, b in zip(before[0], after[0])), value
        assert all(a is b for a, b in zip(before[1] + before[2], after[1] + after[2])), value
        assert all(before[3][name][axis] is after[3][name][axis] for name in before[3] for axis in before[3][name])
        assert len(dialog.xy_axis.images) == 1 and len(dialog.ytof_axis.images) == 1, value
        assert (new_xy is not xy) is ("xy" in recomputed), (value, "xy")
        assert (new_ytof is not ytof) is ("ytof" in recomputed), (value, "ytof")
        xy, ytof = new_xy, new_ytof
    close(dialog)


@pytest.mark.parametrize("background, reason", [
    ([0, 10, 150, 160], "sentinel"), ([121, 130], "four bounds"), ("120, 130", "flat list"), (None, "not set"),
    ([], "four bounds"), ([120, 125, 160], "four bounds"), ([0, 0, 0, 160], "sentinel"), ([0, 0, 0, 0], "sentinel")])
def test_an_unusable_background_is_shown_as_not_set_with_its_reason(background, reason):
    """V14 (B5): an entry the reducer cannot use, or none at all, is "not set" in the status line with the reason. No
    bands are drawn and no values are invented: the spins say "not set". The dialog does not write the row's own
    entry, so it does not hold OK back."""
    dialog = make_dialog(BkgROI=background)
    assert "not set" in dialog.status.text() and reason in dialog.status.text(), dialog.status.text()
    assert dialog.ok_button.isEnabled()
    assert visible_bands(dialog) == []
    assert [spin.value() for spin in dialog.bkg_spins] == [UNSET] * 4
    assert [spin.text() for spin in dialog.bkg_spins] == ["not set"] * 4
    close(dialog)


def test_an_adjacent_background_follows_the_peak_once_it_is_set():
    """V14 (B4, B5): [133, 149, 0, 0] needs the peak's edges. With no peak it is "not set" with that reason; typing the
    peak draws the two bands the reducer will average for it, the note goes, and OK leaves the entry as it is."""
    dialog = make_dialog(RB_Ymin=None, RB_Ymax=None)
    assert "adjacent to the peak" in dialog.status.text() and visible_bands(dialog) == []
    type_into(dialog.peak_spins[0], 136)
    type_into(dialog.peak_spins[1], 146)
    for name, band in zip(("bkg_low", "bkg_high"), roi_estimate.background_bands([133, 149, 0, 0], 136, 146)):
        for axis, vertical in ON[name]:
            assert span_edges(dialog.overlays[name][axis], vertical) == band, (name, axis)
    assert dialog.status.text() == "" and dialog.ok_button.isEnabled()
    assert dialog.changes() == {"RB_Ymin": 136, "RB_Ymax": 146}
    close(dialog)


def test_four_bounds_are_drawn_before_a_peak_is_set():
    """V14 (B4): four explicit bounds do not use the peak, so they are drawn while it is not set."""
    dialog = make_dialog(RB_Ymin=None, BkgROI=[120, 125, 160, 165])
    assert span_edges(dialog.overlays["bkg_low"]["y_axis"], True) == (120, 125)
    assert span_edges(dialog.overlays["bkg_high"]["y_axis"], True) == (160, 165)
    close(dialog)


def empty_column_use_bs():
    """What the dialog is given for a row of a document whose useBS column is [] (the reducer's default, on): the
    reachable route to None."""
    from lr_reduction.settings_document import SettingsDocument

    doc = SettingsDocument()
    doc.add_angle(RB_Ymin=PEAK[0], RB_Ymax=PEAK[1], BkgROI=[133, 149, 0, 0])
    doc.set("useBS", [])
    return doc.angle_row(0)["useBS"]


@pytest.mark.parametrize("use_bs, subtracted", [(False, False), (0, False), (True, True), (None, True), ("[]", True)],
                         ids=["False", "0", "True", "None", "an empty column"])
def test_a_background_the_reducer_does_not_subtract_is_drawn_and_labelled(use_bs, subtracted):
    """V14, V16 (B5, the types table's useBS row): the bands are drawn in every case. Only with useBS off for the row
    (False or 0) does the status line say they are not subtracted. None, and the None an empty column gives, are on,
    as the reducer reads them."""
    if use_bs == "[]":
        use_bs = empty_column_use_bs()
        assert use_bs is None
    dialog = make_dialog(useBS=use_bs)
    assert dialog.overlays["bkg_low"]["y_axis"].get_visible()
    assert ("not subtracted" in dialog.status.text()) is not subtracted, dialog.status.text()
    close(dialog)


@pytest.mark.parametrize("name, value", [("RB_Ymin", 400), ("RB_Ymin", -3), ("RB_Ymin", "136"), ("RB_Ymin", True),
                                         ("RB_Ymin", 136.5), ("RB_Ymax", 304), ("RB_Ymax", float("inf"))])
def test_a_peak_edge_that_is_not_a_row_of_the_detector_is_not_set_with_a_note(name, value):
    """B5, the types table: a peak edge this detector cannot show (past its rows, negative, text, a bool, a fraction)
    starts "not set", and the status line names the value; OK waits for a peak."""
    dialog = make_dialog(**{name: value})
    spin = dialog.peak_spins[0 if name == "RB_Ymin" else 1]
    assert spin.value() == UNSET and spin.text() == "not set"
    assert f"{name} {value!r} is not a row of this detector" in dialog.status.text(), dialog.status.text()
    assert not dialog.ok_button.isEnabled()
    close(dialog)


def test_whole_floats_are_pixels():
    """A JSON 136.0 is row 136, as the data layer's _whole reads it, and [50.0, 200.0] is that x range: nothing is
    noted, and nothing is reported."""
    dialog = make_dialog(RB_Ymin=136.0, data_x_range=[50.0, 200.0])
    assert dialog.peak_spins[0].value() == 136 and [spin.value() for spin in dialog.x_spins] == [50, 200]
    assert dialog.status.text() == ""
    assert dialog.changes() == {}
    close(dialog)


@pytest.mark.parametrize("x_range", [[50], [50, 300], [200, 50], [50.5, 200], "50,200", None, [True, 200]],
                         ids=["one value", "past the detector", "reversed", "a fraction", "text", "None", "a bool"])
def test_a_data_x_range_that_is_not_two_pixels_is_noted_and_written_only_if_changed(x_range):
    """Types table, data_x_range: a value that is not two ascending pixels of this detector is noted in the status
    line; the x spins start at the whole detector, and the x range is reported only once it is changed here."""
    dialog = make_dialog(data_x_range=x_range)
    assert [spin.value() for spin in dialog.x_spins] == [0, N_X - 1]
    assert f"data_x_range {x_range!r}" in dialog.status.text(), dialog.status.text()
    assert dialog.changes() == {}
    type_into(dialog.x_spins[0], 60)
    assert dialog.changes() == {"data_x_range": [60, N_X - 1]}
    close(dialog)


def test_an_unset_peak_draws_no_overlay_until_both_edges_are_set():
    """V15′ (types table, RB_Ymin/RB_Ymax None): with an edge "not set", no peak overlay is drawn on any plot. It
    appears on all three once both edges are set, at the spins' edges."""
    def visible(dialog):
        return sorted(axis for axis, artist in dialog.overlays["peak"].items() if artist.get_visible())

    dialog = make_dialog(RB_Ymin=None, RB_Ymax=155)
    assert visible(dialog) == []
    type_into(dialog.peak_spins[0], 145)
    assert visible(dialog) == ["xy_axis", "y_axis", "ytof_axis"]
    for axis, vertical in ON["peak"]:
        assert span_edges(dialog.overlays["peak"][axis], vertical) == (145, 155), axis
    close(dialog)
    dialog = make_dialog(RB_Ymin=None, RB_Ymax=None)
    type_into(dialog.peak_spins[0], 145)
    assert visible(dialog) == []
    type_into(dialog.peak_spins[1], 155)
    assert visible(dialog) == ["xy_axis", "y_axis", "ytof_axis"]
    close(dialog)


def test_ok_is_unavailable_for_an_inverted_peak_or_a_partial_background():
    """V15 (B10): OK is unavailable while the peak is empty or inverted, or the background is partly set."""
    dialog = make_dialog()
    assert dialog.ok_button.isEnabled()
    type_into(dialog.peak_spins[0], 150)  # above RB_Ymax 146
    assert not dialog.ok_button.isEnabled() and "the peak is reversed" in dialog.status.text()
    type_into(dialog.peak_spins[0], 136)
    assert dialog.ok_button.isEnabled()
    for spin, value in zip(dialog.bkg_spins, (120, 125, 160, 165)):
        type_into(spin, value)
    dialog.bkg_spins[3].setValue(UNSET)  # "not set": a partly set background
    assert not dialog.ok_button.isEnabled()
    close(dialog)
    dialog = make_dialog(RB_Ymin=None)
    assert not dialog.ok_button.isEnabled() and "OK waits for a peak" in dialog.status.text()
    close(dialog)


@pytest.mark.parametrize("bounds, words", [
    ((UNSET, 125, 160, 165), "give all four bounds"),
    ((120, 125, 160, UNSET), "give all four bounds"),
    ((0, 5, 160, 165), "sentinel"),
    ((0, 0, 160, 165), "sentinel"),
    ((125, 120, 160, 165), "must ascend"),
], ids=["first unset", "last unset", "one 0", "two 0s", "out of order"])
def test_a_background_edited_into_one_the_dialog_cannot_write_waits_with_the_reason(bounds, words):
    """V15 (B10): a background edited here is written as four ascending non-zero pixels. Anything else (partly set,
    a 0 bound, the reducer's sentinel, or bounds out of order) draws no bands and says why, OK waits, and it is never
    reported."""
    dialog = make_dialog()
    for spin, value in zip(dialog.bkg_spins, bounds):
        type_into(spin, value)
    assert not dialog.ok_button.isEnabled()
    assert "OK waits" in dialog.status.text() and words in dialog.status.text(), dialog.status.text()
    assert visible_bands(dialog) == []
    assert "BkgROI" not in dialog.changes()
    close(dialog)


def test_a_background_cleared_here_is_the_rows_own_again():
    """V15 (B9, B10): a background edited and then cleared to "not set" in all four is the row's own again: its bands
    are drawn as the reducer will average them, the status line says the row keeps it, OK is available, and nothing
    about the background is reported."""
    dialog = make_dialog()
    for spin, value in zip(dialog.bkg_spins, (120, 125, 160, 165)):
        type_into(spin, value)
    for spin in dialog.bkg_spins:
        type_into(spin, UNSET)
    assert dialog.ok_button.isEnabled()
    assert "the row keeps its own" in dialog.status.text(), dialog.status.text()
    for name, band in zip(("bkg_low", "bkg_high"), roi_estimate.background_bands([133, 149, 0, 0], *PEAK)):
        assert span_edges(dialog.overlays[name]["y_axis"], True) == band, name
    assert dialog.changes() == {}
    close(dialog)


@pytest.mark.parametrize("steps, name", [
    ((("x_spins", 0, 220),), "the x range"),
    ((("tof_spins", 1, 20000), ("tof_spins", 0, 30000)), "the TOF filter"),
], ids=["x", "TOF"])
def test_a_reversed_range_is_named_and_ok_waits(steps, name):
    """B10, L3: a range typed past its other end (typing "190" passes through 1 and 19 on the way) is a step, never
    data: it is named in the status line, OK waits, the plots keep their last state, and nothing raises."""
    dialog = make_dialog()
    for spins, index, value in steps:
        type_into(getattr(dialog, spins)[index], value)
    assert f"{name} is reversed" in dialog.status.text(), dialog.status.text()
    assert "Error" not in dialog.status.text(), dialog.status.text()
    assert not dialog.ok_button.isEnabled()
    overlay = {"the x range": ("x_range", "x_axis"), "the TOF filter": ("tof_filter", "tof_axis")}[name]
    low, high = span_edges(dialog.overlays[overlay[0]][overlay[1]], True)
    assert low <= high
    close(dialog)


@pytest.mark.parametrize("slot", ["_values_changed", "_estimate", "_set_log_scale"])
def test_an_error_inside_a_slot_is_a_status_line_not_an_abort(monkeypatch, slot):
    """L3: an exception out of a PyQt slot reaches qFatal() and aborts the launcher. Each of the dialog's slots shows
    the error in its status line instead; the error is injected into what each one calls."""
    dialog = make_dialog()

    def injected(*_a, **_k):
        raise RuntimeError("injected")

    if slot == "_values_changed":
        monkeypatch.setattr(roi_dialog.roi_estimate, "profile_y", injected)
        type_into(dialog.peak_spins[0], 138)
    elif slot == "_estimate":
        monkeypatch.setattr(roi_dialog.roi_estimate, "estimate_peak_range", injected)
        QTest.mouseClick(dialog.estimate_button, QtCore.Qt.LeftButton)
    else:
        monkeypatch.setattr(roi_dialog, "_plain_log_ticks", injected)
        shown(dialog)
        dialog.log_check.setChecked(False)
        # on the box itself: the check box is as wide as its grid cell, and a click past its label is no click
        indicator = QtCore.QPoint(6, dialog.log_check.height() // 2)
        QTest.mouseClick(dialog.log_check, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier, indicator)
        assert dialog.log_check.isChecked()
    assert "RuntimeError: injected" in dialog.status.text(), dialog.status.text()
    close(dialog)


def test_the_log_toggle_switches_every_profile():
    """V17 (B6): the log check box, clicked, switches the three profiles to linear, and clicked again back to log.
    On log again the plain formatters are applied again, so no label carries mathtext (V11′)."""
    dialog = shown(make_dialog())
    on_the_box = QtCore.QPoint(6, dialog.log_check.height() // 2)  # the check box is as wide as its grid cell
    profiles = (dialog.y_axis, dialog.tof_axis, dialog.x_axis)
    QTest.mouseClick(dialog.log_check, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier, on_the_box)
    dialog.canvas.draw()
    assert [axes.get_yscale() for axes in profiles] == ["linear"] * 3
    QTest.mouseClick(dialog.log_check, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier, on_the_box)
    assert [axes.get_yscale() for axes in profiles] == ["log"] * 3
    assert not mathtext(dialog), mathtext(dialog)[:5]
    close(dialog)


def test_the_profiles_and_the_ytof_image_open_on_their_data():
    """V18 (B3′; advisories A1 and T1): each axes opens on its data:
    - the TOF profile and the Y-TOF image across the TOF edges (a span made at (0, 1) once pulled x = 0 into the
      Y-TOF axes);
    - the X profile across the detector;
    - the Y profile around the ROIs, inside the detector.
    A peak nudge and a filter change leave the TOF and Y-TOF limits where they are."""
    events = make_events()
    dialog = make_dialog(events)
    dialog.canvas.draw()
    edges = roi_estimate.tof_edges(events)
    span = (edges[0], edges[-1])
    assert dialog.tof_axis.get_xlim() == pytest.approx(span) and dialog.ytof_axis.get_xlim() == pytest.approx(span)
    assert dialog.x_axis.get_xlim() == pytest.approx((0, N_X - 1))
    low, high = dialog.y_axis.get_xlim()
    assert 0 <= low <= 133 and 149 <= high <= N_Y - 1, (low, high)
    type_into(dialog.peak_spins[0], 138)
    type_into(dialog.tof_spins[0], 20000)
    dialog.canvas.draw()
    assert dialog.tof_axis.get_xlim() == pytest.approx(span) and dialog.ytof_axis.get_xlim() == pytest.approx(span)
    close(dialog)


def test_the_dialog_carries_no_geometry_literal_and_reads_no_file():
    """Acceptance 4: no geometry, distance, band literal or file reader in the dialog's module."""
    source = inspect.getsource(roi_dialog)
    for literal in ("304", "256", "15.75", "252.7", "h5py", "get_lam_range"):
        assert literal not in source, literal
