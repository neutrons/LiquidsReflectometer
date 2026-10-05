"""The ROI pop-out (roi-popout-dialog): #197's dialog re-seated on ``lr_reduction.roi_estimate``.

Events are built in memory through the real ``RunEvents``, never a mock, and no test reads the data submodule (A4).
Gestures go through QTest: clicks, typed values, and a press-move-release on the canvas at coordinates taken from
the axes' ``transData`` (measured, not reasoned: L2).
"""

import inspect

import numpy as np
import pytest
from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtTest import QTest

from launcher.apps import roi_dialog
from launcher.apps.roi_dialog import ROISelectionDialog
from lr_reduction import roi_estimate
from lr_reduction.roi_estimate import RunEvents

pytestmark = pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")

N_X, N_Y = 256, 304
PEAK = (136, 146)


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
    """V2 (B3, L1): each image's array is the data layer's, unmodified: get_array() is the data, whatever the display
    does. origin is "lower", and the extent puts each pixel's centre on its index."""
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


def drag(dialog, axis, x0, x1):
    """Press at data x0, move to x1, release, at mid-height of `axis`, on the canvas widget.

    The press and release are QTest's. The moves are the QMouseEvents the platform delivers, sent to the canvas:
    QTest.mouseMove delivers no move to an offscreen widget (measured: the canvas saw button_press_event and
    button_release_event, no motion_notify_event, and the selector set nothing). Sent this way, matplotlib sees two
    motion_notify_events and the SpanSelector its drag."""
    canvas = dialog.canvas
    canvas.draw()
    ratio = canvas.devicePixelRatioF() or 1.0
    y_display = axis.transAxes.transform((0, 0.5))[1]

    def point(x):
        x_display = axis.transData.transform((x, 0))[0]
        return QtCore.QPoint(int(round(x_display / ratio)), int(round(canvas.height() - y_display / ratio)))

    QTest.mousePress(canvas, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier, point(x0))
    for x in ((x0 + x1) / 2, x1):
        move = QtGui.QMouseEvent(QtCore.QEvent.MouseMove, QtCore.QPointF(point(x)), QtCore.Qt.NoButton,
                                 QtCore.Qt.LeftButton, QtCore.Qt.NoModifier)
        QtWidgets.QApplication.sendEvent(canvas, move)
    QTest.mouseRelease(canvas, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier, point(x1))
    QtWidgets.QApplication.processEvents()


@pytest.mark.parametrize("mode, spins", [("peak", "peak_spins"), ("left", "bkg_low"), ("right", "bkg_high")])
def test_dragging_on_the_y_profile_sets_the_chosen_range(mode, spins):
    """V4 (B6): with the radio button chosen, a drag on the Y profile sets that range's spins and overlays."""
    dialog = make_dialog()
    dialog.resize(950, 1000)
    dialog.show()
    QTest.qWaitForWindowExposed(dialog)
    button = next(button for button, button_mode in dialog.mode_buttons if button_mode == mode)
    QTest.mouseClick(button, QtCore.Qt.LeftButton)
    target = {"peak_spins": dialog.peak_spins, "bkg_low": dialog.bkg_spins[:2], "bkg_high": dialog.bkg_spins[2:]}[spins]
    drag(dialog, dialog.y_axis, 120, 128)
    assert [spin.value() for spin in target] == [120, 128]
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


def test_cancel_reports_nothing():
    """V6: edited, then cancelled: the dialog reports no change."""
    dialog = make_dialog()
    type_into(dialog.peak_spins[0], 140)
    type_into(dialog.x_spins[1], 180)
    dialog.reject()
    assert dialog.changes() == {}
    close(dialog)


def test_ok_reports_only_what_changed():
    """V7 (B9): OK reports only the fields whose values differ from those the dialog opened with. An untouched
    [a, b, 0, 0] background, or an unset one, is not reported, so the table keeps it byte for byte."""
    for background in ([133, 149, 0, 0], None):
        dialog = make_dialog(BkgROI=background)
        dialog.accept()
        assert dialog.changes() == {}, background
        close(dialog)
    dialog = make_dialog()
    type_into(dialog.peak_spins[0], 140)
    dialog.accept()
    assert dialog.changes() == {"RB_Ymin": 140}
    close(dialog)
    dialog = make_dialog()
    for spin, value in zip(dialog.bkg_spins, (120, 125, 160, 165)):
        type_into(spin, value)
    type_into(dialog.x_spins[0], 60)
    dialog.accept()
    assert dialog.changes() == {"BkgROI": [120, 125, 160, 165], "data_x_range": [60, 200]}
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


def test_file_text_and_log_ticks_never_reach_the_math_parser():
    """V11 (B12, F7, L6): the run title is file text, so it is not parsed as maths; and the log axes and colorbars use
    a plain formatter, so after a draw no tick or colorbar label carries mathtext ("$")."""
    dialog = make_dialog(title="$\\foo$ run")
    assert dialog.xy_axis.title.get_parse_math() is False
    dialog.canvas.draw()
    labels = [label.get_text() for axis in dialog.figure.axes
              for label in axis.get_xticklabels() + axis.get_yticklabels()]
    assert labels and not [text for text in labels if "$" in text], [text for text in labels if "$" in text][:5]
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


def test_a_nudge_moves_artists_and_rebuilds_nothing():
    """V13 (B11, R8): changing a value moves artists; the axes, images, lines and overlays are the same objects
    afterwards, and an image whose inputs did not change keeps its array."""
    dialog = make_dialog()
    objects = (list(dialog.figure.axes), [axis.images[0] for axis in (dialog.xy_axis, dialog.ytof_axis)],
               [axis.lines[0] for axis in (dialog.y_axis, dialog.tof_axis, dialog.x_axis)],
               {name: dict(artists) for name, artists in dialog.overlays.items()})
    xy_array = dialog.xy_axis.images[0].get_array()
    type_into(dialog.peak_spins[0], 138)
    after = (list(dialog.figure.axes), [axis.images[0] for axis in (dialog.xy_axis, dialog.ytof_axis)],
             [axis.lines[0] for axis in (dialog.y_axis, dialog.tof_axis, dialog.x_axis)],
             {name: dict(artists) for name, artists in dialog.overlays.items()})
    assert all(a is b for a, b in zip(objects[0], after[0]))
    assert all(a is b for a, b in zip(objects[1] + objects[2], after[1] + after[2]))
    assert all(objects[3][name][axis] is after[3][name][axis] for name in objects[3] for axis in objects[3][name])
    assert dialog.xy_axis.images[0].get_array() is xy_array
    close(dialog)


@pytest.mark.parametrize("background, reason", [
    ([0, 10, 150, 160], "sentinel"), ([121, 130], "four bounds"), ("120, 130", "flat list"), (None, "not set")])
def test_an_unusable_background_is_shown_as_not_set_with_its_reason(background, reason):
    """V14 (B5): an entry the reducer cannot use, or none at all, is "not set" in the status line with the reason. No
    bands are drawn and no values are invented."""
    dialog = make_dialog(BkgROI=background)
    assert "not set" in dialog.status.text() and reason in dialog.status.text(), dialog.status.text()
    for name in ("bkg_low", "bkg_high"):
        assert not any(artist.get_visible() for artist in dialog.overlays[name].values()), name
    assert [spin.value() for spin in dialog.bkg_spins] == [roi_dialog.UNSET] * 4
    close(dialog)


def test_a_background_the_reducer_does_not_subtract_is_drawn_and_labelled():
    """V14 (B5): with useBS off for the row, the bands are drawn and the status line says they are not subtracted."""
    dialog = make_dialog(useBS=False)
    assert dialog.overlays["bkg_low"]["y_axis"].get_visible()
    assert "not subtracted" in dialog.status.text()
    close(dialog)


def test_ok_is_unavailable_for_an_inverted_peak_or_a_partial_background():
    """V15 (B10): OK is unavailable while the peak is empty or inverted, or the background is partly set."""
    dialog = make_dialog()
    assert dialog.ok_button.isEnabled()
    type_into(dialog.peak_spins[0], 150)  # above RB_Ymax 146
    assert not dialog.ok_button.isEnabled()
    type_into(dialog.peak_spins[0], 136)
    assert dialog.ok_button.isEnabled()
    for spin, value in zip(dialog.bkg_spins, (120, 125, 160, 165)):
        type_into(spin, value)
    dialog.bkg_spins[3].setValue(roi_dialog.UNSET)  # "not set": a partly set background
    assert not dialog.ok_button.isEnabled()
    close(dialog)
    dialog = make_dialog(RB_Ymin=None)
    assert not dialog.ok_button.isEnabled()
    close(dialog)


def test_the_dialog_carries_no_geometry_literal_and_reads_no_file():
    """Acceptance 4: no geometry, distance, band literal or file reader in the dialog's module."""
    source = inspect.getsource(roi_dialog)
    for literal in ("304", "256", "15.75", "252.7", "h5py", "get_lam_range"):
        assert literal not in source, literal
