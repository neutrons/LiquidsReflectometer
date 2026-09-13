import pytest
from qtpy import QtCore
from qtpy.QtWidgets import QMessageBox


def _check(widget, name, state=QtCore.Qt.Checked):
    for i in range(widget.file_list.count()):
        item = widget.file_list.item(i)
        if item.text() == name:
            item.setCheckState(state)


def _states(widget):
    return {
        widget.file_list.item(i).text(): widget.file_list.item(i).checkState()
        for i in range(widget.file_list.count())
    }


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_refresh_button_exists():
    from launcher.apps.overplot import Overplot

    tab = Overplot()
    assert hasattr(tab, "refresh_btn")
    assert tab.refresh_btn.text() == "Refresh"
    assert tab.refresh_btn.toolTip()


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_refresh_no_folder(monkeypatch):
    """Captures the guard's warning — the prior form asserted only the value
    it had just set, and passed even with the guard deleted."""
    from launcher.apps.overplot import Overplot

    captured = {}
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **_k: captured.setdefault("warn", a))
    tab = Overplot()
    tab.folder_edit.setText("")
    tab.refresh()
    assert "warn" in captured
    assert tab.folder_edit.text() == ""


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_refresh_preserves_check_state(tmp_path):
    from launcher.apps.overplot import Overplot

    a = tmp_path / "a.dat"
    a.write_text("0.01 0.5\n")
    b = tmp_path / "b.dat"
    b.write_text("0.01 0.7\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    _check(tab, "a.dat")
    tab.refresh()
    states = _states(tab)
    assert states["a.dat"] == QtCore.Qt.Checked
    assert states["b.dat"] == QtCore.Qt.Unchecked


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_refresh_picks_up_new_files(tmp_path):
    from launcher.apps.overplot import Overplot

    a = tmp_path / "a.dat"
    a.write_text("0.01 0.5\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    (tmp_path / "c.dat").write_text("0.01 0.9\n")
    tab.refresh()
    items = [tab.file_list.item(i).text() for i in range(tab.file_list.count())]
    assert "c.dat" in items


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_refresh_drops_removed(tmp_path, monkeypatch):
    from launcher.apps.overplot import Overplot

    a = tmp_path / "a.dat"
    a.write_text("0.01 0.5\n")
    b = tmp_path / "b.dat"
    b.write_text("0.01 0.7\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    _check(tab, "b.dat")
    b.unlink()
    captured = {}
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **_k: captured.setdefault("info", a))
    tab.refresh()
    items = [tab.file_list.item(i).text() for i in range(tab.file_list.count())]
    assert "b.dat" not in items
    assert "info" in captured


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_refresh_rereads_content(tmp_path):
    from launcher.apps.overplot import Overplot

    a = tmp_path / "a.dat"
    a.write_text("# Q R dR\n0.01 0.5 0.05\n0.02 0.4 0.04\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    _check(tab, "a.dat")
    tab.plot_selected()
    a.write_text("# Q R dR\n0.01 0.9 0.09\n0.02 0.7 0.07\n")
    tab.refresh()
    ax = tab.figure.axes[0]
    line = ax.lines[0]
    ydata = line.get_ydata()
    assert ydata[0] == pytest.approx(0.9)


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_refresh_unmounted_folder(monkeypatch):
    from launcher.apps.overplot import Overplot

    tab = Overplot()
    tab.folder_edit.setText("/nonexistent/path/12345")
    captured = {}
    monkeypatch.setattr(QMessageBox, "critical", lambda *a, **_k: captured.setdefault("crit", a))
    tab.refresh()
    assert "crit" in captured
    assert tab.folder_edit.text() == "/nonexistent/path/12345"


def _indicator_point(view, item):
    """Centre of the check indicator for *item*, in viewport coordinates.

    The style option must be initialized from the view (viewOptions) and carry
    HasCheckIndicator; a bare QStyleOptionViewItem yields an empty rect whose
    centre is (0, 0), which lands outside the indicator and silently tests
    nothing.
    """
    from qtpy.QtWidgets import QStyle, QStyleOptionViewItem

    option = view.viewOptions()
    option.rect = view.visualItemRect(item)
    option.features |= QStyleOptionViewItem.HasCheckIndicator
    indicator = view.style().subElementRect(QStyle.SE_ItemViewItemCheckIndicator, option, view)
    assert not indicator.isEmpty(), "indicator rect not resolved; the click would test nothing"
    return indicator.center()


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_indicator_click_toggles_check(tmp_path):
    """A real click on the indicator, through the view's own delegate.

    v2 connected itemClicked to a toggle slot; the delegate had already
    toggled by then, so the slot toggled it back and the box became
    unclickable. An emit-based test bypasses the delegate and cannot see that,
    so this drives QTest.mouseClick at real coordinates instead.
    """
    from qtpy.QtCore import Qt
    from qtpy.QtTest import QTest

    from launcher.apps.overplot import Overplot

    (tmp_path / "a.dat").write_text("0.01 0.5\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    tab.show()
    try:
        item = tab.file_list.item(0)
        assert item.checkState() == Qt.Unchecked
        QTest.mouseClick(tab.file_list.viewport(), Qt.LeftButton, Qt.NoModifier, _indicator_point(tab.file_list, item))
        assert item.checkState() == Qt.Checked
        QTest.mouseClick(tab.file_list.viewport(), Qt.LeftButton, Qt.NoModifier, _indicator_point(tab.file_list, item))
        assert item.checkState() == Qt.Unchecked
    finally:
        tab.close()


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_text_click_does_not_toggle_check(tmp_path):
    """Matches deployed exp: only the indicator toggles, not the whole row."""
    from qtpy.QtCore import Qt
    from qtpy.QtTest import QTest

    from launcher.apps.overplot import Overplot

    (tmp_path / "a.dat").write_text("0.01 0.5\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    tab.show()
    try:
        item = tab.file_list.item(0)
        rect = tab.file_list.visualItemRect(item)
        text_point = QtCore.QPoint(rect.right() - 4, rect.center().y())
        QTest.mouseClick(tab.file_list.viewport(), Qt.LeftButton, Qt.NoModifier, text_point)
        assert item.checkState() == Qt.Unchecked
    finally:
        tab.close()


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_selection_is_disabled(tmp_path):
    """The B3 halves pinned separately: selection is retired as an input."""
    from qtpy.QtWidgets import QListWidget

    from launcher.apps.overplot import Overplot

    (tmp_path / "a.dat").write_text("0.01 0.5\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    assert tab.file_list.selectionMode() == QListWidget.NoSelection
    tab.file_list.item(0).setSelected(True)
    assert tab.file_list.selectedItems() == []


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_selected_but_unchecked_never_plots_even_if_selection_returns(tmp_path):
    """Regression probe for the other half: if a future edit restores a
    selection mode, plot_selected must still ignore highlight-only rows."""
    from qtpy.QtWidgets import QListWidget

    from launcher.apps.overplot import Overplot

    (tmp_path / "a.dat").write_text("0.01 0.5\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    tab.file_list.setSelectionMode(QListWidget.MultiSelection)
    tab.file_list.item(0).setSelected(True)
    assert tab.file_list.selectedItems() != []
    tab.plot_selected()
    assert tab.figure.axes == []


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_highlight_without_check_does_not_plot(tmp_path):
    """B3: highlight-only rows were a hidden input to plot_selected."""
    from launcher.apps.overplot import Overplot

    (tmp_path / "a.dat").write_text("0.01 0.5\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    item = tab.file_list.item(0)
    item.setSelected(True)
    tab.plot_selected()
    assert tab.figure.axes == []


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_refresh_preserves_checks_across_sort_shift(tmp_path):
    """Kills the named revert: index-keyed check restore."""
    from launcher.apps.overplot import Overplot

    (tmp_path / "b.dat").write_text("0.01 0.5\n")
    (tmp_path / "c.dat").write_text("0.01 0.7\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    _check(tab, "c.dat")
    # sorts ahead of both, shifting every index by one
    (tmp_path / "0aaa.dat").write_text("0.01 0.9\n")
    tab.refresh()
    states = _states(tab)
    assert states["c.dat"] == QtCore.Qt.Checked
    assert states["b.dat"] == QtCore.Qt.Unchecked
    assert states["0aaa.dat"] == QtCore.Qt.Unchecked


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_refresh_button_click_picks_up_new_file(tmp_path):
    """Kills the named revert: refresh_btn.clicked disconnected."""
    from launcher.apps.overplot import Overplot

    (tmp_path / "a.dat").write_text("0.01 0.5\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    (tmp_path / "c.dat").write_text("0.01 0.9\n")
    tab.refresh_btn.click()
    items = [tab.file_list.item(i).text() for i in range(tab.file_list.count())]
    assert "c.dat" in items


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_refresh_clears_canvas_when_everything_vanished(tmp_path, monkeypatch):
    """F2: an all-failed replot must not leave the previous curves on screen
    while the tooltip claims a fresh refresh."""
    from launcher.apps.overplot import Overplot

    a = tmp_path / "a.dat"
    a.write_text("0.01 0.5\n0.02 0.4\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    _check(tab, "a.dat")
    tab.plot_selected()
    assert tab.figure.axes and tab.figure.axes[0].lines

    # The file disappears but stays checked: the rebuild drops the row, so the
    # replot has nothing to draw.
    a.unlink()
    monkeypatch.setattr(QMessageBox, "information", lambda *_a, **_k: None)
    tab.refresh()
    assert all(not ax.lines for ax in tab.figure.axes)
    assert tab.settings.value("overplot_last_refresh")


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_invalid_typed_folder_does_not_replace_stored_one(tmp_path):
    """F3: same data-loss class as B4 — a typed path that does not exist must
    not overwrite the folder the user had stored."""
    from launcher.apps.overplot import Overplot

    good = tmp_path / "good"
    good.mkdir()
    (good / "a.dat").write_text("0.01 0.5\n")
    s = QtCore.QSettings()
    s.setValue("overplot_folder", str(good))
    s.sync()

    tab = Overplot()
    assert tab.folder_edit.text() == str(good)
    tab.folder_edit.setText("/nonexistent/typed/path")
    tab.plot_mode_combo.setCurrentText("Reflectivity")  # routes through save_settings
    assert QtCore.QSettings().value("overplot_folder") == str(good)


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_refresh_rebuild_preserves_active_filter(tmp_path):
    """Pins the filter-honoring rebuild the reviewer named as revertible."""
    from launcher.apps.overplot import Overplot

    (tmp_path / "keep_a.dat").write_text("0.01 0.5\n")
    (tmp_path / "other_b.dat").write_text("0.01 0.7\n")
    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    tab.filter_edit.setText("keep")
    tab.apply_filter("keep")
    tab.refresh()
    shown = [tab.file_list.item(i).text() for i in range(tab.file_list.count())]
    assert shown == ["keep_a.dat"]
