"""View-level tests for the settings-editor tab (T2).

Thin on purpose: nearly all of the editor's behaviour lives in
`SettingsDocument` and is tested without Qt in
`tests/unit/lr_reduction/test_settings_document.py`. What remains here is the
wiring — that a real user gesture reaches the document, and that the table
edits the row it was told to edit rather than the row that happens to be
selected.

Signals are exercised through `QTest` rather than `.emit()`. Calling `.emit()`
proves a slot is connected to a signal you chose to fire; it does not prove Qt
fires that signal for the gesture the user actually makes, which is the
campaign's signature defect class (S2-v2's double-toggle).
"""

import json

import pytest
from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtTest import QTest

from launcher.app_identity import APP_NAME, ORG_NAME
from launcher.apps.settings_editor import SettingsEditorTab
from lr_reduction import field_spec as fs
from lr_reduction.new_reduction_from_file import save_config_json
from lr_reduction.nr_reduction_config import NRReductionConfig
from lr_reduction.settings_document import SettingsDocument

pytestmark = pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")


def test_tab_constructs():
    tab = SettingsEditorTab()
    assert tab.document is not None
    assert tab.angle_table.columnCount() == len(fs.PER_ANGLE_NAMES)


def test_tab_installs_the_shared_identity_when_none_is_set():
    """S3's contract: every settings layer must resolve to one store.

    `ensure_identity()` is non-clobbering, so the fixture's throwaway
    organization survives; clearing it first is what makes the call observable.
    """
    QtCore.QCoreApplication.setOrganizationName("")
    SettingsEditorTab()
    assert QtCore.QCoreApplication.organizationName() == ORG_NAME
    assert QtCore.QCoreApplication.applicationName() == APP_NAME


def test_typing_in_a_scalar_editor_reaches_the_document():
    """A real keystroke gesture, not a synthesised signal."""
    tab = SettingsEditorTab()
    editor = tab.editors["Sname"]
    editor.clear()
    QTest.keyClicks(editor, "typed_name")
    QTest.keyClick(editor, QtCore.Qt.Key_Return)
    assert tab.document.get("Sname") == "typed_name"


def test_toggling_a_checkbox_reaches_the_document():
    """Space, not a click at the widget centre — measured, not assumed.

    A text-less QCheckBox accepts clicks only inside its ~14 px indicator
    (`SE_CheckBoxClickRect`). An un-shown widget carries the default 640x480
    rect, so `rect().center()` is (320, 240) and lands nowhere near it; even
    shown and laid out at 174x15 the centre is still outside. The first version
    of this test clicked the centre and failed for that reason — a geometry
    assumption dressed up as a user gesture, which is the same stated-vs-
    measured trap the plan warns about one level up.

    `Key_Space` is a genuine activation gesture, goes through the same
    `toggled` path, and does not depend on pixel geometry.
    """
    tab = SettingsEditorTab()
    box = tab.editors["Normalize"]
    assert tab.document.get("Normalize") is False
    QTest.keyClick(box, QtCore.Qt.Key_Space)
    assert tab.document.get("Normalize") is True


def test_a_checkbox_is_clickable_across_its_whole_visible_extent():
    """No inert space: the widget must not be wider than it is clickable."""
    tab = SettingsEditorTab()
    tab.show()
    QtWidgets.QApplication.instance().processEvents()
    box = tab.editors["Normalize"]
    option = QtWidgets.QStyleOptionButton()
    box.initStyleOption(option)
    clickable = box.style().subElementRect(
        QtWidgets.QStyle.SE_CheckBoxClickRect, option, box
    )
    assert clickable.contains(box.rect().center())


def test_choosing_an_enumerated_value_reaches_the_document():
    tab = SettingsEditorTab()
    combo = tab.editors["peak_type"]
    combo.setCurrentText("gauss")
    assert tab.document.get("peak_type") == "gauss"


def test_add_angle_button_grows_both_the_table_and_the_document():
    tab = SettingsEditorTab()
    assert tab.angle_table.rowCount() == 0
    QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
    assert tab.angle_table.rowCount() == 1
    assert tab.document.n_angles == 1
    QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
    assert tab.angle_table.rowCount() == 2
    assert tab.document.n_angles == 2


def test_editing_a_cell_updates_the_row_that_was_edited_not_the_selected_one():
    """The active-row-as-hidden-input trap, at the layer where it is born.

    The table is greenfield in this slug, so this bug would be INTRODUCED here
    rather than inherited. Row 2 is selected while row 0 is edited; if the
    handler consulted `currentRow()` the write would land on row 2.
    """
    tab = SettingsEditorTab()
    for _ in range(3):
        QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
    column = fs.PER_ANGLE_NAMES.index("DBname")

    tab.angle_table.setCurrentCell(2, column)
    assert tab.angle_table.currentRow() == 2

    tab.angle_table.item(0, column).setText("row_zero.dat")

    assert tab.document.get("DBname")[0] == "row_zero.dat"
    assert tab.document.get("DBname")[2] is None


def test_remove_angle_removes_the_selected_row_only():
    tab = SettingsEditorTab()
    for name in ("a.dat", "b.dat", "c.dat"):
        QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
        column = fs.PER_ANGLE_NAMES.index("DBname")
        tab.angle_table.item(tab.angle_table.rowCount() - 1, column).setText(name)

    tab.angle_table.setCurrentCell(1, 0)
    QTest.mouseClick(tab.remove_angle_button, QtCore.Qt.LeftButton)

    assert tab.document.get("DBname") == ["a.dat", "c.dat"]
    assert tab.angle_table.rowCount() == 2


def test_loading_a_settings_file_repopulates_the_view(tmp_path, monkeypatch):
    seed = tmp_path / "seed.json"
    seed.write_text(json.dumps({"Sname": "seeded_run", "qmax": 0.33}))
    monkeypatch.setattr(
        QtWidgets.QFileDialog,
        "getOpenFileName",
        staticmethod(lambda *_a, **_k: (str(seed), "")),
    )

    tab = SettingsEditorTab()
    tab.load_settings()

    assert tab.document.get("Sname") == "seeded_run"
    assert tab.editors["Sname"].text() == "seeded_run"


def test_saving_writes_a_file_the_document_can_read_back(tmp_path, monkeypatch):
    target = tmp_path / "written.json"
    monkeypatch.setattr(
        QtWidgets.QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *_a, **_k: (str(target), "")),
    )

    tab = SettingsEditorTab()
    tab.document.set("Sname", "written_run")
    tab.save_settings()

    assert json.loads(target.read_text())["Sname"] == "written_run"


def test_a_field_prompt_is_available_for_every_editor():
    """The "guide the user" ask: every editor carries its FIELD_SPEC help."""
    tab = SettingsEditorTab()
    for name, editor in tab.editors.items():
        assert fs.get(name).help in editor.toolTip()


def test_validation_report_shows_the_validation_message():
    """Assert the validation LINE, and get there by a gesture.

    The first version asserted only that the field NAME appeared in the report.
    It does — in the changed-vs-seed section, which prints every edited field —
    so the test passed with validate() stubbed to return nothing, and even with
    refresh_report() replaced by a hard-coded string. Asserting the message text
    and reaching it through a real edit closes both.
    """
    tab = SettingsEditorTab()
    editor = tab.editors["Qline_threshold"]
    editor.clear()
    QTest.keyClicks(editor, "4.0")
    QTest.keyClick(editor, QtCore.Qt.Key_Return)

    text = tab.report.toPlainText()
    assert "is above 1.0" in text
    for message in tab.document.validate():
        assert message in text


def test_changed_vs_seed_report_lists_edits():
    tab = SettingsEditorTab()
    editor = tab.editors["Sname"]
    editor.clear()
    QTest.keyClicks(editor, "changed_name")
    QTest.keyClick(editor, QtCore.Qt.Key_Return)
    tab.refresh_report()
    assert "changed_name" in tab.report.toPlainText()


def test_the_launcher_carries_the_tab():
    """The tab is reachable from the shipped entry point, not just importable.

    `new_launcher.py` previously carried a commented-out import of
    `JSONSettingsBuilderTab`, a module that never existed — so "the launcher
    mentions a settings builder" was true and meaningless. This asserts the
    live wiring instead.
    """
    from launcher.new_launcher import ReductionInterface

    window = ReductionInterface()
    titles = [window.tabText(i) for i in range(window.count())]
    assert "Settings editor" in titles
    assert isinstance(window.settings_editor_tab, SettingsEditorTab)


# --------------------------------------------------------------------------
# C1 — no gesture may abort the process
# --------------------------------------------------------------------------


def test_editing_a_short_per_angle_column_does_not_raise():
    """An IndexError here reaches qFatal() and takes every tab down with it."""
    tab = SettingsEditorTab()
    QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
    QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
    tab.document.set("DBname", ["only_one.dat"])
    tab.refresh_angles()

    column = fs.PER_ANGLE_NAMES.index("DBname")
    tab.angle_table.item(1, column).setText("second.dat")

    assert tab.document.get("DBname") == ["only_one.dat", "second.dat"]


def test_a_malformed_settings_file_is_reported_not_fatal(tmp_path, monkeypatch):
    """`{"tof_min": 5}` used to raise TypeError out of the slot."""
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"tof_min": 5}))
    monkeypatch.setattr(
        QtWidgets.QFileDialog, "getOpenFileName",
        staticmethod(lambda *_a, **_k: (str(bad), "")),
    )
    tab = SettingsEditorTab()
    tab.load_settings()
    assert "tof_min" in tab.report.toPlainText()


def test_a_slot_that_raises_reports_into_the_panel(monkeypatch):
    tab = SettingsEditorTab()

    def explode(*_a, **_k):
        raise RuntimeError("synthetic slot failure")

    monkeypatch.setattr(tab.document, "add_angle", explode)
    QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
    assert "synthetic slot failure" in tab.report.toPlainText()


# --------------------------------------------------------------------------
# C2 — a cell edit must reach the document as the declared type
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, typed, expected",
    [
        pytest.param("RB_Ymin", "150", 150, id="int-cell"),
        pytest.param("useBS", "False", False, id="bool-cell"),
        pytest.param("tof_min", "1.5", 1.5, id="float-cell"),
    ],
)
def test_a_table_edit_stores_the_declared_type(name, typed, expected):
    tab = SettingsEditorTab()
    QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
    column = fs.PER_ANGLE_NAMES.index(name)
    tab.angle_table.item(0, column).setText(typed)
    stored = tab.document.get(name)[0]
    assert stored == expected
    assert type(stored) is type(expected)


def test_a_scalar_list_edit_stores_a_list():
    """No clear(): the displayed text has to be what coerce can read back.

    The earlier version called editor.clear() before typing, which erased the
    rendered text before it could be parsed — so it passed while construction
    and refresh rendered lists differently and only one round-tripped.
    """
    tab = SettingsEditorTab()
    editor = tab.editors["data_x_range"]
    QTest.keyClicks(editor, "")
    editor.setText("60, 210")
    QTest.keyClick(editor, QtCore.Qt.Key_Return)
    assert tab.document.get("data_x_range") == [60, 210]


def test_a_scalar_list_survives_focus_out_with_no_typing():
    """The corruption needed no edit at all.

    `data_x_range` is the first editor in tab order and has no validator (its
    type is list[int], so the int/float gate misses it), so `editingFinished`
    fires on a bare focus-out. If the rendering does not round-trip, simply
    tabbing through the form rewrites the document — and nothing is reported,
    because a list of strings is still a list.
    """
    tab = SettingsEditorTab()
    before = tab.document.get("data_x_range")
    editor = tab.editors["data_x_range"]
    assert editor.text() == "50, 200"
    editor.editingFinished.emit()
    assert tab.document.get("data_x_range") == before
    assert tab.document.validate() == []


def test_a_per_angle_nested_cell_survives_a_round_trip():
    """BkgROI is list[list[int]]; str() of it is repr, which re-parses wrong."""
    doc = SettingsDocument()
    doc.add_angle(BkgROI=[120, 130, 140, 150])
    tab = SettingsEditorTab(document=doc)
    column = fs.PER_ANGLE_NAMES.index("BkgROI")
    assert tab.angle_table.item(0, column).text() == "120, 130, 140, 150"
    # Re-commit the displayed text, as any edit to the row does.
    tab.angle_table.item(0, column).setText(tab.angle_table.item(0, column).text())
    assert doc.get("BkgROI")[0] == [120, 130, 140, 150]
    assert doc.validate() == []


# --------------------------------------------------------------------------
# C5 / C6 / should-fixes
# --------------------------------------------------------------------------


def test_theta_source_is_a_choice_not_a_checkbox():
    tab = SettingsEditorTab()
    editor = tab.editors["useCalcTheta"]
    assert isinstance(editor, QtWidgets.QComboBox)
    editor.setCurrentText("sample_angle")
    assert tab.document.get("useCalcTheta") == "sample_angle"


def test_a_combo_displays_a_value_outside_its_choices():
    """Driven through set_document — the Load path, which is the one a user takes.

    Passing the document to the constructor exercised `_build_editor`, where
    `_show_in_combo` already ran. The refresh path called a bare
    `setCurrentText`, which is a silent no-op on a non-editable combo, so the
    widget kept displaying the previous document's value.
    """
    tab = SettingsEditorTab()
    tab.set_document(SettingsDocument.from_dict({"peak_type": "sombrero"}))
    assert tab.editors["peak_type"].currentText() == "sombrero"
    assert any("peak_type" in m for m in tab.document.validate())


def test_a_combo_follows_the_document_when_a_field_is_omitted():
    """The fully silent case: widget says sample angle, reduction uses detector.

    Load a file setting `useCalcTheta`, then one that omits it. The combo used
    to keep showing the old value, and re-selecting it emitted nothing because
    the text never changed — so the saved file disagreed with the screen.
    """
    tab = SettingsEditorTab()
    tab.set_document(SettingsDocument.from_dict({"useCalcTheta": "sample_angle"}))
    assert tab.editors["useCalcTheta"].currentText() == "sample_angle"

    tab.set_document(SettingsDocument.from_dict({"Sname": "week2"}))
    assert tab.document.get("useCalcTheta") is False
    assert tab.editors["useCalcTheta"].currentText() == ""


def test_an_injected_document_renders_its_angles():
    """T3's exact path: __init__ used to render zero rows for a passed document."""
    doc = SettingsDocument()
    doc.add_angle(DBname="a.dat")
    doc.add_angle(DBname="b.dat")
    tab = SettingsEditorTab(document=doc)
    assert tab.angle_table.rowCount() == 2
    column = fs.PER_ANGLE_NAMES.index("DBname")
    assert tab.angle_table.item(1, column).text() == "b.dat"


def test_set_document_replaces_everything_shown():
    tab = SettingsEditorTab()
    replacement = SettingsDocument()
    replacement.add_angle(DBname="new.dat")
    replacement.set("Sname", "replaced")
    tab.set_document(replacement)
    assert tab.angle_table.rowCount() == 1
    assert tab.editors["Sname"].text() == "replaced"


def test_table_sorting_is_disabled():
    """One sortItems() decouples visual row order from document index."""
    tab = SettingsEditorTab()
    assert tab.angle_table.isSortingEnabled() is False


def test_the_checkbox_leg_of_show_follows_the_document():
    """`_show`'s three legs each need a pin; only the line-edit and combo had one.

    A checkbox that does not follow the document displays the previous file's
    value after a Load, which is the same silent disagreement the combo had.
    """
    tab = SettingsEditorTab()
    assert tab.editors["Normalize"].isChecked() is False

    tab.set_document(SettingsDocument.from_dict({"Normalize": True}))
    assert tab.editors["Normalize"].isChecked() is True

    tab.set_document(SettingsDocument.from_dict({"Sname": "no_normalize"}))
    assert tab.editors["Normalize"].isChecked() is False


def test_an_enumerated_editor_offers_the_declared_spellings():
    """The combo cannot author a case variant, which is why entry normalises."""
    tab = SettingsEditorTab()
    combo = tab.editors["DetResFn"]
    offered = [combo.itemText(i) for i in range(combo.count())]
    assert offered == list(fs.DET_RES_CHOICES)


# --------------------------------------------------------------------------
# editor-load-fidelity — what the scientist sees for a reducer-written file
# --------------------------------------------------------------------------


def _reducer_written_settings(directory):
    """A settings file shaped and written by the reduction itself.

    Twin of `_reducer_shaped_config` in
    tests/unit/lr_reduction/test_settings_document.py — copied, not imported:
    `test-launcher` runs from the repository root and `test-reduction` from
    tests/, so a cross-suite import would be the fragile part. `useBS` is ints
    (`nr_reduction_calc.py:103`, `new_reduction_from_template.py:182`); the
    runtime record is one scalar per call (`nr_reduction_calc.py:385-391`).
    """
    config = NRReductionConfig()
    config.experiment_id = "IPTS-00000"
    config.RBnum = [201282, 201283, 201284]
    config.DBname = ["db_a.dat", "db_b.dat", "db_c.dat"]
    config.method_per_run = ["meanTheta"]
    config.RB_Ymin = [140, 141, 142]
    config.RB_Ymax = [150, 151, 152]
    config.BkgROI = [[120, 130], [121, 131], [122, 132]]
    config.useBS = [1] * 3
    config.useBS[2] = 0
    config.useCalcTheta = "detector_angle"
    config.LambdaMinUse = 2.95
    config.LambdaMaxUse = 6.1
    path = directory / "run_settings.json"
    save_config_json(path, config)
    return path


def _load(tab, path, monkeypatch):
    """Load through the button's own slot, with the dialog answering `path`."""
    monkeypatch.setattr(
        QtWidgets.QFileDialog, "getOpenFileName", staticmethod(lambda *_a, **_k: (str(path), ""))
    )
    tab.load_settings()


def _column_text(tab, name):
    column = fs.PER_ANGLE_NAMES.index(name)
    return [tab.angle_table.item(row, column).text() for row in range(tab.angle_table.rowCount())]


def test_a_loaded_background_switch_reads_true_or_false(tmp_path, monkeypatch):
    tab = SettingsEditorTab()
    _load(tab, _reducer_written_settings(tmp_path), monkeypatch)
    assert _column_text(tab, "useBS") == ["true", "true", "false"]


def test_a_reducer_written_file_reports_no_problems(tmp_path, monkeypatch):
    tab = SettingsEditorTab()
    _load(tab, _reducer_written_settings(tmp_path), monkeypatch)
    assert tab.report.toPlainText() == "No problems found."


@pytest.mark.parametrize("name", ["LambdaMinUse", "LambdaMaxUse"])
def test_typing_into_the_runtime_record_changes_nothing(tmp_path, monkeypatch, name):
    """Keystrokes, not setText: a read-only QLineEdit refuses the gesture, and
    Return still emits editingFinished, so a live connection would show here."""
    tab = SettingsEditorTab()
    _load(tab, _reducer_written_settings(tmp_path), monkeypatch)
    editor = tab.editors[name]
    shown, held = editor.text(), tab.document.get(name)
    assert editor.isReadOnly()
    QTest.keyClicks(editor, "9.9")
    QTest.keyClick(editor, QtCore.Qt.Key_Return)
    assert editor.text() == shown
    assert repr(tab.document.get(name)) == repr(held)


@pytest.mark.parametrize("name", ["LambdaMinUse", "LambdaMaxUse"])
def test_a_programmatic_edit_signal_on_the_runtime_record_changes_nothing(tmp_path, monkeypatch, name):
    """The pathological leg: a signal no user gesture produced must not reach the
    document either, so the record editors are not connected at all."""
    tab = SettingsEditorTab()
    _load(tab, _reducer_written_settings(tmp_path), monkeypatch)
    held = tab.document.get(name)
    tab.editors[name].editingFinished.emit()
    assert repr(tab.document.get(name)) == repr(held)


def test_the_runtime_record_shows_what_each_loaded_file_recorded(tmp_path, monkeypatch):
    tab = SettingsEditorTab()
    _load(tab, _reducer_written_settings(tmp_path), monkeypatch)
    assert tab.editors["LambdaMinUse"].text() == "2.95"

    listed = tmp_path / "listed.json"
    listed.write_text(json.dumps({"LambdaMinUse": [2.95, 3.1]}))
    _load(tab, listed, monkeypatch)
    assert tab.editors["LambdaMinUse"].text() == "2.95, 3.1"

    absent = tmp_path / "absent.json"
    absent.write_text(json.dumps({"Sname": "no_record"}))
    _load(tab, absent, monkeypatch)
    assert tab.editors["LambdaMinUse"].text() == ""


def test_toggling_a_loaded_switch_saves_ones_and_zeros(tmp_path, monkeypatch):
    tab = SettingsEditorTab()
    _load(tab, _reducer_written_settings(tmp_path), monkeypatch)
    tab.angle_table.item(1, fs.PER_ANGLE_NAMES.index("useBS")).setText("false")
    assert "useBS: [True, True, False] -> [True, False, False]" in tab.report.toPlainText()

    target = tmp_path / "saved.json"
    monkeypatch.setattr(
        QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *_a, **_k: (str(target), ""))
    )
    tab.save_settings()
    saved = json.loads(target.read_text())["useBS"]
    assert [type(v) for v in saved] == [int, int, int]
    assert saved == [1, 0, 0]


def test_an_injected_integer_switch_renders_true_or_false():
    config = NRReductionConfig()
    config.useBS = [1, 1, 0]
    tab = SettingsEditorTab(document=SettingsDocument(config))
    assert _column_text(tab, "useBS") == ["true", "true", "false"]


def test_a_numeric_column_holding_ones_and_zeros_keeps_its_numbers(tmp_path, monkeypatch):
    """Only the boolean column reads true/false.

    The reducer fills an empty `ThetaShift` with `[0] * n` and an empty
    `ScaleFactor` with `[1] * n` (`nr_reduction_calc.py:101,105`) on the config
    it then saves, so a file saved after a run can carry integer 0s and 1s in
    numeric columns too (the 1s survive unless prior-combination replaces them,
    `new_reduction_from_file.py:96`). Added when mutation M22 (the cell-text
    helper applied to every column) survived the battery: no test held a numeric
    1/0, and a scale factor would have read "true".
    """
    path = tmp_path / "after_a_run.json"
    path.write_text(json.dumps({"useBS": [1, 1, 0], "ScaleFactor": [1, 1, 1], "ThetaShift": [0, 0, 0]}))
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    assert _column_text(tab, "useBS") == ["true", "true", "false"]
    assert _column_text(tab, "ScaleFactor") == ["1", "1", "1"]
    assert _column_text(tab, "ThetaShift") == ["0", "0", "0"]


def test_a_hand_written_integer_useGravity_is_reported_and_saved_as_written(tmp_path, monkeypatch):
    """v2 (review 8b62952): the reducer reads useGravity with `is True`, so a file holding 1 reduces
    with gravity correction OFF. Load -> Save must neither hide nor change that."""
    path = tmp_path / "hand_written.json"
    path.write_text(json.dumps({"useGravity": 1}))
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    assert "(useGravity)" in tab.report.toPlainText()

    target = tmp_path / "saved.json"
    monkeypatch.setattr(
        QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *_a, **_k: (str(target), ""))
    )
    tab.save_settings()
    assert '"useGravity": 1' in target.read_text()


# --------------------------------------------------------------------------
# editor-angle-count — the reduction's angle count; surplus rows; Add keeps one index
# --------------------------------------------------------------------------


def _surplus_settings(directory):
    """Three angles by RBnum, useBS four long — IPTS-36119's reduce_settings.json shape."""
    path = directory / "reduce_settings.json"
    path.write_text(json.dumps({
        "RBnum": [201282, 201283, 201284],
        "DBname": ["db_a.dat", "db_b.dat", "db_c.dat"],
        "RB_Ymin": [140, 141, 142],
        "RB_Ymax": [150, 151, 152],
        "BkgROI": [[120, 130], [121, 131], [122, 132]],
        "useBS": [1, 1, 1, 1],
    }))
    return path


def _row_label(tab, row):
    item = tab.angle_table.verticalHeaderItem(row)
    return item.text() if item is not None else ""


def test_a_surplus_file_reads_no_problems_and_shows_its_note(tmp_path, monkeypatch):
    tab = SettingsEditorTab()
    _load(tab, _surplus_settings(tmp_path), monkeypatch)
    text = tab.report.toPlainText()
    assert text.startswith("No problems found.")
    notes = text.split("Notes:", 1)[1]
    assert "(useBS)" in notes and "1 extra" in notes


def test_the_surplus_row_is_marked_and_the_others_are_not(tmp_path, monkeypatch):
    """The queryable property: the row's vertical header says "surplus"."""
    tab = SettingsEditorTab()
    _load(tab, _surplus_settings(tmp_path), monkeypatch)
    assert tab.angle_table.rowCount() == 4
    assert [("surplus" in _row_label(tab, row)) for row in range(4)] == [False, False, False, True]


def test_add_angle_inserts_the_reductions_next_angle_before_the_surplus_rows(tmp_path, monkeypatch):
    """v2 (G6 revised, review 1568397 Q-1): the new angle is row m, directly after the last
    real angle; the surplus row moves down and stays surplus. v1 appended after it, so a
    saved RBnum read [..., null, <new run>]."""
    tab = SettingsEditorTab()
    _load(tab, _surplus_settings(tmp_path), monkeypatch)
    QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
    for name, text in (("DBname", "new.dat"), ("RBnum", "999999")):
        tab.angle_table.item(3, fs.PER_ANGLE_NAMES.index(name)).setText(text)
    doc = tab.document
    assert doc.get("DBname")[3] == "new.dat"
    assert doc.reduction_angles == 4
    assert repr(doc.get("useBS")) == "[True, True, True, None, True]"
    assert [("surplus" in _row_label(tab, row)) for row in range(5)] == [False, False, False, False, True]
    target = tmp_path / "saved.json"
    monkeypatch.setattr(
        QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *_a, **_k: (str(target), ""))
    )
    tab.save_settings()
    assert json.loads(target.read_text())["RBnum"] == [201282, 201283, 201284, 999999]


def test_removing_the_surplus_row_clears_its_note(tmp_path, monkeypatch):
    tab = SettingsEditorTab()
    _load(tab, _surplus_settings(tmp_path), monkeypatch)
    tab.angle_table.setCurrentCell(3, 0)
    QTest.mouseClick(tab.remove_angle_button, QtCore.Qt.LeftButton)
    assert "(useBS)" not in tab.report.toPlainText()
    assert tab.document.get("useBS") == [True, True, True]


def test_a_new_document_writes_the_reducers_defaults_as_empty_lists(tmp_path, monkeypatch):
    """F5: the editor's own file must reduce; an unset default list is written [] (G7)."""
    tab = SettingsEditorTab()
    for _ in range(2):
        QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
    for row in range(2):
        for name, text in (("DBname", f"db_{row}.dat"), ("RB_Ymin", "140"), ("RB_Ymax", "150"), ("BkgROI", "120, 130")):
            tab.angle_table.item(row, fs.PER_ANGLE_NAMES.index(name)).setText(text)
    target = tmp_path / "authored.json"
    monkeypatch.setattr(
        QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *_a, **_k: (str(target), ""))
    )
    tab.save_settings()
    saved = json.loads(target.read_text())
    assert saved["ThetaShift"] == []
    assert saved["method_per_run"] == []


@pytest.mark.parametrize("row, still_surplus", [(0, True), (3, False)], ids=["edit-an-angle", "edit-the-surplus-row"])
def test_the_surplus_marks_follow_an_edit(tmp_path, monkeypatch, row, still_surplus):
    """v2 (G8): marks are re-derived after an edit, as after Load, Add and Remove. Editing an
    angle leaves the surplus row surplus and the panel clean; a DBname typed into the surplus
    row makes it a real angle."""
    tab = SettingsEditorTab()
    _load(tab, _surplus_settings(tmp_path), monkeypatch)
    tab.angle_table.item(row, fs.PER_ANGLE_NAMES.index("DBname")).setText("edited.dat")
    assert ("surplus" in _row_label(tab, 3)) is still_surplus
    if still_surplus:
        text = tab.report.toPlainText()
        assert text.startswith("No problems found.")
        assert "(useBS)" in text.split("Notes:", 1)[1]


# --------------------------------------------------------------------------
# editor-angle-count v3 — G9 at the gesture: an edit of a compact list shows what it wrote; a λ
# typed into a surplus row of a derived list is refused, visibly
# --------------------------------------------------------------------------


def test_an_edit_of_an_empty_default_list_shows_the_values_it_wrote(tmp_path, monkeypatch):
    """G9 writes entries the user did not type (the reducer's own value at the other angles), so the
    column is re-drawn from the document: what is shown is what is held."""
    tab = SettingsEditorTab()
    _load(tab, _surplus_settings(tmp_path), monkeypatch)
    tab.angle_table.item(1, fs.PER_ANGLE_NAMES.index("ThetaShift")).setText("0.01")
    assert repr(tab.document.get("ThetaShift")) == "[0, 0.01, 0]"
    assert _column_text(tab, "ThetaShift") == ["0", "0.01", "0", ""]
    assert tab.report.toPlainText().startswith("No problems found.")


def test_a_lambda_typed_into_a_surplus_row_is_refused_and_the_cell_shows_it(tmp_path, monkeypatch):
    """S3 at the gesture: the document keeps λ derived, the cell goes back to unset, and the panel
    says why."""
    tab = SettingsEditorTab()
    _load(tab, _surplus_settings(tmp_path), monkeypatch)
    tab.angle_table.item(3, fs.PER_ANGLE_NAMES.index("LambdaMin")).setText("3.0")
    assert tab.document.get("LambdaMin") is None
    assert _column_text(tab, "LambdaMin") == ["", "", "", ""]
    problems = tab.report.toPlainText().split("Notes:", 1)[0]
    assert "(LambdaMin)" in problems and "surplus row 4" in problems


def test_emptying_the_lambda_that_alone_reached_the_last_row_leaves_the_table_drawable(tmp_path, monkeypatch):
    """W21 (added when it survived): clearing λ's last value turns it back into None, and if it was the
    only list reaching the last row, the table now shows a row past the document's. Re-drawing the
    column reads that row as empty rather than raising inside the edit."""
    path = tmp_path / "lambda_longest.json"
    path.write_text(json.dumps({
        "RBnum": [201282, 201283, 201284], "DBname": ["db_a.dat", "db_b.dat", "db_c.dat"],
        "RB_Ymin": [140, 141, 142], "RB_Ymax": [150, 151, 152], "BkgROI": [[120, 130], [121, 131], [122, 132]],
        "LambdaMin": [2.5, None, None, None],
    }))
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    tab.angle_table.item(0, fs.PER_ANGLE_NAMES.index("LambdaMin")).setText("")
    assert tab.document.get("LambdaMin") is None
    assert "Could not complete" not in tab.report.toPlainText()
    assert _column_text(tab, "LambdaMin") == ["", "", "", ""]


# --------------------------------------------------------------------------
# editor-combos — the wheel never changes a drop-down (item 2); the Angles table's enumerated columns
# are drop-downs, and a compact column shows the value the reduction uses (item 8; C1-C7)
# --------------------------------------------------------------------------

_THREE_ANGLES = {
    "RBnum": [201282, 201283, 201284],
    "DBname": ["db_a.dat", "db_b.dat", "db_c.dat"],
    "RB_Ymin": [140, 141, 142],
    "RB_Ymax": [150, 151, 152],
    "BkgROI": [[120, 130], [121, 131], [122, 132]],
}
_SCALAR_CHOICES = [f.name for f in fs.FIELD_SPEC if f.allowed and not f.per_angle]


def _wheel(widget, delta=-120):
    """One wheel notch over `widget`, sent to the widget object (not to a coordinate).

    Qt 5 propagates only spontaneous wheel events to the parent (QApplication::notify), so a test that
    sends one can show that the drop-down leaves it unaccepted, not that the list underneath scrolls.
    """
    centre = QtCore.QPointF(widget.rect().center())
    event = QtGui.QWheelEvent(
        centre, QtCore.QPointF(widget.mapToGlobal(widget.rect().center())), QtCore.QPoint(0, 0),
        QtCore.QPoint(0, delta), QtCore.Qt.NoButton, QtCore.Qt.NoModifier, QtCore.Qt.NoScrollPhase, False,
    )
    QtWidgets.QApplication.sendEvent(widget, event)
    return event


def _away(combo):
    """The notch that would move `combo` off its current item: down (the next item) unless it is the last.

    One notch, not a down-and-up pair: two opposite notches cancel, and the test would pass even when
    the wheel moves the selection."""
    return -120 if combo.currentIndex() < combo.count() - 1 else 120


def _settings_file(directory, values):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "settings.json"
    path.write_text(json.dumps(values))
    return path


def _direct_beam_settings(directory, names, extra=None):
    """The three-angle settings with _DBpath_override pointing at a folder holding `names`."""
    folder = directory / "transmission"
    folder.mkdir(parents=True)
    for name in names:
        (folder / name).write_text("")
    return _settings_file(directory, {**_THREE_ANGLES, "_DBpath_override": str(folder), **(extra or {})}), folder


def _open_cell_editor(tab, row, name):
    """Open the cell's editor the way the table's edit triggers do, and return it."""
    index = tab.angle_table.model().index(row, fs.PER_ANGLE_NAMES.index(name))
    tab.angle_table.edit(index)
    QtWidgets.QApplication.processEvents()
    return tab.angle_table.indexWidget(index)


def _leave(editor):
    """Press Return in the cell's editor and let the commit happen: Qt 5 queues a delegate's commit on
    Return (QueuedConnection), so asserting before processing events would see the document unchanged
    whether or not the cell writes."""
    QTest.keyClick(editor, QtCore.Qt.Key_Return)
    QtWidgets.QApplication.processEvents()


def _choose(editor, text):
    """Choose `text` the way a keyboard user does: arrow keys on the closed drop-down, then Return."""
    target = editor.findText(text)
    assert target >= 0, text
    key = QtCore.Qt.Key_Down if target > editor.currentIndex() else QtCore.Qt.Key_Up
    for _ in range(editor.count()):
        if editor.currentIndex() == target:
            break
        QTest.keyClick(editor, key)
    assert editor.currentText() == text
    _leave(editor)


def _items(editor):
    return [editor.itemText(i) for i in range(editor.count())]


def _shown(tab, row, name):
    """(text, italic) as the column's delegate renders the cell: the queryable form of what is displayed,
    and of C7's mark (an implied value is shown in italics)."""
    column = fs.PER_ANGLE_NAMES.index(name)
    delegate = tab.angle_table.itemDelegateForColumn(column)
    if delegate is None:  # a plain text cell: what its item holds (PyQt cannot call the C++ default's hook)
        item = tab.angle_table.item(row, column)
        return (item.text(), item.font().italic()) if item is not None else ("", False)
    option = QtWidgets.QStyleOptionViewItem()
    delegate.initStyleOption(option, tab.angle_table.model().index(row, column))
    return option.text, option.font.italic()


@pytest.mark.parametrize("delta", [-120, 120], ids=["down", "up"])
@pytest.mark.parametrize("name", _SCALAR_CHOICES)
def test_the_wheel_never_changes_a_scalar_drop_down(name, delta):
    """V1 (item 2). Measured at the base: one notch over DetResFn turned rectangular into gaussian, and
    the document with it. The drop-down leaves the event to its parent."""
    tab = SettingsEditorTab()
    combo = tab.editors[name]
    shown = combo.currentText()
    event = _wheel(combo, delta)
    assert combo.currentText() == shown
    assert tab.document.changed_vs_seed() == {}
    assert not event.isAccepted()


def test_the_wheel_never_changes_a_drop_down_that_has_focus():
    """V2. After a click the drop-down keeps focus, so "ignore unless focused" leaves item 2's failure
    reachable: pick a value, keep scrolling with the pointer over it (plan A1)."""
    tab = SettingsEditorTab()
    tab.resize(900, 500)
    tab.show()
    QTest.qWaitForWindowExposed(tab)
    QtWidgets.QApplication.setActiveWindow(tab)
    combo = tab.editors["DetResFn"]
    combo.setFocus()
    QtWidgets.QApplication.processEvents()
    assert combo.hasFocus()
    shown = combo.currentText()
    _wheel(combo, _away(combo))
    assert combo.currentText() == shown
    assert tab.document.changed_vs_seed() == {}
    tab.close()


@pytest.mark.parametrize("row", [0, 3], ids=["angle", "surplus-row"])
@pytest.mark.parametrize("name", ["method_per_run", "useBS", "DBname"])
def test_the_wheel_never_changes_a_table_drop_down(tmp_path, monkeypatch, name, row):
    """V3, V15: the open editor is a closed drop-down until its list is shown."""
    path, _ = _direct_beam_settings(tmp_path, ["db_a.dat", "db_z.dat"], {"useBS": [1, 1, 1, 1]})
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    editor = _open_cell_editor(tab, row, name)
    assert isinstance(editor, QtWidgets.QComboBox)
    assert editor.count() > 1  # something the wheel could move to
    shown = editor.currentText()
    _wheel(editor, _away(editor))
    assert editor.currentText() == shown
    _leave(editor)
    assert tab.document.changed_vs_seed() == {}


def test_the_q_method_cell_offers_exactly_the_declared_methods_and_unset():
    """V4: from the declaration (reduction_domains.METHOD_CHOICES), not a literal list; "" is unset."""
    tab = SettingsEditorTab(SettingsDocument.from_dict(_THREE_ANGLES))
    assert _items(_open_cell_editor(tab, 0, "method_per_run")) == ["", *fs.METHOD_CHOICES]


def test_choosing_a_q_method_writes_the_declared_spelling_into_the_edited_row_not_the_selected_one():
    """V5, C6: row 2 is selected while row 0's drop-down is edited."""
    tab = SettingsEditorTab(SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["meanTheta"] * 3}))
    tab.angle_table.setCurrentCell(2, 0)
    _choose(_open_cell_editor(tab, 0, "method_per_run"), "constantQ")
    assert tab.document.get("method_per_run") == ["constantQ", "meanTheta", "meanTheta"]


def test_the_background_cell_offers_true_false_and_unset_and_stores_a_bool(tmp_path, monkeypatch):
    """V6, C3: the text "false" would be truthy to the reducer (subtracting background the scientist
    switched off), and 0 is not the document's spelling."""
    tab = SettingsEditorTab()
    _load(tab, _reducer_written_settings(tmp_path), monkeypatch)  # useBS [1, 1, 0]
    editor = _open_cell_editor(tab, 1, "useBS")
    assert _items(editor) == ["", "true", "false"]
    _choose(editor, "false")
    assert tab.document.get("useBS")[1] is False


def test_a_drop_down_writes_the_row_it_is_in_after_rows_move():
    """V7, C6: after Add x3 and Remove of row 0, the drop-down now in row 0 writes document index 0."""
    tab = SettingsEditorTab()
    for _ in range(3):
        QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
    for row, name in enumerate(("a.dat", "b.dat", "c.dat")):
        tab.angle_table.item(row, fs.PER_ANGLE_NAMES.index("DBname")).setText(name)
    tab.angle_table.setCurrentCell(0, 0)
    QTest.mouseClick(tab.remove_angle_button, QtCore.Qt.LeftButton)
    tab.angle_table.setCurrentCell(1, 0)
    _choose(_open_cell_editor(tab, 0, "useBS"), "false")
    assert tab.document.get("DBname") == ["b.dat", "c.dat"]
    assert repr(tab.document.get("useBS")) == "[False, True]"


def test_the_direct_beam_cell_offers_the_folders_files_and_completes_typed_text(tmp_path, monkeypatch):
    """V8, C4: *.txt and *.dat in the resolved folder, sorted; notes.md is not offered."""
    path, _ = _direct_beam_settings(tmp_path, ["db_b.dat", "db_a.txt", "notes.md"])
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    editor = _open_cell_editor(tab, 0, "DBname")
    assert editor.isEditable()
    assert _items(editor) == ["db_a.txt", "db_b.dat"]
    completer = editor.completer()
    completer.setCompletionPrefix("db_b")
    assert completer.currentCompletion() == "db_b.dat"


def test_a_direct_beam_name_not_in_the_folder_can_be_typed_and_is_stored_verbatim(tmp_path, monkeypatch):
    """V9, C4: the list is a help, not a constraint."""
    path, _ = _direct_beam_settings(tmp_path, ["db_a.dat"])
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    editor = _open_cell_editor(tab, 1, "DBname")
    editor.lineEdit().selectAll()
    QTest.keyClicks(editor.lineEdit(), "elsewhere 1.dat")
    _leave(editor)
    assert tab.document.get("DBname") == ["db_a.dat", "elsewhere 1.dat", "db_c.dat"]


def test_an_empty_or_missing_direct_beam_folder_offers_nothing_and_typing_still_works(tmp_path, monkeypatch):
    """C4: the folder does not exist (an IPTS whose shared/transmission is absent): an empty list, no
    exception, and a typed name is stored."""
    path = _settings_file(tmp_path, {**_THREE_ANGLES, "_DBpath_override": str(tmp_path / "absent")})
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    editor = _open_cell_editor(tab, 0, "DBname")
    assert _items(editor) == []
    editor.lineEdit().selectAll()
    QTest.keyClicks(editor.lineEdit(), "typed.dat")
    _leave(editor)
    assert tab.document.get("DBname")[0] == "typed.dat"


def _offered(tab):
    editor = _open_cell_editor(tab, 0, "DBname")
    names = _items(editor)
    QTest.keyClick(editor, QtCore.Qt.Key_Escape)
    QtWidgets.QApplication.processEvents()
    return names


@pytest.mark.parametrize("how", ["edit-the-folder", "load-another-file"])
def test_the_direct_beam_list_follows_the_folder(tmp_path, monkeypatch, how):
    """V10, C5: the next time a cell offers names, they come from the folder the document resolves now."""
    first, _ = _direct_beam_settings(tmp_path / "one", ["first.dat"])
    second, second_folder = _direct_beam_settings(tmp_path / "two", ["second.dat"])
    tab = SettingsEditorTab()
    _load(tab, first, monkeypatch)
    assert _offered(tab) == ["first.dat"]
    if how == "edit-the-folder":
        editor = tab.editors["_DBpath_override"]
        editor.clear()
        QTest.keyClicks(editor, str(second_folder))
        QTest.keyClick(editor, QtCore.Qt.Key_Return)
    else:
        _load(tab, second, monkeypatch)
    assert _offered(tab) == ["second.dat"]


def test_a_loaded_case_variant_shows_as_its_choice_and_is_kept_until_one_is_chosen(tmp_path, monkeypatch):
    """V11, C2: the reducer lower-cases (nr_reduction_calc.py:82), so 'meantheta' is meanTheta. Displaying
    it as that is not an edit: leaving the cell without choosing keeps the file's spelling."""
    path = _settings_file(tmp_path, {**_THREE_ANGLES, "method_per_run": ["meantheta"] * 3})
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    assert _shown(tab, 0, "method_per_run") == ("meanTheta", False)
    editor = _open_cell_editor(tab, 0, "method_per_run")
    assert editor.currentText() == "meanTheta"
    _leave(editor)
    assert tab.document.get("method_per_run") == ["meantheta"] * 3


@pytest.mark.parametrize("name, held, text", [("method_per_run", "sombrero", "sombrero"), ("useBS", 2, "2")])
def test_an_out_of_domain_cell_value_is_shown_as_itself_and_reported(tmp_path, monkeypatch, name, held, text):
    """V12, C2: never silently replaced by the first item, which a save would then write back."""
    path = _settings_file(tmp_path, {**_THREE_ANGLES, name: [held] * 3})
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    assert _shown(tab, 0, name)[0] == text
    editor = _open_cell_editor(tab, 0, name)
    assert editor.currentText() == text
    _leave(editor)
    assert tab.document.get(name) == [held] * 3
    assert any(f"({name})" in line for line in tab.document.validate())


def test_a_compact_column_shows_the_value_the_reduction_uses_marked_as_implied(tmp_path, monkeypatch):
    """V13, C7: with a one-entry method_per_run and an empty useBS, the reduction uses constantQ and
    background subtraction on at every angle; an empty cell said otherwise (F8). Showing it writes
    nothing; choosing in such a cell writes the list out as editor-angle-count's G9 defines."""
    path = _settings_file(tmp_path, {**_THREE_ANGLES, "method_per_run": ["constantQ"], "useBS": []})
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    assert [_shown(tab, row, "method_per_run") for row in range(3)] == [
        ("constantQ", False), ("constantQ", True), ("constantQ", True)]
    assert [_shown(tab, row, "useBS") for row in range(3)] == [("true", True)] * 3
    assert tab.document.changed_vs_seed() == {}
    _choose(_open_cell_editor(tab, 1, "method_per_run"), "meanTheta")
    assert tab.document.get("method_per_run") == ["constantQ", "meanTheta", "constantQ"]
    assert [_shown(tab, row, "method_per_run") for row in range(3)] == [
        ("constantQ", False), ("meanTheta", False), ("constantQ", False)]


def test_a_surplus_row_shows_what_it_holds_never_an_implied_value(tmp_path, monkeypatch):
    """V13's surplus leg, C7: the reduction never reads a surplus row, so nothing is implied there."""
    tab = SettingsEditorTab()
    _load(tab, _surplus_settings(tmp_path), monkeypatch)  # useBS x4, method_per_run []
    assert _shown(tab, 3, "useBS") == ("true", False)
    assert _shown(tab, 3, "method_per_run") == ("", False)
    assert _shown(tab, 2, "method_per_run") == ("meanTheta", True)


def test_displaying_implied_values_writes_nothing_and_a_save_keeps_the_lists_compact(tmp_path, monkeypatch):
    """V14 (Save): the implied values are shown, not held; the file still says "use your default"."""
    tab = SettingsEditorTab()
    _load(tab, _settings_file(tmp_path, {**_THREE_ANGLES, "method_per_run": [], "useBS": []}), monkeypatch)
    assert _shown(tab, 0, "useBS") == ("true", True)
    target = tmp_path / "saved.json"
    monkeypatch.setattr(
        QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *_a, **_k: (str(target), ""))
    )
    tab.save_settings()
    saved = json.loads(target.read_text())
    assert saved["method_per_run"] == [] and saved["useBS"] == []


def test_add_angle_keeps_compact_columns_compact_and_its_row_shows_the_implied_values():
    """V14 (Add): the new angle at index m shows the values the reduction will use there; choosing in it
    writes the list out (G9)."""
    tab = SettingsEditorTab(SettingsDocument.from_dict(_THREE_ANGLES))
    QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
    assert tab.document.get("method_per_run") == [] and tab.document.get("useBS") == []
    assert _shown(tab, 3, "method_per_run") == ("meanTheta", True)
    assert _shown(tab, 3, "useBS") == ("true", True)
    _choose(_open_cell_editor(tab, 3, "method_per_run"), "constantQ")
    assert tab.document.get("method_per_run") == ["meanTheta"] * 3 + ["constantQ"]
    assert not any("(method_per_run)" in line for line in tab.document.validate())


def test_choosing_unset_in_an_implied_cell_writes_nothing():
    """V14 (choose "unset", compact): the list was compact, so the cell stays implied."""
    tab = SettingsEditorTab(SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["constantQ"]}))
    _choose(_open_cell_editor(tab, 1, "method_per_run"), "")
    assert tab.document.get("method_per_run") == ["constantQ"]
    assert tab.document.changed_vs_seed() == {}
    assert _shown(tab, 1, "method_per_run") == ("constantQ", True)


def test_choosing_unset_in_a_held_cell_unsets_that_entry():
    """V14 (choose "unset", full list): the entry becomes None, and the some-unset rule names the angle."""
    tab = SettingsEditorTab(SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["constantQ"] * 3}))
    _choose(_open_cell_editor(tab, 1, "method_per_run"), "")
    assert tab.document.get("method_per_run") == ["constantQ", None, "constantQ"]
    assert any("(method_per_run)" in line and "angles [1]" in line for line in tab.document.validate())


def test_removing_rows_redraws_the_drop_downs_from_the_document(tmp_path, monkeypatch):
    """V14 (Remove): a surplus row, then a real one; what is shown follows the document."""
    tab = SettingsEditorTab()
    _load(tab, _surplus_settings(tmp_path), monkeypatch)  # useBS x4, method_per_run []
    tab.angle_table.setCurrentCell(3, 0)
    QTest.mouseClick(tab.remove_angle_button, QtCore.Qt.LeftButton)
    assert "(useBS)" not in tab.report.toPlainText()
    tab.angle_table.setCurrentCell(0, 0)
    QTest.mouseClick(tab.remove_angle_button, QtCore.Qt.LeftButton)
    assert tab.document.get("DBname") == ["db_b.dat", "db_c.dat"]
    assert [_shown(tab, row, "method_per_run") for row in range(2)] == [("meanTheta", True)] * 2
    assert [_shown(tab, row, "useBS") for row in range(2)] == [("true", False)] * 2


def test_a_choice_in_a_surplus_row_is_a_surplus_value(tmp_path, monkeypatch):
    """V14 (choose, surplus row): G9 fills the angles with meanTheta; the note names the surplus entry
    and the row stays marked."""
    tab = SettingsEditorTab()
    _load(tab, _surplus_settings(tmp_path), monkeypatch)  # m = 3, n = 4, method_per_run []
    _choose(_open_cell_editor(tab, 3, "method_per_run"), "constantQ")
    assert tab.document.get("method_per_run") == ["meanTheta"] * 3 + ["constantQ"]
    notes = tab.report.toPlainText().split("Notes:", 1)[1]
    assert "(method_per_run)" in notes and "1 extra" in notes
    assert "surplus" in _row_label(tab, 3)


def test_a_capped_direct_beam_list_says_so(tmp_path, monkeypatch):
    """Added at GREEN, for a branch the RED set did not construct: a folder holding more names than the
    cap offers the first MAX_CANDIDATES, and the drop-down's tooltip says how many there were."""
    from lr_reduction.settings_document import MAX_CANDIDATES

    names = [f"db_{k:05d}.dat" for k in range(MAX_CANDIDATES + 1)]
    path, _ = _direct_beam_settings(tmp_path, names)
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    editor = _open_cell_editor(tab, 0, "DBname")
    assert editor.count() == MAX_CANDIDATES
    assert f"first {MAX_CANDIDATES} of {MAX_CANDIDATES + 1}" in editor.toolTip()
