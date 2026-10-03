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
from qtpy import QtCore, QtWidgets
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
