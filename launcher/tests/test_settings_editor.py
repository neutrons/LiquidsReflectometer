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
from qtpy import QtCore, QtGui, QtWidgets, sip
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


def test_the_theta_control_is_a_choice_not_a_checkbox():
    """C5, rewritten for editor-defaults-and-theta: its entries are labels now, and the document gets the stored
    value the label stands for."""
    tab = SettingsEditorTab()
    editor = tab.editors["useCalcTheta"]
    assert isinstance(editor, QtWidgets.QComboBox)
    editor.setCurrentText("trust sample angle")
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
    assert tab.editors["useCalcTheta"].currentText() == "trust sample angle"

    tab.set_document(SettingsDocument.from_dict({"Sname": "week2"}))
    assert tab.document.get("useCalcTheta") is False
    assert tab.editors["useCalcTheta"].currentText() == "False"  # V5 (editor-defaults-and-theta): no blank entry


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


def _settle():
    """Let deferred work run: the cell opens on a zero-delay timer after a click, its list opens after the
    editor is shown, and Qt 5 queues a delegate's commit and close."""
    for _ in range(3):
        QtWidgets.QApplication.processEvents()


def _open_cell_editor(tab, row, name):
    """Open the cell's editor through the view's own entry point (what every edit trigger calls), and return it."""
    index = tab.angle_table.model().index(row, fs.PER_ANGLE_NAMES.index(name))
    tab.angle_table.edit(index)
    _settle()
    return tab.angle_table.indexWidget(index)


def _list_shown(combo):
    """Is the drop-down's list open?"""
    return combo.view().isVisible()


def _commit(editor):
    """Press Return in a cell's editor and let the commit happen (typed text in the direct-beam cell)."""
    QTest.keyClick(editor, QtCore.Qt.Key_Return)
    _settle()


def _dismiss(editor):
    """Leave a cell without choosing: Escape closes its open list, and Escape on the closed drop-down closes
    the editor and returns to the grid (APG)."""
    if editor.view().isVisible():
        QTest.keyClick(editor.view(), QtCore.Qt.Key_Escape)
        _settle()
    if not sip.isdeleted(editor) and not editor.isHidden():  # the view hides an editor it closes
        QTest.keyClick(editor, QtCore.Qt.Key_Escape)
        _settle()


def _choose(combo, text):
    """Choose `text` the way a user does: open the list if it is not open (one click on the drop-down), move
    to the item with the arrow keys inside the list, and press Return there. With the list closed, the arrow
    keys do not change the value at all (C10)."""
    target = combo.findText(text)
    assert target >= 0, text
    _choose_at(combo, target)


def _choose_at(combo, target):
    """`_choose` by position: two entries can show the same text (a raw string "True" beside the entry True)."""
    if not combo.view().isVisible():
        QTest.mouseClick(combo, QtCore.Qt.LeftButton)
        _settle()
    view = combo.view()
    assert view.isVisible()
    for _ in range(combo.count() + 1):
        row = view.currentIndex().row()
        if row == target:
            break
        QTest.keyClick(view, QtCore.Qt.Key_Down if target > row else QtCore.Qt.Key_Up)
    assert view.currentIndex().row() == target
    QTest.keyClick(view, QtCore.Qt.Key_Return)
    _settle()


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
    """V3', V15 (v3, T-1): a cell opens with its list shown, and an open list ignores the wheel whatever
    the combo does. So the list is closed first (Escape in it) to leave the closed drop-down C10 leaves,
    and the wheel goes to that."""
    path, _ = _direct_beam_settings(tmp_path, ["db_a.dat", "db_z.dat"], {"useBS": [1, 1, 1, 1]})
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    editor = _open_cell_editor(tab, row, name)
    assert isinstance(editor, QtWidgets.QComboBox)
    assert editor.count() > 1  # something the wheel could move to
    QTest.keyClick(editor.view(), QtCore.Qt.Key_Escape)
    _settle()
    assert not editor.view().isVisible()
    assert tab.angle_table.indexWidget(tab.angle_table.model().index(row, fs.PER_ANGLE_NAMES.index(name))) is editor
    shown = editor.currentText()
    _wheel(editor, _away(editor))
    assert editor.currentText() == shown
    _dismiss(editor)
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
    _commit(editor)
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
    _commit(editor)
    assert tab.document.get("DBname")[0] == "typed.dat"


def _offered(tab):
    editor = _open_cell_editor(tab, 0, "DBname")
    names = _items(editor)
    _dismiss(editor)
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


def test_a_loaded_case_variant_shows_as_the_file_spells_it_and_is_kept(tmp_path, monkeypatch):
    """V11, as v2's C9 rewrites it. The reducer lower-cases (nr_reduction_calc.py:82), so 'meantheta' is
    meanTheta, but the file's spelling is the file's: the cell and its list read as the file does, and
    leaving the cell without choosing writes nothing. v1 displayed the declared spelling ('meanTheta'), and
    re-choosing it rewrote the file (PR #36, finding 2)."""
    path = _settings_file(tmp_path, {**_THREE_ANGLES, "method_per_run": ["meantheta"] * 3})
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    assert _shown(tab, 0, "method_per_run") == ("meantheta", False)
    editor = _open_cell_editor(tab, 0, "method_per_run")
    assert editor.currentText() == "meantheta"
    _dismiss(editor)
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
    _dismiss(editor)
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
    editor = _open_cell_editor(tab, 1, "method_per_run")
    assert editor.currentText() == "constantQ"  # the editor opens on what the cell shows
    _choose(editor, "meanTheta")
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



# Added when frame rows of the battery survived GREEN (ledger scripts/mutations-editor-combos.py).


def test_an_implied_value_is_drawn_in_the_placeholder_colour():
    """F12: italics and the placeholder colour together mark a value the cell does not hold."""
    tab = SettingsEditorTab(SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["constantQ"]}))
    column = fs.PER_ANGLE_NAMES.index("method_per_run")
    delegate = tab.angle_table.itemDelegateForColumn(column)

    def colours(row):
        option = QtWidgets.QStyleOptionViewItem()
        delegate.initStyleOption(option, tab.angle_table.model().index(row, column))
        return option.palette.color(QtGui.QPalette.Text), option.palette.color(QtGui.QPalette.PlaceholderText)

    text, placeholder = colours(1)
    assert text == placeholder
    text, placeholder = colours(0)
    assert text != placeholder


def test_an_implied_cell_says_where_its_value_comes_from():
    """F13: the tooltip a hover shows, through the view's own tooltip path (a help event to the viewport)."""
    tab = SettingsEditorTab(SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["constantQ"]}))
    tab.resize(1100, 700)
    tab.show()
    QTest.qWaitForWindowExposed(tab)
    index = tab.angle_table.model().index(1, fs.PER_ANGLE_NAMES.index("method_per_run"))
    tab.angle_table.scrollTo(index)
    centre = tab.angle_table.visualRect(index).center()
    viewport = tab.angle_table.viewport()
    QtWidgets.QApplication.sendEvent(
        viewport, QtGui.QHelpEvent(QtCore.QEvent.ToolTip, centre, viewport.mapToGlobal(centre)))
    assert "the reduction uses constantQ" in QtWidgets.QToolTip.text()
    QtWidgets.QToolTip.hideText()
    tab.close()


def test_a_drop_down_does_not_take_focus_from_the_wheel():
    """F4: QComboBox defaults to WheelFocus, so a list scrolled past one would leave the focus, and the
    keyboard's arrow keys, on it. Qt only gives focus for a real (spontaneous) wheel event, which a test
    cannot send, so this asserts the policy that decides it."""
    tab = SettingsEditorTab(SettingsDocument.from_dict(_THREE_ANGLES))
    for name in _SCALAR_CHOICES:
        assert tab.editors[name].focusPolicy() == QtCore.Qt.StrongFocus, name
    for name in ("method_per_run", "useBS", "DBname"):
        editor = _open_cell_editor(tab, 0, name)
        assert editor.focusPolicy() == QtCore.Qt.StrongFocus, name
        QTest.keyClick(editor, QtCore.Qt.Key_Escape)
        QtWidgets.QApplication.processEvents()


def test_a_direct_beam_list_at_the_cap_does_not_say_it_was_cut(tmp_path, monkeypatch):
    """F16: exactly MAX_CANDIDATES names is the whole folder."""
    from lr_reduction.settings_document import MAX_CANDIDATES

    path, _ = _direct_beam_settings(tmp_path, [f"db_{k:05d}.dat" for k in range(MAX_CANDIDATES)])
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    editor = _open_cell_editor(tab, 0, "DBname")
    assert editor.count() == MAX_CANDIDATES
    assert editor.toolTip() == ""


def test_a_change_in_another_column_repaints_the_implied_cells(tmp_path, monkeypatch):
    """F23: a DBname typed into a surplus row makes it an angle, so its Q-method cell now implies meanTheta,
    but nothing in that column changed. While other rows stay surplus, the row header keeps its width, so
    nothing else repaints that cell: refresh_marks must."""
    tab = SettingsEditorTab()
    _load(tab, _settings_file(tmp_path, {**_THREE_ANGLES, "useBS": [1] * 6}), monkeypatch)
    tab.resize(1100, 700)
    tab.show()
    QTest.qWaitForWindowExposed(tab)
    method = tab.angle_table.model().index(3, fs.PER_ANGLE_NAMES.index("method_per_run"))
    tab.angle_table.scrollTo(method)
    QtWidgets.QApplication.processEvents()
    assert _shown(tab, 3, "method_per_run") == ("", False)
    painted = []

    class Recorder(QtCore.QObject):
        def eventFilter(self, _watched, event):
            if event.type() == QtCore.QEvent.Paint:
                painted.append(QtGui.QRegion(event.region()))
            return False

    recorder = Recorder()
    tab.angle_table.viewport().installEventFilter(recorder)
    tab.angle_table.item(3, fs.PER_ANGLE_NAMES.index("DBname")).setText("d.dat")
    QtWidgets.QApplication.processEvents()
    assert _shown(tab, 3, "method_per_run") == ("meanTheta", True)
    assert any(region.contains(tab.angle_table.visualRect(method)) for region in painted)
    tab.close()


# --------------------------------------------------------------------------
# editor-combos v2 — the human's gate (PR #36 comment, 2026-10-04 00:30Z): C8 a table drop-down is visible
# and opens on one gesture; C9 its choices read as the file does and re-choosing the held value writes
# nothing; C10 a drop-down lets go of focus after a choice, and only a choice in its open list changes it.
# --------------------------------------------------------------------------

_SIX_ANGLES = {
    "RBnum": [201282 + k for k in range(6)],
    "DBname": [f"db_{c}.dat" for c in "abcdef"],
    "RB_Ymin": [140 + k for k in range(6)],
    "RB_Ymax": [150 + k for k in range(6)],
    "BkgROI": [[120, 130]] * 6,
}
_DROP_DOWN_COLUMNS = ["method_per_run", "useBS", "DBname"]


def _shown_tab(tab):
    tab.resize(1100, 700)
    tab.show()
    QTest.qWaitForWindowExposed(tab)
    QtWidgets.QApplication.setActiveWindow(tab)
    return tab


@pytest.mark.parametrize("name", _DROP_DOWN_COLUMNS)
def test_one_click_opens_a_table_drop_down_and_its_list(tmp_path, monkeypatch, name):
    """V16, C8 (finding 1): one click on the cell, no double-click, opens the drop-down with its list shown."""
    path, _ = _direct_beam_settings(tmp_path, ["db_a.dat", "db_b.dat"])
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    _shown_tab(tab)
    index = tab.angle_table.model().index(1, fs.PER_ANGLE_NAMES.index(name))
    tab.angle_table.scrollTo(index)
    _settle()
    QTest.mouseClick(tab.angle_table.viewport(), QtCore.Qt.LeftButton, pos=tab.angle_table.visualRect(index).center())
    _settle()
    editor = tab.angle_table.indexWidget(index)
    assert isinstance(editor, QtWidgets.QComboBox)
    assert _list_shown(editor)
    _dismiss(editor)
    assert tab.document.changed_vs_seed() == {}
    tab.close()


@pytest.mark.parametrize(
    "key, modifier",
    [(QtCore.Qt.Key_Return, QtCore.Qt.NoModifier), (QtCore.Qt.Key_Enter, QtCore.Qt.KeypadModifier),
     (QtCore.Qt.Key_Space, QtCore.Qt.NoModifier), (QtCore.Qt.Key_F2, QtCore.Qt.NoModifier),
     (QtCore.Qt.Key_Down, QtCore.Qt.AltModifier)],
    ids=["Return", "Enter", "Space", "F2", "Alt+Down"],
)
def test_a_key_on_the_focused_cell_opens_its_drop_down(key, modifier):
    """V16, C8: the cell's widget is entered from the grid with Enter or F2 (APG grid pattern), and the list
    is opened with Space or Alt+Down (APG combobox/menu button)."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["meanTheta"] * 3})))
    column = fs.PER_ANGLE_NAMES.index("method_per_run")
    tab.angle_table.setFocus()
    tab.angle_table.setCurrentCell(1, column)
    QTest.keyClick(tab.angle_table, key, modifier)
    _settle()
    editor = tab.angle_table.indexWidget(tab.angle_table.model().index(1, column))
    assert isinstance(editor, QtWidgets.QComboBox)
    assert _list_shown(editor)
    _dismiss(editor)
    tab.close()


def test_down_in_the_table_moves_to_the_next_row():
    """C8 names Down among the keys that open a drop-down, after the menu-button example. In a table, APG's
    grid pattern gives Down to the grid ("Moves focus one cell down"), and a drop-down column that opened
    on Down would leave no keyboard way down it. Alt+Down opens (APG combobox); Down moves. Deliberate,
    and flagged in the commit."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["meanTheta"] * 3})))
    column = fs.PER_ANGLE_NAMES.index("method_per_run")
    tab.angle_table.setFocus()
    tab.angle_table.setCurrentCell(0, column)
    QTest.keyClick(tab.angle_table, QtCore.Qt.Key_Down)
    _settle()
    assert tab.angle_table.currentRow() == 1
    assert tab.angle_table.indexWidget(tab.angle_table.model().index(1, column)) is None
    tab.close()


class _ArrowRecorder(QtWidgets.QProxyStyle):
    """Records what the style is asked to draw in the table: drop-down arrows, item backgrounds (panels) and
    item bodies (text). The painted control's state, queryable."""

    def __init__(self):
        super().__init__()
        self.arrows, self.panels, self.items = [], [], []

    def drawPrimitive(self, element, option, painter, widget=None):
        if element == QtWidgets.QStyle.PE_IndicatorArrowDown:
            self.arrows.append(QtCore.QRect(option.rect))
        elif element == QtWidgets.QStyle.PE_PanelItemViewItem:
            self.panels.append(QtCore.QRect(option.rect))
        super().drawPrimitive(element, option, painter, widget)

    def drawControl(self, element, option, painter, widget=None):
        if element == QtWidgets.QStyle.CE_ItemViewItem:
            self.items.append(QtCore.QRect(option.rect))
        super().drawControl(element, option, painter, widget)


def test_a_table_drop_down_shows_its_arrow_at_rest():
    """V16, C8: every drop-down cell draws its arrow beside its value without being opened; a text cell
    draws none."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({**_THREE_ANGLES, "useBS": [1, 1, 0]})))
    recorder = _ArrowRecorder()
    tab.angle_table.setStyle(recorder)
    tab.angle_table.viewport().grab()
    model = tab.angle_table.model()
    for name in _DROP_DOWN_COLUMNS:
        cell = tab.angle_table.visualRect(model.index(0, fs.PER_ANGLE_NAMES.index(name)))
        arrows = [arrow for arrow in recorder.arrows if cell.contains(arrow)]
        assert arrows, name
        # the whole cell keeps its item background (selection), and the value's text stays clear of the arrow
        assert cell in recorder.panels, name
        bodies = [item for item in recorder.items if cell.contains(item)]
        assert bodies and not any(body.intersects(arrows[0]) for body in bodies), name
    text_cell = tab.angle_table.visualRect(model.index(0, fs.PER_ANGLE_NAMES.index("ThetaShift")))
    assert not any(text_cell.intersects(arrow) for arrow in recorder.arrows)
    tab.close()


@pytest.mark.parametrize("case", ["lower", "upper"])
def test_the_q_method_choices_read_as_the_file_does_and_re_choosing_writes_nothing(tmp_path, monkeypatch, case):
    """V17, C9 (finding 2, the human's reproduction): a reducer-written file holds 'meantheta'. Its column
    offers the methods in that spelling; re-choosing the held value in row 2 is the identity (v1 wrote
    'meanTheta'); a new choice is written in the file's spelling."""
    spell = str.lower if case == "lower" else str.upper
    path = _settings_file(tmp_path, {**_SIX_ANGLES, "method_per_run": [spell("meanTheta")] * 6})
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    editor = _open_cell_editor(tab, 2, "method_per_run")
    assert _items(editor) == ["", *(spell(choice) for choice in fs.METHOD_CHOICES)]
    _choose(editor, spell("meanTheta"))
    assert repr(tab.document.get("method_per_run")) == repr([spell("meanTheta")] * 6)
    assert tab.document.changed_vs_seed() == {}
    _choose(_open_cell_editor(tab, 2, "method_per_run"), spell("constantQ"))
    assert tab.document.get("method_per_run")[2] == spell("constantQ")


def test_a_fresh_declared_or_mixed_column_offers_the_declared_spellings(tmp_path, monkeypatch):
    """V17, C9: no case variant (a fresh document, declared spellings) or mixed casing -> the declared
    spellings; a held variant is still one of the items (never substituted, C2), so re-choosing it writes
    nothing."""
    fresh = SettingsEditorTab(SettingsDocument.from_dict(_THREE_ANGLES))
    editor = _open_cell_editor(fresh, 0, "method_per_run")
    assert _items(editor) == ["", *fs.METHOD_CHOICES]
    _dismiss(editor)
    mixed = SettingsEditorTab(SettingsDocument.from_dict(
        {**_THREE_ANGLES, "method_per_run": ["meantheta", "constantQ", "meantheta"]}))
    editor = _open_cell_editor(mixed, 1, "method_per_run")
    assert _items(editor) == ["", *fs.METHOD_CHOICES]
    _dismiss(editor)
    editor = _open_cell_editor(mixed, 0, "method_per_run")
    assert set(_items(editor)) == {"", *fs.METHOD_CHOICES, "meantheta"}
    assert editor.currentText() == "meantheta"
    _choose(editor, "meantheta")
    assert mixed.document.changed_vs_seed() == {}


@pytest.mark.parametrize("name, choice, stored", [
    ("method_per_run", "constantQ", "constantQ"), ("useBS", "false", False), ("DBname", "db_c.dat", "db_c.dat")])
def test_a_cell_drop_down_lets_go_of_focus_after_a_choice(tmp_path, monkeypatch, name, choice, stored):
    """V18, C10 (finding 3; widened in v3 for T-2). After a choice in the open list of any cell column,
    including a direct-beam name picked from the folder's list, the choice is written, the drop-down is
    gone, and the table has the focus, so a later arrow key moves in the grid and changes no value."""
    path, _ = _direct_beam_settings(tmp_path, ["db_a.dat", "db_b.dat", "db_c.dat"],
                                    {"method_per_run": ["meanTheta"] * 3, "useBS": [1, 1, 1]})
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    _shown_tab(tab)
    _choose(_open_cell_editor(tab, 1, name), choice)
    assert repr(tab.document.get(name)[1]) == repr(stored)
    assert QtWidgets.QApplication.focusWidget() is tab.angle_table
    QTest.keyClick(tab.angle_table, QtCore.Qt.Key_Down)
    _settle()
    assert repr(tab.document.get(name)[1]) == repr(stored)
    tab.close()


def test_a_scalar_drop_down_lets_go_of_focus_after_a_choice_and_keys_never_step_it():
    """V18, C10: a choice in its list moves the focus to the panel. A closed drop-down reached with Tab opens
    on Down or Space and changes only by a choice in its open list (APG select-only combobox)."""
    tab = _shown_tab(SettingsEditorTab())
    combo = tab.editors["DetResFn"]
    other = next(combo.itemText(i) for i in range(combo.count()) if combo.itemText(i) != combo.currentText())
    _choose(combo, other)
    assert tab.document.get("DetResFn") == other
    assert QtWidgets.QApplication.focusWidget() in (tab.scalar_panel, tab.scalar_panel.focusProxy())
    combo.setFocus()
    _settle()
    letter = next(combo.itemText(i)[0] for i in range(combo.count()) if combo.itemText(i)[0] != other[0])
    keys = (QtCore.Qt.Key_Down, QtCore.Qt.Key_Up, QtCore.Qt.Key_Space, QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter,
            getattr(QtCore.Qt, "Key_" + letter.upper()))
    for key in keys:
        QTest.keyClick(combo, key)
        _settle()
        assert combo.view().isVisible()
        assert tab.document.get("DetResFn") == other
        QTest.keyClick(combo.view(), QtCore.Qt.Key_Escape)
        _settle()
    tab.close()


def test_typing_with_the_direct_beam_list_open_starts_a_new_name(tmp_path, monkeypatch):
    """C8 with C4: the direct-beam cell opens with its list shown, and the list has the keyboard. A typed
    character closes the list and starts a new name in the line edit, as typing into the selected text would;
    Return stores it, and the focus goes back to the table (C10)."""
    path, _ = _direct_beam_settings(tmp_path, ["db_a.dat", "db_b.dat"])
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    _shown_tab(tab)
    editor = _open_cell_editor(tab, 1, "DBname")
    assert editor.view().isVisible()
    assert editor.currentText() == "db_b.dat"  # it opens on the held name
    QTest.keyClicks(editor.view(), "new 1.dat")
    _settle()
    assert not editor.view().isVisible()
    assert editor.currentText() == "new 1.dat"
    _commit(editor)
    assert tab.document.get("DBname") == ["db_a.dat", "new 1.dat", "db_c.dat"]
    assert not isinstance(QtWidgets.QApplication.focusWidget(), QtWidgets.QComboBox)
    tab.close()



def test_leaving_a_cell_without_a_choice_writes_nothing():
    """C2, C7: moving the focus off an open cell commits it (Qt's delegate does so on focus-out), but no
    choice was made, so nothing is written. Here an implied value would otherwise be written out."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["constantQ"]})))
    index = tab.angle_table.model().index(1, fs.PER_ANGLE_NAMES.index("method_per_run"))
    editor = _open_cell_editor(tab, 1, "method_per_run")
    QTest.keyClick(editor.view(), QtCore.Qt.Key_Escape)
    _settle()
    tab.angle_table.setFocus()
    _settle()
    assert tab.angle_table.indexWidget(index) is None
    assert tab.document.changed_vs_seed() == {}
    tab.close()


def test_return_on_a_text_cell_opens_no_drop_down():
    """C8's keys open a drop-down cell. Return on a plain text cell is not one of them."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict(_THREE_ANGLES)))
    column = fs.PER_ANGLE_NAMES.index("ThetaShift")
    tab.angle_table.setFocus()
    tab.angle_table.setCurrentCell(1, column)
    QTest.keyClick(tab.angle_table, QtCore.Qt.Key_Return)
    _settle()
    assert tab.angle_table.indexWidget(tab.angle_table.model().index(1, column)) is None
    tab.close()



# --------------------------------------------------------------------------
# editor-combos v3 — review/editor-combos @ d3ee364: C11 only a deliberate choice writes (U-1: click + Return
# on a direct-beam cell holding a name outside the listed folder wrote the folder's first file); C8' the
# grid convention in the table; C9' re-choosing the value a cell shows is the identity, implied included.
# --------------------------------------------------------------------------

_V19_HELD = {
    # column: {held state: (settings over _THREE_ANGLES, the row, Add first?)}
    "method_per_run": {
        "listed": ({"method_per_run": ["meanTheta"] * 3}, 1, False),
        "empty": ({"method_per_run": ["meanTheta", None, "meanTheta"]}, 1, False),
        "new-row": ({"method_per_run": ["meanTheta"] * 3}, 3, True),
        "implied": ({"method_per_run": ["constantQ"]}, 1, False),
    },
    "useBS": {
        "listed": ({"useBS": [1, 1, 0]}, 1, False),
        "empty": ({"useBS": [1, None, 1]}, 1, False),
        "new-row": ({"useBS": [1, 1, 0]}, 3, True),
        "implied": ({"useBS": []}, 1, False),
    },
    "DBname": {
        "listed": ({}, 1, False),  # db_b.dat, in the folder
        "outside": ({"DBname": ["A2_div10_Cd.txt", "A2_div10_Cd.txt", "db_c.dat"]}, 1, False),
        "empty": ({"DBname": ["db_a.dat", None, "db_c.dat"]}, 1, False),
        "new-row": ({}, 3, True),
    },
}
_V19_SHOWS = {
    # what each held state's cell opens showing: the value as held (C2), the implied one (C7), or "" for unset
    "method_per_run": {"listed": "meanTheta", "empty": "", "new-row": "", "implied": "constantQ"},
    "useBS": {"listed": "true", "empty": "", "new-row": "", "implied": "true"},
    "DBname": {"listed": "db_b.dat", "outside": "A2_div10_Cd.txt", "empty": "", "new-row": ""},
}
_V19_GESTURES = ["return", "escape", "tab", "click-away", "arrow-return", "click-item"]
_V19_CASES = [
    pytest.param(name, held, gesture, id=f"{name}-{held}-{gesture}")
    for name, states in _V19_HELD.items() for held in states for gesture in _V19_GESTURES
]


def _saved_text(tab, path):
    tab.document.save(path)
    return path.read_bytes()


@pytest.mark.parametrize("name, held, gesture", _V19_CASES)
def test_only_a_deliberate_choice_in_a_cell_writes(tmp_path, monkeypatch, name, held, gesture):
    """V19, C11: three columns x held state x gesture, each cell opened with one click on a shown tab. Only
    a deliberate choice writes: a move to another item with the arrow keys and then Return, or a click on an
    item. Return, Escape, Tab or a click elsewhere, with no move, leave the cell exactly as held. The "outside"
    state is the reducer-written norm for the direct-beam column: a held name the listed folder does not
    contain. For the file, the names live in shared/transmission/Aug2026/ while the editor lists
    shared/transmission; here the folder simply lacks the held name."""
    values, row, add = _V19_HELD[name][held]
    path, _ = _direct_beam_settings(tmp_path, ["db_a.dat", "db_b.dat", "db_c.dat"], values)
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    _shown_tab(tab)
    if add:
        QTest.mouseClick(tab.add_angle_button, QtCore.Qt.LeftButton)
        _settle()
    before_changes = tab.document.changed_vs_seed()
    before_text = _saved_text(tab, tmp_path / "before.json")
    index = tab.angle_table.model().index(row, fs.PER_ANGLE_NAMES.index(name))
    tab.angle_table.scrollTo(index)
    _settle()
    QTest.mouseClick(tab.angle_table.viewport(), QtCore.Qt.LeftButton, pos=tab.angle_table.visualRect(index).center())
    _settle()
    editor = tab.angle_table.indexWidget(index)
    assert isinstance(editor, QtWidgets.QComboBox) and editor.view().isVisible()
    view = editor.view()
    shown = editor.currentText()
    # The cell opens on what it shows, and its list on that item, or on none when the held name is not listed.
    # Qt makes row 0 current there when the list takes the focus, and Return then chose it (U-1).
    assert shown == _V19_SHOWS[name][held]
    assert view.currentIndex().row() == editor.findText(shown)
    chosen = None
    if gesture == "return":
        QTest.keyClick(view, QtCore.Qt.Key_Return)
    elif gesture == "escape":
        QTest.keyClick(view, QtCore.Qt.Key_Escape)
    elif gesture == "tab":
        QTest.keyClick(view, QtCore.Qt.Key_Tab)
    elif gesture == "click-away":
        QTest.mouseClick(tab.report, QtCore.Qt.LeftButton)
    elif gesture == "arrow-return":
        QTest.keyClick(view, QtCore.Qt.Key_Down)
        chosen = editor.itemText(view.currentIndex().row())
        if chosen in ("", shown):  # move once more to a different, non-empty item
            QTest.keyClick(view, QtCore.Qt.Key_Down)
            chosen = editor.itemText(view.currentIndex().row())
        QTest.keyClick(view, QtCore.Qt.Key_Return)
    else:
        target = next(i for i in range(editor.count()) if editor.itemText(i) not in ("", shown))
        chosen = editor.itemText(target)
        # The item's rect can be wider than the list's viewport (measured: 217 px in a 101 px list), so the
        # click goes to the viewport's horizontal centre on the item's row, inside both.
        row_rect = view.visualRect(view.model().index(target, 0))
        QTest.mouseClick(view.viewport(), QtCore.Qt.LeftButton,
                         pos=QtCore.QPoint(view.viewport().rect().center().x(), row_rect.center().y()))
    _settle()
    if gesture in ("return", "arrow-return", "click-item"):
        # Return, and a choice, close the cell and give the focus back to the table (C10), whether or not the list
        # had a row to choose.
        assert tab.angle_table.indexWidget(index) is None
        assert QtWidgets.QApplication.focusWidget() is tab.angle_table
    if chosen is None:
        # Then the user goes on to another cell. Qt's delegate commits an editor that is still open when it loses
        # the focus, and that commit must write nothing either.
        if not sip.isdeleted(editor) and view.isVisible():
            QTest.keyClick(view, QtCore.Qt.Key_Escape)
        tab.angle_table.setFocus()
        _settle()
        assert tab.document.changed_vs_seed() == before_changes
        assert _saved_text(tab, tmp_path / "after.json") == before_text
    else:
        assert chosen not in ("", shown)
        assert repr(tab.document.get(name)[row]) == repr(fs.get(name).coerce_element(chosen))
    tab.close()


@pytest.mark.parametrize("name, letter", [("method_per_run", "c"), ("useBS", "f")])
def test_typing_to_another_item_in_a_cells_list_and_pressing_return_writes_it(name, letter):
    """C11: in a list that is not editable, typing moves to the item it names, a deliberate move as an arrow key
    is. Return then chooses that item, and the focus goes back to the table (C10)."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict(
        {**_THREE_ANGLES, "method_per_run": ["meanTheta"] * 3, "useBS": [1, 1, 1]})))
    editor = _open_cell_editor(tab, 1, name)
    shown = editor.currentText()
    QTest.keyClick(editor.view(), letter)
    chosen = editor.itemText(editor.view().currentIndex().row())
    assert chosen not in ("", shown) and chosen.startswith(letter)
    QTest.keyClick(editor.view(), QtCore.Qt.Key_Return)
    _settle()
    assert repr(tab.document.get(name)[1]) == repr(fs.get(name).coerce_element(chosen))
    assert QtWidgets.QApplication.focusWidget() is tab.angle_table
    tab.close()


def test_choosing_the_implied_value_a_cell_shows_writes_nothing():
    """V21, C9' (the Integrator's D-b, decided): in a compact column, choosing the implied value a cell
    shows is the identity, and the list stays compact. Choosing a different value writes the list out
    (G9)."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict(
        {**_THREE_ANGLES, "method_per_run": ["constantQ"], "useBS": []})))
    for name, implied in (("method_per_run", "constantQ"), ("useBS", "true")):
        editor = _open_cell_editor(tab, 1, name)
        assert editor.currentText() == implied
        QTest.keyClick(editor.view(), QtCore.Qt.Key_Down)
        QTest.keyClick(editor.view(), QtCore.Qt.Key_Up)  # a deliberate move, back to the shown value
        QTest.keyClick(editor.view(), QtCore.Qt.Key_Return)
        _settle()
        assert tab.document.changed_vs_seed() == {}, name
    assert tab.document.get("method_per_run") == ["constantQ"] and tab.document.get("useBS") == []
    _choose(_open_cell_editor(tab, 1, "method_per_run"), "meanTheta")
    assert tab.document.get("method_per_run") == ["constantQ", "meanTheta", "constantQ"]
    tab.close()


@pytest.mark.parametrize("name, held", [("DBname", "outside"), ("DBname", "listed"), ("method_per_run", "listed"),
                                        ("useBS", "listed"), ("method_per_run", "implied")])
def test_hovering_over_another_item_and_pressing_return_writes_nothing(tmp_path, monkeypatch, name, held):
    """C11: a hover moves the list's current row, and is not a choice. Return after it, with no arrow-key
    move, leaves the cell as held. QComboBox's own list would select the hovered row (measured on the v2
    code)."""
    values, row, _ = _V19_HELD[name][held]
    path, _ = _direct_beam_settings(tmp_path, ["db_a.dat", "db_b.dat", "db_c.dat"], values)
    tab = SettingsEditorTab()
    _load(tab, path, monkeypatch)
    _shown_tab(tab)
    before = _saved_text(tab, tmp_path / "before.json")
    editor = _open_cell_editor(tab, row, name)
    view = editor.view()
    other = next(i for i in range(editor.count()) if editor.itemText(i) not in ("", editor.currentText()))
    x = view.viewport().rect().center().x()
    # Two moves: QTest's mouse position is global, so a move to where an earlier test left the pointer would
    # send no event and hover nothing (measured: the full suite left row 1 current).
    QTest.mouseMove(view.viewport(), QtCore.QPoint(x, view.visualRect(view.currentIndex()).center().y()))
    QTest.mouseMove(view.viewport(), QtCore.QPoint(x, view.visualRect(view.model().index(other, 0)).center().y()))
    _settle()
    assert view.currentIndex().row() == other  # the hover moved the current row
    QTest.keyClick(view, QtCore.Qt.Key_Return)
    _settle()
    assert tab.document.changed_vs_seed() == {}
    assert _saved_text(tab, tmp_path / "after.json") == before
    tab.close()



# --------------------------------------------------------------------------
# editor-defaults-and-theta — a new file starts at gaussian / 1.0; "Apply theta calculation" offers False / True /
# trust sample angle and stores the canonical value
# --------------------------------------------------------------------------

_THETA_ENTRIES = ["False", "True", "trust sample angle"]


def _scalar_label(tab, name):
    editor = tab.editors[name]
    return editor.parentWidget().layout().labelForField(editor).text()


def test_a_new_file_in_the_tab_starts_at_gaussian_and_1():
    """V1, D1: the tab opened with no document. The starting values are the seed, so nothing shows as changed."""
    tab = SettingsEditorTab()
    assert tab.editors["DetResFn"].currentText() == "gaussian"
    assert tab.editors["DetSigma"].text() == "1.0"
    assert tab.document.changed_vs_seed() == {}


def test_a_new_file_saved_from_the_tab_states_the_resolution_pair_and_false(tmp_path, monkeypatch):
    """§5, common: open the tab, add an angle, save."""
    tab = SettingsEditorTab()
    tab.add_angle_button.click()
    target = tmp_path / "new.json"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *_a, **_k: (str(target), "")))
    tab.save_settings()
    saved = json.loads(target.read_text())
    assert (saved["DetResFn"], saved["DetSigma"]) == ("gaussian", 1.0)
    assert saved["useCalcTheta"] is False


def test_the_theta_control_reads_apply_theta_calculation_and_offers_three_entries():
    """V2, D3: exactly the scientists' three entries, in their order, and no blank one."""
    tab = SettingsEditorTab()
    assert _scalar_label(tab, "useCalcTheta") == "Apply theta calculation"
    assert _items(tab.editors["useCalcTheta"]) == _THETA_ENTRIES


@pytest.mark.parametrize("entry, stored", [("True", "detector_angle"), ("trust sample angle", "sample_angle"),
                                           ("False", False)])
def test_choosing_a_theta_entry_stores_its_canonical_value(entry, stored):
    """V3, D4: from the keyboard, through the combo's own list. The stored value is the canonical one, never the
    label text, and the off entry stores False itself (not 0, not None)."""
    held = "sample_angle" if stored is False else False
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"useCalcTheta": held})))
    combo = tab.editors["useCalcTheta"]
    combo.setFocus()
    QTest.keyClick(combo, QtCore.Qt.Key_Space)  # the menu-button key opens the list
    _settle()
    _choose(combo, entry)
    value = tab.document.get("useCalcTheta")
    assert value == stored and type(value) is type(stored)
    tab.close()


@pytest.mark.parametrize("loaded, entry", [
    (True, "True"), ("detector_angle", "True"), ("Detector_Angle", "True"),
    ("sample_angle", "trust sample angle"), ("SAMPLE_ANGLE", "trust sample angle"),
    (False, "False"), (None, "False"), (0, "False"), ("", "False"),
])
def test_loading_an_accepted_theta_spelling_selects_its_entry(loaded, entry):
    """V4, D5: through set_document, the Load path; the panel is quiet about it."""
    tab = SettingsEditorTab()
    tab.set_document(SettingsDocument.from_dict({"useCalcTheta": loaded}))
    assert tab.editors["useCalcTheta"].currentText() == entry
    assert not [message for message in tab.document.validate() if "useCalcTheta" in message]


@pytest.mark.parametrize("loaded", ["true", "TRUE", "True", "False", "trust sample angle", 1, "detector",
                                    ["detector_angle"]])
def test_a_theta_value_the_reducer_rejects_is_shown_as_itself_kept_and_reported(loaded):
    """V6, D6: the reducer raises on these, so the editor does not guess what they meant. A string that spells
    an entry ("True", "False", "trust sample angle") is an entry of its own, not that entry: it is found by the
    value it holds, so the selection is asserted by position and type (advisory A1), not by its text."""
    tab = SettingsEditorTab()
    tab.set_document(SettingsDocument.from_dict({"useCalcTheta": loaded}))
    combo = tab.editors["useCalcTheta"]
    assert combo.currentIndex() == 3
    assert combo.currentData() == loaded and type(combo.currentData()) is type(loaded)
    assert combo.currentText() == str(loaded)
    assert _items(combo) == [*_THETA_ENTRIES, str(loaded)]
    held = tab.document.get("useCalcTheta")
    assert held == loaded and type(held) is type(loaded)
    assert "useCalcTheta" in tab.report.toPlainText()


def test_a_document_given_to_the_tab_keeps_its_own_resolution_values():
    """V7, D2: the tab does not re-initialise a document it was given."""
    tab = SettingsEditorTab(SettingsDocument())
    assert (tab.document.get("DetResFn"), tab.document.get("DetSigma")) == ("rectangular", 0.8)
    assert tab.editors["DetResFn"].currentText() == "rectangular" and tab.editors["DetSigma"].text() == "0.8"


# held state: (settings, the position of the entry shown, the position of another entry, the value it stores)
_THETA_HELD = {
    "false": ({}, 0, 1, "detector_angle"),
    "detector_angle": ({"useCalcTheta": "detector_angle"}, 1, 2, "sample_angle"),
    "sample_angle": ({"useCalcTheta": "sample_angle"}, 2, 0, False),
    "rejected": ({"useCalcTheta": "true"}, 3, 1, "detector_angle"),
    "rejected-True": ({"useCalcTheta": "True"}, 3, 1, "detector_angle"),  # a raw entry that spells the entry True
}
_THETA_REJECTED = ("rejected", "rejected-True")


def _theta_state(tab):
    combo = tab.editors["useCalcTheta"]
    held = tab.document.get("useCalcTheta")
    return (held, type(held), combo.currentIndex(), _items(combo), tab.document.changed_vs_seed(),
            tab.report.toPlainText())


@pytest.mark.parametrize("operation", ["choose-shown", "choose-another", "escape"])
@pytest.mark.parametrize("held", list(_THETA_HELD))
def test_each_gesture_on_the_theta_control_in_each_held_state(held, operation):
    """§3's operation x state table for useCalcTheta, the rows with the list open. Re-choosing the shown entry is
    the identity: a rejected value, its entry and its problem line stay. Another entry writes its canonical value
    as the one change and clears a problem line; once the choice has returned, exactly the three entries are left
    (D3′, V10). Escape in the open list writes nothing. The closed-combo rows (Tab out, click-away, wheel, Up,
    Down) are V8, test_the_closed_theta_control_ignores_every_gesture_in_each_held_state; the second-Load rows
    are V9's two tests."""
    values, shown, other, stored = _THETA_HELD[held]
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict(values)))
    combo = tab.editors["useCalcTheta"]
    before = _theta_state(tab)
    assert combo.currentIndex() == shown
    if operation == "escape":
        QTest.mouseClick(combo, QtCore.Qt.LeftButton)
        _settle()
        assert combo.view().isVisible()
        QTest.keyClick(combo.view(), QtCore.Qt.Key_Escape)
        _settle()
    else:
        _choose_at(combo, shown if operation == "choose-shown" else other)
    after = tab.document.get("useCalcTheta")
    reported = any("useCalcTheta" in message for message in tab.document.validate())
    if operation == "choose-another":
        assert after == stored and type(after) is type(stored)
        assert list(tab.document.changed_vs_seed()) == ["useCalcTheta"]
        assert not reported
        assert _items(combo) == _THETA_ENTRIES and combo.currentIndex() == other
    else:
        assert _theta_state(tab) == before
        assert reported is (held in _THETA_REJECTED)
    tab.close()


@pytest.mark.parametrize("gesture", ["tab-out", "click-away", "wheel", "up", "down"])
@pytest.mark.parametrize("held", list(_THETA_HELD))
def test_the_closed_theta_control_ignores_every_gesture_in_each_held_state(held, gesture):
    """V8, §3 (T-1): with the combo closed and focused, and no deliberate choice, Tab out, a click elsewhere,
    the wheel, Up and Down change nothing. The held value and its type, the entry shown, the item list, "Changed
    from the seed" and the report all stay. Up and Down open the list, the scalars' menu-button convention, and
    do not step the value; the list is then closed with Escape."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict(_THETA_HELD[held][0])))
    combo = tab.editors["useCalcTheta"]
    combo.setFocus()
    _settle()
    assert combo.hasFocus() and not combo.view().isVisible()
    before = _theta_state(tab)
    if gesture == "tab-out":
        QTest.keyClick(combo, QtCore.Qt.Key_Tab)
        _settle()
        assert not combo.hasFocus()
    elif gesture == "click-away":
        QTest.mouseClick(tab.report, QtCore.Qt.LeftButton)
        _settle()
        assert not combo.hasFocus()
    elif gesture == "wheel":
        _wheel(combo, _away(combo))
        _settle()
    else:
        QTest.keyClick(combo, QtCore.Qt.Key_Up if gesture == "up" else QtCore.Qt.Key_Down)
        _settle()
        assert combo.view().isVisible()
        QTest.keyClick(combo.view(), QtCore.Qt.Key_Escape)
        _settle()
    assert _theta_state(tab) == before
    tab.close()


@pytest.mark.parametrize("held", list(_THETA_HELD))
def test_loading_a_second_file_that_omits_the_key_leaves_exactly_the_three_entries(held):
    """V9, the rejection's D-1 reproduced from every held state: the second file holds no useCalcTheta, so the
    document holds False, the control shows False, and the list is exactly the three entries (D3′)."""
    tab = SettingsEditorTab()
    tab.set_document(SettingsDocument.from_dict(_THETA_HELD[held][0]))
    tab.set_document(SettingsDocument.from_dict({"Sname": "second"}))
    combo = tab.editors["useCalcTheta"]
    assert tab.document.get("useCalcTheta") is False
    assert combo.currentIndex() == 0 and _items(combo) == _THETA_ENTRIES
    assert not [message for message in tab.document.validate() if "useCalcTheta" in message]


@pytest.mark.parametrize("second, index, raw", [
    ("detector_angle", 1, None), (True, 1, None), ("sample_angle", 2, None), ("yes", 3, "yes"), ("True", 3, "True"),
])
@pytest.mark.parametrize("held", list(_THETA_HELD))
def test_loading_a_second_file_shows_its_value_with_only_its_own_raw_entry(held, second, index, raw):
    """V9, D3′: after a second Load the list is the three entries plus, at most, the raw form of the value now
    held, found by type and value. From the string "True", a JSON true is the entry True (position 1), and no
    fourth entry is left."""
    tab = SettingsEditorTab()
    tab.set_document(SettingsDocument.from_dict(_THETA_HELD[held][0]))
    tab.set_document(SettingsDocument.from_dict({"useCalcTheta": second}))
    combo = tab.editors["useCalcTheta"]
    assert _items(combo) == _THETA_ENTRIES + ([raw] if raw is not None else [])
    assert combo.currentIndex() == index
    if raw is not None:
        assert combo.currentData() == second and type(combo.currentData()) is type(second)


@pytest.mark.parametrize("held", _THETA_REJECTED)
def test_choosing_true_over_a_rejected_value_leaves_exactly_the_three_entries(held):
    """V10, D3′: the choice replaces the raw value. Once it has returned, the raw entry is gone, the problem line
    with it, and "Changed from the seed" has the one line. From the string "True" the entry True differs only by
    position, and the choice is still made."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict(_THETA_HELD[held][0])))
    combo = tab.editors["useCalcTheta"]
    _choose_at(combo, 1)
    assert tab.document.get("useCalcTheta") == "detector_angle"
    assert _items(combo) == _THETA_ENTRIES and combo.currentIndex() == 1
    assert not [message for message in tab.document.validate() if "useCalcTheta" in message]
    assert list(tab.document.changed_vs_seed()) == ["useCalcTheta"]
    tab.close()


def test_a_raw_entry_in_a_plain_drop_down_goes_when_its_value_does():
    """V11, D3′'s plain path, on DetResFn: its list is the declared spellings, and a held value outside them is
    an entry of its own only while it is held."""
    tab = SettingsEditorTab()
    tab.set_document(SettingsDocument.from_dict({"DetResFn": "foo"}))
    combo = tab.editors["DetResFn"]
    assert combo.currentText() == "foo" and _items(combo) == [*fs.DET_RES_CHOICES, "foo"]
    assert any("DetResFn" in message for message in tab.document.validate())
    tab.set_document(SettingsDocument.from_dict({"DetResFn": "gaussian"}))
    assert _items(combo) == list(fs.DET_RES_CHOICES) and combo.currentText() == "gaussian"


# --------------------------------------------------------------------------
# editor-paths-header — IPTS and the two input paths at the top of the tab; derived paths are shown, never written
# --------------------------------------------------------------------------

_HEADER_FIELDS = ("experiment_id", "_NEXUSpathRB_override", "_DBpath_override")
_PATH_FIELDS = ("_NEXUSpathRB_override", "_DBpath_override")
_PATH_TAIL = {"_NEXUSpathRB_override": "nexus", "_DBpath_override": "shared/transmission"}
_NO_IPTS = "set an IPTS or type a path"


def _type_into(editor, text):
    """Replace a line edit's text as a user does: focus it, select all, type (or Delete), Return."""
    editor.setFocus()
    _settle()
    editor.selectAll()
    if text:
        QTest.keyClicks(editor, text)
    else:
        QTest.keyClick(editor, QtCore.Qt.Key_Delete)
    QTest.keyClick(editor, QtCore.Qt.Key_Return)
    _settle()


def _focus_through(tab, editor):
    """Focus a control and leave it without typing."""
    editor.setFocus()
    _settle()
    tab.angle_table.setFocus()
    _settle()


def test_the_header_holds_the_ipts_and_the_two_input_paths_and_the_list_does_not():
    """V1, P1: moved, not duplicated. Each name is one widget, in the header above the Angles table, and the
    scrolling list has no editor for it."""
    tab = SettingsEditorTab()
    header = tab.paths_header
    for name in _HEADER_FIELDS:
        assert header.isAncestorOf(tab.editors[name]), name
        assert not tab.scalar_panel.isAncestorOf(tab.editors[name]), name
    listed = {widget.toolTip().split(" — ")[0] for widget in tab.scalar_panel.findChildren(QtWidgets.QWidget)}
    assert not listed & set(_HEADER_FIELDS)
    layout = tab.layout()
    order = [layout.itemAt(i).widget() for i in range(layout.count())]
    assert order.index(header) < order.index(tab.findChild(QtWidgets.QSplitter))


def test_typing_an_ipts_number_updates_both_derived_paths_and_writes_no_override():
    """V2, P2, P3, P5: the document gets IPTS-36119; both paths display the derived folders as placeholders; both
    overrides stay None."""
    tab = _shown_tab(SettingsEditorTab())
    _type_into(tab.editors["experiment_id"], "36119")
    assert tab.document.get("experiment_id") == "IPTS-36119"
    assert tab.editors["experiment_id"].text() == "IPTS-36119"
    for name in _PATH_FIELDS:
        assert tab.document.get(name) is None
        assert tab.editors[name].text() == ""
        assert tab.editors[name].placeholderText() == f"/SNS/REF_L/IPTS-36119/{_PATH_TAIL[name]}"
    tab.close()


def test_a_save_after_an_ipts_writes_null_for_both_overrides(tmp_path, monkeypatch):
    """V3, Q5: a written override freezes an absolute path into the file, so a derived one is saved as null."""
    tab = _shown_tab(SettingsEditorTab())
    _type_into(tab.editors["experiment_id"], "36119")
    target = tmp_path / "new.json"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *_a, **_k: (str(target), "")))
    tab.save_settings()
    saved = json.loads(target.read_text())
    assert saved["experiment_id"] == "IPTS-36119"
    assert saved["_NEXUSpathRB_override"] is None and saved["_DBpath_override"] is None
    tab.close()


@pytest.mark.parametrize("name", _PATH_FIELDS)
def test_a_typed_path_is_the_override_and_does_not_follow_the_ipts(name):
    """V4, P3: typing sets the override; a later IPTS leaves it, while the other path follows."""
    tab = _shown_tab(SettingsEditorTab())
    _type_into(tab.editors[name], "/data/typed")
    assert tab.document.get(name) == "/data/typed"
    _type_into(tab.editors["experiment_id"], "36119")
    assert tab.document.get(name) == "/data/typed" and tab.editors[name].text() == "/data/typed"
    other = next(field for field in _PATH_FIELDS if field != name)
    assert tab.document.get(other) is None
    assert tab.editors[other].placeholderText() == f"/SNS/REF_L/IPTS-36119/{_PATH_TAIL[other]}"
    tab.close()


@pytest.mark.parametrize("name", _PATH_FIELDS)
def test_clearing_a_path_returns_it_to_derived(name):
    """V5, P4: None, not "" (an override of "" is Path("") to the reduction, the current directory)."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "IPTS-1", name: "/data/typed"})))
    _type_into(tab.editors[name], "")
    assert tab.document.get(name) is None
    assert tab.editors[name].text() == "" and tab.editors[name].placeholderText() == f"/SNS/REF_L/IPTS-1/{_PATH_TAIL[name]}"
    tab.close()


@pytest.mark.parametrize("held", [None, "/data/loaded"], ids=["derived", "loaded-override"])
@pytest.mark.parametrize("name", _PATH_FIELDS)
def test_focusing_through_a_path_without_typing_writes_nothing(name, held):
    """V6, P3: the real focus path on a shown tab. A line edit reports editingFinished on every focus-out; with no
    typing nothing is written (the base's data_x_range focus-out corruption is the precedent)."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "IPTS-1", name: held})))
    _focus_through(tab, tab.editors[name])
    _focus_through(tab, tab.editors["experiment_id"])
    assert tab.document.changed_vs_seed() == {}
    assert tab.document.get(name) == held
    tab.close()


def test_focusing_through_an_ipts_the_file_spelled_as_a_number_writes_nothing():
    """V6's IPTS leg (battery F2): a file may hold a bare number. Normalising is for what the user types, so a focus
    change with no typing leaves the file's spelling as loaded."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "36119"})))
    _focus_through(tab, tab.editors["experiment_id"])
    assert tab.document.get("experiment_id") == "36119"
    assert tab.document.changed_vs_seed() == {}
    tab.close()


@pytest.mark.parametrize("ipts", ["../x", "/abs"])
def test_an_ipts_the_panel_reports_is_not_shown_as_the_reductions_folder(ipts):
    """§5, pathological (battery F6): a traversed or absolute IPTS is kept as loaded and reported, and neither path
    presents a folder derived from it; each says the IPTS is not a folder name, not that none is set."""
    tab = SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": ipts}))
    assert tab.document.get("experiment_id") == ipts
    assert any("experiment_id" in message for message in tab.document.validate())
    for name in _PATH_FIELDS:
        assert tab.editors[name].text() == ""
        assert tab.editors[name].placeholderText() == "the IPTS is not a folder name; type a path"


@pytest.mark.parametrize("name", _PATH_FIELDS)
def test_a_second_load_shows_the_second_files_paths_and_nothing_of_the_first(name):
    """V7, P6: file A has an override, file B none; the control shows B's derived path, not A's text."""
    tab = SettingsEditorTab()
    tab.set_document(SettingsDocument.from_dict({"experiment_id": "IPTS-1", name: "/data/a"}))
    assert tab.editors[name].text() == "/data/a"
    tab.set_document(SettingsDocument.from_dict({"experiment_id": "IPTS-2"}))
    assert tab.editors[name].text() == ""
    assert tab.editors[name].placeholderText() == f"/SNS/REF_L/IPTS-2/{_PATH_TAIL[name]}"


@pytest.mark.parametrize("name", _PATH_FIELDS)
def test_browse_sets_the_override_and_a_cancelled_browse_writes_nothing(name, monkeypatch):
    """V8, P3: a folder chosen with the control's Browse button is an explicit edit; a cancelled dialog ("") is
    not."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "IPTS-1"})))
    monkeypatch.setattr(QtWidgets.QFileDialog, "getExistingDirectory", staticmethod(lambda *_a, **_k: ""))
    QTest.mouseClick(tab.path_browse[name], QtCore.Qt.LeftButton)
    _settle()
    assert tab.document.changed_vs_seed() == {}
    monkeypatch.setattr(QtWidgets.QFileDialog, "getExistingDirectory",
                        staticmethod(lambda *_a, **_k: "/data/browsed"))
    QTest.mouseClick(tab.path_browse[name], QtCore.Qt.LeftButton)
    _settle()
    assert tab.document.get(name) == "/data/browsed" and tab.editors[name].text() == "/data/browsed"
    tab.close()


def test_a_derived_path_is_shown_as_a_placeholder_never_as_text():
    """V9, P2: the queryable property is text() == "" with placeholderText() set to the derived path. A derived
    value drawn as typed text would be an input the document does not hold."""
    tab = SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "IPTS-1"}))
    for name in _PATH_FIELDS:
        editor = tab.editors[name]
        assert editor.text() == "" and editor.placeholderText() == f"/SNS/REF_L/IPTS-1/{_PATH_TAIL[name]}"


def test_clearing_the_ipts_stores_an_empty_name_never_none(tmp_path, monkeypatch):
    """V10, P5, F7: "" and a str, not None (with None the path properties raise TypeError); both paths say no IPTS is
    set; the direct-beam listing answers ([], 0); a save writes "experiment_id": ""."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "IPTS-1"})))
    _type_into(tab.editors["experiment_id"], "")
    value = tab.document.get("experiment_id")
    assert value == "" and type(value) is str
    assert tab._last_error is None  # V10' (v2): a swallowed TypeError would leave the value check passing
    for name in _PATH_FIELDS:
        assert tab.editors[name].placeholderText() == _NO_IPTS
    assert tab.document.candidates("DBname") == ([], 0)
    target = tmp_path / "out.json"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *_a, **_k: (str(target), "")))
    tab.save_settings()
    assert json.loads(target.read_text())["experiment_id"] == ""
    tab.close()


def test_the_direct_beam_list_follows_the_path_typed_in_the_header(tmp_path):
    """V11, F8: the cell's listing reads the effective DBpath (override if set, else derived), so a folder typed in
    the header is what the next cell open lists, and clearing it lists the derived folder again (absent here)."""
    folder = tmp_path / "db"
    folder.mkdir()
    (folder / "a.txt").write_text("")
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "IPTS-0"})))
    _type_into(tab.editors["_DBpath_override"], str(folder))
    assert tab.document.candidates("DBname") == (["a.txt"], 1)
    _type_into(tab.editors["_DBpath_override"], "")
    assert tab.document.candidates("DBname") == ([], 0)
    tab.close()


# Held state of a path control (plan §3): D an unset override with an IPTS, N an unset override with no IPTS, S a string
# override, X a non-string override from a malformed file.
_PATH_STATES = {"D": None, "N": None, "S": "/data/held", "X": 5}  # the override each state holds


def _path_settings(state, name):
    settings = {"experiment_id": "" if state == "N" else "IPTS-1"}
    if _PATH_STATES[state] is not None:
        settings[name] = _PATH_STATES[state]
    return settings

_PATH_OPERATIONS = ["display", "type-ipts", "clear-ipts", "type-path", "type-whitespace", "clear-path",
                    "focus-through", "browse", "browse-derived", "browse-cancelled", "save", "load-second"]


def _panel_changed(tab):
    """The fields the panel's "Changed from the seed" section names: tab.report itself, not the model."""
    text = tab.report.toPlainText()
    if "Changed from the seed:" not in text:
        return []
    section = text.split("Changed from the seed:", 1)[1]
    return [line.strip()[2:].split(":", 1)[0] for line in section.splitlines() if line.strip().startswith("- ")]


def _panel_reports(tab, name):
    """Whether the panel's problems name the field (a problem line carries "(<field>)")."""
    return f"({name})" in tab.report.toPlainText().split("Changed from the seed:", 1)[0]


def _expected_path_cell(state, operation, name):
    """What §3's table requires after `operation` in `state`: (the value held, the control's text, its placeholder,
    the fields changed from the seed)."""
    held = _PATH_STATES[state]
    ipts = "" if state == "N" else "IPTS-1"

    def derived(ipts):
        return f"/SNS/REF_L/{ipts}/{_PATH_TAIL[name]}" if ipts else _NO_IPTS

    text = "" if held is None else str(held)
    to_derived = [] if held is None else [name]  # S and X change when they return to derived; D and N do not
    if operation in ("display", "focus-through", "browse-cancelled", "save"):
        return held, text, derived(ipts), []
    if operation == "type-ipts":
        return held, text, derived("IPTS-36119"), ["experiment_id"]
    if operation == "clear-ipts":
        return held, text, _NO_IPTS, [] if state == "N" else ["experiment_id"]
    if operation == "type-path":
        return "/data/typed", "/data/typed", derived(ipts), [name]
    if operation in ("type-whitespace", "clear-path"):
        return None, "", derived(ipts), to_derived
    if operation == "browse":
        return "/data/browsed", "/data/browsed", derived(ipts), [name]
    if operation == "browse-derived":
        if state == "N":  # no derived folder exists: what the dialog returns is an override (P7)
            return "/data/browsed", "/data/browsed", _NO_IPTS, [name]
        return None, "", derived(ipts), to_derived
    return None, "", derived("IPTS-2"), []  # load-second: the second file's state only


@pytest.mark.parametrize("operation", _PATH_OPERATIONS)
@pytest.mark.parametrize("state", list(_PATH_STATES))
@pytest.mark.parametrize("name", _PATH_FIELDS)
def test_each_operation_on_a_path_control_in_each_held_state(name, state, operation, tmp_path, monkeypatch):
    """V12: plan §3's operation × state table, every cell, for both path controls. Each asserts:
    - the document value and its type;
    - the control's text and its placeholder;
    - the fields changed from the seed, both in the model and in the panel's "Changed from the seed" section
      (tab.report);
    - whether the panel's problems name the field (only a non-string override is reported);
    - in a cell that changes nothing, a panel text equal to the text before.
    v1's docstring said it checked "the report"; it read only the model."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict(_path_settings(state, name))))
    editor = tab.editors[name]
    panel_before = tab.report.toPlainText()
    if operation == "type-ipts":
        _type_into(tab.editors["experiment_id"], "36119")
    elif operation == "clear-ipts":
        _type_into(tab.editors["experiment_id"], "")
    elif operation == "type-path":
        _type_into(editor, "/data/typed")
    elif operation == "type-whitespace":
        _type_into(editor, "   ")
    elif operation == "clear-path":
        _type_into(editor, "")
    elif operation == "focus-through":
        _focus_through(tab, editor)
    elif operation in ("browse", "browse-derived", "browse-cancelled"):
        chosen = {"browse": "/data/browsed", "browse-cancelled": ""}.get(
            operation, "/data/browsed" if state == "N" else f"/SNS/REF_L/IPTS-1/{_PATH_TAIL[name]}")
        monkeypatch.setattr(QtWidgets.QFileDialog, "getExistingDirectory", staticmethod(lambda *_a, **_k: chosen))
        QTest.mouseClick(tab.path_browse[name], QtCore.Qt.LeftButton)
        _settle()
    elif operation == "save":
        target = tmp_path / "out.json"
        monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName",
                            staticmethod(lambda *_a, **_k: (str(target), "")))
        tab.save_settings()
        saved = json.loads(target.read_text())[name]
        held = tab.document.get(name)
        assert saved == held and type(saved) is type(held)
    elif operation == "load-second":
        tab.set_document(SettingsDocument.from_dict({"experiment_id": "IPTS-2"}))
    value, text, placeholder, changed = _expected_path_cell(state, operation, name)
    held = tab.document.get(name)
    assert held == value and type(held) is type(value)
    assert editor.text() == text and editor.placeholderText() == placeholder
    assert sorted(tab.document.changed_vs_seed()) == sorted(changed)
    assert sorted(_panel_changed(tab)) == sorted(changed)
    assert _panel_reports(tab, name) is (type(held) not in (str, type(None)))
    if not changed and operation != "load-second":
        assert tab.report.toPlainText() == panel_before
    tab.close()


@pytest.mark.parametrize("slot", ["type-ipts", "type-path", "clear-path", "browse"])
def test_the_panel_follows_each_header_write_slot(slot, monkeypatch):
    """V13 (v2, B-1): each of the header's three write slots refreshes the report panel itself. Starting from a
    malformed direct-beam override (5, reported), the panel's "Changed from the seed" names the field written.
    After a path write (typed, cleared, browsed) the panel's problem line for it is gone. Removing
    refresh_report() from any one slot fails its legs."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "IPTS-1", "_DBpath_override": 5})))
    assert _panel_reports(tab, "_DBpath_override")
    field = "_DBpath_override"
    if slot == "type-ipts":
        _type_into(tab.editors["experiment_id"], "36119")
        field = "experiment_id"
    elif slot == "type-path":
        _type_into(tab.editors["_DBpath_override"], "/data/typed")
    elif slot == "clear-path":
        _type_into(tab.editors["_DBpath_override"], "")
    else:
        monkeypatch.setattr(QtWidgets.QFileDialog, "getExistingDirectory",
                            staticmethod(lambda *_a, **_k: "/data/browsed"))
        QTest.mouseClick(tab.path_browse["_DBpath_override"], QtCore.Qt.LeftButton)
        _settle()
    assert field in _panel_changed(tab)
    if slot != "type-ipts":
        assert not _panel_reports(tab, "_DBpath_override")
    tab.close()


def test_the_integrators_reproduction_leaves_no_problem_and_both_changes_on_the_panel():
    """V13 (v2): the rejection's reproduction of B-1, verbatim."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "IPTS-1", "_DBpath_override": 5})))
    _type_into(tab.editors["_DBpath_override"], "/data/typed")
    _type_into(tab.editors["experiment_id"], "36119")
    panel = tab.report.toPlainText()
    assert "No problems found." in panel
    assert sorted(_panel_changed(tab)) == ["_DBpath_override", "experiment_id"]
    tab.close()


@pytest.mark.parametrize("state", ["D", "S", "X"])
@pytest.mark.parametrize("name", _PATH_FIELDS)
def test_browse_never_writes_the_derived_folder(name, state, monkeypatch):
    """V14 (v2, U-1, P7), through the Browse button on a shown tab.
    - From D the dialog opens at the derived folder, and Choose without navigating returns it: nothing is written,
      with the value, text, placeholder, "Changed from the seed" and the panel as before.
    - From S and X the dialog opens at the held override; choosing the derived folder holds no override, as
      clearing does. The field is changed, the placeholder is the derived folder, and X's problem line is gone."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict(_path_settings(state, name))))
    derived = f"/SNS/REF_L/IPTS-1/{_PATH_TAIL[name]}"
    opened = []

    def dialog(_parent, _caption, start, *_args, **_kwargs):
        opened.append(start)
        return start if state == "D" else derived

    monkeypatch.setattr(QtWidgets.QFileDialog, "getExistingDirectory", staticmethod(dialog))
    panel = tab.report.toPlainText()
    QTest.mouseClick(tab.path_browse[name], QtCore.Qt.LeftButton)
    _settle()
    editor = tab.editors[name]
    assert opened == [derived if state == "D" else str(_PATH_STATES[state])]
    assert tab.document.get(name) is None
    assert editor.text() == "" and editor.placeholderText() == derived
    if state == "D":
        assert tab.document.changed_vs_seed() == {} and tab.report.toPlainText() == panel
    else:
        assert list(tab.document.changed_vs_seed()) == [name] and _panel_changed(tab) == [name]
        assert not _panel_reports(tab, name)
    tab.close()


@pytest.mark.parametrize("tail, written", [("/", False), ("/.", False), ("/../other", True)],
                         ids=["trailing-separator", "dot-segment", "sibling-folder"])
@pytest.mark.parametrize("name", _PATH_FIELDS)
def test_browse_compares_the_chosen_folder_as_a_path(name, tail, written, monkeypatch):
    """V14 (v2, P7): the derived folder is recognised as a path (os.path.normpath on both sides), so a trailing
    separator or a "." segment is still the derived folder and writes nothing. A sibling folder is another folder,
    and is written as typed by the dialog."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "IPTS-1"})))
    chosen = f"/SNS/REF_L/IPTS-1/{_PATH_TAIL[name]}{tail}"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getExistingDirectory", staticmethod(lambda *_a, **_k: chosen))
    QTest.mouseClick(tab.path_browse[name], QtCore.Qt.LeftButton)
    _settle()
    assert tab.document.get(name) == (chosen if written else None)
    tab.close()


@pytest.mark.parametrize("name", _PATH_FIELDS)
def test_a_typed_path_is_stored_stripped_and_whitespace_alone_is_derived(name):
    """V15 (v2, P8): "   " is no override (None, the derived placeholder shows); "/some/dir " is "/some/dir"."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "IPTS-1"})))
    editor = tab.editors[name]
    _type_into(editor, "   ")
    assert tab.document.get(name) is None
    assert editor.text() == "" and editor.placeholderText() == f"/SNS/REF_L/IPTS-1/{_PATH_TAIL[name]}"
    _type_into(editor, "/some/dir ")
    assert tab.document.get(name) == "/some/dir" and editor.text() == "/some/dir"
    tab.close()


@pytest.mark.parametrize("name", _PATH_FIELDS)
def test_a_cleared_path_is_saved_as_null(name, tmp_path, monkeypatch):
    """V5' (v2, test advisory A1): clearing a held path, then Save, writes null for it."""
    tab = _shown_tab(SettingsEditorTab(SettingsDocument.from_dict({"experiment_id": "IPTS-1", name: "/data/held"})))
    _type_into(tab.editors[name], "")
    target = tmp_path / "out.json"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *_a, **_k: (str(target), "")))
    tab.save_settings()
    assert json.loads(target.read_text())[name] is None
    tab.close()


@pytest.mark.parametrize("name", _PATH_FIELDS)
def test_a_second_load_with_an_override_shows_it_as_text(name):
    """V7' (v2, test advisory A2): file A has no override, file B one; the control shows B's text as a value."""
    tab = SettingsEditorTab()
    tab.set_document(SettingsDocument.from_dict({"experiment_id": "IPTS-1"}))
    tab.set_document(SettingsDocument.from_dict({"experiment_id": "IPTS-2", name: "/data/b"}))
    assert tab.editors[name].text() == "/data/b"
    assert tab.editors["experiment_id"].text() == "IPTS-2"
