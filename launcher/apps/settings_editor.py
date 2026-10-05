#!/usr/bin/python3
"""Settings-editor tab: author a reduction settings JSON with guidance.

A thin view over :class:`~lr_reduction.settings_document.SettingsDocument`.
Every widget here is built from :mod:`lr_reduction.field_spec`, so the fields
the editor offers, the prompts it shows and the values it accepts all come from
one table rather than from hand-written widget code that drifts from the config
class.

Supersedes the never-existent ``JSONSettingsBuilderTab`` that
``new_launcher.py`` carried as a commented-out import.

**The row-index rule.** Every per-angle write goes through
``SettingsDocument.set_angle_field(row, name, value)`` with the row Qt reports
for the edited cell. Nothing here consults ``currentRow()`` to decide *what* to
edit — the selection is used only to choose which row the Remove button
deletes, where it is the actual input rather than a hidden one. Reading config
from the selected row instead of the acted-on row is a known reduction-GUI bug
class, and this table is new code, so the trap would be introduced here.

**Drop-downs.** Every enumerated field is a drop-down, and only a choice in its
open list changes it. The mouse wheel never does, and the keys that would step a
closed one open its list instead (``NoWheelComboBox``). After a choice, the focus
goes to the drop-down's container. In the Angles table they are item delegates
rather than per-cell widgets: the cells stay text items and show the arrow
painted at rest. One click, or Enter, Space or Alt+Down on the current cell,
opens a cell's list, and Down still moves down the grid (APG grid pattern). A
choice is written into the item, so ``_on_cell_changed`` remains the one write
path. The choices read as the file spells them (``_in_file_spelling``).
"""

import contextlib
import functools
import os
import traceback
from pathlib import Path

from qtpy import QtCore, QtGui, QtWidgets

from launcher.app_identity import ensure_identity
from launcher.apps import roi_dialog
from launcher.apps.roi_dialog import ROISelectionDialog
from lr_reduction import field_spec as fs
from lr_reduction import roi_estimate
from lr_reduction.settings_document import (
    SettingsDocument,
    file_spelling,
    load_start_folder,
    normalise_experiment_id,
    settings_folders,
)

#: Above this, populating the table freezes the GUI thread for seconds and
#: costs ~1100x the file size in memory. A settings file with more angles than
#: this is a mistake, not a workload.
MAX_TABLE_ROWS = 500

#: Events the ROI pop-out reads from a run (roi-popout-data's stride sampling over the whole run, never its first
#: N): a choice of pixel ranges needs a sample, and a long run then opens quickly. #197's value (A7).
MAX_ROI_EVENTS = 2_000_000


#: Editor property set when the user chooses an item in a table drop-down.
_CHOSEN = "chosen"

#: Editor property: the text the cell showed when it opened, held or implied. A choice equal to it is the
#: identity (C9').
_SHOWN = "shown"

#: Keys that move an open list's current item: a deliberate move (C11).
_MOVE_KEYS = {
    QtCore.Qt.Key_Up, QtCore.Qt.Key_Down, QtCore.Qt.Key_PageUp, QtCore.Qt.Key_PageDown,
    QtCore.Qt.Key_Home, QtCore.Qt.Key_End,
}

#: The keys that open a closed drop-down instead of stepping its value (C10): with the list shut, the
#: value changes only by a choice in the open list (APG select-only combobox).
_LIST_KEYS = {
    QtCore.Qt.Key_Up, QtCore.Qt.Key_Down, QtCore.Qt.Key_PageUp, QtCore.Qt.Key_PageDown,
    QtCore.Qt.Key_Home, QtCore.Qt.Key_End, QtCore.Qt.Key_F4,
}

def section_state_key(name):
    """The QSettings key a list section's state is stored under: the section's declared name, never its position."""
    return f"settings_editor/sections/{name}"


def _stored_expanded(value):
    """A stored section state, read by its meaning. A new process reads the strings "true"/"false" back from INI
    (measured, Qt 5.15), and bool("false") is True. Anything else, or nothing stored, is expanded."""
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().lower() in ("true", "false"):
        return value.strip().lower() == "true"
    return True


class _SectionHeading(QtWidgets.QToolButton):
    """A list section's heading: a checkable button, checked while the section is expanded. A click or Space toggles
    it, as any button does; Return and Enter do too, for a heading reached with Tab."""

    def keyPressEvent(self, event):
        modifiers = event.modifiers() & ~QtCore.Qt.KeypadModifier
        if event.key() in (QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter) and modifiers == QtCore.Qt.NoModifier:
            self.click()
            return
        super().keyPressEvent(event)


class _Section(QtWidgets.QWidget):
    """One section of the editor's list: a heading that collapses and expands the fields under it.

    Collapsed, the fields are hidden: they take no space, and Tab passes from the heading to the next section's.
    Nothing else changes. They are still the tab's editors, so a Load refreshes them, and their values are still in
    the document, in validate() and in a save.
    """

    def __init__(self, title, expanded, parent=None):
        super().__init__(parent)
        self.heading = _SectionHeading()
        self.heading.setText(title)
        self.heading.setCheckable(True)
        self.heading.setChecked(expanded)
        self.heading.setAutoRaise(True)
        self.heading.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.heading.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        self.heading.setToolTip("Collapse or expand this section. The state is remembered for you.")
        font = self.heading.font()
        font.setBold(True)
        self.heading.setFont(font)
        self.body = QtWidgets.QWidget()
        self.form = QtWidgets.QFormLayout()
        self.body.setLayout(self.form)
        layout = QtWidgets.QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.heading)
        layout.addWidget(self.body)
        self.setLayout(layout)
        self.heading.toggled.connect(self._show_body)
        self._show_body(expanded)

    def _show_body(self, expanded):
        self.heading.setArrowType(QtCore.Qt.DownArrow if expanded else QtCore.Qt.RightArrow)
        self.body.setVisible(expanded)


#: The header's path overrides (fs.HEADER_NAMES without the IPTS), and what a path shows when nothing is derived.
_HEADER_PATHS = tuple(name for name in fs.HEADER_NAMES if name != "experiment_id")
_NO_IPTS = "set an IPTS or type a path"
_NOT_A_FOLDER = "the IPTS is not a folder name; type a path"


def _later(owner, slot):
    """Run `slot` once on the next pass of the event loop, on a timer owned by `owner`.

    The timer dies with its owner, so a slot meant for an editor that has been closed in the meantime never
    runs on a deleted widget. PyQt5 sends an exception from a timer slot to qFatal.
    """
    timer = QtCore.QTimer(owner)
    timer.setSingleShot(True)
    timer.timeout.connect(slot)
    timer.timeout.connect(timer.deleteLater)
    timer.start(0)


def _in_file_spelling(field, held):
    """The choices of an enumerated per-angle field, in the column's own spelling (C9).

    The reducer lower-cases ``method_per_run`` (``nr_reduction_calc.py:82``), and its own files hold
    'meantheta'. When every case variant of a choice that the column holds shares one spelling (all lower
    case, or all upper case), the choices are offered in that spelling. The list then reads as the file
    does, and a new choice is written in the file's convention. With no variant (a fresh document, or
    declared spellings only) or mixed casing, the declared spellings are offered. A held value is always
    one of the items: ``_ChoiceDelegate.setEditorData`` adds one the list lacks.
    """
    declared = [str(choice) for choice in field.allowed]
    entries = held if isinstance(held, (list, tuple)) else ()
    variants = {entry for entry in entries if isinstance(entry, str) and field.canonical(entry) in declared}
    if variants and not all(variant == field.canonical(variant) for variant in variants):
        for spell in (str.lower, str.upper):
            if all(variant == spell(variant) for variant in variants):
                return [spell(choice) for choice in declared]
    return declared


class NoWheelComboBox(QtWidgets.QComboBox):
    """A drop-down that only a choice in its open list changes (the scientists' item 2; C10).

    ``QComboBox`` steps through its items on a wheel event, focused or not. So
    scrolling the settings list with the pointer crossing one changed a setting,
    and the document with it (measured: ``DetResFn`` rectangular -> gaussian).
    The event is ignored instead, which leaves it to the parent. Inferred, not
    measured: Qt passes a real (spontaneous) wheel event on to the parent, so
    the list scrolls (``QApplication::notify``). A test cannot send a
    spontaneous event, which is why the tests assert the event left unaccepted
    instead. Also inferred: an open pop-up list still scrolls, since the wheel
    then reaches its list view and not the combo.

    The keys that step a closed ``QComboBox`` open its list here instead
    (APG's select-only combobox: "the only way users can set its value is by
    selecting a value in the popup"). That covers the arrows, Page, Home, End
    and F4. On a non-editable one it also covers Space, Return, Enter and a
    typed letter. An editable one's typing still reaches its line edit. The
    wheel is ignored even with focus: after a click the combo keeps focus, and
    the reported failure was a scroll that followed a click (plan A1).
    ``StrongFocus`` means the wheel does not take focus either.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)

    def wheelEvent(self, event):
        event.ignore()

    def keyPressEvent(self, event):
        key = event.key()
        opens = key in _LIST_KEYS
        if not self.isEditable():
            opens = opens or key in (QtCore.Qt.Key_Space, QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter)
            opens = opens or (bool(event.text()) and event.text().isprintable())
        if opens:
            self.showPopup()
            event.accept()
            return
        super().keyPressEvent(event)


class _DropDownDelegate(QtWidgets.QStyledItemDelegate):
    """What every table drop-down shares (C8): it shows its arrow at rest and opens on one click.

    The arrow is drawn beside the cell's value whether or not the cell is open, so a drop-down looks like one.
    A left click opens the cell's editor, which opens its list. That is the menu-button pattern the human
    cites, without the double-click Qt's default edit trigger needs. The editor opens on the next pass of the
    event loop, on a persistent index, so a click that also changes the selection lands first. Nothing is
    built per cell at rest: the arrow is painted, and only visible cells are painted.
    """

    def __init__(self, tab, field):
        super().__init__(tab.angle_table)
        self._tab = tab
        self._field = field

    def paint(self, painter, option, index):
        widget = option.widget
        style = widget.style() if widget is not None else QtWidgets.QApplication.style()
        side = max(8, min(option.rect.height() - 4, 14))
        background = QtWidgets.QStyleOptionViewItem(option)
        self.initStyleOption(background, index)
        background.text = ""
        style.drawPrimitive(QtWidgets.QStyle.PE_PanelItemViewItem, background, painter, widget)
        content = QtWidgets.QStyleOptionViewItem(option)
        content.rect = option.rect.adjusted(0, 0, -(side + 6), 0)
        super().paint(painter, content, index)
        arrow = QtWidgets.QStyleOption()
        arrow.rect = QtCore.QRect(option.rect.right() - side - 3, option.rect.center().y() - side // 2, side, side)
        arrow.palette = option.palette
        arrow.state = QtWidgets.QStyle.State_Enabled
        style.drawPrimitive(QtWidgets.QStyle.PE_IndicatorArrowDown, arrow, painter, widget)

    def editorEvent(self, event, model, option, index):
        if event.type() == QtCore.QEvent.MouseButtonRelease and event.button() == QtCore.Qt.LeftButton:
            view = self.parent()
            cell = QtCore.QPersistentModelIndex(index)
            _later(view, lambda: view.edit(QtCore.QModelIndex(cell)) if cell.isValid() else None)
            return True
        return super().editorEvent(event, model, option, index)

    def _choose(self, editor):
        """A choice was made in the editor's list: write it unless it is what the cell showed, close the editor,
        and let go of focus (C10).

        Choosing the value a cell shows is the identity (C9'), an implied one included: in a compact column the
        reduction already uses that value at every angle, so writing it would only turn the list explicit.
        """
        gestures = getattr(editor, "gestures", None)
        if gestures is not None and gestures.deliberate and editor.currentText() != editor.property(_SHOWN):
            editor.setProperty(_CHOSEN, True)
            self.commitData.emit(editor)
        self._leave(editor)

    def _leave(self, editor):
        """Close the editor; the view takes the focus back."""
        self.closeEditor.emit(editor, QtWidgets.QAbstractItemDelegate.NoHint)


class _ChoiceDelegate(_DropDownDelegate):
    """The drop-down of an enumerated per-angle column: its choices, plus "" for unset.

    It opens with its list shown. A choice in the list writes it into the cell's
    item, so ``_on_cell_changed`` and its row-from-the-signal rule stay the one
    write path (C6). Then the editor closes and the focus goes back to the table
    (C10). Closing the list without a choice (Escape, a click elsewhere) writes
    nothing and leaves the closed drop-down, whose keys reopen the list;
    another Escape returns to the grid. The choices are in the column's own spelling
    (``_in_file_spelling``, C9), and the held value is always one of them, so
    choosing the value a cell holds changes no data and writes nothing.

    What it shows:
    * a held value as it is held: the file's spelling, and an out-of-domain
      value as itself, never the first item (C2);
    * an empty cell as the value the reduction uses there, in italics, when the
      document implies one (``SettingsDocument.implied_entry``, C7).

    The implied value is asked of the document when a cell is drawn or opened,
    not stored on the items. Only visible cells are drawn, so filling a 500-row
    table does no implied work. The value also cannot go stale when another
    edit moves the angle count; ``refresh_marks`` repaints after every change.
    """

    def _choices(self):
        if self._field.element_type == "bool":
            return ["true", "false"]
        return _in_file_spelling(self._field, self._tab.document.get(self._field.name))

    def _implied(self, index):
        """The text of the value the reduction uses in this cell when it holds none, or ``""``."""
        value = self._tab.document.implied_entry(index.row(), self._field.name)
        return "" if value is None else self._tab._cell_text(self._field, value)

    def createEditor(self, parent, _option, _index):
        editor = NoWheelComboBox(parent)
        editor.addItems(["", *self._choices()])
        editor.setProperty(_CHOSEN, False)
        editor.gestures = _ListKeys(editor, self._leave)
        editor.activated.connect(lambda _index, editor=editor: self._choose(editor))
        _later(editor, editor.gestures.open_list)
        return editor

    def setEditorData(self, editor, index):
        shown = index.data(QtCore.Qt.DisplayRole) or self._implied(index)
        at = editor.findText(shown)
        if at < 0:
            editor.addItem(shown)
            at = editor.count() - 1
        editor.setCurrentIndex(at)
        editor.setProperty(_SHOWN, shown)

    def setModelData(self, editor, model, index):
        if editor.property(_CHOSEN):
            model.setData(index, editor.currentText(), QtCore.Qt.EditRole)

    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        if index.data(QtCore.Qt.DisplayRole):
            return
        implied = self._implied(index)
        if implied:
            option.text = implied
            font = QtGui.QFont(option.font)
            font.setItalic(True)
            option.font = font
            palette = QtGui.QPalette(option.palette)
            palette.setColor(QtGui.QPalette.Text, palette.color(QtGui.QPalette.PlaceholderText))
            option.palette = palette

    def helpEvent(self, event, view, option, index):
        """An implied cell's tooltip says where its value comes from."""
        implied = self._implied(index) if event.type() == QtCore.QEvent.ToolTip else ""
        if implied and not index.data(QtCore.Qt.DisplayRole):
            QtWidgets.QToolTip.showText(
                event.globalPos(),
                f"Not set here: the reduction uses {implied} for this angle. Choose a value to set it.",
                view,
            )
            return True
        return super().helpEvent(event, view, option, index)


class _ListKeys(QtCore.QObject):
    """The gestures in a table drop-down's open list: only a deliberate choice writes (C11).

    Qt makes row 0 current when the list takes the focus although the user moved
    nowhere, and the held value may not be listed at all (a direct-beam name
    outside the folder, the reducer-written norm). Return there chose that row
    and wrote the folder's first file (review d3ee364, U-1). QComboBox's list
    chooses its current row on the ShortcutOverride event that comes before the
    key press (Qt 5.15, traced), so a key-press filter cannot take Return over.
    Instead:

    * the list opens on the combo's own current item, or on none when the held
      value is not listed (``open_list``). Return with no move then chooses the
      value the cell shows, which is the identity (C9'), or no row at all;
    * ``deliberate`` records a deliberate move: an arrow key, Page, Home or End,
      type-ahead in a list that is not editable, or a mouse press on an item. A
      hover moves the current row too, and Return chooses the hovered row, so
      ``_choose`` writes nothing without ``deliberate``;
    * a Return that reaches this filter with the list still open found no row to
      choose: it closes the list and leaves the cell as held. After a choice the
      cell is already closed, and closing it again does nothing.

    In an editable list, a typed character closes the list and starts a new name
    in the line edit, since an open list takes the keystrokes (measured). Owned by
    the combo, so it goes when the editor goes.
    """

    def __init__(self, combo, leave):
        super().__init__(combo)
        self._combo = combo
        self._leave = leave
        self.deliberate = False
        combo.view().installEventFilter(self)
        combo.view().viewport().installEventFilter(self)

    def open_list(self):
        """Show the list on the combo's own current item, or on none when the held value is not listed."""
        combo = self._combo
        combo.showPopup()
        combo.view().setCurrentIndex(combo.model().index(combo.currentIndex(), combo.modelColumn()))

    def eventFilter(self, watched, event):
        if event.type() == QtCore.QEvent.MouseButtonPress and watched is self._combo.view().viewport():
            self.deliberate = True  # a press on an item; its release is the choice (QComboBox)
            return False
        if event.type() != QtCore.QEvent.KeyPress:
            return False
        key = event.key()
        if key in _MOVE_KEYS:
            self.deliberate = True
            return False
        if key in (QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter):
            self._combo.hidePopup()
            self._leave(self._combo)
            return True
        if event.text() and event.text().isprintable():
            if not self._combo.isEditable():
                self.deliberate = True  # type-ahead moves to an item
                return False
            if self._combo.view().isVisible():
                # The character that closes the list starts a new name, as typing into the selected text
                # the cell opened with would; opening the list drops that selection (measured: appended).
                self._combo.hidePopup()
                self._combo.lineEdit().selectAll()
            QtWidgets.QApplication.sendEvent(
                self._combo.lineEdit(),
                QtGui.QKeyEvent(QtCore.QEvent.KeyPress, event.key(), event.modifiers(), event.text()),
            )
            return True
        return False


class _CandidatesDelegate(_DropDownDelegate):
    """An editable drop-down of the file names in the folder a column's field declares (``DBname``, C4).

    It opens with the folder's names listed. A character typed then closes the
    list and goes to the name (``_ListKeys``). Typing completes
    inline from the names. A name picked from the list is written, the editor
    closes, and the focus goes back to the table (C10). A typed name that is not
    in the folder is stored as typed on Return.

    The completer stays inline on purpose. A pop-up completer is a parentless
    top-level window that QCompleter owns through a raw pointer
    (``qcompleter.cpp``, Qt 5.15). Anything that deletes top-level windows then
    frees it twice: the launcher tests' teardown did, and aborted
    (``QCompleter::~QCompleter``, gdb).

    The names come from ``SettingsDocument.candidates``, asked for each time a
    cell opens, so they follow ``experiment_id`` and the direct-beam path
    wherever those change (C5). The folder is listed only when a cell asks:
    never per keystroke, never in a refresh. When the folder holds more names
    than the cap, the tooltip says so.
    """

    def createEditor(self, parent, _option, _index):
        editor = NoWheelComboBox(parent)
        editor.setEditable(True)
        editor.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        editor.completer().setCompletionMode(QtWidgets.QCompleter.InlineCompletion)
        editor.completer().setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        names, total = self._tab.document.candidates(self._field.name)
        editor.addItems(names)
        if total > len(names):
            editor.setToolTip(
                f"Showing the first {len(names)} of {total} files in the folder; type a name to use another."
            )
        editor.setProperty(_CHOSEN, False)
        editor.gestures = _ListKeys(editor, self._leave)
        editor.activated.connect(lambda _index, editor=editor: self._choose(editor))
        _later(editor, editor.gestures.open_list)
        return editor

    def setEditorData(self, editor, index):
        # Select the held name when the folder has it. Opening the list takes the focus from the line edit,
        # and an editable QComboBox then reads a name that matches an item other than its current one as a
        # new choice: it emits activated, which closed the cell the moment it opened (traced).
        text = index.data(QtCore.Qt.DisplayRole) or ""
        editor.setCurrentIndex(editor.findText(text))
        editor.setEditText(text)
        editor.lineEdit().selectAll()
        editor.setProperty(_SHOWN, text)

    def setModelData(self, editor, model, index):
        # A name picked deliberately from the list, or typed by the user; never a text Qt put in the line edit
        # (a row the list made current on its own, applied when it closed).
        if editor.property(_CHOSEN) or editor.lineEdit().isModified():
            model.setData(index, editor.currentText(), QtCore.Qt.EditRole)


class _FileDialogSidebar(QtCore.QObject):
    """Gives the next file dialog shown the IPTS's settings folders as its sidebar (editor-ipts-inference, I5).

    The static ``QFileDialog`` calls take no sidebar, and the tests' autouse net (``conftest.no_qfiledialog``)
    stubs exactly those calls so that no test can block on a modal dialog. So the dialogs stay static, and this
    filter, installed on the application for the one call, sets the sidebar of the dialog the call builds when it is
    shown. Measured offscreen on Qt 5.15 (the ledger's ``editor-ipts-inference-probes.py``):
    the filter sees the static dialog's Show event and the sidebar holds. Qt's own dialog only
    (``DontUseNativeDialog``): a native dialog builds no sidebar.

    When the dialog hides, it gets its own sidebar back. Qt saves a dialog's sidebar ("shortcuts") to the user's
    QtProject.conf when the dialog is destroyed, and every later Qt 5 file dialog of the user's starts from it:
    left in place, the IPTS's folders would replace the user's own sidebar in every Qt application.
    """

    def __init__(self, folders):
        super().__init__()
        self._urls = [QtCore.QUrl.fromLocalFile(folder) for folder in folders]
        self._own = None

    def eventFilter(self, watched, event):  # noqa: N802 -- Qt's name
        # Never raise here: an exception out of a PyQt virtual reaches qFatal() and aborts the launcher.
        try:
            if isinstance(watched, QtWidgets.QFileDialog):
                if event.type() == QtCore.QEvent.Show:
                    self._own = watched.sidebarUrls()
                    watched.setSidebarUrls(self._urls)
                elif event.type() == QtCore.QEvent.Hide and self._own is not None:
                    watched.setSidebarUrls(self._own)
                    self._own = None
        except Exception:  # noqa: BLE001
            traceback.print_exc()
        return False


@contextlib.contextmanager
def _file_dialog_sidebar(folders):
    """For one static file-dialog call: the dialog it shows gets ``folders`` as its sidebar (none: unchanged)."""
    application = QtWidgets.QApplication.instance()
    if not folders or application is None:
        yield
        return
    sidebar = _FileDialogSidebar(folders)
    application.installEventFilter(sidebar)
    try:
        yield
    finally:
        application.removeEventFilter(sidebar)


def guarded(method):
    """Report an exception into the panel instead of letting it leave the slot.

    An unhandled exception in a Qt slot under PyQt5 reaches ``qFatal()``, which
    calls ``abort()``: the whole launcher dies and every other tab loses its
    unsaved state. A settings editor reads files it did not write, so a
    malformed one must be a message, never a process death.
    """

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        try:
            return method(self, *args, **kwargs)
        except Exception as exc:  # noqa: BLE001 -- the point is to catch everything
            self.report_problem(exc)
            return None

    return wrapper


class SettingsEditorTab(QtWidgets.QWidget):
    """Editor for one :class:`SettingsDocument`."""

    def __init__(self, document=None, parent=None):
        # First, before any QSettings is constructed: QSettings derives its
        # path from the application identity, so a tab that builds one before
        # the identity is installed would resolve to a different store than the
        # rest of the launcher (S3's adoption contract).
        ensure_identity()
        super().__init__(parent)

        # A tab opened with no document is a new file: it starts at the instrument's current operation
        # (SettingsDocument.for_new_file). A document given to the tab is shown as it holds.
        self.document = document if document is not None else SettingsDocument.for_new_file()
        self.settings = QtCore.QSettings()
        self.editors = {}
        # Guards the table's cellChanged signal while the view writes into it,
        # so repopulating from the document does not echo back as user edits.
        self._populating = False
        self._rows_hidden = 0
        self._last_error = None

        layout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)
        layout.addWidget(self._build_toolbar())
        layout.addWidget(self._build_paths_header())

        splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        splitter.addWidget(self._build_angle_panel())
        splitter.addWidget(self._build_scalar_panel())
        splitter.addWidget(self._build_report_panel())
        layout.addWidget(splitter)

        # Render whatever the document already holds. __init__ used to call only
        # refresh_report(), so an injected document — the exact path a
        # resolution layer uses — displayed zero angle rows.
        self.set_document(self.document)

    # -- construction ------------------------------------------------------

    def _build_toolbar(self):
        bar = QtWidgets.QWidget()
        row = QtWidgets.QHBoxLayout()
        bar.setLayout(row)

        self.load_button = QtWidgets.QPushButton("Load settings...")
        self.load_button.setToolTip(
            "Seed from a settings JSON, or from the header of a pre-reduced .dat"
        )
        self.load_button.clicked.connect(lambda _checked=False: self.load_settings())
        row.addWidget(self.load_button)

        self.save_button = QtWidgets.QPushButton("Save settings...")
        self.save_button.clicked.connect(lambda _checked=False: self.save_settings())
        row.addWidget(self.save_button)

        row.addStretch(1)
        return bar

    def _build_paths_header(self):
        """IPTS and the two input paths it roots, above the angles: these fields' only editors (fs.HEADER_NAMES).

        A path override is written only by an explicit edit: typing in its control, or a folder from its Browse
        button. While it is unset, the control shows the folder the reduction derives from the IPTS as its
        placeholder (_show_derived_paths), so the derived path is visible but never held: a written override would
        freeze an absolute path into the file, where an unset one is derived again from the IPTS on load.
        """
        box = QtWidgets.QGroupBox("Experiment")
        self.paths_header = box
        form = QtWidgets.QFormLayout()
        box.setLayout(form)
        field = fs.get("experiment_id")
        ipts = QtWidgets.QLineEdit()
        ipts.setToolTip(f"{field.name} — {field.help} A number is stored as IPTS-<number>.")
        ipts.editingFinished.connect(lambda widget=ipts: self._on_ipts_edited(widget))
        self.editors[field.name] = ipts
        form.addRow(field.label, ipts)
        self.path_browse = {}
        for name in _HEADER_PATHS:
            field = fs.get(name)
            edit = QtWidgets.QLineEdit()
            edit.setToolTip(f"{field.name} — {field.help} Clear it to use the derived folder again.")
            edit.editingFinished.connect(lambda name=name, widget=edit: self._on_path_edited(name, widget))
            browse = QtWidgets.QPushButton("Browse...")
            browse.clicked.connect(lambda _checked=False, name=name: self._browse_path(name))
            row = QtWidgets.QHBoxLayout()
            row.addWidget(edit, 1)
            row.addWidget(browse)
            self.editors[name] = edit
            self.path_browse[name] = browse
            form.addRow(field.label, row)
        return box

    def _build_angle_panel(self):
        panel = QtWidgets.QGroupBox("Angles")
        box = QtWidgets.QVBoxLayout()
        panel.setLayout(box)

        self.angle_table = QtWidgets.QTableWidget(0, len(fs.PER_ANGLE_NAMES))
        self.angle_table.setHorizontalHeaderLabels(
            [fs.get(name).label for name in fs.PER_ANGLE_NAMES]
        )
        for column, name in enumerate(fs.PER_ANGLE_NAMES):
            self.angle_table.horizontalHeaderItem(column).setToolTip(fs.get(name).help)
        # Explicitly off, and asserted in the row-isolation test. The header is
        # clickable, and a single sortItems() would decouple visual row order
        # from document index — re-introducing the active-row bug this slug is
        # built to avoid, through the back door.
        self.angle_table.setSortingEnabled(False)
        # Drop-downs for the enumerated columns and the direct-beam names. Kept
        # here as well as on the table, so PyQt does not collect them.
        self._cell_delegates = {}
        for column, name in enumerate(fs.PER_ANGLE_NAMES):
            field = fs.get(name)
            if field.candidates_folder:
                delegate = _CandidatesDelegate(self, field)
            elif field.allowed or field.element_type == "bool":
                delegate = _ChoiceDelegate(self, field)
            else:
                continue
            self.angle_table.setItemDelegateForColumn(column, delegate)
            self._cell_delegates[name] = delegate
        self.angle_table.cellChanged.connect(self._on_cell_changed)
        self.angle_table.installEventFilter(self)
        box.addWidget(self.angle_table)

        buttons = QtWidgets.QHBoxLayout()
        self.add_angle_button = QtWidgets.QPushButton("Add angle")
        self.add_angle_button.clicked.connect(lambda _checked=False: self.add_angle())
        buttons.addWidget(self.add_angle_button)

        self.remove_angle_button = QtWidgets.QPushButton("Remove angle")
        self.remove_angle_button.clicked.connect(lambda _checked=False: self.remove_selected_angle())
        buttons.addWidget(self.remove_angle_button)

        # roi-popout-dialog (B1): the selected row's run on its detector images and profiles, its ROIs adjustable.
        self.select_roi_button = QtWidgets.QPushButton("Select ROI")
        self.select_roi_button.clicked.connect(lambda _checked=False: self.select_roi())
        buttons.addWidget(self.select_roi_button)
        self.angle_table.currentCellChanged.connect(lambda *_cells: self._update_select_roi_button())
        self._update_select_roi_button()
        buttons.addStretch(1)
        box.addLayout(buttons)
        return panel

    def _update_select_roi_button(self):
        """B1: enabled exactly while the Angles table has a current row, and the ROI plots can be drawn."""
        if roi_dialog.Figure is None:
            self.select_roi_button.setEnabled(False)
            self.select_roi_button.setToolTip("Unavailable: matplotlib's Qt backend could not be imported")
            return
        self.select_roi_button.setToolTip("Show the selected angle's run on the detector and adjust its ROIs")
        self.select_roi_button.setEnabled(self.angle_table.currentRow() >= 0)

    def _build_scalar_panel(self):
        scroll = QtWidgets.QScrollArea()
        # Where a scalar drop-down hands the focus after a choice (C10).
        self.scalar_panel = scroll
        scroll.setWidgetResizable(True)
        inner = QtWidgets.QWidget()
        column = QtWidgets.QVBoxLayout()
        inner.setLayout(column)

        # One collapsible section per declared group, in the scientists' order (fs.SECTION_ORDER), each opened as
        # this user left it. Keyed by the section's name, so a reordering cannot hand one section another's state.
        self.sections = {}
        for group in fs.SECTION_ORDER:
            # The header's fields have their editors there, and only there (fs.HEADER_NAMES).
            scalars = [f for f in fs.fields_in(group) if not f.per_angle and f.name not in fs.HEADER_NAMES]
            if not scalars:
                continue
            section = _Section(group, _stored_expanded(self.settings.value(section_state_key(group))))
            section.heading.toggled.connect(lambda expanded, group=group: self._on_section_toggled(group, expanded))
            for field in scalars:
                editor = self._build_editor(field)
                editor.setToolTip(f"{field.name} — {field.help}")
                self.editors[field.name] = editor
                section.form.addRow(field.label, editor)
            self.sections[group] = section
            column.addWidget(section)

        column.addStretch(1)
        scroll.setWidget(inner)
        return scroll

    def _build_editor(self, field):
        """One widget per field, chosen from the declared type and allowed set."""
        value = self.document.get(field.name)

        if field.type == "bool":
            editor = QtWidgets.QCheckBox()
            # A text-less QCheckBox responds to clicks only within its ~14 px
            # indicator (SE_CheckBoxClickRect), but a form layout will happily
            # stretch the widget to the column width. That leaves most of a
            # visibly-wide control inert, which reads as a broken checkbox.
            # Measured: stretched to 174 px, a click at the widget centre does
            # nothing. Fixing the size to the hint makes the clickable area and
            # the visible extent the same thing.
            editor.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
            editor.toggled.connect(
                lambda checked, name=field.name: self._set_scalar(name, bool(checked))
            )
            self._show(field, editor, value)
            return editor

        if field.allowed:
            editor = NoWheelComboBox()
            if field.choice_labels:
                # Entries in words (Field.choice_labels): each item carries the value it stores.
                for choice, text in field.choice_labels:
                    editor.addItem(text, choice)
            else:
                # A blank first entry for the tri-state fields, where a falsy value
                # means "off" and is the class default.
                if field.falsy_means_off:
                    editor.addItem("")
                editor.addItems([str(a) for a in field.allowed])
            # A choice in its list is the last thing it does: the focus goes to the panel, so a later arrow
            # key or wheel changes nothing (C10).
            editor.activated.connect(lambda _index: self.scalar_panel.setFocus(QtCore.Qt.OtherFocusReason))
            editor.currentIndexChanged.connect(lambda _index, name=field.name, editor=editor: self._on_choice(name, editor))
            self._show(field, editor, value)
            return editor

        editor = QtWidgets.QLineEdit()
        if field.runtime_owned:
            # The reduction's record of what it used (LambdaMinUse/LambdaMaxUse):
            # shown, never edited. It is deliberately not connected to
            # anything, so no signal, typed or programmatic, reaches the
            # document. Measured before this change: a bare editingFinished
            # turned the recorded 2.95 into [2.95], changing its shape with no
            # keystroke at all.
            editor.setReadOnly(True)
            self._show(field, editor, value)
            return editor
        if field.type in ("int", "float"):
            validator = (
                QtGui.QIntValidator() if field.type == "int" else QtGui.QDoubleValidator()
            )
            editor.setValidator(validator)
        # editingFinished, not textChanged: a partially typed number ("0.", "-")
        # is not a value to store, and writing on every keystroke would put the
        # document through states the user never asked for.
        editor.editingFinished.connect(
            lambda name=field.name, widget=editor: self._on_scalar_edited(name, widget)
        )
        self._show(field, editor, value)
        return editor

    @staticmethod
    def _show(field, editor, value):
        """Put `value` into `editor`. The ONLY place a value becomes widget state.

        Construction and refresh each used to implement this, and they
        disagreed: construction rendered a list through `_as_text`
        ("50, 200"), the refresh through `str()` ("[50, 200]"), and only the
        first survives being read back by `Field.coerce`. Since `__init__` now
        routes through `set_document` -> `refresh_scalars`, the divergent one
        ran on every tab open — so `data_x_range`, the first editor in tab
        order, was corrupted by a bare focus-out with no typing at all.

        Signals are blocked throughout: displaying a value is not an edit, and
        letting it echo back would rewrite the document from its own rendering.
        """
        was = editor.blockSignals(True)
        try:
            if isinstance(editor, QtWidgets.QCheckBox):
                editor.setChecked(bool(value))
            elif isinstance(editor, QtWidgets.QComboBox):
                SettingsEditorTab._show_in_combo(editor, field, value)
            else:
                editor.setText(SettingsEditorTab._as_text(value))
        finally:
            editor.blockSignals(was)

    def _on_choice(self, name, editor):
        """Store a choice made in a scalar combo, then show the document's value again once the signal has returned.

        The second step drops a raw entry the choice replaced (D3′, ``_show_in_combo``). It waits for the signal to
        return because removing items inside the combo's own signal changes the model under the handler, the known
        Qt trap. The write follows ``currentIndexChanged``: the position is what tells two entries with the same
        text apart (the string "True" beside the entry True). Qt 5.15 emits ``currentTextChanged`` on every index
        change of a non-editable combo as well (measured), so the two agree there. The position is the one that does
        not depend on that. Re-choosing the entry shown changes no position, so it writes nothing (C9′).
        """
        field = fs.get(name)
        self._set_scalar(name, self._chosen_value(field, editor, editor.currentText()))
        _later(editor, lambda: self._show(field, editor, self.document.get(name)))

    @staticmethod
    def _chosen_value(field, editor, text):
        """What a choice in an enumerated combo stores: the item's value where the entries are words
        (``Field.choice_labels``), else the text coerced to the field's type, "" being off."""
        if field.choice_labels:
            return editor.currentData()
        return field.coerce(text) if text else False

    @staticmethod
    def _show_in_combo(editor, field, value):
        """Display `value`, even when it is not one of the offered choices.

        A combo asked to show an unknown value silently displays its first item
        instead, so a file holding an out-of-set value would look like a valid
        one — and saving would then write the substituted value back. Adding the
        stray value as an entry keeps what is shown equal to what is held;
        validate() is what reports it as a problem.

        Where the entries are words, a value the reducer accepts shows its word
        (``Field.label_for``: ``True`` and any case of a name included), and any
        other value is an entry of its own.

        A raw entry lives exactly as long as the value it shows is held (D3′). Each
        refresh removes the raw entries and adds back one only for the value held
        now, so a Load, a new document, or a choice that replaced a raw value leaves
        the offered entries plus at most one raw entry: nothing in the list writes a
        value no loaded file holds. One rule for every scalar combo, labelled or not.
        Whether a value has an offered entry is decided by type and value
        (``_offered_index``), never by its text: the string ``"True"`` is not the
        entry ``True``.
        """
        at = SettingsEditorTab._offered_index(field, value)
        offered = SettingsEditorTab._offered_count(field)
        while editor.count() > offered:
            editor.removeItem(editor.count() - 1)
        if at is None:
            editor.addItem(str(value), value)
            at = editor.count() - 1
        editor.setCurrentIndex(at)

    @staticmethod
    def _offered_count(field):
        """The entries a scalar combo offers before any raw one: its words, or its choices after a blank first entry
        where falsy means off."""
        if field.choice_labels:
            return len(field.choice_labels)
        return len(field.allowed) + (1 if field.falsy_means_off else 0)

    @staticmethod
    def _offered_index(field, value):
        """The position of the offered entry that shows `value`: ``-1`` for none (an unset plain choice), ``None``
        when `value` needs a raw entry of its own. Matched by type and value."""
        if field.choice_labels:
            label = field.label_for(value)
            return None if label is None else [text for _, text in field.choice_labels].index(label)
        if value is None or value is False or value == "":
            return 0 if field.falsy_means_off else -1
        blank = 1 if field.falsy_means_off else 0
        for position, choice in enumerate(field.allowed):
            if type(choice) is type(value) and choice == value:
                return blank + position
        return None

    @staticmethod
    def _as_text(value):
        """Render a stored value for a single-line editor."""
        if value is None:
            return ""
        if isinstance(value, (list, tuple)):
            return ", ".join("" if v is None else str(v) for v in value)
        return str(value)

    @staticmethod
    def _row_header(row, surplus):
        """A row beyond the reduction's angle count says so, and says why.

        The row stays visible and editable: hiding it would hide a held value
        and let Add reuse its slot. Remove angle on it drops the surplus entries.
        """
        if not surplus:
            return QtWidgets.QTableWidgetItem(str(row + 1))
        item = QtWidgets.QTableWidgetItem(f"{row + 1} (surplus)")
        item.setToolTip(
            "Beyond the angles the reduction uses: it never reads these entries. "
            "Remove angle on this row drops them."
        )
        return item

    @staticmethod
    def _cell_text(field, value):
        """Render one Angles-table cell.

        A boolean column shows ``true``/``false``, the scientists' spelling and
        the one a two-item drop-down can take over, for ``True``/``1`` and
        ``False``/``0`` alike, because the reducer writes the integers. Anything
        else renders as :meth:`_as_text` does, so a stray ``2`` or ``"0"``
        stays visible as itself beside the problem ``validate()`` reports. Both
        spellings read back through ``Field.coerce_element``.
        """
        if field.element_type == "bool":
            boolean = fs.as_boolean(value)
            if boolean is not None:
                return "true" if boolean else "false"
        return SettingsEditorTab._as_text(value)

    def _build_report_panel(self):
        panel = QtWidgets.QGroupBox("Validation and changes")
        box = QtWidgets.QVBoxLayout()
        panel.setLayout(box)
        self.report = QtWidgets.QPlainTextEdit()
        self.report.setReadOnly(True)
        box.addWidget(self.report)
        return panel

    # -- editing -----------------------------------------------------------

    def _set_scalar(self, name, value):
        """Store a scalar and refresh the panel, reporting rather than aborting."""
        try:
            self.document.set(name, value)
            self.refresh_report()
        except Exception as exc:  # noqa: BLE001
            self.report_problem(exc)

    @guarded
    def _on_section_toggled(self, name, expanded):
        """Remember a section's state for this user, under the section's name. The section has already shown or
        hidden its fields. A store that cannot be written costs only the memory of the state for the next session,
        so the failure is printed and the slot returns normally. Written as "true"/"false", the text a new process
        reads back from INI either way."""
        try:
            self.settings.setValue(section_state_key(name), "true" if expanded else "false")
        except Exception:  # noqa: BLE001
            traceback.print_exc()

    @guarded
    def _on_ipts_edited(self, widget):
        """An IPTS typed in the header: stored as its directory name (normalise_experiment_id), and both derived
        paths follow it. A line edit reports editingFinished on every focus-out; with no typing, nothing is
        written."""
        if not widget.isModified():
            return
        widget.setModified(False)
        self.document.set("experiment_id", normalise_experiment_id(widget.text()))
        self._show(fs.get("experiment_id"), widget, self.document.get("experiment_id"))
        self._show_derived_paths()
        self.refresh_report()

    @guarded
    def _on_path_edited(self, name, widget):
        """A path typed in the header is that override, stripped. Emptied, or only spaces, it is ``None`` again and
        derived, never ``""``: an override of "" or "   " is ``Path("")`` / ``Path("   ")`` to the reduction, not
        the folder it derives. No typing, no write."""
        if not widget.isModified():
            return
        widget.setModified(False)
        self.document.set(name, widget.text().strip() or None)
        self._show(fs.get(name), widget, self.document.get(name))
        self.refresh_report()

    @guarded
    def _browse_path(self, name):
        """A folder chosen with a header path's Browse button is that override, unless it is the folder the
        reduction derives now. A cancelled dialog writes nothing.

        The derived folder is never written as an override: written, it would freeze an absolute path into the
        file, and a file reused for another experiment would then read this one's folders. So choosing the derived
        folder holds no override. That writes nothing when there is none (the dialog opens there, and Choose
        without navigating re-chooses what is shown), and returns a held override to derived, as clearing does.
        Folders are compared with os.path.normpath, never resolve(), which follows symlinks on a facility mount and
        differs between machines.
        """
        editor = self.editors[name]
        derived = self.document.derived_path(name)
        start = editor.text() or derived or ""
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, f"Choose the {fs.get(name).label.lower()}", start)
        if not folder:
            return
        if derived is not None and os.path.normpath(folder) == os.path.normpath(derived):
            if self.document.get(name) is None:
                return
            folder = None
        self.document.set(name, folder)
        self._show(fs.get(name), editor, folder)
        self.refresh_report()

    def _show_derived_paths(self):
        """Each header path's placeholder: the folder the reduction derives while its override is unset, or why
        there is none. A placeholder is never the control's text, so a derived path cannot be read back, or
        saved, as a typed one."""
        for name in _HEADER_PATHS:
            derived = self.document.derived_path(name)
            if derived is None:
                derived = _NOT_A_FOLDER if self.document.get("experiment_id") else _NO_IPTS
            self.editors[name].setPlaceholderText(derived)

    @guarded
    def _on_scalar_edited(self, name, widget):
        self.document.set(name, fs.get(name).coerce(widget.text()))
        self.refresh_report()

    def eventFilter(self, watched, event):
        """Open a drop-down cell from the keyboard (C8).

        Return, Enter or Space on the current cell enters it, as do Alt+Down and
        F4: APG's grid pattern (Enter enters a cell's widget) and its combobox
        (Alt+Down opens the list). F2 already opens any cell (the view's edit
        key). Down stays the grid's ("Moves focus one cell down"), so a drop-down
        column can still be walked with the keyboard.
        """
        if watched is self.angle_table and event.type() == QtCore.QEvent.KeyPress and self._opens_a_cell(event):
            index = self.angle_table.currentIndex()
            if (index.isValid() and fs.PER_ANGLE_NAMES[index.column()] in self._cell_delegates
                    and self.angle_table.indexWidget(index) is None):
                self.angle_table.edit(index)
                return True
        return super().eventFilter(watched, event)

    @staticmethod
    def _opens_a_cell(event):
        key = event.key()
        modifiers = event.modifiers() & ~QtCore.Qt.KeypadModifier
        if key in (QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter, QtCore.Qt.Key_Space, QtCore.Qt.Key_F4):
            return modifiers == QtCore.Qt.NoModifier
        return key == QtCore.Qt.Key_Down and modifiers == QtCore.Qt.AltModifier

    @guarded
    def _on_cell_changed(self, row, column):
        """Write one per-angle value, to the row Qt says was edited.

        `row` comes from the signal — the cell that actually changed. Using
        `self.angle_table.currentRow()` here instead would be the active-row
        bug: with row 2 selected, editing row 0 would write to row 2.

        Coercion goes through the field's ELEMENT type. Every per-angle field is
        a `list[...]`, so a coercer that only understood "int" and "float" left
        every cell as text — `useBS` holding the string "False", which is truthy,
        subtracts background the scientist switched off.
        """
        if self._populating:
            return
        name = fs.PER_ANGLE_NAMES[column]
        item = self.angle_table.item(row, column)
        text = item.text() if item is not None else ""
        field = fs.get(name)
        value = field.coerce_element(text)
        # An enumerated column's drop-down offers its choices in the file's own
        # spelling (C9), and a choice is stored as offered. coerce_element would
        # turn 'constantq' into the declared 'constantQ', rewriting a file
        # whose convention is lower case (PR #36, finding 2).
        if field.allowed and isinstance(value, str) and value != text and value.lower() == text.lower():
            value = text
        self.document.set_angle_field(row, name, value)
        self.refresh_column(column)
        self.refresh_marks()
        self.refresh_report()

    @guarded
    def add_angle(self):
        self.document.add_angle()
        self.refresh_angles()
        self.refresh_report()

    @guarded
    def remove_selected_angle(self):
        """Remove the selected row.

        Here the selection IS the input — the user is saying "this one" — which
        is different from consulting it to decide where an edit lands.
        """
        row = self.angle_table.currentRow()
        if row < 0:
            return
        self.document.remove_angle(row)
        self.refresh_angles()
        self.refresh_report()

    @guarded
    def select_roi(self):
        """B1, B2, B9: the ROI pop-out for the selected row, and what it reports written to that row.

        The row is read once, at the gesture, and passed on: nothing later consults the selection, so the row the
        dialog was opened for is the row written (the active-row trap; E2). Only the fields the dialog reports
        changed are written, through the document: RB_Ymin, RB_Ymax and BkgROI with ``set_angle_field``, which pads a
        short column only as far as the row (review 1568397; E2's short leg); data_x_range, shared by every angle,
        with ``set`` (E8). The dialog writes nothing. The slot writes no settings file and no data file: the one thing
        it records is the folder a chosen run came from, in the launcher's QSettings (``roi_nexus_dir``), so that the
        next file dialog opens there (B2; E4 watches every file and key).
        """
        row = self.angle_table.currentRow()
        if row < 0:
            return
        loaded = self._events_for_row(row)
        if loaded is None:  # the file dialog was cancelled
            return
        events, title, band = loaded
        dialog = ROISelectionDialog(events, self._roi_values(row), title=title, tof_band=band, parent=self)
        try:
            if dialog.exec_() != QtWidgets.QDialog.Accepted:
                return
            changes = dialog.changes()
        finally:
            dialog.deleteLater()  # released, never destroy() (E7): exec_ has hidden it, and Qt frees it in the loop
        if row >= self.document.n_angles:
            raise IndexError(f"angle {row + 1} no longer exists, so nothing was written")
        for name in ("RB_Ymin", "RB_Ymax", "BkgROI"):
            if name in changes:
                self.document.set_angle_field(row, name, changes[name])
        if "data_x_range" in changes:
            self.document.set("data_x_range", changes["data_x_range"])
        self.refresh_angles()
        self.refresh_scalars()
        self.refresh_report()

    def _roi_values(self, row):
        """The row's values the pop-out shows, as the document holds them (a short column reads None)."""
        angle = self.document.angle_row(row)
        values = {name: angle.get(name) for name in ("RB_Ymin", "RB_Ymax", "BkgROI", "tof_min", "tof_max", "useBS")}
        values["data_x_range"] = self.document.get("data_x_range")
        return values

    def _events_for_row(self, row):
        """B2: the row's run, as ``(events, title, tof_band)``, or None when the user cancels choosing a file.

        The file is the reducer's own name for the row's run, ``NEXUSpathRB / REF_L_<RBnum>.nxs.h5``
        (``nr_reduction_calc.py:325``), when RBnum is set and the file exists. Otherwise (an authored file holds no
        RBnum) a file dialog asks, starting in that NeXus folder when it exists, else where it was last. The view
        filter starts at the run's chopper band when it has a chopper log, else the full TOF span (B8). A read
        failure raises, and ``@guarded`` reports it in the panel.
        """
        runs = self.document.get("RBnum")
        run = runs[row] if isinstance(runs, list) and row < len(runs) else None
        try:
            folder = Path(self.document.config.NEXUSpathRB)
        except TypeError:  # an experiment_id of None makes the path property raise
            folder = None
        path = None
        if run is not None and folder is not None:
            path = folder / f"REF_L_{run}.nxs.h5"  # the reducer's own name for the run (nr_reduction_calc.py:325)
            if not path.is_file():
                path = None
        if path is None:
            start = str(folder) if folder is not None and folder.is_dir() else self.settings.value("roi_nexus_dir", "")
            chosen, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, f"The NeXus file of angle {row + 1}", start, "NeXus (*.nxs.h5);;All files (*)"
            )
            if not chosen:
                return None
            path = Path(chosen)
            self.settings.setValue("roi_nexus_dir", str(path.parent))
        events = roi_estimate.load_event_pixels(path, max_events=MAX_ROI_EVENTS)
        try:
            meta = roi_estimate.read_nexus_metadata(path)
            title = f"{meta['title']} (run {meta['run_number']})"
        except (OSError, KeyError, ValueError):
            meta, title = None, path.name
        try:
            band = roi_estimate.lambda_to_tof(roi_estimate.chopper_lambda_range(path), meta["start_time"])
        except (OSError, KeyError, ValueError, TypeError):
            band = None  # no chopper log: the view filter starts at the full span
        return events, title, band

    # -- refresh -----------------------------------------------------------

    @guarded
    def set_document(self, document):
        """Adopt a document and render all of it.

        The single entry point a resolution layer uses: replacing the document
        without the three refreshes leaves the view showing the previous one.
        A document adopted here is shown as it holds: its IPTS is resolved by a
        Load only (``load_settings``; editor-ipts-inference v2, design A3).
        """
        self.document = document
        self.refresh_angles()
        self.refresh_scalars()
        self.refresh_report()

    def report_problem(self, exc):
        """Show a failure in the panel instead of letting it kill the process."""
        self._last_error = f"{type(exc).__name__}: {exc}"
        traceback.print_exc()
        self.report.setPlainText(
            f"Could not complete that action.\n\n  {self._last_error}\n\n"
            f"The settings in this tab are unchanged."
        )

    def refresh_angles(self):
        self._populating = True
        try:
            shown = min(self.document.n_angles, MAX_TABLE_ROWS)
            self._rows_hidden = self.document.n_angles - shown
            self.angle_table.setRowCount(shown)
            for row in range(shown):
                values = self.document.angle_row(row)
                for column, name in enumerate(fs.PER_ANGLE_NAMES):
                    value = values[name]
                    # Through _as_text, not str(): repr of a nested list
                    # ("[120, 130]") is re-parsed by coerce_element into
                    # ['[120', '130]'], so BkgROI was corrupted by any edit to
                    # its row.
                    self.angle_table.setItem(
                        row, column, QtWidgets.QTableWidgetItem(self._cell_text(fs.get(name), value))
                    )
            self.refresh_marks()
        finally:
            self._populating = False

    def refresh_column(self, column):
        """Re-draw one column from the document, after an edit to it.

        An edit can change cells the user did not type in: writing out a list
        the reducer filled itself puts its value at the other angles, and a
        refused λ leaves the edited cell unset. The text is set on the existing
        items, never by replacing them, because this runs inside the table's own
        cellChanged signal for one of them.
        """
        field = fs.get(fs.PER_ANGLE_NAMES[column])
        populating, self._populating = self._populating, True
        try:
            for row in range(self.angle_table.rowCount()):
                # A row past n_angles is left from before an edit that emptied the
                # last λ entry back to None; it holds nothing.
                value = self.document.angle_row(row)[field.name] if row < self.document.n_angles else None
                text = self._cell_text(field, value)
                item = self.angle_table.item(row, column)
                if item is None:
                    self.angle_table.setItem(row, column, QtWidgets.QTableWidgetItem(text))
                elif item.text() != text:
                    item.setText(text)
        finally:
            self._populating = populating

    def refresh_marks(self):
        """Mark the rows beyond the reduction's angle count, from the document as it is now, and repaint the cells.

        Run after every change that can move the count or a list's state: Load,
        Add, Remove, and a cell edit (an angle-defining value typed into a
        surplus row makes it an angle; a choice writes a compact list out). The
        only place a mark is decided. The drop-down columns draw an implied value
        from the document as they paint (``_ChoiceDelegate``). A change elsewhere
        can alter it without touching their items, so the viewport is repainted
        here.
        """
        count = self.document.reduction_angles
        for row in range(self.angle_table.rowCount()):
            self.angle_table.setVerticalHeaderItem(row, self._row_header(row, surplus=row >= count))
        self.angle_table.viewport().update()

    def refresh_scalars(self):
        for name, editor in self.editors.items():
            self._show(fs.get(name), editor, self.document.get(name))
        self._show_derived_paths()

    def refresh_report(self):
        lines = []
        if self._rows_hidden:
            lines.append(
                f"Showing the first {MAX_TABLE_ROWS} of {self.document.n_angles} angles; "
                f"{self._rows_hidden} are not displayed."
            )
            lines.append("")
        problems = self.document.validate()
        if problems:
            lines.append("Problems:")
            lines.extend(f"  - {message}" for message in problems)
        else:
            lines.append("No problems found.")

        # Notes are true of a file that reduces (surplus entries; a default the
        # reducer will fill), so they sit in their own section, apart from the
        # problems.
        notes = self.document.notes()
        if notes:
            lines.append("")
            lines.append("Notes:")
            lines.extend(f"  - {note}" for note in notes)

        # A changed value is spelled as the file holds it (file_spelling, the rules save() applies): one
        # value, one spelling between this report and the saved file. The Angles-table cell keeps the
        # scientists' true/false (_cell_text).
        changed = self.document.changed_vs_seed()
        if changed:
            lines.append("")
            lines.append("Changed from the seed:")
            count = self.document.reduction_angles
            for name in sorted(changed):
                before, after = changed[name]
                field = fs.BY_NAME.get(name)
                lines.append(f"  - {name}: {file_spelling(field, before)} -> {file_spelling(field, after, count)}")
        self.report.setPlainText("\n".join(lines))

    # -- files -------------------------------------------------------------

    def _settings_dialog_folders(self):
        """Where the Load and Save dialogs open, and their sidebar (I5, A5): the IPTS's shared folder and its
        settings folders, unless the remembered folder is already under that IPTS (``load_start_folder``)."""
        ipts = self.document.get("experiment_id")
        remembered = self.settings.value("settings_editor_dir", "")
        return load_start_folder(ipts, remembered), settings_folders(ipts)

    @guarded
    def load_settings(self):
        start, sidebar = self._settings_dialog_folders()
        with _file_dialog_sidebar(sidebar):
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self,
                "Load reduction settings",
                start,
                "Settings (*.json *.dat);;All files (*)",
                "",
                QtWidgets.QFileDialog.DontUseNativeDialog,
            )
        if not path:
            return
        # The refreshes are INSIDE the try. They were outside it, and the catch
        # was only (ValueError, OSError), so a file whose per-angle value is not
        # a sequence ({"tof_min": 5}) raised TypeError out of the slot and
        # aborted the launcher.
        try:
            document = SettingsDocument.from_file(path)
            # The Load's IPTS (editor-ipts-inference, I1), with the IPTS the header holds now as the field's: the
            # file's own; else its runs'; else the header's; else the folder a file without runs came from. Here
            # only (v2, A3): a document injected or adopted again is not resolved, so nothing infers over an IPTS
            # the user cleared or typed (I6).
            document.resolve_ipts(normalise_experiment_id(self.editors["experiment_id"].text()))
            self.set_document(document)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "Could not load settings", str(exc))
            self.report_problem(exc)
            return
        self.settings.setValue("settings_editor_dir", str(Path(path).parent))

    @guarded
    def save_settings(self):
        start, sidebar = self._settings_dialog_folders()
        with _file_dialog_sidebar(sidebar):
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self,
                "Save reduction settings",
                start,
                "Settings (*.json);;All files (*)",
                "",
                QtWidgets.QFileDialog.DontUseNativeDialog,
            )
        if not path:
            return
        # load_from_file dispatches on the suffix, so a name saved without a
        # recognised one cannot be reloaded. Checked against the accepted set
        # rather than "has a suffix": "settings_0.5deg" has suffix ".5deg",
        # which is not a suffix anyone meant.
        if Path(path).suffix.lower() not in (".json", ".dat"):
            path = path + ".json"
        try:
            self.document.save(path)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "Could not save settings", str(exc))
            return
        self.settings.setValue("settings_editor_dir", str(Path(path).parent))
