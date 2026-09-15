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
"""

import functools
import traceback
from pathlib import Path

from qtpy import QtCore, QtGui, QtWidgets

from launcher.app_identity import ensure_identity
from lr_reduction import field_spec as fs
from lr_reduction.settings_document import SettingsDocument

#: Above this, populating the table freezes the GUI thread for seconds and
#: costs ~1100x the file size in memory. A settings file with more angles than
#: this is a mistake, not a workload.
MAX_TABLE_ROWS = 500


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

        self.document = document if document is not None else SettingsDocument()
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
        self.angle_table.cellChanged.connect(self._on_cell_changed)
        box.addWidget(self.angle_table)

        buttons = QtWidgets.QHBoxLayout()
        self.add_angle_button = QtWidgets.QPushButton("Add angle")
        self.add_angle_button.clicked.connect(lambda _checked=False: self.add_angle())
        buttons.addWidget(self.add_angle_button)

        self.remove_angle_button = QtWidgets.QPushButton("Remove angle")
        self.remove_angle_button.clicked.connect(lambda _checked=False: self.remove_selected_angle())
        buttons.addWidget(self.remove_angle_button)
        buttons.addStretch(1)
        box.addLayout(buttons)
        return panel

    def _build_scalar_panel(self):
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        inner = QtWidgets.QWidget()
        column = QtWidgets.QVBoxLayout()
        inner.setLayout(column)

        for group in fs.GROUPS:
            scalars = [f for f in fs.fields_in(group) if not f.per_angle]
            if not scalars:
                continue
            box = QtWidgets.QGroupBox(group)
            grid = QtWidgets.QFormLayout()
            box.setLayout(grid)
            for field in scalars:
                editor = self._build_editor(field)
                editor.setToolTip(f"{field.name} — {field.help}")
                self.editors[field.name] = editor
                grid.addRow(field.label, editor)
            column.addWidget(box)

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
            editor = QtWidgets.QComboBox()
            # A blank first entry for the tri-state fields, where a falsy value
            # means "off" and is the class default.
            if field.falsy_means_off:
                editor.addItem("")
            editor.addItems([str(a) for a in field.allowed])
            editor.currentTextChanged.connect(
                lambda text, name=field.name: self._set_scalar(
                    name, fs.get(name).coerce(text) if text else False
                )
            )
            self._show(field, editor, value)
            return editor

        editor = QtWidgets.QLineEdit()
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

    @staticmethod
    def _show_in_combo(editor, field, value):
        """Display `value`, even when it is not one of the offered choices.

        A combo asked to show an unknown value silently displays its first item
        instead, so a file holding an out-of-set value would look like a valid
        one — and saving would then write the substituted value back. Adding the
        stray value as an entry keeps what is shown equal to what is held;
        validate() is what reports it as a problem.
        """
        if value is None or value is False or value == "":
            editor.setCurrentIndex(0 if field.falsy_means_off else -1)
            return
        text = str(value)
        if editor.findText(text) < 0:
            editor.addItem(text)
        editor.setCurrentText(text)

    @staticmethod
    def _as_text(value):
        """Render a stored value for a single-line editor."""
        if value is None:
            return ""
        if isinstance(value, (list, tuple)):
            return ", ".join("" if v is None else str(v) for v in value)
        return str(value)

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
    def _on_scalar_edited(self, name, widget):
        self.document.set(name, fs.get(name).coerce(widget.text()))
        self.refresh_report()

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
        value = fs.get(name).coerce_element(item.text() if item is not None else "")
        self.document.set_angle_field(row, name, value)
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

    # -- refresh -----------------------------------------------------------

    @guarded
    def set_document(self, document):
        """Adopt a document and render all of it.

        The single entry point a resolution layer uses: replacing the document
        without the three refreshes leaves the view showing the previous one.
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
                    # _as_text, not str(): repr of a nested list ("[120, 130]")
                    # is re-parsed by coerce_element into ['[120', '130]'], so
                    # BkgROI was corrupted by any edit to its row.
                    self.angle_table.setItem(
                        row, column, QtWidgets.QTableWidgetItem(self._as_text(value))
                    )
        finally:
            self._populating = False

    def refresh_scalars(self):
        for name, editor in self.editors.items():
            self._show(fs.get(name), editor, self.document.get(name))

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

        changed = self.document.changed_vs_seed()
        if changed:
            lines.append("")
            lines.append("Changed from the seed:")
            for name in sorted(changed):
                before, after = changed[name]
                lines.append(f"  - {name}: {before!r} -> {after!r}")
        self.report.setPlainText("\n".join(lines))

    # -- files -------------------------------------------------------------

    @guarded
    def load_settings(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Load reduction settings",
            self.settings.value("settings_editor_dir", ""),
            "Settings (*.json *.dat);;All files (*)",
        )
        if not path:
            return
        # The refreshes are INSIDE the try. They were outside it, and the catch
        # was only (ValueError, OSError), so a file whose per-angle value is not
        # a sequence ({"tof_min": 5}) raised TypeError out of the slot and
        # aborted the launcher.
        try:
            self.set_document(SettingsDocument.from_file(path))
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "Could not load settings", str(exc))
            self.report_problem(exc)
            return
        self.settings.setValue("settings_editor_dir", str(Path(path).parent))

    @guarded
    def save_settings(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save reduction settings",
            self.settings.value("settings_editor_dir", ""),
            "Settings (*.json);;All files (*)",
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
