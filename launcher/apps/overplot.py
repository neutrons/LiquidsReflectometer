#!/usr/bin/python3
import os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from qtpy import QtCore
from qtpy.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

# Try a few possible Qt backends for matplotlib (Qt5, QtAgg, Qt4)
FigureCanvas = None
NavigationToolbar = None
try:
    # Preferred for matplotlib <=3.x with Qt5
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
except Exception:
    try:
        # Matplotlib 3.5+ provides a qt abstraction backend
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
    except Exception:
        try:
            # Older installations might have qt4agg
            from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
            from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
        except Exception:
            FigureCanvas = None
            NavigationToolbar = None

FOLDER_DIRECTIVE = "Click to choose a folder containing .dat files"


# Two reduction writers produce reflectivity: output.py emits
# "# Q [1/Angstrom] R dR dQ [FWHM]" and save_reduced_data.py emits
# "# columns = Q, R, dR, dQ (sigma)[, ...]". Both sit after a metadata
# preamble — real files carry the marker at lines 13-20 — so the scan must
# follow the comment block, not a fixed line count.
REFLECTIVITY_MARKERS = ("q [1/angstrom]", "columns = q, r, dr, dq")
DIRECT_BEAM_HEADER_MARKER = "lambda intensity error"


def classify_file(path):
    """Classify a plot input by its leading comment block.

    Returns one of "reflectivity", "direct_beam", or "unknown". Scans the
    leading `#` lines (capped at 200) rather than a fixed window: the
    facility writers put the column marker after a variable-length preamble,
    one line per stitched angle.
    """
    try:
        # Both writers emit non-ASCII (Å) in the preamble above the marker, so
        # the decode must not depend on the ambient locale: a UnicodeDecodeError
        # here would reproduce B1's symptom (everything classifies "unknown").
        with open(path, encoding="utf-8", errors="replace") as fh:
            for _ in range(200):
                line = fh.readline()
                if not line:
                    break
                if line.strip() and not line.startswith("#"):
                    # First real data row: the header block is over.
                    break
                low = line.lower()
                if any(marker in low for marker in REFLECTIVITY_MARKERS):
                    return "reflectivity"
                if DIRECT_BEAM_HEADER_MARKER in low:
                    return "direct_beam"
    except (OSError, UnicodeDecodeError):
        return "unknown"
    return "unknown"


def _axis_labels(mode, transform):
    """Return (xlabel, ylabel) for a resolved selection mode."""
    if mode == "reflectivity":
        return "Q [1/Å]", ("R · Q⁴" if transform == "R*Q^4" else "R")
    if mode == "direct_beam":
        return "λ [Å]", "I"
    return "x", "y"


class Overplot(QWidget):
    """Tab for overplotting multiple .dat files.

    Expects files to be whitespace-delimited with 3 or 4 columns:
    x, y, err_y, err_x (err_x optional). If 3 columns are present assume no x-error.
    """

    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle("Overplot")

        self.settings = QtCore.QSettings()
        self._files = []  # full list of .dat filenames in chosen folder
        self._last_sel_mode = None

        # Main layout: two columns (controls | canvas)
        main_layout = QHBoxLayout()
        self.setLayout(main_layout)

        # Left: controls
        controls = QVBoxLayout()

        # Folder chooser row
        folder_row = QHBoxLayout()
        self.choose_folder_btn = QPushButton("Choose folder")
        folder_row.addWidget(self.choose_folder_btn)
        self.folder_edit = QLineEdit(self)
        self.folder_edit.setPlaceholderText(FOLDER_DIRECTIVE)
        folder_row.addWidget(self.folder_edit)
        controls.addLayout(folder_row)

        # Plot mode (Auto / Reflectivity / Direct Beam)
        plot_mode_row = QHBoxLayout()
        self.plot_mode_label = QLabel("Plot mode")
        plot_mode_row.addWidget(self.plot_mode_label)
        self.plot_mode_combo = QComboBox(self)
        self.plot_mode_combo.addItems(["Auto", "Reflectivity", "Direct Beam"])
        plot_mode_row.addWidget(self.plot_mode_combo)
        controls.addLayout(plot_mode_row)

        # Filter box
        self.filter_edit = QLineEdit(self)
        self.filter_edit.setPlaceholderText("Filter filenames (substring)")
        controls.addWidget(self.filter_edit)

        # File list
        self.file_list = QListWidget(self)
        # The check box is the single notion of "chosen": highlight-only rows
        # used to plot as a hidden input that refresh() could not see, so a
        # refresh reported success while leaving a stale plot on screen.
        # NoSelection retires that second notion; the view's own delegate still
        # toggles the indicator on click, so no click handler is needed here —
        # adding one toggles the state straight back and the box stops working.
        self.file_list.setSelectionMode(QListWidget.NoSelection)
        controls.addWidget(self.file_list)

        # Select all / Deselect all
        sel_row = QHBoxLayout()
        self.select_all_btn = QPushButton("Select all")
        self.deselect_all_btn = QPushButton("Deselect all")
        sel_row.addWidget(self.select_all_btn)
        sel_row.addWidget(self.deselect_all_btn)
        controls.addLayout(sel_row)

        # X scale option
        xscale_row = QHBoxLayout()
        self.xscale_label = QLabel("X scale")
        xscale_row.addWidget(self.xscale_label)
        self.xscale_combo = QComboBox(self)
        self.xscale_combo.addItems(["linear", "log"])
        xscale_row.addWidget(self.xscale_combo)
        controls.addLayout(xscale_row)

        # Y transform option (None, R*Q^4) — R*Q^4 only meaningful for reflectivity
        ytransform_row = QHBoxLayout()
        self.ytransform_label = QLabel("Y transform")
        ytransform_row.addWidget(self.ytransform_label)
        self.ytransform_combo = QComboBox(self)
        self.ytransform_combo.addItems(["None", "R*Q^4"])
        ytransform_row.addWidget(self.ytransform_combo)
        controls.addLayout(ytransform_row)

        # Y scale info
        yscale_row = QHBoxLayout()
        self.yscale_label = QLabel("Y scale (fixed)")
        yscale_row.addWidget(self.yscale_label)
        self.yscale_val = QLabel("log")
        yscale_row.addWidget(self.yscale_val)
        controls.addLayout(yscale_row)

        # Plot button
        self.plot_btn = QPushButton("Plot selected")
        controls.addWidget(self.plot_btn)

        # Clear plot button
        self.clear_btn = QPushButton("Clear plot")
        controls.addWidget(self.clear_btn)

        # Refresh button — rescan the folder and replot checked files
        self.refresh_btn = QPushButton("Refresh")
        self._refresh_tooltip_base = (
            "Re-read this folder. Keeps current selection where files still exist; replots checked files."
        )
        self._update_refresh_tooltip(self._read_last_refresh())
        controls.addWidget(self.refresh_btn)

        # Stretch to push controls to top
        controls.addStretch()

        main_layout.addLayout(controls, 1)

        # Right: matplotlib canvas + toolbar (if available). If no Qt backend is
        # available for matplotlib, fall back to external pyplot windows.
        canvas_layout = QVBoxLayout()
        if FigureCanvas is not None:
            self.figure = Figure(figsize=(5, 4))
            self.canvas = FigureCanvas(self.figure)
            self.toolbar = NavigationToolbar(self.canvas, self)
            canvas_layout.addWidget(self.toolbar)
            canvas_layout.addWidget(self.canvas)
        else:
            # Placeholder label explaining fallback
            self.figure = None
            self.canvas = None
            self.toolbar = None
            placeholder = QLabel(
                "Embedded plotting unavailable: matplotlib Qt backend not found.\nPlots will open in external windows."
            )
            placeholder.setWordWrap(True)
            canvas_layout.addWidget(placeholder)

        main_layout.addLayout(canvas_layout, 3)

        # Connections
        self.choose_folder_btn.clicked.connect(self.choose_folder)
        self.plot_btn.clicked.connect(self.plot_selected)
        self.clear_btn.clicked.connect(self.clear_plot)
        self.refresh_btn.clicked.connect(self.refresh)
        self.filter_edit.textChanged.connect(self.apply_filter)
        self.select_all_btn.clicked.connect(self.select_all)
        self.deselect_all_btn.clicked.connect(self.deselect_all)
        self.folder_edit.editingFinished.connect(self.folder_changed)
        self.plot_mode_combo.currentTextChanged.connect(self._on_plot_mode_changed)

        # Populate from previous session
        self.read_settings()
        # After loading settings, apply the user's chosen plot mode to the
        # transform combo's enabled state (Auto leaves it enabled until a plot
        # resolves the per-file modes).
        self._on_plot_mode_changed(self.plot_mode_combo.currentText())

    def read_settings(self):
        _folder = self.settings.value("overplot_folder", "")
        if _folder and os.path.isdir(_folder):
            self.folder_edit.setText(_folder)
            self.populate_file_list(_folder)

        _xscale = self.settings.value("overplot_xscale", "linear")
        if _xscale in ("linear", "log"):
            self.xscale_combo.setCurrentText(_xscale)
        _ytransform = self.settings.value("overplot_ytransform", "None")
        if _ytransform in ("None", "R*Q^4"):
            self.ytransform_combo.setCurrentText(_ytransform)
        _mode = self.settings.value("overplot_mode", "Auto")
        if _mode in ("Auto", "Reflectivity", "Direct Beam"):
            self.plot_mode_combo.setCurrentText(_mode)

    def save_settings(self):
        # Only persist a folder we actually have. read_settings leaves the
        # field empty when the stored path is not a directory (an unmounted
        # /SNS at launch), and __init__ reaches save_settings via
        # _on_plot_mode_changed — writing "" there would destroy the user's
        # stored path just for opening the app before the mount came up.
        folder = self.folder_edit.text().strip()
        if folder and os.path.isdir(folder):
            self.settings.setValue("overplot_folder", folder)
        self.settings.setValue("overplot_xscale", self.xscale_combo.currentText())
        self.settings.setValue("overplot_mode", self.plot_mode_combo.currentText())

    def choose_folder(self):
        _dir = QFileDialog.getExistingDirectory(None, "Select a folder:", os.path.expanduser("~"), QFileDialog.ShowDirsOnly)
        if os.path.isdir(_dir):
            self.folder_edit.setText(_dir)
            self.populate_file_list(_dir)
            self.save_settings()

    def populate_file_list(self, folder):
        self.file_list.clear()
        try:
            self._files = self._scan_folder(folder)
        except OSError as e:
            QMessageBox.critical(self, "Error", f"Could not list directory: {e}")
            self._files = []
            return

        for f in self._files:
            item = QListWidgetItem(f)
            item.setCheckState(QtCore.Qt.Unchecked)
            self.file_list.addItem(item)

    def _scan_folder(self, folder):
        """List `.dat`/`.txt` filenames in *folder*, sorted. Raises OSError."""
        files = sorted(os.listdir(folder))
        return [f for f in files if f.lower().endswith(".dat") or f.lower().endswith(".txt")]

    def _read_last_refresh(self):
        v = self.settings.value("overplot_last_refresh", "")
        return v if isinstance(v, str) else ""

    def _update_refresh_tooltip(self, last):
        suffix = f" Last refreshed: {last or 'never'}."
        self.refresh_btn.setToolTip(self._refresh_tooltip_base + suffix)

    def refresh(self):
        """Two-phase refresh: rescan folder preserving check state, then
        replot whatever is still checked."""
        folder = self.folder_edit.text().strip()
        if not folder or folder == FOLDER_DIRECTIVE:
            QMessageBox.warning(self, "Refresh", "Choose a folder first.")
            return

        pre_checked = {
            self.file_list.item(i).text()
            for i in range(self.file_list.count())
            if self.file_list.item(i).checkState() == QtCore.Qt.Checked
        }

        try:
            on_disk = self._scan_folder(folder)
        except OSError as e:
            QMessageBox.critical(self, "Refresh failed", f"Could not list {folder}: {e}")
            return

        self._files = on_disk
        # Rebuild list applying any currently-typed filter so the visible rows
        # stay consistent with the filter box; apply_filter reads self._files.
        self.apply_filter(self.filter_edit.text())

        # Restore check state from pre_checked for files still on disk.
        for i in range(self.file_list.count()):
            item = self.file_list.item(i)
            if item.text() in pre_checked:
                item.setCheckState(QtCore.Qt.Checked)

        still_checked = pre_checked & set(on_disk)
        dropped = pre_checked - still_checked
        if still_checked:
            self.plot_selected()
        elif pre_checked:
            # Everything that was plotted has gone from disk. Leaving those
            # curves up while stamping "last refreshed" is the freshness lie
            # this button exists to prevent.
            self.clear_plot()

        # Stamped after the replot, so the tooltip cannot claim a refresh the
        # canvas has not caught up with.
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.settings.setValue("overplot_last_refresh", timestamp)
        self._update_refresh_tooltip(timestamp)
        if dropped:
            QMessageBox.information(
                self,
                "Refresh",
                f"{len(dropped)} previously-plotted file(s) no longer present on disk:\n"
                + "\n".join(sorted(dropped)),
            )

    def apply_filter(self, text):
        # preserve checked state
        checked = {self.file_list.item(i).text() for i in range(self.file_list.count()) if self.file_list.item(i).checkState() == QtCore.Qt.Checked}
        self.file_list.clear()
        text = text.strip().lower()
        for f in self._files:
            if text == "" or text in f.lower():
                item = QListWidgetItem(f)
                item.setCheckState(QtCore.Qt.Checked if f in checked else QtCore.Qt.Unchecked)
                self.file_list.addItem(item)

    def select_all(self):
        for i in range(self.file_list.count()):
            item = self.file_list.item(i)
            item.setCheckState(QtCore.Qt.Checked)

    def deselect_all(self):
        for i in range(self.file_list.count()):
            item = self.file_list.item(i)
            item.setCheckState(QtCore.Qt.Unchecked)

    def folder_changed(self):
        folder = self.folder_edit.text().strip()

        if not folder:
            return
        if os.path.isdir(folder):
            self.populate_file_list(folder)
            self.save_settings()
        else:
            QMessageBox.warning(self, "Invalid folder", f"Folder does not exist:\n{folder}")


    def _set_rq4_enabled(self, enabled):
        """Toggle the Qt::ItemIsEnabled bit on the R*Q^4 combo entry."""
        idx = self.ytransform_combo.findText("R*Q^4")
        if idx < 0:
            return
        item = self.ytransform_combo.model().item(idx)
        if item is None:
            return
        flags = item.flags()
        if enabled:
            item.setFlags(flags | QtCore.Qt.ItemIsEnabled)
        else:
            item.setFlags(flags & ~QtCore.Qt.ItemIsEnabled)
            if self.ytransform_combo.currentText() == "R*Q^4":
                self.ytransform_combo.setCurrentText("None")

    def _on_plot_mode_changed(self, mode):
        """User picked an explicit Plot mode; gate the transform combo."""
        if mode == "Direct Beam":
            self._set_rq4_enabled(False)
        else:
            # Reflectivity or Auto: re-enable; Auto's per-file gating
            # happens at plot_selected time.
            self._set_rq4_enabled(True)
        self.save_settings()

    def _resolve_modes(self, paths):
        """Return (per_file_modes, selection_mode) for the user's combo choice."""
        combo_mode = self.plot_mode_combo.currentText()
        if combo_mode == "Reflectivity":
            return ["reflectivity"] * len(paths), "reflectivity"
        if combo_mode == "Direct Beam":
            return ["direct_beam"] * len(paths), "direct_beam"
        # Auto. "unknown" is a kind: a selection mixing recognized and
        # unrecognized files is heterogeneous, and must fall back to x/y
        # rather than let one recognized file label all the others.
        per = [classify_file(p) for p in paths]
        kinds = set(per)
        if not kinds or kinds == {"unknown"}:
            return per, "unknown"
        if len(kinds) == 1:
            return per, next(iter(kinds))
        return per, "mixed"

    def _prepare_data(self, path, transform):
        try:
            data = np.loadtxt(path)
        except Exception as e:
            raise RuntimeError(f"Failed to load {path}: {e}")

        if data is None:
            raise RuntimeError(f"No data in {path}")

        if data.ndim == 1 and data.size == 0:
            raise RuntimeError(f"Empty data in {path}")

        if data.ndim == 1:
            data = data[np.newaxis, :]

        if data.shape[1] < 2:
            raise RuntimeError(f"File does not have at least 2 columns: {path}")

        x = data[:, 0].astype(float)
        y = data[:, 1].astype(float)
        ey = None
        ex = None
        if data.shape[1] >= 3:
            ey = data[:, 2].astype(float)
        if data.shape[1] >= 4:
            ex = data[:, 3].astype(float)

        # Apply transform if requested
        if transform == "R*Q^4":
            # y -> y * x**4
            # error propagation: sigma_f^2 = (Q^4)^2 * sigma_R^2 + (4*R*Q^3)^2 * sigma_Q^2
            x4 = x ** 4
            y_trans = y * x4

            sigma_R = ey if ey is not None else None
            sigma_Q = ex if ex is not None else None

            if sigma_R is None and sigma_Q is None:
                ey_trans = None
            else:
                t1 = (x4 ** 2) * (sigma_R ** 2) if sigma_R is not None else 0
                t2 = ((4 * y * (x ** 3)) ** 2) * (sigma_Q ** 2) if sigma_Q is not None else 0
                ey_trans = np.sqrt(t1 + t2)

            y = y_trans
            ey = ey_trans

        # For log y plotting we must have y > 0 and finite
        mask = np.isfinite(x) & np.isfinite(y) & (y > 0)
        if not np.any(mask):
            raise RuntimeError(f"No positive/finite y values to plot in {path}")

        x_out = x[mask]
        y_out = y[mask]
        ey_out = ey[mask] if (ey is not None) else None
        ex_out = ex[mask] if (ex is not None) else None
        return x_out, y_out, ey_out, ex_out

    def plot_selected(self):
        # collect selected items (checked or selected)
        items = []
        for i in range(self.file_list.count()):
            item = self.file_list.item(i)
            if item.checkState() == QtCore.Qt.Checked:
                items.append(item.text())

        if len(items) == 0:
            QMessageBox.warning(self, "No files", "No files selected to plot")
            return

        folder = self.folder_edit.text().strip()
        if not os.path.isdir(folder):
            QMessageBox.critical(self, "Invalid folder", "The selected folder is not valid")
            return

        paths = [os.path.join(folder, fname) for fname in items]
        per_file_modes, sel_mode = self._resolve_modes(paths)
        self._last_sel_mode = sel_mode

        # Honest detection: if the user-overrode mode disagrees with any
        # file's header, list the mismatched filenames so they know.
        combo_mode = self.plot_mode_combo.currentText()
        if combo_mode in ("Reflectivity", "Direct Beam"):
            override_kind = "reflectivity" if combo_mode == "Reflectivity" else "direct_beam"
            detected = [classify_file(p) for p in paths]
            mismatched = [items[i] for i, kind in enumerate(detected) if kind not in (override_kind, "unknown")]
            if mismatched:
                QMessageBox.warning(
                    self,
                    "Plot mode override",
                    f"{len(mismatched)} file(s) have headers inconsistent with {combo_mode}: {', '.join(mismatched)}",
                )

        # R*Q^4 only meaningful for a homogeneous reflectivity selection.
        self._set_rq4_enabled(sel_mode == "reflectivity")
        transform = self.ytransform_combo.currentText()
        if sel_mode != "reflectivity":
            transform = "None"

        xlabel, ylabel = _axis_labels(sel_mode, transform)

        if self.canvas is not None:
            # embedded canvas plotting
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            any_plotted = False
            for fname, fmode in zip(items, per_file_modes):
                path = os.path.join(folder, fname)
                try:
                    x, y, ey, ex = self._prepare_data(path, transform)
                except Exception as e:
                    QMessageBox.warning(self, "Load/format error", str(e))
                    continue

                # Annotate mixed-selection labels with detected kind so the
                # plot is honest about its contents.
                if sel_mode == "mixed":
                    label = f"{os.path.basename(fname)} [{fmode}]"
                else:
                    label = os.path.basename(fname)
                try:
                    if ey is not None:
                        ax.errorbar(x, y, yerr=ey, label=label, fmt='-o')
                    else:
                        ax.plot(x, y, '-o', label=label)
                    any_plotted = True
                except Exception as e:
                    QMessageBox.warning(self, "Plot error", f"Failed to plot {fname}: {e}")

            if not any_plotted:
                # The figure was cleared above; push that to the widget so the
                # canvas cannot keep showing the previous plot while the dialog
                # says nothing was drawn.
                if self.canvas is not None:
                    self.canvas.draw()
                QMessageBox.information(self, "No data", "No data was plotted")
                return

            ax.set_yscale('log')
            if self.xscale_combo.currentText() == 'log':
                ax.set_xscale('log')
            else:
                ax.set_xscale('linear')

            ax.legend()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            self.figure.tight_layout()
            self.canvas.draw()
            self.save_settings()
        else:
            # external pyplot fallback
            plt.figure()
            ax = plt.gca()
            any_plotted = False
            for fname, fmode in zip(items, per_file_modes):
                path = os.path.join(folder, fname)
                try:
                    x, y, ey, ex = self._prepare_data(path, transform)
                except Exception as e:
                    QMessageBox.warning(self, "Load/format error", str(e))
                    continue

                # Annotate mixed-selection labels with detected kind so the
                # plot is honest about its contents.
                if sel_mode == "mixed":
                    label = f"{os.path.basename(fname)} [{fmode}]"
                else:
                    label = os.path.basename(fname)
                try:
                    if ey is not None:
                        ax.errorbar(x, y, yerr=ey, label=label, fmt='-o')
                    else:
                        ax.plot(x, y, '-o', label=label)
                    any_plotted = True
                except Exception as e:
                    QMessageBox.warning(self, "Plot error", f"Failed to plot {fname}: {e}")

            if not any_plotted:
                QMessageBox.information(self, "No data", "No data was plotted")
                return

            ax.set_yscale('log')
            if self.xscale_combo.currentText() == 'log':
                ax.set_xscale('log')
            else:
                ax.set_xscale('linear')

            ax.legend()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            plt.tight_layout()
            plt.show()

    def clear_plot(self):
        """Clear the embedded canvas or close external matplotlib windows."""
        if self.canvas is not None and self.figure is not None:
            try:
                self.figure.clear()
                self.canvas.draw()
            except Exception:
                # If something goes wrong, fallback to closing all pyplot figures
                plt.close('all')
        else:
            plt.close('all')
