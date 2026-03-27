#!/usr/bin/python3
import os
from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import (
    QFileDialog,
    QGridLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QComboBox,
    QWidget,
    QLineEdit,
    QMessageBox,
    QHBoxLayout,
    QVBoxLayout,
)

import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

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

        # Main layout: two columns (controls | canvas)
        main_layout = QHBoxLayout()
        self.setLayout(main_layout)

        # Left: controls
        controls = QVBoxLayout()

        # Folder chooser row
        folder_row = QHBoxLayout()
        self.choose_folder_btn = QPushButton("Choose folder")
        folder_row.addWidget(self.choose_folder_btn)
        self.folder_label = QLabel(self)
        self.folder_label.setText(FOLDER_DIRECTIVE)
        folder_row.addWidget(self.folder_label)
        controls.addLayout(folder_row)

        # Filter box
        self.filter_edit = QLineEdit(self)
        self.filter_edit.setPlaceholderText("Filter filenames (substring)")
        controls.addWidget(self.filter_edit)

        # File list
        self.file_list = QListWidget(self)
        self.file_list.setSelectionMode(QListWidget.MultiSelection)
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

        # Y transform option (None, R*Q^4)
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
        self.filter_edit.textChanged.connect(self.apply_filter)
        self.select_all_btn.clicked.connect(self.select_all)
        self.deselect_all_btn.clicked.connect(self.deselect_all)

        # Populate from previous session
        self.read_settings()

    def read_settings(self):
        _folder = self.settings.value("overplot_folder", "")
        if _folder and os.path.isdir(_folder):
            self.folder_label.setText(_folder)
            self.populate_file_list(_folder)

        _xscale = self.settings.value("overplot_xscale", "linear")
        if _xscale in ("linear", "log"):
            self.xscale_combo.setCurrentText(_xscale)
        _ytransform = self.settings.value("overplot_ytransform", "None")
        if _ytransform in ("None", "R*Q^4"):
            self.ytransform_combo.setCurrentText(_ytransform)

    def save_settings(self):
        self.settings.setValue("overplot_folder", self.folder_label.text())
        self.settings.setValue("overplot_xscale", self.xscale_combo.currentText())

    def choose_folder(self):
        _dir = QFileDialog.getExistingDirectory(None, "Select a folder:", os.path.expanduser("~"), QFileDialog.ShowDirsOnly)
        if os.path.isdir(_dir):
            self.folder_label.setText(_dir)
            self.populate_file_list(_dir)
            self.save_settings()

    def populate_file_list(self, folder):
        self.file_list.clear()
        self._files = []
        try:
            files = sorted(os.listdir(folder))
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Could not list directory: {e}")
            return
        for f in files:
            if f.lower().endswith(".dat"):
                self._files.append(f)

        # show all initially
        for f in self._files:
            item = QListWidgetItem(f)
            item.setCheckState(QtCore.Qt.Unchecked)
            self.file_list.addItem(item)

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
            item.setSelected(True)

    def deselect_all(self):
        for i in range(self.file_list.count()):
            item = self.file_list.item(i)
            item.setCheckState(QtCore.Qt.Unchecked)
            item.setSelected(False)
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
            if item.checkState() == QtCore.Qt.Checked or item.isSelected():
                items.append(item.text())

        if len(items) == 0:
            QMessageBox.warning(self, "No files", "No files selected to plot")
            return

        folder = self.folder_label.text()
        if not os.path.isdir(folder):
            QMessageBox.critical(self, "Invalid folder", "The selected folder is not valid")
            return

        transform = self.ytransform_combo.currentText()

        if self.canvas is not None:
            # embedded canvas plotting
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            any_plotted = False
            for fname in items:
                path = os.path.join(folder, fname)
                try:
                    x, y, ey, ex = self._prepare_data(path, transform)
                except Exception as e:
                    QMessageBox.warning(self, "Load/format error", str(e))
                    continue

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
            ax.set_xlabel('Q')
            ax.set_ylabel('R' if transform == "None" else 'R*Q^4')
            self.figure.tight_layout()
            self.canvas.draw()
            self.save_settings()
        else:
            # external pyplot fallback
            plt.figure()
            ax = plt.gca()
            any_plotted = False
            for fname in items:
                path = os.path.join(folder, fname)
                try:
                    x, y, ey, ex = self._prepare_data(path, transform)
                except Exception as e:
                    QMessageBox.warning(self, "Load/format error", str(e))
                    continue

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
            ax.set_xlabel('x')
            ax.set_ylabel('y')
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
