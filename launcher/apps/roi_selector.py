#!/usr/bin/python3
import os

import h5py
import numpy as np

# Try to import scipy curve fit for better peak finding (optional)
try:
    from scipy.optimize import curve_fit
    HAS_SCIPY = True
except Exception:
    curve_fit = None
    HAS_SCIPY = False
# Qt and Matplotlib imports (some are optional depending on environment)
try:
    from qtpy import QtCore
    from qtpy.QtGui import QBrush, QColor, QFont
    from qtpy.QtWidgets import (
        QCheckBox,
        QComboBox,
        QDialog,
        QDialogButtonBox,
        QFileDialog,
        QFormLayout,
        QGridLayout,
        QGroupBox,
        QHBoxLayout,
        QLabel,
        QLineEdit,
        QListWidget,
        QMessageBox,
        QPushButton,
        QSizePolicy,
        QSpinBox,
        QTableWidget,
        QTableWidgetItem,
        QVBoxLayout,
        QWidget,
    )
except Exception:
    # allow headless import failures; variables below will be None/missing and
    # the rest of the code should guard against unavailable GUI elements.
    try:
        import qtpy.QtCore as QtCore
    except Exception:
        QtCore = None

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.colors import LogNorm
    from matplotlib.figure import Figure
except Exception:
    FigureCanvas = None
    NavigationToolbar = None
    Figure = None
    plt = None
    LogNorm = None

class DBPerRunDialog(QDialog):
    """Simple dialog to edit per-run DB_file and q_method lists."""
    def __init__(self, parent=None, runs=None, initial=None, q_initial=None):
        super().__init__(parent)
        self.setWindowTitle("Edit per-run settings")
        self.runs = runs or []
        self.initial = list(initial) if isinstance(initial, (list, tuple)) else ([""] * len(self.runs))
        self.q_initial = list(q_initial) if isinstance(q_initial, (list, tuple)) else ([""] * len(self.runs))

        layout = QVBoxLayout()
        self.table = QTableWidget(self)
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(["Run", "DB_file", "q_method"])
        self.table.setRowCount(len(self.runs))
        q_options = ["meanTheta", "constantTOF", "constantQ"]
        for i, r in enumerate(self.runs):
            it = QTableWidgetItem(str(r))
            it.setFlags(it.flags() & ~QtCore.Qt.ItemIsEditable)
            self.table.setItem(i, 0, it)
            dbv = self.initial[i] if i < len(self.initial) else ''
            self.table.setItem(i, 1, QTableWidgetItem(dbv or ''))
            combo = QComboBox(self.table)
            combo.addItems(q_options)
            qv = self.q_initial[i] if i < len(self.q_initial) else ''
            if qv in q_options:
                combo.setCurrentIndex(q_options.index(qv))
            self.table.setCellWidget(i, 2, combo)

        layout.addWidget(self.table)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.setLayout(layout)

    def get_values(self):
        dbs = []
        qs = []
        for i in range(self.table.rowCount()):
            try:
                db_item = self.table.item(i, 1)
                dbs.append(db_item.text().strip() if db_item is not None else '')
            except Exception:
                dbs.append('')
            try:
                widget = self.table.cellWidget(i, 2)
                qs.append(widget.currentText() if widget is not None else '')
            except Exception:
                qs.append('')
        return {'db_files': dbs, 'q_methods': qs}

class SaveTemplateDialog(QDialog):
    def __init__(self, parent=None, defaults=None, runs=None):
        super().__init__(parent)
        self.defaults = defaults or {}
        self.runs = runs or []
        self.per_run_db = None
        layout = QVBoxLayout()
        form = QFormLayout()
        # Save directory + filename (split fields)
        h_dir = QHBoxLayout()
        self.dir_edit = QLineEdit(self)
        self.dir_browse_btn = QPushButton("...")
        h_dir.addWidget(self.dir_edit)
        h_dir.addWidget(self.dir_browse_btn)
        form.addRow("Save folder:", h_dir)

        self.filename_edit = QLineEdit(self)
        form.addRow("Filename:", self.filename_edit)
        # Note: per-run settings table moved to main ROI selector UI.
        # Save dialog will read per-run settings from its parent when saving.

        self.autoscale_cb = QCheckBox(self)
        self.autoscale_cb.setChecked(self.defaults.get("autoscale", True))
        form.addRow("autoscale:", self.autoscale_cb)

        self.use_calc_theta_cb = QCheckBox(self)
        self.use_calc_theta_cb.setChecked(self.defaults.get("use_calc_theta", True))
        form.addRow("use_calc_theta:", self.use_calc_theta_cb)

        # set defaults if provided
        try:
            default_dir = self.defaults.get('default_dir')
            default_name = self.defaults.get('default_name')
            if default_dir:
                self.dir_edit.setText(str(default_dir))
            if default_name:
                self.filename_edit.setText(str(default_name))
        except Exception:
            pass

        layout.addLayout(form)

        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        # connect browse (directory chooser)
        self.dir_browse_btn.clicked.connect(self._browse_directory)

        self.setLayout(layout)

    def get_values(self):
        # join directory and filename into a single path
        dirv = self.dir_edit.text().strip()
        fname = self.filename_edit.text().strip()
        # always ensure .xml suffix on filename
        if fname and not fname.lower().endswith('.xml'):
            fname = fname + '.xml'
        full = fname if not dirv else os.path.join(dirv, fname)
        # Read per-run settings from parent ROISelector when available
        parent = self.parent()
        db_files = None
        q_methods = None
        try:
            if parent is not None:
                db_files = getattr(parent, 'per_run_db', None)
                q_methods = getattr(parent, 'per_run_qmethods', None)
        except Exception:
            db_files = getattr(self, 'per_run_db', None)
            q_methods = getattr(self, 'per_run_qmethods', None)

        # Maintain legacy DB_file single-string key
        db_single = ''
        if isinstance(db_files, (list, tuple)) and db_files:
            try:
                db_single = ','.join([str(x) for x in db_files])
            except Exception:
                db_single = ''

        return {
            "path": full,
            "dir": dirv,
            "filename": fname,
            "DB_file": db_single,
            "DB_files": db_files,
            "q_methods": q_methods,
            "autoscale": self.autoscale_cb.isChecked(),
            "use_calc_theta": self.use_calc_theta_cb.isChecked(),
        }

    def _browse_directory(self):
        d = QFileDialog.getExistingDirectory(self, "Select folder to save template", os.getcwd())
        if d:
            self.dir_edit.setText(d)

    def _open_db_per_run(self):
        # open the per-run editor dialog with current runs and split values
        # Prefer any previously set per-run values stored in this dialog.
        if isinstance(self.per_run_db, (list, tuple)) and self.per_run_db:
            cur_list = list(self.per_run_db)
        else:
            cur = self.defaults.get("DB_file", "")
            cur_list = [s.strip() for s in cur.split(',') if s.strip()] if cur else []
        # prepare initial q_method list from parent per_run_qmethods if available
        q_init = []
        parent = self.parent()
        try:
            runs = self.runs or []
            if parent and hasattr(parent, 'per_run_qmethods') and isinstance(parent.per_run_qmethods, (list, tuple)) and len(parent.per_run_qmethods) == len(runs):
                q_init = list(parent.per_run_qmethods)
            else:
                for run in runs:
                    q_init.append('')
        except Exception:
            q_init = []

        dlg = DBPerRunDialog(self, runs=self.runs, initial=cur_list, q_initial=q_init)
        if dlg.exec_() == QDialog.Accepted:
            vals = dlg.get_values()
            # vals is a dict with 'db_files' and 'q_methods'
            if isinstance(vals, dict):
                self.per_run_db = vals.get('db_files', [])
                self.per_run_qmethods = vals.get('q_methods', [])
                # notify parent (ROISelector) to refresh its per-run table if present
                try:
                    if parent is not None and hasattr(parent, '_refresh_per_run_table'):
                        parent._refresh_per_run_table()
                except Exception:
                    pass
            else:
                # fallback to legacy list
                self.per_run_db = vals
                self.per_run_qmethods = None
                try:
                    if parent is not None and hasattr(parent, '_refresh_per_run_table'):
                        parent._refresh_per_run_table()
                except Exception:
                    pass


class ROISelector(QWidget):
    """ROI selector tab: load nexus event file and interactively pick ROIs.

    Shows 5 plots simultaneously and overlays ROI dashed lines. Allows saving
    out an XML-like template with ROI and a few extra settings.
    """

    def __init__(self):
        super().__init__()
        self.setWindowTitle("ROI selector")
        self.settings = QtCore.QSettings()

        # Defaults (detector geometry)
        self.n_y = 304
        self.n_x = 256
        self.tof_bin = 50

        # Layout
        main_layout = QGridLayout()
        self.setLayout(main_layout)
        # make left column narrower on initial load (left:1, right:5)
        main_layout.setColumnStretch(0, 1)
        main_layout.setColumnStretch(1, 5)

        # Left panel: controls
        control_box = QGroupBox("Controls")
        c_layout = QFormLayout()
        control_box.setLayout(c_layout)

        self.ipts_edit = QLineEdit(self)
        c_layout.addRow("IPTS number:", self.ipts_edit)

        # template controls moved up so they are near the run list
        self.template_cb = QCheckBox("Start from template", self)
        c_layout.addRow(self.template_cb)

        # template path with browse (kept compact)
        h_tpath = QHBoxLayout()
        self.template_path_edit = QLineEdit(self)
        self.template_path_edit.setMinimumWidth(300)
        self.template_browse_btn = QPushButton("...")
        h_tpath.addWidget(self.template_path_edit)
        h_tpath.addWidget(self.template_browse_btn)
        c_layout.addRow("Template path:", h_tpath)

        # internal run edit (not shown) used to track current run when loading from run list
        self.run_edit = QLineEdit(self)
        # keep the run_edit hidden (we don't show a standalone run number field)
        self.run_edit.hide()

        # multi-run input: comma-separated run list and a widget to navigate/select
        h_runs = QHBoxLayout()
        self.run_list_edit = QLineEdit(self)
        self.run_list_edit.setPlaceholderText("e.g. 1001,1002,1005")
        self.load_runs_btn = QPushButton("Load runs")
        h_runs.addWidget(self.run_list_edit)
        h_runs.addWidget(self.load_runs_btn)
        c_layout.addRow("Run list:", h_runs)

        # list widget showing runs loaded
        # QListWidget already imported at module top
        self.run_list_widget = QListWidget(self)
        self.run_list_widget.setMaximumHeight(100)
        c_layout.addRow(self.run_list_widget)

        # Per-run settings table
        self.per_run_table = QTableWidget(self)
        # Columns: Seq, Run, DB_file, q_method
        self.per_run_table.setColumnCount(4)
        self.per_run_table.setHorizontalHeaderLabels(["Seq", "Run", "DB_file", "q_method"])
        try:
            self.per_run_table.setMinimumHeight(120)
            self.per_run_table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        except Exception:
            try:
                self.per_run_table.setMaximumHeight(150)
            except Exception:
                pass
        # allow editing DB_file cells; q_method is provided as a combobox per row
        self.per_run_table.itemChanged.connect(self._on_per_run_item_changed)
        # place the label on the line above the table to save horizontal space
        per_run_label = QLabel("Per-run settings:")
        c_layout.addRow(per_run_label)
        c_layout.addRow(self.per_run_table)


        # buttons to store ROIs for the current run and to save combined template
        self.store_run_rois_btn = QPushButton("Store ROIs for run")
        self.save_combined_btn = QPushButton("Save combined template")
        self.reset_btn = QPushButton("Reset")
        c_layout.addRow(self.store_run_rois_btn)
        c_layout.addRow(self.save_combined_btn)
        c_layout.addRow(self.reset_btn)

        # optional log color scale for heatmaps
        self.log_color_cb = QCheckBox("Log color scale", self)
        # default to log color scale as requested
        self.log_color_cb.setChecked(True)
        c_layout.addRow(self.log_color_cb)

        # ROI fields group (two-column grid to avoid squishing)
        roi_group = QGroupBox("ROIs")
        roi_layout = QGridLayout()
        roi_group.setLayout(roi_layout)

        # y range
        self.ymin_spin = QSpinBox(self)
        self.ymin_spin.setRange(0, self.n_y - 1)
        self.ymax_spin = QSpinBox(self)
        self.ymax_spin.setRange(0, self.n_y - 1)
        roi_layout.addWidget(QLabel("Y min:"), 0, 0)
        roi_layout.addWidget(self.ymin_spin, 0, 1)
        roi_layout.addWidget(QLabel("Y max:"), 0, 2)
        roi_layout.addWidget(self.ymax_spin, 0, 3)

        # x range
        self.xmin_spin = QSpinBox(self)
        self.xmin_spin.setRange(0, self.n_x - 1)
        self.xmax_spin = QSpinBox(self)
        self.xmax_spin.setRange(0, self.n_x - 1)
        roi_layout.addWidget(QLabel("X min:"), 1, 0)
        roi_layout.addWidget(self.xmin_spin, 1, 1)
        roi_layout.addWidget(QLabel("X max:"), 1, 2)
        roi_layout.addWidget(self.xmax_spin, 1, 3)

        # tof range
        self.tofmin_spin = QSpinBox(self)
        self.tofmin_spin.setRange(0, 10 ** 9)
        self.tofmax_spin = QSpinBox(self)
        self.tofmax_spin.setRange(0, 10 ** 9)
        roi_layout.addWidget(QLabel("TOF min:"), 2, 0)
        roi_layout.addWidget(self.tofmin_spin, 2, 1)
        roi_layout.addWidget(QLabel("TOF max:"), 2, 2)
        roi_layout.addWidget(self.tofmax_spin, 2, 3)

        # background roi (4 ints)
        self.bkg1_min = QSpinBox(self)
        self.bkg1_max = QSpinBox(self)
        self.bkg2_min = QSpinBox(self)
        self.bkg2_max = QSpinBox(self)
        for s in (self.bkg1_min, self.bkg1_max, self.bkg2_min, self.bkg2_max):
            s.setRange(0, self.n_y - 1)
        roi_layout.addWidget(QLabel("Bkg1 min:"), 3, 0)
        roi_layout.addWidget(self.bkg1_min, 3, 1)
        roi_layout.addWidget(QLabel("Bkg1 max:"), 3, 2)
        roi_layout.addWidget(self.bkg1_max, 3, 3)
        roi_layout.addWidget(QLabel("Bkg2 min:"), 4, 0)
        roi_layout.addWidget(self.bkg2_min, 4, 1)
        roi_layout.addWidget(QLabel("Bkg2 max:"), 4, 2)
        roi_layout.addWidget(self.bkg2_max, 4, 3)

        # background mode selector: 'Use 2 background' or 'No background'
        from qtpy.QtWidgets import QComboBox
        self.bkg_mode_cb = QComboBox(self)
        # three options: two, one, or none
        self.bkg_mode_cb.addItems(["Use 2 background", "Use 1 background", "No background"])
        self.bkg_mode_cb.setCurrentIndex(0)
        roi_layout.addWidget(QLabel("Background mode:"), 5, 0)
        roi_layout.addWidget(self.bkg_mode_cb, 5, 1, 1, 3)
        # connect change to handler
        self.bkg_mode_cb.currentIndexChanged.connect(self._set_bkg_mode)

        # add ROI group into the control layout
        c_layout.addRow(roi_group)

        # connect spinboxes to update plots when editing finished (press Return)
        for spin in (self.ymin_spin, self.ymax_spin, self.xmin_spin, self.xmax_spin,
                     self.tofmin_spin, self.tofmax_spin, self.bkg1_min, self.bkg1_max,
                     self.bkg2_min, self.bkg2_max):
            try:
                spin.editingFinished.connect(self.update_rois)
            except Exception:
                try:
                    spin.valueChanged.connect(self.update_rois)
                except Exception:
                    pass

        # Update and save buttons
        self.update_btn = QPushButton("Update ROIs")
        self.update_btn.setToolTip("Redraw plots with current ROI settings")
        # single Save button removed; use 'Save combined template' which handles single or multiple runs
        c_layout.addRow(self.update_btn)

        # place the control box into the main layout (left column)
        main_layout.addWidget(control_box, 0, 0)

        # Right: matplotlib canvas showing 5 plots
        plots_box = QGroupBox("Plots")
        p_layout = QVBoxLayout()
        plots_box.setLayout(p_layout)

        # initialize data holders
        self.counts_xy = None
        self.y_vs_tof = None
        self.tof_edges = None
        self.seq_num = None
        self._chopper_val = None
        self._drag_cids = []
        # rectangle selectors and zoom state for per-plot zoom
        self._rect_selectors = {}
        self._zoom_active_axes = set()
        # per-run stored ROI dict: runnum -> roi dict
        self.per_run_rois = {}
        # per-run sequence mapping (runnum -> seq index from template)
        self.per_run_seq = {}

        # Title label for run
        self.title_label = QLabel("")
        self.title_label.setToolTip("Run title (Base Title) read from the NeXus file if available")
        p_layout.addWidget(self.title_label)

        if FigureCanvas is not None:
            self.figure = Figure(figsize=(10, 8))
            self.canvas = FigureCanvas(self.figure)
            self.toolbar = NavigationToolbar(self.canvas, self)
            p_layout.addWidget(self.toolbar)
            p_layout.addWidget(self.canvas)
        else:
            self.figure = None
            self.canvas = None
            self.toolbar = None
            p_layout.addWidget(QLabel("Embedded plotting unavailable. Install a Qt matplotlib backend."))

        main_layout.addWidget(plots_box, 0, 1)

        # Status line (sequence/run)
        self.status_label = QLabel("")
        self.status_label.setToolTip("Shows loaded run and sequence number")
        main_layout.addWidget(self.status_label, 1, 0, 1, 2)

        # Connect signals
        # run list / multi-run buttons
        self.load_runs_btn.clicked.connect(self.load_runs)
        self.run_list_widget.itemClicked.connect(self._on_run_selected)
        self.store_run_rois_btn.clicked.connect(self.save_rois_for_current_run)
        self.save_combined_btn.clicked.connect(self.save_combined_template)
        self.reset_btn.clicked.connect(self.reset_session)
        # update and browse
        self.update_btn.clicked.connect(self.update_rois)
        self.template_browse_btn.clicked.connect(self._browse_template)
        # ensure the per-run table shows headers immediately
        try:
            self._refresh_per_run_table()
        except Exception:
            pass
    # update when log-color checkbox toggled (created below)

    def _set_bkg_mode(self, index: int):
        """Set background mode based on combo box index.

        Index meanings:
          0 -> Use 2 background (enable both bkg1 and bkg2)
          1 -> Use 1 background (enable only bkg1)
          2 -> No background (disable both and collapse values to Y ROI)
        """
        try:
            mode = self.bkg_mode_cb.currentText()
        except Exception:
            mode = None
        if mode is None:
            return
        m = str(mode).lower()
        if m.startswith('use 2') or '2' in m:
            # enable both
            try:
                for s in (self.bkg1_min, self.bkg1_max, self.bkg2_min, self.bkg2_max):
                    s.setEnabled(True)
            except Exception:
                pass
        elif '1' in m or 'one' in m:
            # enable only bkg1, collapse bkg2
            try:
                self.bkg1_min.setEnabled(True)
                self.bkg1_max.setEnabled(True)
                self.bkg2_min.setEnabled(False)
                self.bkg2_max.setEnabled(False)
                # collapse bkg2 to Y ROI
                y0 = self.ymin_spin.value()
                y1 = self.ymax_spin.value()
                self.bkg2_min.setValue(y0)
                self.bkg2_max.setValue(y1)
            except Exception:
                pass
        else:
            # No background
            try:
                for s in (self.bkg1_min, self.bkg1_max, self.bkg2_min, self.bkg2_max):
                    s.setEnabled(False)
                y0 = self.ymin_spin.value()
                y1 = self.ymax_spin.value()
                self.bkg1_min.setValue(y0)
                self.bkg1_max.setValue(y1)
                self.bkg2_min.setValue(y0)
                self.bkg2_max.setValue(y1)
            except Exception:
                pass
        return

    def _refresh_per_run_table(self):
        """Populate the per-run settings table from self.runs, self.per_run_db, and self.per_run_qmethods."""
        try:
            runs = getattr(self, 'runs', []) or []
            self.per_run_table.blockSignals(True)
            self.per_run_table.setRowCount(len(runs))
            q_options = ["meanTheta", "constantTOF", "constantQ"]
            for i, run in enumerate(runs):
                # Seq (editable) - used to map to template entries
                seq_val = ''
                try:
                    if isinstance(self.per_run_seq, dict) and run in self.per_run_seq:
                        seq_val = str(self.per_run_seq.get(run, '') or '')
                except Exception:
                    seq_val = ''
                item_seq = QTableWidgetItem(seq_val)
                self.per_run_table.setItem(i, 0, item_seq)

                # Run (read-only)
                item_run = QTableWidgetItem(str(run))
                item_run.setFlags(item_run.flags() & ~QtCore.Qt.ItemIsEditable)
                self.per_run_table.setItem(i, 1, item_run)

                # DB_file cell
                dbv = ''
                try:
                    if isinstance(self.per_run_db, (list, tuple)) and i < len(self.per_run_db):
                        dbv = str(self.per_run_db[i]) or ''
                except Exception:
                    dbv = ''
                item_db = QTableWidgetItem(dbv)
                self.per_run_table.setItem(i, 2, item_db)

                # q_method combobox
                q_combo = QComboBox(self.per_run_table)
                q_combo.addItems(q_options)
                qv = ''
                try:
                    if isinstance(self.per_run_qmethods, (list, tuple)) and i < len(self.per_run_qmethods):
                        qv = self.per_run_qmethods[i] or ''
                except Exception:
                    qv = ''
                if qv in q_options:
                    q_combo.setCurrentIndex(q_options.index(qv))
                self.per_run_table.setCellWidget(i, 3, q_combo)
                # connect per-row combo handler
                try:
                    q_combo.currentIndexChanged.connect(self._on_qcombo_changed)
                except Exception:
                    pass
            self.per_run_table.blockSignals(False)
        except Exception:
            pass

    def _on_per_run_item_changed(self, item: QTableWidgetItem):
        # update self.per_run_db when a DB_file cell is edited
        try:
            col = item.column()
            row = item.row()
            runs = getattr(self, 'runs', []) or []
            if col == 2:
                # DB_file edited
                val = item.text().strip() if item.text() is not None else ''
                if not isinstance(self.per_run_db, list) or len(self.per_run_db) != len(runs):
                    self.per_run_db = [''] * len(runs)
                if 0 <= row < len(self.per_run_db):
                    self.per_run_db[row] = val
            elif col == 0:
                # Seq edited - store into per_run_seq[run] if possible
                val = item.text().strip() if item.text() is not None else ''
                try:
                    seq_val = int(val) if val else None
                except Exception:
                    seq_val = None
                try:
                    if not isinstance(self.per_run_seq, dict):
                        self.per_run_seq = {}
                    if 0 <= row < len(runs):
                        runnum = runs[row]
                        self.per_run_seq[runnum] = seq_val
                except Exception:
                    pass
        except Exception:
            pass

    def _on_qcombo_changed(self, _index):
        # find which row emitted the signal and update per_run_qmethods
        try:
            sender = self.sender()
            for r in range(self.per_run_table.rowCount()):
                if self.per_run_table.cellWidget(r, 3) is sender:
                    runs = getattr(self, 'runs', []) or []
                    if not isinstance(self.per_run_qmethods, list) or len(self.per_run_qmethods) != len(runs):
                        self.per_run_qmethods = [''] * len(runs)
                    try:
                        self.per_run_qmethods[r] = sender.currentText()
                    except Exception:
                        self.per_run_qmethods[r] = ''
                    break
        except Exception:
            pass

    def _file_path_from_fields(self):
        ipts = self.ipts_edit.text().strip()
        run = self.run_edit.text().strip()
        if not ipts or not run:
            return None
        return f"/SNS/REF_L/IPTS-{ipts}/nexus/REF_L_{run}.nxs.h5"

    def _browse_template(self):
        # Default to IPTS folder if IPTS provided, otherwise use cwd
        ipts = self.ipts_edit.text().strip() if hasattr(self, 'ipts_edit') else ''
        if ipts:
            start_dir = f"/SNS/REF_L/IPTS-{ipts}"
            if not os.path.isdir(start_dir):
                start_dir = os.getcwd()
        else:
            start_dir = os.getcwd()
        fname, _ = QFileDialog.getOpenFileName(self, "Open template", start_dir, "XML files (*.xml);;All files (*)")
        if fname:
            self.template_path_edit.setText(fname)

    def _get_chopper_estimate(self, h5file):
        """Attempt to compute a TOF min/max from chopper log values in the
        provided NeXus file. Uses the following algorithm (per request):

            scaled_width = 3.5
            chopper_lam = f['entry/DASlogs/BL4B:Det:TH:BL:Lambda/value'][0]
            chopper_speed = f['entry/DASlogs/BL4B:Det:TH:BL:Frequency/value'][0]
            wl_min = chopper_lam - (scaled_width / 2) * 60.0 / chopper_speed
            wl_max = chopper_lam + (scaled_width / 2) * 60.0 / chopper_speed

        The wavelengths are converted to TOF using distance = 1830 mm and the
        relation: t = lambda * distance / 3956

        Returns
        -------
        tuple(int_tof_min, int_tof_max) or None
            Integer TOF min/max (same units used by event_time_offset in the
            NeXus file) or None if the required log fields are not present.
        """
        scaled_width = 3.5
        try:
            # read chopper wavelength and frequency from expected paths
            # these may be scalars or small arrays; pick first element if needed
            try:
                chopper_lam = h5file['entry/DASlogs/BL4B:Det:TH:BL:Lambda/value'][0]
            except Exception:
                chopper_lam = None
            try:
                chopper_speed = h5file['entry/DASlogs/BL4B:Det:TH:BL:Frequency/value'][0]
            except Exception:
                chopper_speed = None

            if chopper_lam is None or chopper_speed is None:
                return None

            # ensure numeric
            chopper_lam = float(chopper_lam)
            chopper_speed = float(chopper_speed)
            if chopper_speed == 0:
                return None

            wl_min = chopper_lam - (scaled_width / 2.0) * 60.0 / chopper_speed
            wl_max = chopper_lam + (scaled_width / 2.0) * 60.0 / chopper_speed

            # convert wavelength (Angstrom) to TOF using distance 1830 mm
            dist_m = 15.75
            tof_min = int(252.78 * wl_min * dist_m)
            tof_max = int(252.78 * wl_max * dist_m)
            print(f"Chopper estimate: lam={chopper_lam:.3f}A, speed={chopper_speed:.1f}Hz -> "
                  f"wl_min={wl_min:.3f}A, wl_max={wl_max:.3f}A -> "
                  f"tof_min={tof_min}, tof_max={tof_max}")

            # guard: ensure min < max
            if tof_min >= tof_max:
                return None

            return tof_min, tof_max
        except Exception:
            return None

    def load_file(self):
        path = self._file_path_from_fields()
        if not path or not os.path.isfile(path):
            QMessageBox.critical(self, "File error", f"File not found: {path}")
            return

        try:
            with h5py.File(path, "r") as f:
                e_offset = np.array(f['entry/bank1_events/event_time_offset'][:])
                event_id = np.array(f['entry/bank1_events/event_id'][:])
                # try to read sequence number for template indexing
                try:
                    seq_num = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Num/value'][0]
                    self.seq_num = int(seq_num)
                except Exception:
                    self.seq_num = None
                # try to read title/base title
                try:
                    base_title = None
                    try:
                        base_title = f['entry/base_title'][0]
                    except Exception:
                        try:
                            base_title = f['entry/title'][0]
                        except Exception:
                            base_title = None
                    if base_title is not None:
                        # base_title may be bytes in some files; decode if necessary
                        try:
                            if isinstance(base_title, (bytes, bytearray)):
                                b = base_title.decode('utf-8', errors='ignore')
                                self.title_label.setText(b)
                            else:
                                self.title_label.setText(str(base_title))
                        except Exception:
                            self.title_label.setText(str(base_title))
                except Exception:
                    pass
                # Try to estimate a chopper-related value from logs/datasets to help pick TOF ranges
                try:
                    ch_val = self._get_chopper_estimate(f)
                    self._chopper_val = ch_val
                except Exception:
                    self._chopper_val = None
                # chopper-based TOF estimate (if available) is computed by
                # _get_chopper_estimate() above and stored in self._chopper_val.
                # No additional heuristic scanning is required here.
        except Exception as e:
            QMessageBox.critical(self, "Read error", f"Failed to read event data: {e}")
            return

        self.event_id = event_id
        self.e_offset = e_offset

        # Build counts per pixel
        # Map event ids to x,y using detector geometry: event_id = x * n_y + y
        xvals = event_id // self.n_y
        yvals = event_id % self.n_y
        valid = (xvals >= 0) & (xvals < self.n_x) & (yvals >= 0) & (yvals < self.n_y)
        x_good = xvals[valid]
        y_good = yvals[valid]

        counts = np.zeros((self.n_y, self.n_x), dtype=int)
        np.add.at(counts, (y_good, x_good), 1)
        self.counts_xy = counts

        # useful projections for ROI guessing and background selection
        y_sum = counts.sum(axis=1)
        x_sum = counts.sum(axis=0)
        self._y_sum = y_sum
        self._x_sum = x_sum

        # title already read from file (if available)

        # create tof bins
        # Prefer chopper-derived TOF range when available; otherwise default
        # to the min/max TOF present in the file.
        ch = getattr(self, '_chopper_val', None)
        e_min = int(np.nanmin(e_offset))
        e_max = int(np.nanmax(e_offset))
        if ch and isinstance(ch, tuple) and len(ch) == 2:
            try:
                tof_min, tof_max = int(ch[0]), int(ch[1])
            except Exception:
                tof_min, tof_max = e_min, e_max
            # Validate and expand tiny or out-of-range windows so they cover
            # a sensible range relative to the event TOF span.
            if tof_max <= tof_min:
                tof_min, tof_max = e_min, e_max
            # If the chopper window is completely outside event offsets, ignore it
            if tof_max < e_min or tof_min > e_max:
                tof_min, tof_max = e_min, e_max
            # Ensure minimum width: at least two tof bins (so plotting has a visible range)
            min_width = max(self.tof_bin * 2, int(0.02 * (e_max - e_min)))
            if (tof_max - tof_min) < min_width:
                center = (tof_min + tof_max) // 2
                half = max(min_width // 2, 1)
                tof_min = max(e_min, center - half)
                tof_max = min(e_max, center + half)
        else:
            tof_min, tof_max = e_min, e_max
        tof_edges = np.arange(tof_min, tof_max + self.tof_bin, self.tof_bin)
        self.tof_edges = tof_edges

        # y vs tof collapsed along x
        # compute bin indices for each event
        tof_idx = np.digitize(e_offset[valid], tof_edges) - 1
        n_tof = len(tof_edges) - 1
        if n_tof <= 0:
            # no bins available, create a single-bin array and clamp indices
            n_tof = 1
            tof_edges = np.array([int(np.nanmin(e_offset)), int(np.nanmax(e_offset)) + 1])
            tof_idx = np.zeros_like(tof_idx)
        y_vs_tof = np.zeros((self.n_y, n_tof), dtype=int)
        # avoid indexing errors by masking out-of-range bin indices
        valid_tof_mask = (tof_idx >= 0) & (tof_idx < n_tof)
        if np.any(valid_tof_mask):
            np.add.at(y_vs_tof, (y_good[valid_tof_mask], tof_idx[valid_tof_mask]), 1)
        self.y_vs_tof = y_vs_tof

        # Set initial ROI values either from template or from improved peak detection
        if self.template_cb.isChecked() and os.path.isfile(self.template_path_edit.text()):
            templ_path = self.template_path_edit.text()
            applied = False
            try:
                # Prefer the project's reduction_template_reader if available
                try:
                    from lr_reduction import new_reduction_template_reader as rtr
                except Exception:
                    rtr = None

                if rtr is not None:
                    try:
                        with open(templ_path, 'r') as tf:
                            txt = tf.read()
                        entries = rtr.from_xml(txt)
                    except Exception:
                        entries = []

                    # Build a sequence-index -> (DB_file, q_method) mapping from entries
                    seq_map = {}
                    try:
                        for xi, ent in enumerate(entries):
                            try:
                                seq_map[xi + 1] = (getattr(ent, 'DB_file', None), getattr(ent, 'q_method', None))
                            except Exception:
                                seq_map[xi + 1] = (None, None)
                    except Exception:
                        seq_map = {}

                    # Prefer the sequence index read from the NeXus file
                    seq = getattr(self, 'seq_num', None)
                    chosen = None
                    if seq is not None and isinstance(seq, int) and seq >= 1:
                        if 1 <= seq <= len(entries):
                            chosen = entries[seq - 1]
                        else:
                            QMessageBox.warning(self, "Template", f"Template has {len(entries)} entries but sequence {seq} requested; using defaults.")

                    # If not found by seq, try to match by run number in data_files
                    if chosen is None:
                        try:
                            runnum = int(self.run_edit.text().strip())
                        except Exception:
                            runnum = None
                        if runnum is not None:
                            for ent in entries:
                                df = getattr(ent, 'data_files', None)
                                try:
                                    if df and int(runnum) in df:
                                        chosen = ent
                                        break
                                except Exception:
                                    continue

                    # apply chosen entry if found (same as before)
                    if chosen is not None:
                        try:
                            if getattr(chosen, 'data_peak_range', None):
                                self.ymin_spin.setValue(int(chosen.data_peak_range[0]))
                                self.ymax_spin.setValue(int(chosen.data_peak_range[1]))
                            if getattr(chosen, 'data_x_range', None):
                                self.xmin_spin.setValue(int(chosen.data_x_range[0]))
                                self.xmax_spin.setValue(int(chosen.data_x_range[1]))
                            if getattr(chosen, 'tof_range', None):
                                self.tofmin_spin.setValue(int(chosen.tof_range[0]))
                                self.tofmax_spin.setValue(int(chosen.tof_range[1]))
                            if getattr(chosen, 'background_roi', None):
                                b = chosen.background_roi
                                if len(b) >= 4:
                                    self.bkg1_min.setValue(int(b[0])); self.bkg1_max.setValue(int(b[1]))
                                    self.bkg2_min.setValue(int(b[2])); self.bkg2_max.setValue(int(b[3]))
                            # set background mode according to template flags (subtract_background, two_backgrounds)
                            try:
                                sb = getattr(chosen, 'subtract_background', True)
                                two = getattr(chosen, 'two_backgrounds', True)
                                if not sb:
                                    # No background
                                    self.bkg_mode_cb.setCurrentText("No background")
                                else:
                                    # Use two backgrounds when template requests it; otherwise default to two
                                    self.bkg_mode_cb.setCurrentText("Use 2 background")
                            except Exception:
                                pass
                            applied = True
                        except Exception:
                            applied = False

                    # Now populate per-run DB/q lists using seq_map and known seq numbers per run
                    try:
                        runs_loaded = getattr(self, 'runs', None) or []
                    except Exception:
                        runs_loaded = []
                    if runs_loaded:
                        per_db = [''] * len(runs_loaded)
                        per_q = [''] * len(runs_loaded)
                        # fill existing entries if present
                        try:
                            if isinstance(self.per_run_db, (list, tuple)) and len(self.per_run_db) == len(runs_loaded):
                                per_db = list(self.per_run_db)
                        except Exception:
                            pass
                        try:
                            if isinstance(self.per_run_qmethods, (list, tuple)) and len(self.per_run_qmethods) == len(runs_loaded):
                                per_q = list(self.per_run_qmethods)
                        except Exception:
                            pass

                        # determine seq for each run (prefer stored per_run_rois[run]['seq']; fall back to current seq for current run)
                        for i, r in enumerate(runs_loaded):
                            seq_for_run = None
                            try:
                                if isinstance(self.per_run_seq, dict) and r in self.per_run_seq:
                                    seq_for_run = self.per_run_seq.get(r)
                            except Exception:
                                seq_for_run = None
                            # if this row corresponds to the current run being loaded, use self.seq_num
                            try:
                                cur_run = int(self.run_edit.text().strip())
                            except Exception:
                                cur_run = None
                            if cur_run == r and getattr(self, 'seq_num', None) is not None:
                                seq_for_run = getattr(self, 'seq_num', None)
                            if seq_for_run is not None and seq_for_run in seq_map:
                                per_db[i] = seq_map[seq_for_run][0] or ''
                                per_q[i] = seq_map[seq_for_run][1] or ''
                                # record seq into per_run_seq for that run so Seq column shows it
                                try:
                                    if not isinstance(self.per_run_seq, dict):
                                        self.per_run_seq = {}
                                    self.per_run_seq[r] = seq_for_run
                                except Exception:
                                    pass
                        self.per_run_db = per_db
                        self.per_run_qmethods = per_q
                        try:
                            self._refresh_per_run_table()
                        except Exception:
                            pass
                    else:
                        # single-run fallback
                        try:
                            if seq is not None and seq in seq_map:
                                self.per_run_db = [seq_map[seq][0]]
                                self.per_run_qmethods = [seq_map[seq][1]]
                                try:
                                    # also record onto per_run_seq for current run if possible
                                    cur_run = int(self.run_edit.text().strip())
                                    if cur_run:
                                        if not isinstance(self.per_run_seq, dict):
                                            self.per_run_seq = {}
                                        self.per_run_seq[cur_run] = seq
                                except Exception:
                                    pass
                            else:
                                self.per_run_db = []
                                self.per_run_qmethods = []
                        except Exception:
                            self.per_run_db = []
                            self.per_run_qmethods = []

                # Fallback: simple XML parsing for older/simple templates
                if not applied:
                    try:
                        import xml.dom.minidom as md
                        dom = md.parse(templ_path)
                        def get(tag):
                            el = dom.getElementsByTagName(tag)
                            if el and el[0].firstChild:
                                return el[0].firstChild.data.strip()
                            return None
                        data_peak_from = get('from_peak_pixels')
                        data_peak_to = get('to_peak_pixels')
                        if data_peak_from and data_peak_to:
                            self.ymin_spin.setValue(int(data_peak_from))
                            self.ymax_spin.setValue(int(data_peak_to))
                            applied = True
                        tof_from = get('from_tof_range')
                        tof_to = get('to_tof_range')
                        if tof_from and tof_to:
                            self.tofmin_spin.setValue(int(float(tof_from)))
                            self.tofmax_spin.setValue(int(float(tof_to)))
                            applied = True
                        # read older-style background tags if present
                        try:
                            back1_from = get('back_roi1_from')
                            back1_to = get('back_roi1_to')
                            back2_from = get('back_roi2_from')
                            back2_to = get('back_roi2_to')
                            if back1_from is not None and back1_to is not None:
                                self.bkg1_min.setValue(int(back1_from))
                                self.bkg1_max.setValue(int(back1_to))
                            if back2_from is not None and back2_to is not None:
                                self.bkg2_min.setValue(int(back2_from))
                                self.bkg2_max.setValue(int(back2_to))
                            # background flag
                            background_flag = get('background_flag')
                            two_bg = get('two_backgrounds')
                            if background_flag is not None and background_flag.lower() in ('false', '0'):
                                try:
                                    self.bkg_mode_cb.setCurrentText("No background")
                                except Exception:
                                    pass
                            else:
                                try:
                                    self.bkg_mode_cb.setCurrentText("Use 2 background")
                                except Exception:
                                    pass
                        except Exception:
                            pass
                        # read DB_file and q_method if present in the simple XML
                        try:
                            db_simple = get('DB_file') or get('DBFile') or None
                            q_simple = get('q_method') or get('qMethod') or None
                            if db_simple is not None or q_simple is not None:
                                try:
                                    runs_loaded = getattr(self, 'runs', None) or []
                                except Exception:
                                    runs_loaded = []

                                # prepare per-run lists preserving existing values where possible
                                per_db = []
                                per_q = []
                                if isinstance(self.per_run_db, (list, tuple)) and len(self.per_run_db) == len(runs_loaded):
                                    per_db = list(self.per_run_db)
                                else:
                                    per_db = [None] * len(runs_loaded)
                                if isinstance(self.per_run_qmethods, (list, tuple)) and len(self.per_run_qmethods) == len(runs_loaded):
                                    per_q = list(self.per_run_qmethods)
                                else:
                                    per_q = [None] * len(runs_loaded)

                                # prioritize sequence-based assignment (like other template fields)
                                try:
                                    seq = getattr(self, 'seq_num', None)
                                except Exception:
                                    seq = None
                                try:
                                    cur_run = int(self.run_edit.text().strip())
                                except Exception:
                                    cur_run = None

                                # if sequence is available and RefLData nodes exist, take the db/q from the seq-th entry
                                try:
                                    all_ds = dom.getElementsByTagName('RefLData')
                                except Exception:
                                    all_ds = []

                                chosen_db = None
                                chosen_q = None
                                if seq is not None and seq >= 1 and all_ds and seq - 1 < len(all_ds):
                                    ds = all_ds[seq - 1]
                                    el_db = ds.getElementsByTagName('DB_file')
                                    if el_db and el_db[0].firstChild:
                                        chosen_db = el_db[0].firstChild.data.strip()
                                    el_q = ds.getElementsByTagName('q_method')
                                    if el_q and el_q[0].firstChild:
                                        chosen_q = el_q[0].firstChild.data.strip()
                                else:
                                    # fallback to the simple single-value tags
                                    chosen_db = db_simple
                                    chosen_q = q_simple

                                if runs_loaded:
                                    if cur_run is not None and cur_run in runs_loaded:
                                        idx = runs_loaded.index(cur_run)
                                        per_db[idx] = chosen_db
                                        per_q[idx] = chosen_q
                                    elif seq is not None and seq - 1 < len(runs_loaded):
                                        per_db[seq - 1] = chosen_db
                                        per_q[seq - 1] = chosen_q
                                    else:
                                        # fallback: set first run
                                        per_db[0] = chosen_db
                                        per_q[0] = chosen_q
                                    # finalize lists
                                    per_db = [x if x is not None else '' for x in per_db]
                                    per_q = [x if x is not None else '' for x in per_q]
                                    self.per_run_db = per_db
                                    self.per_run_qmethods = per_q
                                else:
                                    self.per_run_db = [chosen_db] if chosen_db is not None else []
                                    self.per_run_qmethods = [chosen_q] if chosen_q is not None else []

                                try:
                                    self._refresh_per_run_table()
                                except Exception:
                                    pass
                        except Exception:
                            pass
                    except Exception:
                        applied = False

            except Exception:
                applied = False

            if not applied:
                # If template couldn't be applied for this run/seq, warn and fall back
                QMessageBox.warning(self, "Template", "Could not apply a template entry for this run/sequence; using automatic defaults")
                self._guess_rois(counts, y_vs_tof, tof_edges)
        else:
            # improved automatic selection using peak and FWHM estimation
            self._guess_rois(counts, y_vs_tof, tof_edges)

        # background ROIs: choose defaults relative to detected Y ROI (offsets +/ -)
        ymin = self.ymin_spin.value()
        ymax = self.ymax_spin.value()
        # Bkg defaults depend on the selected background mode
        mode = None
        try:
            mode = self.bkg_mode_cb.currentText()
        except Exception:
            mode = None
        use_background = mode is not None and not mode.lower().startswith('no')
        # Bkg1 default: just below the Y ROI by 5 pixels
        if use_background:
            b1min = max(0, ymin - 5)
            #b1max = max(b1min, ymin - 1)
            b2min = max(0, b1min - 8)
            # Bkg2 default: above the Y ROI by 8 pixels
            #b2min = min(self.n_y - 1, ymax + 1)
            b1max = min(self.n_y - 1, ymax + 5)
            b2max = min(self.n_y - 1, b1max + 8)
        else:
            b1min = ymin
            b1max = ymin
            b2min = ymax
            b2max = ymax

        self.bkg1_min.setValue(int(b1min))
        self.bkg1_max.setValue(int(b1max))
        self.bkg2_min.setValue(int(b2min))
        self.bkg2_max.setValue(int(b2max))
        # Ensure spinboxes enabled state matches the selected mode
        try:
            self._set_bkg_mode(self.bkg_mode_cb.currentIndex())
        except Exception:
            pass

        # Update status label with run/sequence info
        runnum = self.run_edit.text().strip()
        seq = getattr(self, 'seq_num', None)
        if seq is not None:
            self.status_label.setText(f"Run: {runnum}    Seq: {seq}")
        else:
            self.status_label.setText(f"Run: {runnum}")
        # Plot everything
        self._draw_plots()

    def _draw_plots(self):
        if self.figure is None:
            # fallback: open external pyplot windows
            fig = plt.figure(figsize=(12, 8))
            ax1 = fig.add_subplot(231)
            ax2 = fig.add_subplot(232)
            ax3 = fig.add_subplot(233)
            ax4 = fig.add_subplot(234)
            ax5 = fig.add_subplot(235)
        else:
            self.figure.clear()
            axs = []
            axs.append(self.figure.add_subplot(231))
            axs.append(self.figure.add_subplot(232))
            axs.append(self.figure.add_subplot(233))
            axs.append(self.figure.add_subplot(234))
            axs.append(self.figure.add_subplot(235))
            ax1, ax2, ax3, ax4, ax5 = axs
            fig = self.figure

        counts = self.counts_xy
        # 1) x vs y pixel heatmap
        norm1 = None
        try:
            if getattr(self, 'log_color_cb', None) is not None and self.log_color_cb.isChecked():
                pos = counts[counts > 0]
                vmin = float(np.min(pos)) if pos.size else 1.0
                vmax = float(np.max(counts)) if counts.size else 1.0
                norm1 = LogNorm(vmin=max(1e-3, vmin), vmax=max(1.0, vmax))
        except Exception:
            norm1 = None
        im1 = ax1.imshow(counts, aspect='auto', origin='lower', norm=norm1)
        fig.colorbar(im1, ax=ax1)

        # 2) y pixel vs TOF heatmap (collapsed along x)
        if self.y_vs_tof is not None:
            norm2 = None
            try:
                if getattr(self, 'log_color_cb', None) is not None and self.log_color_cb.isChecked():
                    pos2 = self.y_vs_tof[self.y_vs_tof > 0]
                    vmin2 = float(np.min(pos2)) if pos2.size else 1.0
                    vmax2 = float(np.max(self.y_vs_tof)) if self.y_vs_tof.size else 1.0
                    norm2 = LogNorm(vmin=max(1e-3, vmin2), vmax=max(1.0, vmax2))
            except Exception:
                norm2 = None
            im2 = ax2.imshow(self.y_vs_tof, aspect='auto', origin='lower', extent=(self.tof_edges[0], self.tof_edges[-1], 0, self.n_y), norm=norm2)
            # axis labels
            ax2.set_xlabel('TOF')
            ax2.set_ylabel('Y pixel')
            fig.colorbar(im2, ax=ax2)

        # 3) total counts vs y-pixel (we plot below with masking for log scale)

        # 4) total counts vs x-pixel
        x_total = counts.sum(axis=0)
        ax4.plot(np.arange(self.n_x), x_total)
        ax4.set_xlabel('X pixel')
        ax4.set_ylabel('Counts')

        # 5) TOF vs total counts (integrated)
        tof_total = self.y_vs_tof.sum(axis=0) if self.y_vs_tof is not None else np.zeros(len(self.tof_edges) - 1)
        ax5.plot((self.tof_edges[:-1] + self.tof_edges[1:]) / 2, tof_total)
        ax5.set_xlabel('TOF')
        ax5.set_ylabel('Counts')

        # axis labels for heatmap
        ax1.set_xlabel('X pixel')
        ax1.set_ylabel('Y pixel')

        # Overlay ROI dashed lines
        ymin = self.ymin_spin.value()
        ymax = self.ymax_spin.value()
        xmin = self.xmin_spin.value()
        xmax = self.xmax_spin.value()
        tofmin = self.tofmin_spin.value()
        tofmax = self.tofmax_spin.value()

        # pickable overlay lines on heatmaps and their color (green normally, black for log scale)
        try:
            use_log = getattr(self, 'log_color_cb', None) is not None and self.log_color_cb.isChecked()
        except Exception:
            use_log = False
        overlay_color = 'black' if use_log else 'green'

        heat_artists = []
        try:
            h1_ymin = ax1.axhline(y=ymin, color=overlay_color, linestyle='--')
            h1_ymax = ax1.axhline(y=ymax, color=overlay_color, linestyle='--')
            h1_xmin = ax1.axvline(x=xmin, color=overlay_color, linestyle='--')
            h1_xmax = ax1.axvline(x=xmax, color=overlay_color, linestyle='--')
            heat_artists.extend([h1_ymin, h1_ymax, h1_xmin, h1_xmax])
            # make heatmap overlays pickable so users can drag directly on heatmaps
            try:
                for ha in (h1_ymin, h1_ymax, h1_xmin, h1_xmax):
                    ha.set_picker(7)
                    ha.set_pickradius(7)
            except Exception:
                pass
        except Exception:
            heat_artists = []

        if self.y_vs_tof is not None:
            try:
                h2_ymin = ax2.axhline(y=ymin, color=overlay_color, linestyle='--')
                h2_ymax = ax2.axhline(y=ymax, color=overlay_color, linestyle='--')
                h2_tofmin = ax2.axvline(x=tofmin, color=overlay_color, linestyle='--')
                h2_tofmax = ax2.axvline(x=tofmax, color=overlay_color, linestyle='--')
                heat_artists.extend([h2_ymin, h2_ymax, h2_tofmin, h2_tofmax])
                try:
                    for ha in (h2_ymin, h2_ymax, h2_tofmin, h2_tofmax):
                        ha.set_picker(7)
                        ha.set_pickradius(7)
                except Exception:
                    pass
            except Exception:
                pass

        # zoom color plots around ROI +/- 75 pixels where possible (do NOT change TOF x-limits)
        try:
            xcenter = (xmin + xmax) / 2.0
            ycenter = (ymin + ymax) / 2.0
            xlo = max(0, int(round(xcenter - 75)))
            xhi = min(self.n_x - 1, int(round(xcenter + 75)))
            ylo = max(0, int(round(ycenter - 75)))
            yhi = min(self.n_y - 1, int(round(ycenter + 75)))
            ax1.set_xlim(xlo, xhi)
            ax1.set_ylim(ylo, yhi)
            if self.y_vs_tof is not None:
                ax2.set_ylim(ylo, yhi)
        except Exception:
            pass

        # (Rectangle zoom selectors removed - zoom checkboxes were removed from the UI)

        # y-total: vertical lines show Y ROI (use log scale)
        # avoid zeros for log scale by masking
        y_total = counts.sum(axis=1)
        y_masked = np.ma.masked_less_equal(y_total, 0)
        ax3.plot(np.arange(self.n_y), y_masked)
        ax3.set_yscale('log')
        ax3.set_ylim(bottom=max(1e-3, np.nanmin(y_total[y_total>0]) if np.any(y_total>0) else 1e-3))
        # label y-total plot and default x-limits around peak center (+/-40)
        try:
            ax3.set_xlabel('Y pixel')
            ax3.set_ylabel('Counts')
            if y_total.size:
                y_peak = int(np.argmax(y_total))
                ax3.set_xlim(max(0, y_peak - 40), min(self.n_y - 1, y_peak + 40))
        except Exception:
            pass
        # create ROI lines on the y-total and keep references for dragging
        line_ymin = ax3.axvline(x=ymin, color='k', linestyle='--')
        line_ymax = ax3.axvline(x=ymax, color='k', linestyle='--')
        try:
            for ln in (line_ymin, line_ymax):
                ln.set_picker(7)
                ln.set_pickradius(7)
        except Exception:
            pass

        # create small marker handles on ax3 to make grabbing vertical lines easier
        try:
            top3 = ax3.get_ylim()[1]
            handle_ymin = ax3.plot([ymin], [top3], marker='s', color='k', markersize=6)[0]
            handle_ymax = ax3.plot([ymax], [top3], marker='s', color='k', markersize=6)[0]
            for h in (handle_ymin, handle_ymax):
                try:
                    h.set_picker(7)
                except Exception:
                    pass
        except Exception:
            handle_ymin = None
            handle_ymax = None

        # x-total: vertical lines show X ROI
        line_xmin = ax4.axvline(x=xmin, color='k', linestyle='--')
        line_xmax = ax4.axvline(x=xmax, color='k', linestyle='--')
        try:
            for ln in (line_xmin, line_xmax):
                ln.set_picker(7)
                ln.set_pickradius(7)
        except Exception:
            pass
        try:
            top4 = ax4.get_ylim()[1]
            handle_xmin = ax4.plot([xmin], [top4], marker='s', color='k', markersize=6)[0]
            handle_xmax = ax4.plot([xmax], [top4], marker='s', color='k', markersize=6)[0]
            for h in (handle_xmin, handle_xmax):
                try:
                    h.set_picker(7)
                except Exception:
                    pass
        except Exception:
            handle_xmin = None
            handle_xmax = None

        # tof total: vertical lines show TOF ROI
        line_tofmin = ax5.axvline(x=tofmin, color='k', linestyle='--')
        line_tofmax = ax5.axvline(x=tofmax, color='k', linestyle='--')
        try:
            for ln in (line_tofmin, line_tofmax):
                ln.set_picker(7)
                ln.set_pickradius(7)
        except Exception:
            pass
        try:
            top5 = ax5.get_ylim()[1]
            handle_tofmin = ax5.plot([tofmin], [top5], marker='s', color='k', markersize=6)[0]
            handle_tofmax = ax5.plot([tofmax], [top5], marker='s', color='k', markersize=6)[0]
            for h in (handle_tofmin, handle_tofmax):
                try:
                    h.set_picker(7)
                except Exception:
                    pass
        except Exception:
            handle_tofmin = None
            handle_tofmax = None

        # background ranges on y are highlighted as red dashed lines (they refer to Y-pixel ranges)
        # Compute background ranges according to selected mode
        mode_txt = ''
        try:
            mode_txt = self.bkg_mode_cb.currentText().lower()
        except Exception:
            mode_txt = ''
        # default: Use 2 background (explicit ranges)
        if mode_txt.startswith('use 2'):
            b1min = self.bkg1_min.value()
            b1max = self.bkg1_max.value()
            b2min = self.bkg2_min.value()
            b2max = self.bkg2_max.value()
        elif '1' in mode_txt or 'one' in mode_txt:
            # Use 1 background: bkg1 fields define outer left and outer right limits.
            # Background regions are then [b1min, ymin-1] and [ymax+1, b1max].
            outer_left = self.bkg1_min.value()
            outer_right = self.bkg1_max.value()
            # clamp relative to main ROI
            #b1min = outer_left
            #b1max = max(0, ymin - 1)
            #b2min = min(self.n_y - 1, ymax + 1)
            #b2max = outer_right
            b1min = outer_left
            b1max = outer_right
            b2min = ymin
            b2max = ymax
        else:
            # No background: collapse to main ROI so no background pixels selected
            b1min = ymin
            b1max = ymax
            b2min = ymin
            b2max = ymax
        # mark on heatmaps (visual only)
        for ax in (ax1, ax2):
            ax.axhline(y=b1min, color='r', linestyle='--')
            ax.axhline(y=b1max, color='r', linestyle='--')
            ax.axhline(y=b2min, color='r', linestyle='--')
            ax.axhline(y=b2max, color='r', linestyle='--')
        # mark on y-total and keep pickable references so backgrounds can be dragged
        b1min_line = ax3.axvline(x=b1min, color='r', linestyle='--')
        b1max_line = ax3.axvline(x=b1max, color='r', linestyle='--')
        b2min_line = ax3.axvline(x=b2min, color='r', linestyle='--')
        b2max_line = ax3.axvline(x=b2max, color='r', linestyle='--')
        try:
            for ln in (b1min_line, b1max_line, b2min_line, b2max_line):
                ln.set_picker(7)
                ln.set_pickradius(7)
        except Exception:
            pass

        # create small marker handles for background lines to make grabbing easier
        try:
            top3 = ax3.get_ylim()[1]
            handle_b1min = ax3.plot([b1min], [top3], marker='s', color='r', markersize=6)[0]
            handle_b1max = ax3.plot([b1max], [top3], marker='s', color='r', markersize=6)[0]
            handle_b2min = ax3.plot([b2min], [top3], marker='s', color='r', markersize=6)[0]
            handle_b2max = ax3.plot([b2max], [top3], marker='s', color='r', markersize=6)[0]
            for h in (handle_b1min, handle_b1max, handle_b2min, handle_b2max):
                try:
                    h.set_picker(7)
                except Exception:
                    pass
        except Exception:
            handle_b1min = handle_b1max = handle_b2min = handle_b2max = None

        # Set up draggable behavior (only for embedded canvas)
        if self.figure is not None and self.canvas is not None:
            # disconnect previous handlers
            for cid in getattr(self, '_drag_cids', []):
                try:
                    self.canvas.mpl_disconnect(cid)
                except Exception:
                    pass
            self._drag_cids = []

            # map names to line artists and their clamps
            # include heatmap overlay lines as extras so picks on heatmaps also trigger drags
            draggable_map = {
                'ymin': (line_ymin, handle_ymin, (h1_ymin, h2_ymin), 0, self.n_y - 1, self.ymin_spin),
                'ymax': (line_ymax, handle_ymax, (h1_ymax, h2_ymax), 0, self.n_y - 1, self.ymax_spin),
                'xmin': (line_xmin, handle_xmin, (h1_xmin,), 0, self.n_x - 1, self.xmin_spin),
                'xmax': (line_xmax, handle_xmax, (h1_xmax,), 0, self.n_x - 1, self.xmax_spin),
                'tofmin': (line_tofmin, handle_tofmin, (h2_tofmin,), int(self.tof_edges[0]) if self.tof_edges is not None else 0, int(self.tof_edges[-1]) if self.tof_edges is not None else 1, self.tofmin_spin),
                'tofmax': (line_tofmax, handle_tofmax, (h2_tofmax,), int(self.tof_edges[0]) if self.tof_edges is not None else 0, int(self.tof_edges[-1]) if self.tof_edges is not None else 1, self.tofmax_spin),
                'b1min': (b1min_line, handle_b1min, (), 0, self.n_y - 1, self.bkg1_min),
                'b1max': (b1max_line, handle_b1max, (), 0, self.n_y - 1, self.bkg1_max),
                'b2min': (b2min_line, handle_b2min, (), 0, self.n_y - 1, self.bkg2_min),
                'b2max': (b2max_line, handle_b2max, (), 0, self.n_y - 1, self.bkg2_max),
            }

            # no extra artist name mapping needed; heatmap overlay lines are non-pickable
            artist_name_map = {}

            self._active_line = None

            def on_pick(event):
                art = getattr(event, 'artist', None)
                # choose the first matching artist (line or handle)
                for name, tup in draggable_map.items():
                    try:
                        ln = tup[0]
                        handle = tup[1] if len(tup) > 1 else None
                    except Exception:
                        ln = tup[0]
                        handle = None
                    if art is ln or art is handle:
                        self._active_line = name
                        break

            def on_motion(event):
                # if no active drag, switch cursor when hovering near a draggable line
                if self._active_line is None:
                    try:
                        if event.x is None:
                            return
                        # default arrow
                        cur = QtCore.Qt.ArrowCursor
                        ex = event.x
                        for name, tup2 in draggable_map.items():
                            try:
                                ln = tup2[0]
                                xdata = ln.get_xdata()
                                if not xdata:
                                    continue
                                xval = float(xdata[0])
                                # convert data x to pixel coords (take first coord)
                                xpix = ln.axes.transData.transform((xval, 0))[0]
                                if abs(ex - xpix) < 8:
                                    cur = QtCore.Qt.SizeHorCursor
                                    break
                            except Exception:
                                continue
                        self.canvas.setCursor(cur)
                    except Exception:
                        pass
                    return
                if event.inaxes is None or event.xdata is None:
                    return
                name = self._active_line
                tup = draggable_map[name]
                # unpack possibly (line, handle, lo, hi, sp)
                ln = tup[0]
                handle = tup[1] if len(tup) > 1 else None
                lo = tup[-3]
                hi = tup[-2]
                sp = tup[-1]
                x = int(round(event.xdata))
                x = max(lo, min(hi, x))
                # update the line position
                try:
                    ln.set_xdata([x, x])
                    # move handle if present
                    try:
                        if handle is not None:
                            # place handle at top of axis
                            axhl = ln.axes
                            topy = axhl.get_ylim()[1]
                            handle.set_data([x], [topy])
                    except Exception:
                        pass
                except Exception:
                    pass
                # update corresponding spinbox without triggering redraw of whole UI
                try:
                    sp.setValue(int(x))
                except Exception:
                    pass
                self.canvas.draw_idle()

            def on_release(event):
                self._active_line = None

            self._drag_cids.append(self.canvas.mpl_connect('pick_event', on_pick))
            self._drag_cids.append(self.canvas.mpl_connect('motion_notify_event', on_motion))
            self._drag_cids.append(self.canvas.mpl_connect('button_release_event', on_release))

        if self.figure is None:
            plt.tight_layout()
            plt.show()
        else:
            self.figure.tight_layout()
            self.canvas.draw()

    def update_rois(self):
        # redraw plots with updated ROI overlays
        if self.counts_xy is None:
            QMessageBox.warning(self, "No data", "No file loaded")
            return
        self._draw_plots()

    def _guess_rois(self, counts, y_vs_tof, tof_edges):
        """Estimate ROIs by finding peak indices and FWHM-like widths in x and y totals."""
        try:
            y_sum = counts.sum(axis=1)
            x_sum = counts.sum(axis=0)

            # Prefer SciPy-based gaussian fits (Option B) if available; fallback to half-max method
            def _fit_peak(arr, max_width_default=10, nmax=None):
                """Return (min_idx, max_idx) around a fitted peak or fallback window."""
                n = len(arr)
                if nmax is None:
                    nmax = n
                if np.all(arr == 0):
                    # empty signal
                    return 0, min(n - 1, max_width_default)
                peak = int(np.argmax(arr))
                # Try SciPy fit
                if HAS_SCIPY:
                    try:
                        x = np.arange(n)
                        # gaussian with offset: a*exp(-(x-x0)^2/(2*sigma^2)) + c
                        def gauss(x, a, x0, sigma, c):
                            return a * np.exp(-((x - x0) ** 2) / (2.0 * sigma ** 2)) + c

                        p0 = [float(arr.max()), float(peak), max(1.0, float(n) / 10.0), float(np.min(arr))]
                        popt, _ = curve_fit(gauss, x, arr, p0=p0, maxfev=20000)
                        a, x0, sigma, c = popt
                        left = int(max(0, np.floor(x0 - 4 * abs(sigma))))
                        right = int(min(n - 1, np.ceil(x0 + 4 * abs(sigma))))
                        if left >= right:
                            left = max(0, peak - max_width_default)
                            right = min(n - 1, peak + max_width_default)
                        return left, right
                    except Exception:
                        pass
                # Fallback half-max approach
                try:
                    half = float(arr.max()) / 2.0
                    idx = np.where(arr >= half)[0]
                    if len(idx) > 0:
                        return int(idx[0]), int(idx[-1])
                except Exception:
                    pass
                # Final fallback: fixed window around peak
                return max(0, peak - max_width_default), min(n - 1, peak + max_width_default)

            ymin, ymax = _fit_peak(y_sum, max_width_default=10)
            xmin, xmax = _fit_peak(x_sum, max_width_default=10)

            # TOF estimation: prefer chopper-derived window if available, otherwise use TOF histogram peak
            tofmin = None
            tofmax = None
            if getattr(self, '_chopper_val', None) is not None and self.e_offset is not None:
                try:
                    center = int(np.median(self.e_offset))
                    # heuristic window scaling based on chopper value magnitude
                    window = int(max(2000, 2000 + 200 * abs(float(self._chopper_val))))
                    tofmin = max(int(tof_edges[0]) if tof_edges is not None else 0, center - window)
                    tofmax = min(int(tof_edges[-1]) if tof_edges is not None else int(np.nanmax(self.e_offset)), center + window)
                except Exception:
                    tofmin = None
                    tofmax = None

            if (tofmin is None or tofmax is None) and y_vs_tof is not None and tof_edges is not None:
                tof_sum = y_vs_tof.sum(axis=0)
                if tof_sum.max() > 0:
                    try:
                        tpeak = int(np.argmax(tof_sum))
                        halft = float(tof_sum.max()) / 2.0
                        t_idx = np.where(tof_sum >= halft)[0]
                        if len(t_idx) > 0:
                            tofmin = int(tof_edges[t_idx[0]])
                            tofmax = int(tof_edges[min(t_idx[-1] + 1, len(tof_edges) - 1)])
                        else:
                            tofmin = int(max(tof_edges[0], tof_edges[tpeak] - 2000))
                            tofmax = int(min(tof_edges[-1], tof_edges[tpeak] + 2000))
                    except Exception:
                        tofmin = int(max(tof_edges[0], (tof_edges[0] + tof_edges[-1]) // 2 - 2000))
                        tofmax = int(min(tof_edges[-1], (tof_edges[0] + tof_edges[-1]) // 2 + 2000))

            # Final fallback if no TOF info
            if tofmin is None or tofmax is None:
                try:
                    tofmin = 0
                    tofmax = min(10000, int(np.nanmax(self.e_offset)))
                except Exception:
                    tofmin = 0
                    tofmax = 20000

            # If the user did not supply a template, use requested default X and TOF ranges
            use_template = (getattr(self, 'template_cb', None) is not None and self.template_cb.isChecked() and os.path.isfile(self.template_path_edit.text()))
            if not use_template:
                # default X-range while user works out better fitting
                xmin, xmax = 75, 190
                xmin = max(0, min(self.n_x - 1, int(xmin)))
                xmax = max(0, min(self.n_x - 1, int(xmax)))
                # default TOF window: try to compute from chopper lambda request and chopper speed
                try:
                    lam_min = None
                    lam_max = None
                    # try to use nr_tools.get_lam_range if available
                    try:
                        from lr_reduction import nr_tools
                    except Exception:
                        nr_tools = None
                    if nr_tools is not None and getattr(self, '_chopper_lambda', None) is not None and getattr(self, '_chopper_speed', None) is not None:
                        try:
                            lr = nr_tools.get_lam_range(self._chopper_lambda, self._chopper_speed)
                            lam_min, lam_max = float(lr[0]), float(lr[1])
                        except Exception:
                            lam_min = lam_max = None

                    if lam_min is not None and lam_max is not None:
                        # need flight path length to convert lambda->TOF; try discovered value else default to 20.0 m
                        L = getattr(self, '_flight_path', None) or 20.0
                        # conversion factor TOF(µs) = (m_n/h)*1e-10*1e6 * L * lambda(Å) ≈ 252.78 * L * lambda
                        factor = 252.78
                        tofmin = int(max(int(tof_edges[0]) if tof_edges is not None else 0, factor * L * lam_min))
                        tofmax = int(min(int(tof_edges[-1]) if tof_edges is not None else int(np.nanmax(self.e_offset)) if self.e_offset is not None else 100000, factor * L * lam_max))
                    else:
                        try:
                            te0 = int(tof_edges[0]) if tof_edges is not None else 0
                            teN = int(tof_edges[-1]) if tof_edges is not None else int(np.nanmax(self.e_offset)) if self.e_offset is not None else 100000
                        except Exception:
                            te0 = 0
                            teN = 100000
                        tofmin = max(te0, 100)
                        tofmax = min(teN, 100000)
                except Exception:
                    lam_min = lam_max = None

            # set UI values
            self.ymin_spin.setValue(int(ymin))
            self.ymax_spin.setValue(int(ymax))
            self.xmin_spin.setValue(int(xmin))
            self.xmax_spin.setValue(int(xmax))
            self.tofmin_spin.setValue(int(tofmin))
            self.tofmax_spin.setValue(int(tofmax))

        except Exception:
            # fallback simple center window
            y_peak = int(np.argmax(counts.sum(axis=1)))
            self.ymin_spin.setValue(max(0, y_peak - 10))
            self.ymax_spin.setValue(min(self.n_y - 1, y_peak + 10))
            x_peak = int(np.argmax(counts.sum(axis=0)))
            self.xmin_spin.setValue(max(0, x_peak - 10))
            self.xmax_spin.setValue(min(self.n_x - 1, x_peak + 10))
            self.tofmin_spin.setValue(0)
            self.tofmax_spin.setValue(20000)

    # --- Multi-run helpers -------------------------------------------------
    def load_runs(self):
        """Parse the run list edit field and populate the run list widget.
        Supports comma or whitespace separated integers and simple ranges like 100-105.
        """
        import re
        text = self.run_list_edit.text().strip()
        if not text:
            QMessageBox.warning(self, "No runs", "Please enter a comma-separated list of runs")
            return
        tokens = re.split(r"[,\s]+", text)
        runs = []
        for t in tokens:
            if not t:
                continue
            if '-' in t:
                try:
                    a, b = t.split('-', 1)
                    a = int(a); b = int(b)
                    if b >= a:
                        runs.extend(list(range(a, b + 1)))
                except Exception:
                    continue
            else:
                try:
                    runs.append(int(t))
                except Exception:
                    continue

        if not runs:
            QMessageBox.warning(self, "No valid runs", "No valid run numbers found in input")
            return

        # store the loaded runs for later validation and DB editing
        self.runs = runs

        self.run_list_widget.clear()
        for r in runs:
            self.run_list_widget.addItem(str(r))

        # select first and load it
        self.run_list_widget.setCurrentRow(0)
        self.run_edit.setText(str(runs[0]))
        self.load_file()
        # refresh per-run settings table to reflect any loaded template values
        try:
            self._refresh_per_run_table()
        except Exception:
            pass
        # refresh visual marks for stored ROIs (none initially)
        try:
            self._refresh_run_list_marks()
        except Exception:
            pass

    def _on_run_selected(self, item):
        try:
            runnum = int(item.text())
        except Exception:
            return
        self.run_edit.setText(str(runnum))
        self.load_file()
        # if ROIs stored for this run, apply them
        if runnum in self.per_run_rois:
            self._apply_stored_rois(runnum)

    def _refresh_run_list_marks(self):
        """Visually mark runs in the run list which have stored ROIs.

        Stored runs are shown with a green foreground and a bold font.
        """
        try:
            font_stored = QFont()
            font_stored.setBold(True)
            for i in range(self.run_list_widget.count()):
                it = self.run_list_widget.item(i)
                try:
                    rnum = int(it.text())
                except Exception:
                    # try strip trailing markers like ' (stored)'
                    try:
                        rnum = int(it.text().split()[0])
                    except Exception:
                        rnum = None
                if rnum is not None and rnum in self.per_run_rois:
                    it.setForeground(QBrush(QColor('darkgreen')))
                    it.setFont(font_stored)
                else:
                    it.setForeground(QBrush(QColor('black')))
                    f = it.font()
                    f.setBold(False)
                    it.setFont(f)
        except Exception:
            pass

    def save_rois_for_current_run(self):
        try:
            runnum = int(self.run_edit.text().strip())
        except Exception:
            QMessageBox.warning(self, "Run error", "Current run number is invalid")
            return
        roi = {
            'ymin': int(self.ymin_spin.value()),
            'ymax': int(self.ymax_spin.value()),
            'xmin': int(self.xmin_spin.value()),
            'xmax': int(self.xmax_spin.value()),
            'tofmin': int(self.tofmin_spin.value()),
            'tofmax': int(self.tofmax_spin.value()),
            'bkg1': (int(self.bkg1_min.value()), int(self.bkg1_max.value())),
            'bkg2': (int(self.bkg2_min.value()), int(self.bkg2_max.value())),
            'template_path': self.template_path_edit.text().strip(),
            'bkg_mode': self.bkg_mode_cb.currentText(),
            'seq': getattr(self, 'seq_num', None),
        }
        self.per_run_rois[runnum] = roi
        QMessageBox.information(self, "Stored", f"Stored ROIs for run {runnum}")
        try:
            self._refresh_run_list_marks()
        except Exception:
            pass

    def _apply_stored_rois(self, runnum):
        r = self.per_run_rois.get(runnum)
        if not r:
            return
        try:
            self.ymin_spin.setValue(r['ymin'])
            self.ymax_spin.setValue(r['ymax'])
            self.xmin_spin.setValue(r['xmin'])
            self.xmax_spin.setValue(r['xmax'])
            self.tofmin_spin.setValue(r['tofmin'])
            self.tofmax_spin.setValue(r['tofmax'])
            self.bkg1_min.setValue(r['bkg1'][0])
            self.bkg1_max.setValue(r['bkg1'][1])
            self.bkg2_min.setValue(r['bkg2'][0])
            self.bkg2_max.setValue(r['bkg2'][1])
        except Exception:
            pass
        # redraw
        try:
            self.update_rois()
        except Exception:
            pass

    def reset_session(self):
        """Clear all stored inputs and start a fresh session.

        This clears stored per-run ROIs, clears the template path, resets spinboxes
        to sensible defaults, and refreshes the run list visual markers.
        """
        try:
            # clear stored ROI data
            self.per_run_rois = {}
            # clear loaded runs and run list widgets
            try:
                self.runs = []
            except Exception:
                pass
            try:
                self.run_list_widget.clear()
            except Exception:
                pass
            try:
                self.run_list_edit.setText("")
            except Exception:
                pass
            try:
                self.run_edit.setText("")
            except Exception:
                pass

            # clear template path and sequence/title/status markers
            try:
                self.template_path_edit.setText("")
            except Exception:
                pass
            try:
                self.seq_num = None
            except Exception:
                pass
            try:
                self.title_label.setText("")
            except Exception:
                pass
            try:
                self.status_label.setText("")
            except Exception:
                pass

            # reset ROI spinboxes to defaults
            try:
                self.ymin_spin.setValue(0)
                self.ymax_spin.setValue(min(self.n_y - 1, 20))
                self.xmin_spin.setValue(0)
                self.xmax_spin.setValue(min(self.n_x - 1, 100))
                self.tofmin_spin.setValue(0)
                self.tofmax_spin.setValue(20000)
                self.bkg1_min.setValue(0)
                self.bkg1_max.setValue(0)
                self.bkg2_min.setValue(0)
                self.bkg2_max.setValue(0)
            except Exception:
                pass

            # clear any loaded data arrays used for plotting
            try:
                self.counts_xy = None
            except Exception:
                pass
            try:
                self.counts = None
            except Exception:
                pass
            try:
                self.y_vs_tof = None
            except Exception:
                pass
            try:
                self.tof_edges = None
            except Exception:
                pass

            # clear figure/canvas
            try:
                if getattr(self, 'figure', None) is not None:
                    self.figure.clf()
                if getattr(self, 'canvas', None) is not None:
                    self.canvas.draw_idle()
            except Exception:
                pass

            # refresh visual marks
            try:
                self._refresh_run_list_marks()
            except Exception:
                pass

        except Exception:
            QMessageBox.warning(self, "Reset failed", "Failed to reset session")

    def save_combined_template(self):
        if not self.per_run_rois:
            QMessageBox.warning(self, "No ROIs", "No per-run ROIs stored. Use 'Store ROIs for run' first.")
            return
        # compute sensible defaults for save dialog: IPTS-aware folder and filename based on lowest run
        defaults = {"q_method": "meanTheta", "DB_file": "", "autoscale": True, "use_calc_theta": True}
        try:
            # prefer last saved folder in this session if available
            last = getattr(self, '_last_saved_template_dir', None)
            if last:
                defaults['default_dir'] = last
            else:
                ipts = self.ipts_edit.text().strip()
                if ipts:
                    defaults['default_dir'] = f"/SNS/REF_L/IPTS-{ipts}/shared"
        except Exception:
            pass
        try:
            runs = sorted(list(self.per_run_rois.keys()))
            if runs:
                defaults['default_name'] = f"REF_L_{runs[0]}_template.xml"
        except Exception:
            pass
        dlg = SaveTemplateDialog(self, defaults=defaults)
        # ensure we pass the loaded runs so the DB editor can show the correct rows
        try:
            loaded_runs = getattr(self, 'runs', None) or []
        except Exception:
            loaded_runs = []
        dlg = SaveTemplateDialog(self, defaults=defaults, runs=loaded_runs)
        # validate that all loaded runs have stored ROIs
        if loaded_runs:
            missing = [r for r in loaded_runs if r not in self.per_run_rois]
            if missing:
                QMessageBox.warning(self, "Missing ROIs", f"Not all runs have stored ROIs. Missing: {missing}")
                return

        if dlg.exec_() != QDialog.Accepted:
            return
        vals = dlg.get_values()
        print('stored vals', vals)
        path = vals.get('path')
        # Persist per-run q_method into per_run_rois for this session if provided
        try:
            q_methods = vals.get('q_methods') if isinstance(vals, dict) else None
            if q_methods:
                for idx, runnum in enumerate(runs_sorted):
                    try:
                        if idx < len(q_methods) and q_methods[idx]:
                            self.per_run_rois[runnum]['q_method'] = q_methods[idx]
                    except Exception:
                        pass
        except Exception:
            pass
        # DB_file validation: allow comma-separated list matching the number of stored runs
        db_raw = vals.get('DB_file', '') or ''
        db_list = [s.strip() for s in db_raw.split(',') if s.strip()] if db_raw else []
        runs_sorted = sorted(list(self.per_run_rois.keys()))
        if db_list and len(db_list) != len(runs_sorted):
            QMessageBox.warning(self, "DB_file mismatch", f"DB_file must contain {len(runs_sorted)} comma-separated entries (one per run).")
            return
        if not path:
            QMessageBox.warning(self, "No path", "Please provide a save path")
            return

        # Try to use reduction_template_reader if available
        try:
            from lr_reduction import new_reduction_template_reader as rtr
        except Exception:
            rtr = None

        if rtr is not None:
            # If user provided a prior template path, attempt to update it at saved sequence indices
            prior = self.template_path_edit.text().strip()
            data_sets = []
            try:
                if prior and os.path.isfile(prior):
                    with open(prior, 'r') as pf:
                        prior_xml = pf.read()
                    data_sets = rtr.from_xml(prior_xml)
                else:
                    data_sets = []
            except Exception:
                data_sets = []

            # Build entries and place into data_sets at stored seq indices when available
            # iterate with index so we can pull per-run DB_file if provided
            for idx, runnum in enumerate(runs_sorted):
                r = self.per_run_rois[runnum]
                try:
                    ent = rtr.ReductionParameters()
                    ent.data_peak_range = [r['ymin'], r['ymax']]
                    # compute background_roi according to the per-run background mode
                    try:
                        mode_val = (r.get('bkg_mode') or '').lower()
                    except Exception:
                        mode_val = ''
                    if mode_val.startswith('no'):
                        # no background: collapse to peak so subtract_background=False
                        ent.background_roi = [r['ymin'], r['ymin'], r['ymax'], r['ymax']]
                        ent.subtract_background = False
                        ent.two_backgrounds = False
                    elif '1' in mode_val or 'one' in mode_val:
                        # one background: treat bkg1 as outer limits; background regions are
                        # [outer_left, ymin-1] and [ymax+1, outer_right]
                        outer_left = r['bkg1'][0]
                        outer_right = r['bkg1'][1]
                        n_y = getattr(self, 'n_y', None)
                        b1min = outer_left
                        b1max = max(0, r['ymin'])
                        if n_y is not None:
                            b2min = min(n_y - 1, r['ymax'])
                        else:
                            b2min = r['ymax']
                        b2max = outer_right
                        ent.background_roi = [b1min, b1max, b2min, b2max]
                        ent.subtract_background = True
                        ent.two_backgrounds = False
                    else:
                        # use 2 backgrounds: explicit ranges
                        ent.background_roi = [r['bkg1'][0], r['bkg1'][1], r['bkg2'][0], r['bkg2'][1]]
                    ent.tof_range = [r['tofmin'], r['tofmax']]
                    ent.data_x_range = [r['xmin'], r['xmax']]
                    ent.data_files = [int(runnum)]
                    # prefer per-run q_method if provided via the per-run editor
                    q_methods = vals.get('q_methods') if isinstance(vals, dict) else None
                    try:
                        if q_methods and len(q_methods) > idx and q_methods[idx]:
                            ent.q_method = q_methods[idx]
                        else:
                            ent.q_method = vals.get('q_method')
                    except Exception:
                        ent.q_method = vals.get('q_method')
                    # assign per-run DB_file if the user supplied a comma-separated list
                    try:
                        if db_list:
                            ent.DB_file = db_list[idx]
                        else:
                            ent.DB_file = vals.get('DB_file')
                    except Exception:
                        ent.DB_file = vals.get('DB_file')
                    ent.autoscale = vals.get('autoscale')
                    ent.use_calc_theta = vals.get('use_calc_theta')
                    # set background flags from UI mode
                    # set background flags from per-run stored mode
                    try:
                        mode = r.get('bkg_mode') or self.bkg_mode_cb.currentText()
                        if mode is not None and str(mode).lower().startswith('no'):
                            ent.subtract_background = False
                            ent.two_backgrounds = False
                        else:
                            # if 'Use 1 background' selected, set two_backgrounds False
                            if '1' in str(mode) or 'one' in str(mode).lower():
                                ent.subtract_background = True
                                ent.two_backgrounds = False
                            else:
                                ent.subtract_background = True
                                ent.two_backgrounds = True
                    except Exception:
                        pass
                    seq = r.get('seq')
                    if seq is not None and isinstance(seq, int) and seq >= 1:
                        # ensure list long enough
                        while len(data_sets) < seq:
                            data_sets.append(rtr.ReductionParameters())
                        data_sets[seq - 1] = ent
                    else:
                        data_sets.append(ent)
                except Exception:
                    continue

            try:
                out_xml = rtr.to_xml(data_sets)
                with open(path, 'w') as fh:
                    fh.write(out_xml)
                # remember last saved folder for this session
                try:
                    self._last_saved_template_dir = os.path.dirname(path)
                except Exception:
                    pass
                QMessageBox.information(self, "Saved", f"Combined template saved to {path}")
                return
            except Exception as e:
                QMessageBox.warning(self, "Save failed", f"Failed to write combined template via template reader: {e}\nFalling back to simple XML")

        # Fallback: write simple XML with multiple Reduction entries
        try:
            xml = '<Reductions>\n'
            for idx, runnum in enumerate(runs_sorted):
                r = self.per_run_rois[runnum]
                # determine background flags from per-run stored mode (best-effort)
                try:
                    mode = r.get('bkg_mode') or self.bkg_mode_cb.currentText()
                except Exception:
                    mode = None
                m = str(mode).lower() if mode is not None else ''
                bg_flag = 'True' if (m and not m.startswith('no')) or mode is None else 'False'
                # two_backgrounds is False for 'Use 1 background'
                if '1' in m or 'one' in m:
                    two_bg = 'False'
                else:
                    two_bg = 'True' if bg_flag == 'True' else 'False'
                xml += f'  <Reduction run="{runnum}">\n'
                xml += f'    <from_peak_pixels>{r["ymin"]}</from_peak_pixels>\n'
                xml += f'    <to_peak_pixels>{r["ymax"]}</to_peak_pixels>\n'
                xml += f'    <from_pixel>{r["xmin"]}</from_pixel>\n'
                xml += f'    <to_pixel>{r["xmax"]}</to_pixel>\n'
                xml += f'    <from_tof_range>{r["tofmin"]}</from_tof_range>\n'
                xml += f'    <to_tof_range>{r["tofmax"]}</to_tof_range>\n'
                xml += f'    <background_flag>{bg_flag}</background_flag>\n'
                xml += f'    <two_backgrounds>{two_bg}</two_backgrounds>\n'
                # compute back ROI values according to background mode semantics
                if m.startswith('no'):
                    b1min = r['ymin']
                    b1max = r['ymin']
                    b2min = r['ymax']
                    b2max = r['ymax']
                elif '1' in m or 'one' in m:
                    outer_left = r['bkg1'][0]
                    outer_right = r['bkg1'][1]
                    n_y = getattr(self, 'n_y', None)
                    b1min = outer_left
                    b1max = max(0, r['ymin'] - 1)
                    if n_y is not None:
                        b2min = min(n_y - 1, r['ymax'] + 1)
                    else:
                        b2min = r['ymax'] + 1
                    b2max = outer_right
                else:
                    b1min = r['bkg1'][0]
                    b1max = r['bkg1'][1]
                    b2min = r['bkg2'][0]
                    b2max = r['bkg2'][1]
                xml += f'    <back_roi1_from>{b1min}</back_roi1_from>\n'
                xml += f'    <back_roi1_to>{b1max}</back_roi1_to>\n'
                xml += f'    <back_roi2_from>{b2min}</back_roi2_from>\n'
                xml += f'    <back_roi2_to>{b2max}</back_roi2_to>\n'
                # per-run DB_file if supplied as comma-separated list
                try:
                    db_val = db_list[idx] if db_list else vals.get('DB_file', '')
                except Exception:
                    db_val = vals.get('DB_file', '')
                # include per-run q_method if available
                try:
                    q_methods = vals.get('q_methods') if isinstance(vals, dict) else None
                    if q_methods and len(q_methods) > idx and q_methods[idx]:
                        q_val = q_methods[idx]
                    else:
                        q_val = vals.get('q_method')
                except Exception:
                    q_val = vals.get('q_method')
                xml += f'    <q_method>{q_val}</q_method>\n'
                xml += f'    <DB_file>{db_val}</DB_file>\n'
                xml += '  </Reduction>\n'
            xml += '</Reductions>\n'
            with open(path, 'w') as fh:
                fh.write(xml)
            # remember last saved folder for this session
            try:
                self._last_saved_template_dir = os.path.dirname(path)
            except Exception:
                pass
            QMessageBox.information(self, "Saved", f"Combined template saved to {path}")
        except Exception as e:
            QMessageBox.critical(self, "Save error", f"Failed to save combined template: {e}")



