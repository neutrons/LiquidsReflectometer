"""Template-based reduction tab for the reflectometry launcher."""

import logging
import os
import sys

from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QGridLayout, QGroupBox, QHBoxLayout,
    QLabel, QLineEdit, QMessageBox, QProgressBar, QPushButton, QScrollArea,
    QSizePolicy, QVBoxLayout, QWidget,
)

logger = logging.getLogger(__name__)

# Map display names to the lowercase values that NRReductionConfig expects
_METHOD_DISPLAY_TO_INTERNAL = {
    "meanTheta": "meantheta",
    "constantQ": "constantq",
    "constantTOF": "constanttof",
}


def _settings_bool(value):
    """Convert a QSettings value (may be bool, str, or None) to bool."""
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() == "true"
    return False


def _safe_float(text):
    """Convert text to float, returning None if the text is not a valid number."""
    text = text.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _safe_int(text):
    """Convert text to int, returning None if the text is not a valid integer."""
    text = text.strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        # Handle intermediate validator states like "1." or "-"
        f = _safe_float(text)
        if f is not None:
            return int(f)
        return None


TEMPLATE_DIRECTIVE = "Click to choose a template file (*.xml)"
OUTPUT_DIR_DIRECTIVE = "Click to choose an output directory"
DATA_PATH_DIRECTIVE = "Click to choose a NEXUS data directory"
DB_PATH_DIRECTIVE = "Click to choose a DB directory"
DB_FILE_DIRECTIVE = "Click to choose a DB file (*.dat)"
TEMPLATE_SAVE_DIRECTIVE = "Click to choose a template save directory"


def _parse_run_numbers(text):
    """Parse run number text into a list of integers.

    Supports single numbers, comma-separated lists, and ranges (e.g. "211029-211031").
    Mixed formats like "211029,211033-211035" are also supported.

    Raises ValueError for empty or invalid input.
    """
    text = text.strip()
    if not text:
        raise ValueError("Run number text is empty")

    result = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            tokens = part.split("-")
            if len(tokens) != 2:
                raise ValueError(f"Invalid range: {part}")
            try:
                start = int(tokens[0].strip())
                end = int(tokens[1].strip())
            except ValueError:
                raise ValueError(f"Invalid range: {part}")
            if end < start:
                raise ValueError(f"Invalid range (end < start): {part}")
            result.extend(range(start, end + 1))
        else:
            try:
                result.append(int(part))
            except ValueError:
                raise ValueError(f"Invalid run number: {part}")

    if not result:
        raise ValueError("No valid run numbers found")
    return result


class TemplateReduce(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle("Template reduction")
        self.settings = QtCore.QSettings()
        self._worker = None

        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        # === Top section: file/path inputs ===
        top_grid = QGridLayout()
        main_layout.addLayout(top_grid)

        # Row 0: Run numbers
        self.run_numbers_edit = QLineEdit()
        top_grid.addWidget(self.run_numbers_edit, 0, 0)
        run_label = QLabel("Run number(s) — comma-separated or range (e.g. 211029,211030 or 211029-211031)")
        top_grid.addWidget(run_label, 0, 1)

        # Row 1: Template file
        self.choose_template_btn = QPushButton("Template file")
        top_grid.addWidget(self.choose_template_btn, 1, 0)
        self.template_path_label = QLabel(TEMPLATE_DIRECTIVE)
        top_grid.addWidget(self.template_path_label, 1, 1)

        # Row 2: Experiment ID
        self.experiment_id_edit = QLineEdit()
        top_grid.addWidget(self.experiment_id_edit, 2, 0)
        exp_label = QLabel("Experiment ID (e.g. IPTS-36119)")
        top_grid.addWidget(exp_label, 2, 1)

        # Row 3: Output directory
        self.choose_output_btn = QPushButton("Output directory")
        top_grid.addWidget(self.choose_output_btn, 3, 0)
        self.output_dir_label = QLabel(OUTPUT_DIR_DIRECTIVE)
        top_grid.addWidget(self.output_dir_label, 3, 1)

        # Row 4: Data path (NEXUS)
        self.choose_data_path_btn = QPushButton("Data path")
        top_grid.addWidget(self.choose_data_path_btn, 4, 0)
        self.data_path_label = QLabel(DATA_PATH_DIRECTIVE)
        top_grid.addWidget(self.data_path_label, 4, 1)

        # Row 5: DB path
        self.choose_db_path_btn = QPushButton("DB path")
        top_grid.addWidget(self.choose_db_path_btn, 5, 0)
        self.db_path_label = QLabel(DB_PATH_DIRECTIVE)
        top_grid.addWidget(self.db_path_label, 5, 1)

        # Row 6: DB file name
        self.choose_db_file_btn = QPushButton("DB file name")
        top_grid.addWidget(self.choose_db_file_btn, 6, 0)
        self.db_file_label = QLabel(DB_FILE_DIRECTIVE)
        top_grid.addWidget(self.db_file_label, 6, 1)

        # Row 7: Template save path
        self.choose_template_save_btn = QPushButton("Template save path")
        top_grid.addWidget(self.choose_template_save_btn, 7, 0)
        self.template_save_label = QLabel(TEMPLATE_SAVE_DIRECTIVE)
        top_grid.addWidget(self.template_save_label, 7, 1)

        # === Middle section: scrollable options ===
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        scroll.setWidget(scroll_widget)
        main_layout.addWidget(scroll)

        # Processing Flags group
        flags_group = QGroupBox("Processing Flags")
        flags_layout = QGridLayout()
        flags_group.setLayout(flags_layout)

        flags_layout.addWidget(QLabel("Method:"), 0, 0)
        self.method_combo = QComboBox()
        self.method_combo.addItems(["meanTheta", "constantQ", "constantTOF"])
        flags_layout.addWidget(self.method_combo, 0, 1)

        self.normalize_check = QCheckBox("Normalize")
        flags_layout.addWidget(self.normalize_check, 1, 0)
        self.autoscale_check = QCheckBox("AutoScale")
        flags_layout.addWidget(self.autoscale_check, 1, 1)
        self.use_calc_theta_check = QCheckBox("useCalcTheta")
        flags_layout.addWidget(self.use_calc_theta_check, 1, 2)

        scroll_layout.addWidget(flags_group)

        # Plot Options group
        plot_group = QGroupBox("Plot Options")
        plot_layout = QGridLayout()
        plot_group.setLayout(plot_layout)

        self.save_plots_check = QCheckBox("Save plots")
        plot_layout.addWidget(self.save_plots_check, 0, 0)
        self.interact_plots_check = QCheckBox("Interact with plots")
        plot_layout.addWidget(self.interact_plots_check, 0, 1)
        self.plot_rq4_check = QCheckBox("Plot as RQ4")
        plot_layout.addWidget(self.plot_rq4_check, 0, 2)

        scroll_layout.addWidget(plot_group)

        # Q-space Parameters group
        qspace_group = QGroupBox("Q-space Parameters")
        qspace_layout = QGridLayout()
        qspace_group.setLayout(qspace_layout)

        qspace_layout.addWidget(QLabel("qmin:"), 0, 0)
        self.qmin_edit = QLineEdit()
        self.qmin_edit.setValidator(QtGui.QDoubleValidator())
        qspace_layout.addWidget(self.qmin_edit, 0, 1)

        qspace_layout.addWidget(QLabel("qmax:"), 0, 2)
        self.qmax_edit = QLineEdit()
        self.qmax_edit.setValidator(QtGui.QDoubleValidator())
        qspace_layout.addWidget(self.qmax_edit, 0, 3)

        qspace_layout.addWidget(QLabel("dqbin:"), 1, 0)
        self.dqbin_edit = QLineEdit()
        self.dqbin_edit.setValidator(QtGui.QDoubleValidator())
        qspace_layout.addWidget(self.dqbin_edit, 1, 1)

        qspace_layout.addWidget(QLabel("Qline_threshold:"), 1, 2)
        self.qline_threshold_edit = QLineEdit()
        self.qline_threshold_edit.setValidator(QtGui.QDoubleValidator())
        qspace_layout.addWidget(self.qline_threshold_edit, 1, 3)

        scroll_layout.addWidget(qspace_group)

        # Dead Time group
        dt_group = QGroupBox("Dead Time")
        dt_layout = QGridLayout()
        dt_group.setLayout(dt_layout)

        dt_layout.addWidget(QLabel("dead_time:"), 0, 0)
        self.dead_time_edit = QLineEdit()
        self.dead_time_edit.setValidator(QtGui.QDoubleValidator())
        dt_layout.addWidget(self.dead_time_edit, 0, 1)

        dt_layout.addWidget(QLabel("dead_time_tof_step:"), 0, 2)
        self.dead_time_tof_step_edit = QLineEdit()
        self.dead_time_tof_step_edit.setValidator(QtGui.QDoubleValidator())
        dt_layout.addWidget(self.dead_time_tof_step_edit, 0, 3)

        scroll_layout.addWidget(dt_group)

        # Detector / Peak group
        det_group = QGroupBox("Detector / Peak")
        det_layout = QGridLayout()
        det_group.setLayout(det_layout)

        det_layout.addWidget(QLabel("DetResFn:"), 0, 0)
        self.det_res_fn_combo = QComboBox()
        self.det_res_fn_combo.addItems(["rectangular", "gaussian"])
        det_layout.addWidget(self.det_res_fn_combo, 0, 1)

        det_layout.addWidget(QLabel("DetSigma:"), 0, 2)
        self.det_sigma_edit = QLineEdit()
        self.det_sigma_edit.setValidator(QtGui.QDoubleValidator())
        det_layout.addWidget(self.det_sigma_edit, 0, 3)

        det_layout.addWidget(QLabel("peak_pad:"), 1, 0)
        self.peak_pad_edit = QLineEdit()
        self.peak_pad_edit.setValidator(QtGui.QIntValidator())
        det_layout.addWidget(self.peak_pad_edit, 1, 1)

        det_layout.addWidget(QLabel("peak_type:"), 1, 2)
        self.peak_type_combo = QComboBox()
        self.peak_type_combo.addItems(["supergauss", "gauss"])
        det_layout.addWidget(self.peak_type_combo, 1, 3)

        scroll_layout.addWidget(det_group)

        # Emission Time group
        emission_group = QGroupBox("Emission Time")
        emission_layout = QGridLayout()
        emission_group.setLayout(emission_layout)

        self.use_emission_time_check = QCheckBox("use_emission_time")
        self.use_emission_time_check.setChecked(True)
        emission_layout.addWidget(self.use_emission_time_check, 0, 0)

        emission_layout.addWidget(QLabel("IncidentTheta:"), 0, 1)
        self.incident_theta_edit = QLineEdit()
        self.incident_theta_edit.setValidator(QtGui.QDoubleValidator())
        emission_layout.addWidget(self.incident_theta_edit, 0, 2)

        scroll_layout.addWidget(emission_group)

        # Instrument Geometry group
        geom_group = QGroupBox("Instrument Geometry (blank = use defaults from settings.json)")
        geom_layout = QGridLayout()
        geom_group.setLayout(geom_layout)

        geom_fields = [
            ("mmpix", 0, 0), ("dSampDet", 0, 2),
            ("ny", 1, 0), ("nx", 1, 2),
            ("dMod", 2, 0), ("xi_ref", 2, 2),
            ("dS1Samp", 3, 0),
        ]
        self.geom_edits = {}
        for name, row, col in geom_fields:
            geom_layout.addWidget(QLabel(f"{name}:"), row, col)
            edit = QLineEdit()
            edit.setValidator(QtGui.QDoubleValidator())
            geom_layout.addWidget(edit, row, col + 1)
            self.geom_edits[name] = edit

        scroll_layout.addWidget(geom_group)

        # Spacer at bottom of scroll area
        scroll_layout.addStretch()

        # === Bottom section: reduce button, progress ===
        bottom_layout = QHBoxLayout()
        main_layout.addLayout(bottom_layout)

        self.reduce_btn = QPushButton("REDUCE")
        self.reduce_btn.setStyleSheet("background-color: green; font-weight: bold;")
        self.reduce_btn.setMinimumHeight(40)
        bottom_layout.addWidget(self.reduce_btn)

        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        bottom_layout.addWidget(self.cancel_btn)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        main_layout.addWidget(self.progress_bar)

        self.status_label = QLabel("Ready")
        main_layout.addWidget(self.status_label)

        # === Log to stdout checkbox ===
        self.log_stdout_check = QCheckBox("Log to stdout")
        main_layout.addWidget(self.log_stdout_check)

        # === Connections ===
        self.choose_template_btn.clicked.connect(self.template_selection)
        self.choose_output_btn.clicked.connect(self.output_dir_selection)
        self.choose_data_path_btn.clicked.connect(self.data_path_selection)
        self.choose_db_path_btn.clicked.connect(self.db_path_selection)
        self.choose_db_file_btn.clicked.connect(self.db_file_selection)
        self.choose_template_save_btn.clicked.connect(self.template_save_selection)
        self.reduce_btn.clicked.connect(self.reduce)
        self.cancel_btn.clicked.connect(self.on_cancel)

        # Populate from previous session
        self.read_settings()

    # --- File/folder selection dialogs ---
    def template_selection(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select template file", self.template_path_label.text(), "Template file (*.xml)")
        if path and os.path.isfile(path):
            self.template_path_label.setText(path)

    def output_dir_selection(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select output directory", self.output_dir_label.text(), QFileDialog.ShowDirsOnly)
        if path and os.path.isdir(path):
            self.output_dir_label.setText(path)

    def data_path_selection(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select NEXUS data directory", self.data_path_label.text(), QFileDialog.ShowDirsOnly)
        if path and os.path.isdir(path):
            self.data_path_label.setText(path)

    def db_path_selection(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select DB directory", self.db_path_label.text(), QFileDialog.ShowDirsOnly)
        if path and os.path.isdir(path):
            self.db_path_label.setText(path)

    def db_file_selection(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select DB file", self.db_file_label.text(), "DB file (*.dat)")
        if path and os.path.isfile(path):
            self.db_file_label.setText(path)

    def template_save_selection(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select template save directory", self.template_save_label.text(), QFileDialog.ShowDirsOnly)
        if path and os.path.isdir(path):
            self.template_save_label.setText(path)

    # --- Settings persistence ---
    def read_settings(self):
        self.run_numbers_edit.setText(self.settings.value("tmpl_run_numbers", ""))
        self.template_path_label.setText(self.settings.value("tmpl_template_file", TEMPLATE_DIRECTIVE))
        self.experiment_id_edit.setText(self.settings.value("tmpl_experiment_id", ""))
        self.output_dir_label.setText(self.settings.value("tmpl_output_dir", OUTPUT_DIR_DIRECTIVE))
        self.data_path_label.setText(self.settings.value("tmpl_data_path", DATA_PATH_DIRECTIVE))
        self.db_path_label.setText(self.settings.value("tmpl_db_path", DB_PATH_DIRECTIVE))
        self.db_file_label.setText(self.settings.value("tmpl_db_file", DB_FILE_DIRECTIVE))
        self.template_save_label.setText(self.settings.value("tmpl_template_save", TEMPLATE_SAVE_DIRECTIVE))

        self.method_combo.setCurrentText(self.settings.value("tmpl_method", "meanTheta"))
        self.normalize_check.setChecked(_settings_bool(self.settings.value("tmpl_normalize", False)))
        self.autoscale_check.setChecked(_settings_bool(self.settings.value("tmpl_autoscale", False)))
        self.use_calc_theta_check.setChecked(_settings_bool(self.settings.value("tmpl_use_calc_theta", False)))

        self.save_plots_check.setChecked(_settings_bool(self.settings.value("tmpl_save_plots", False)))
        self.interact_plots_check.setChecked(_settings_bool(self.settings.value("tmpl_interact_plots", False)))
        self.plot_rq4_check.setChecked(_settings_bool(self.settings.value("tmpl_plot_rq4", False)))

        self.qmin_edit.setText(self.settings.value("tmpl_qmin", ""))
        self.qmax_edit.setText(self.settings.value("tmpl_qmax", ""))
        self.dqbin_edit.setText(self.settings.value("tmpl_dqbin", ""))
        self.qline_threshold_edit.setText(self.settings.value("tmpl_qline_threshold", ""))

        self.dead_time_edit.setText(self.settings.value("tmpl_dead_time", ""))
        self.dead_time_tof_step_edit.setText(self.settings.value("tmpl_dead_time_tof_step", ""))

        self.det_res_fn_combo.setCurrentText(self.settings.value("tmpl_det_res_fn", "rectangular"))
        self.det_sigma_edit.setText(self.settings.value("tmpl_det_sigma", ""))
        self.peak_pad_edit.setText(self.settings.value("tmpl_peak_pad", ""))
        self.peak_type_combo.setCurrentText(self.settings.value("tmpl_peak_type", "supergauss"))

        self.use_emission_time_check.setChecked(_settings_bool(self.settings.value("tmpl_use_emission_time", True)))
        self.incident_theta_edit.setText(self.settings.value("tmpl_incident_theta", ""))

        self.log_stdout_check.setChecked(_settings_bool(self.settings.value("tmpl_log_stdout", False)))

        for name, edit in self.geom_edits.items():
            edit.setText(self.settings.value(f"tmpl_geom_{name}", ""))

    def save_settings(self):
        self.settings.setValue("tmpl_run_numbers", self.run_numbers_edit.text())
        self.settings.setValue("tmpl_template_file", self.template_path_label.text())
        self.settings.setValue("tmpl_experiment_id", self.experiment_id_edit.text())
        self.settings.setValue("tmpl_output_dir", self.output_dir_label.text())
        self.settings.setValue("tmpl_data_path", self.data_path_label.text())
        self.settings.setValue("tmpl_db_path", self.db_path_label.text())
        self.settings.setValue("tmpl_db_file", self.db_file_label.text())
        self.settings.setValue("tmpl_template_save", self.template_save_label.text())

        self.settings.setValue("tmpl_method", self.method_combo.currentText())
        self.settings.setValue("tmpl_normalize", self.normalize_check.isChecked())
        self.settings.setValue("tmpl_autoscale", self.autoscale_check.isChecked())
        self.settings.setValue("tmpl_use_calc_theta", self.use_calc_theta_check.isChecked())

        self.settings.setValue("tmpl_save_plots", self.save_plots_check.isChecked())
        self.settings.setValue("tmpl_interact_plots", self.interact_plots_check.isChecked())
        self.settings.setValue("tmpl_plot_rq4", self.plot_rq4_check.isChecked())

        self.settings.setValue("tmpl_qmin", self.qmin_edit.text())
        self.settings.setValue("tmpl_qmax", self.qmax_edit.text())
        self.settings.setValue("tmpl_dqbin", self.dqbin_edit.text())
        self.settings.setValue("tmpl_qline_threshold", self.qline_threshold_edit.text())

        self.settings.setValue("tmpl_dead_time", self.dead_time_edit.text())
        self.settings.setValue("tmpl_dead_time_tof_step", self.dead_time_tof_step_edit.text())

        self.settings.setValue("tmpl_det_res_fn", self.det_res_fn_combo.currentText())
        self.settings.setValue("tmpl_det_sigma", self.det_sigma_edit.text())
        self.settings.setValue("tmpl_peak_pad", self.peak_pad_edit.text())
        self.settings.setValue("tmpl_peak_type", self.peak_type_combo.currentText())

        self.settings.setValue("tmpl_use_emission_time", self.use_emission_time_check.isChecked())
        self.settings.setValue("tmpl_incident_theta", self.incident_theta_edit.text())

        self.settings.setValue("tmpl_log_stdout", self.log_stdout_check.isChecked())

        for name, edit in self.geom_edits.items():
            self.settings.setValue(f"tmpl_geom_{name}", edit.text())

    # --- Input validation ---
    def check_inputs(self):
        error = None

        # Validate run numbers
        try:
            _parse_run_numbers(self.run_numbers_edit.text())
        except ValueError as e:
            error = f"Invalid run numbers: {e}"

        # Check template file
        if error is None and not os.path.isfile(self.template_path_label.text()):
            error = "The chosen template file could not be found"

        # Check experiment ID
        if error is None and not self.experiment_id_edit.text().strip():
            error = "Experiment ID is required"

        # Check output directory
        if error is None:
            out_dir = self.output_dir_label.text()
            if out_dir == OUTPUT_DIR_DIRECTIVE or not os.path.isdir(out_dir):
                error = "Please select a valid output directory"

        # Validate numeric fields that have text in them
        if error is None:
            error = self._validate_numeric_fields()

        if error:
            self.show_dialog(error)
            return False
        return True

    def _validate_numeric_fields(self):
        """Validate all numeric text fields. Returns error string or None."""
        float_fields = [
            (self.qmin_edit, "qmin"),
            (self.qmax_edit, "qmax"),
            (self.dqbin_edit, "dqbin"),
            (self.qline_threshold_edit, "Qline_threshold"),
            (self.dead_time_edit, "dead_time"),
            (self.dead_time_tof_step_edit, "dead_time_tof_step"),
            (self.det_sigma_edit, "DetSigma"),
            (self.incident_theta_edit, "IncidentTheta"),
        ]
        for edit, name in float_fields:
            text = edit.text().strip()
            if text and _safe_float(text) is None:
                return f"Invalid value for {name}: '{text}' is not a valid number"

        # Integer field
        text = self.peak_pad_edit.text().strip()
        if text and _safe_int(text) is None:
            return f"Invalid value for peak_pad: '{text}' is not a valid integer"

        # Geometry fields
        for name, edit in self.geom_edits.items():
            text = edit.text().strip()
            if text and _safe_float(text) is None:
                return f"Invalid value for {name}: '{text}' is not a valid number"

        return None

    def show_dialog(self, text, title="Invalid inputs", icon=None):
        if icon is None:
            icon = QMessageBox.Critical
        msgBox = QMessageBox()
        msgBox.setIcon(icon)
        msgBox.setText(text)
        msgBox.setWindowTitle(title)
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec_()

    # --- Build override params dict ---
    def _build_override_params(self):
        params = {}

        # Method: convert display name to lowercase internal value
        display_method = self.method_combo.currentText()
        params["method"] = _METHOD_DISPLAY_TO_INTERNAL.get(display_method, display_method.lower())

        # Boolean flags: always send current state to override template values
        params["Normalize"] = self.normalize_check.isChecked()
        params["AutoScale"] = self.autoscale_check.isChecked()
        params["useCalcTheta"] = self.use_calc_theta_check.isChecked()
        params["plotQ4"] = self.plot_rq4_check.isChecked()

        # Q-space params (only override if non-empty and valid)
        for attr, edit in [("qmin", self.qmin_edit), ("qmax", self.qmax_edit),
                           ("dqbin", self.dqbin_edit), ("Qline_threshold", self.qline_threshold_edit)]:
            val = _safe_float(edit.text())
            if val is not None:
                params[attr] = val

        # Dead time
        val = _safe_float(self.dead_time_edit.text())
        if val is not None:
            params["dead_time"] = val
        val = _safe_float(self.dead_time_tof_step_edit.text())
        if val is not None:
            params["dead_time_tof_step"] = val

        # Detector/Peak
        params["DetResFn"] = self.det_res_fn_combo.currentText()
        val = _safe_float(self.det_sigma_edit.text())
        if val is not None:
            params["DetSigma"] = val
        val = _safe_int(self.peak_pad_edit.text())
        if val is not None:
            params["peak_pad"] = val
        params["peak_type"] = self.peak_type_combo.currentText()

        # Emission time
        params["use_emission_time"] = self.use_emission_time_check.isChecked()
        val = _safe_float(self.incident_theta_edit.text())
        if val is not None:
            params["IncidentTheta"] = val

        # Instrument geometry (only override if non-empty and valid)
        for name, edit in self.geom_edits.items():
            val = _safe_float(edit.text())
            if val is not None:
                params[name] = val

        # Output path
        out_dir = self.output_dir_label.text()
        if out_dir != OUTPUT_DIR_DIRECTIVE and os.path.isdir(out_dir):
            params["Spath"] = out_dir

        # DB path and file
        db_path = self.db_path_label.text()
        if db_path != DB_PATH_DIRECTIVE and os.path.isdir(db_path):
            params["DBpath"] = db_path

        db_file = self.db_file_label.text()
        if db_file != DB_FILE_DIRECTIVE and os.path.isfile(db_file):
            params["DBname"] = [os.path.basename(db_file)]

        return params

    # --- Reduction ---
    def reduce(self):
        # Guard against double-click while already running
        if self._worker is not None and self._worker.isRunning():
            return

        if not self.check_inputs():
            return

        self.save_settings()

        run_numbers = _parse_run_numbers(self.run_numbers_edit.text())
        template_file = self.template_path_label.text()
        experiment_id = self.experiment_id_edit.text().strip()

        datapath = self.data_path_label.text()
        if datapath == DATA_PATH_DIRECTIVE or not os.path.isdir(datapath):
            datapath = None

        template_path = self.template_save_label.text()
        if template_path == TEMPLATE_SAVE_DIRECTIVE or not os.path.isdir(template_path):
            template_path = None

        override_params = self._build_override_params()

        save_plots = self.save_plots_check.isChecked()
        plot_dir = self.output_dir_label.text() if save_plots else None

        log_to_stdout = self.log_stdout_check.isChecked()

        # Log the parameters being sent
        if log_to_stdout:
            print("--- Template Reduction ---")
            print(f"  Run numbers: {run_numbers}")
            print(f"  Template: {template_file}")
            print(f"  Experiment: {experiment_id}")
            print(f"  Data path: {datapath}")
            print(f"  Template save: {template_path}")
            print(f"  Save plots: {save_plots}  Plot dir: {plot_dir}")
            print(f"  Override params: {override_params}")
            print("--------------------------")

        # Disable controls during reduction
        self.reduce_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.progress_bar.setValue(0)
        self.status_label.setText("Starting reduction...")

        from apps.reduction_worker import ReductionWorker

        self._worker = ReductionWorker(
            run_numbers=run_numbers,
            template_file=template_file,
            experiment_id=experiment_id,
            datapath=datapath,
            template_path=template_path,
            override_params=override_params,
            plot=False,
            save_plots=save_plots,
            plot_dir=plot_dir,
            log_to_stdout=log_to_stdout,
        )
        self._worker.progress.connect(self.on_progress)
        self._worker.finished.connect(self.on_finished)
        self._worker.error.connect(self.on_error)
        self._worker.start()

    def on_cancel(self):
        if self._worker:
            self._worker.cancel()
            self.status_label.setText("Cancelling...")
            self.cancel_btn.setEnabled(False)

    def on_progress(self, pct, total, desc):
        self.progress_bar.setValue(pct)
        self.status_label.setText(desc)

    def on_finished(self, results):
        self.reduce_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.progress_bar.setValue(100)

        # Check if reduction was cancelled
        if self._worker and self._worker._cancelled:
            self.status_label.setText("Cancelled")
            self.show_dialog("Reduction was cancelled.", "Cancelled", icon=QMessageBox.Warning)
            return

        self.status_label.setText("Complete")

        # Show interactive plots from main thread if requested
        if self.interact_plots_check.isChecked() and results:
            try:
                from lr_reduction.new_reduction_from_template import plot_reflectivity
                for result in results:
                    if result is not None:
                        plot_reflectivity([result], RQ4=self.plot_rq4_check.isChecked(), show=True)
            except Exception as e:
                self.show_dialog(f"Error showing interactive plots: {e}", "Plot Error")

        out_dir = self.output_dir_label.text()
        msg = f"Reduction completed successfully!\n\nOutput directory: {out_dir}"
        self.show_dialog(msg, "Task completed", icon=QMessageBox.Information)

    def on_error(self, error_msg):
        self.reduce_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.status_label.setText("Error")
        if self.log_stdout_check.isChecked():
            print(f"REDUCTION ERROR:\n{error_msg}", file=sys.stderr)
        self.show_dialog(f"Reduction failed:\n{error_msg}", "Reduction Error")
