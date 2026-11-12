#!/usr/bin/python3

from pathlib import Path

from matplotlib.figure import Figure
from qtpy import QtCore
from qtpy.QtWidgets import (
    QComboBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

try:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
except ImportError:
    try:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
    except ImportError:
        FigureCanvas = None
        NavigationToolbar = None

from launcher.app_identity import ensure_identity
from launcher.apps.settings_editor import guarded
from lr_reduction.new_reduction_time_resolved import reduce_time_list, reduce_time_slices


class TimeResolvedTab(QWidget):
    """Simple tab for time-resolved reflectivity reduction.

    As the launcher's other tabs do, it reports in its own panel (the log area): every input it cannot use, and every
    error a reduction raises (the slots are guarded), never in a modal box. Reduce is disabled while a reduction runs.
    """

    def __init__(self):
        # Before any QSettings is constructed: QSettings derives its path from the application identity (as the
        # settings editor does).
        ensure_identity()
        super().__init__()
        self.setWindowTitle("Time-resolved reduction")
        self.settings = QtCore.QSettings()

        layout = QVBoxLayout(self)

        form = QFormLayout()

        self.experiment_edit = QLineEdit(self)
        self.experiment_edit.setPlaceholderText("IPTS-36776")
        form.addRow("Experiment ID:", self.experiment_edit)

        self.run_edit = QLineEdit(self)
        self.run_edit.setPlaceholderText("e.g. 227164")
        form.addRow("Run number:", self.run_edit)

        self.settings_edit = QLineEdit(self)
        self.settings_edit.setPlaceholderText("/SNS/.../REFL_227158_settings.json")
        settings_row = QHBoxLayout()
        settings_row.addWidget(self.settings_edit)
        settings_btn = QPushButton("Browse")
        settings_btn.clicked.connect(self._browse_settings)
        settings_row.addWidget(settings_btn)
        form.addRow("Settings file:", settings_row)

        self.output_edit = QLineEdit(self)
        self.output_edit.setPlaceholderText("Optional output directory")
        output_row = QHBoxLayout()
        output_row.addWidget(self.output_edit)
        output_btn = QPushButton("Browse")
        output_btn.clicked.connect(self._browse_output)
        output_row.addWidget(output_btn)
        form.addRow("Output dir:", output_row)

        self.subname_edit = QLineEdit(self)
        self.subname_edit.setPlaceholderText("optional label")
        form.addRow("Subname:", self.subname_edit)

        self.mode_combo = QComboBox(self)
        self.mode_combo.addItems(["Number of slices", "Time values"])
        form.addRow("Mode:", self.mode_combo)

        self.num_slices_spin = QSpinBox(self)
        self.num_slices_spin.setMinimum(2)
        self.num_slices_spin.setValue(2)
        form.addRow("Number of slices:", self.num_slices_spin)

        self.start_times_edit = QLineEdit(self)
        self.start_times_edit.setPlaceholderText("0, 80, 120")
        form.addRow("Start times:", self.start_times_edit)

        self.end_times_edit = QLineEdit(self)
        self.end_times_edit.setPlaceholderText("80, 120, 136")
        form.addRow("End times:", self.end_times_edit)

        layout.addLayout(form)

        self.mode_combo.currentIndexChanged.connect(self._update_mode_controls)
        self._update_mode_controls()

        buttons = QHBoxLayout()
        self.process_btn = QPushButton("Reduce")
        self.process_btn.clicked.connect(lambda _checked=False: self._run_reduction())
        buttons.addWidget(self.process_btn)

        self.save_fig_btn = QPushButton("Save figure")
        self.save_fig_btn.clicked.connect(lambda _checked=False: self._save_figure())
        buttons.addWidget(self.save_fig_btn)
        layout.addLayout(buttons)

        self.log_edit = QTextEdit(self)
        self.log_edit.setReadOnly(True)
        self.log_edit.setMinimumHeight(100)
        layout.addWidget(self.log_edit)

        plot_layout = QVBoxLayout()
        if FigureCanvas is not None:
            self.figure = Figure(figsize=(6, 4))
            self.canvas = FigureCanvas(self.figure)
            self.toolbar = NavigationToolbar(self.canvas, self)

            plot_layout.addWidget(self.toolbar)
            plot_layout.addWidget(self.canvas)
        else:
            self.figure = None
            self.canvas = None
            self.toolbar = None
            plot_layout.addWidget(QLabel("Embedded plotting unavailable"))
        layout.addLayout(plot_layout)

        self.read_settings()

    def read_settings(self):
        self.experiment_edit.setText(self.settings.value("time_resolved_experiment_id", ""))
        self.run_edit.setText(self.settings.value("time_resolved_run_number", ""))
        self.settings_edit.setText(self.settings.value("time_resolved_settings_file", ""))
        self.output_edit.setText(self.settings.value("time_resolved_output_dir", ""))
        self.subname_edit.setText(self.settings.value("time_resolved_subname", ""))

        mode = self.settings.value("time_resolved_mode", "Number of slices")
        idx = self.mode_combo.findText(str(mode))
        if idx >= 0:
            self.mode_combo.setCurrentIndex(idx)

        try:
            self.num_slices_spin.setValue(int(self.settings.value("time_resolved_num_slices", self.num_slices_spin.value())))
        except (TypeError, ValueError):  # a stored value that is not a number: keep the default
            pass

        self.start_times_edit.setText(self.settings.value("time_resolved_start_times", ""))
        self.end_times_edit.setText(self.settings.value("time_resolved_end_times", ""))
        self._update_mode_controls()

    def save_settings(self):
        self.settings.setValue("time_resolved_experiment_id", self.experiment_edit.text())
        self.settings.setValue("time_resolved_run_number", self.run_edit.text())
        self.settings.setValue("time_resolved_settings_file", self.settings_edit.text())
        self.settings.setValue("time_resolved_output_dir", self.output_edit.text())
        self.settings.setValue("time_resolved_subname", self.subname_edit.text())
        self.settings.setValue("time_resolved_mode", self.mode_combo.currentText())
        self.settings.setValue("time_resolved_num_slices", int(self.num_slices_spin.value()))
        self.settings.setValue("time_resolved_start_times", self.start_times_edit.text())
        self.settings.setValue("time_resolved_end_times", self.end_times_edit.text())

    def _update_mode_controls(self):
        use_slices = self.mode_combo.currentText() == "Number of slices"
        self.num_slices_spin.setEnabled(use_slices)
        self.start_times_edit.setEnabled(not use_slices)
        self.end_times_edit.setEnabled(not use_slices)

    def _browse_settings(self):
        start = str(Path.home())
        path, _ = QFileDialog.getOpenFileName(self, "Select settings file", start, "JSON files (*.json)")
        if path:
            self.settings_edit.setText(path)

    def _browse_output(self):
        start = str(Path.home())
        path = QFileDialog.getExistingDirectory(self, "Select output directory", start)
        if path:
            self.output_edit.setText(path)

    def report_problem(self, exc):
        """Show a failure in the panel instead of letting it leave the slot (launcher.apps.settings_editor.guarded)."""
        self.log_edit.append(f"ERROR: {type(exc).__name__}: {exc}")

    @guarded
    def _save_figure(self):
        if self.figure is None or not self.figure.axes:
            self.log_edit.append("There is no figure to save yet.")
            return

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save figure",
            str(Path.home()),
            "PNG files (*.png);;PDF files (*.pdf);;SVG files (*.svg)",
        )
        if path:
            self.figure.savefig(path, dpi=200)
            self.log_edit.append(f"Saved figure to {path}")

    def _parse_float_list(self, text):
        if text is None or not str(text).strip():
            return []
        parts = [p.strip() for p in str(text).split(",") if p.strip()]
        return [float(p) for p in parts]

    def _request(self):
        """The reduction the inputs ask for, or None once the panel has said what is missing or malformed."""
        run_text = self.run_edit.text().strip()
        if not run_text:
            self.log_edit.append("Please enter a run number.")
            return None

        try:
            run_number = int(float(run_text))
        except ValueError:
            self.log_edit.append("Run number must be an integer.")
            return None

        experiment_id = self.experiment_edit.text().strip()
        if not experiment_id:
            self.log_edit.append("Please enter an experiment ID, e.g. IPTS-36776.")
            return None

        settings_file = self.settings_edit.text().strip()
        if not settings_file:
            self.log_edit.append("Please choose a settings JSON file.")
            return None

        request = {"run": run_number, "experiment_id": experiment_id, "settings_file": settings_file,
                   "mode": self.mode_combo.currentText()}
        if request["mode"] != "Number of slices":
            starts = self._parse_float_list(self.start_times_edit.text())
            ends = self._parse_float_list(self.end_times_edit.text())
            if len(starts) == 0 or len(ends) == 0:
                self.log_edit.append("Please provide at least one start and end time.")
                return None
            if len(starts) != len(ends):
                self.log_edit.append("Start times and end times must have the same length.")
                return None
            request["starts"], request["ends"] = starts, ends
        return request

    @guarded
    def _run_reduction(self):
        self.log_edit.clear()
        request = self._request()
        if request is None:
            return

        savepath = self.output_edit.text().strip() or None
        if savepath is not None:
            savepath = Path(savepath)

        subname = self.subname_edit.text().strip() or None

        self.save_settings()

        run_number, experiment_id, settings_file = request["run"], request["experiment_id"], request["settings_file"]
        self.log_edit.append(f"Starting time-resolved reduction for run {run_number} in {experiment_id}")

        # Synchronous on the GUI thread (§10 A4): Reduce is unavailable until it returns or raises.
        self.process_btn.setEnabled(False)
        try:
            if request["mode"] == "Number of slices":
                n = self.num_slices_spin.value()
                self.log_edit.append(f"Mode: split into {n} equal time slices")
                output, plots = reduce_time_slices(
                    run_number,
                    settings_file,
                    experiment_id,
                    n,
                    savepath=savepath,
                    plot_time=True,
                    plot_ref=False,
                    subname_input=subname,
                    show_plots=False,
                    figure=self.figure,
                )
                self.log_edit.append(f"Slices processed: {len(output)}")
            else:
                starts, ends = request["starts"], request["ends"]
                self.log_edit.append(f"Mode: custom time ranges ({len(starts)} slice(s))")
                for s, e in zip(starts, ends):
                    self.log_edit.append(f"  -> processing time window: {s} to {e} s")
                output, plots = reduce_time_list(
                    run_number,
                    settings_file,
                    experiment_id,
                    starts,
                    ends,
                    savepath=savepath,
                    plot_time=True,
                    plot_ref=False,
                    subname_input=subname,
                    show_plots=False,
                    figure=self.figure,
                )
                self.log_edit.append(f"Custom windows processed: {len(output)}")
        finally:
            self.process_btn.setEnabled(True)

        # The kinetic plot is drawn into this tab's own figure (figure=self.figure), on its own canvas: only a redraw
        if plots is not None and self.canvas is not None:
            self.canvas.draw_idle()
            self.log_edit.append("Updated embedded plot view")

        self.log_edit.append(f"Completed reduction: {len(output)} output slice(s)")


__all__ = ["TimeResolvedTab"]
