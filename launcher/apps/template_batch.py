#!/usr/bin/python3
import os
import threading
from pathlib import Path

from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import (
    QFileDialog,
    QGridLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QCheckBox,
    QLineEdit,
    QDialog,
    QWidget,
)

try:
    import new_reduction_from_template as new_template
except Exception:
    # fallback if package is installed as lr_reduction
    try:
        from lr_reduction import new_reduction_from_template as new_template
    except Exception:
        new_template = None


def parse_run_list(text):
    """Parse text like '1234-1236,1240,1250:1252' -> list of ints."""
    runs = set()
    if text is None:
        return []
    for part in str(text).split(","):
        p = part.strip()
        if p == "":
            continue
        if "-" in p:
            a, b = p.split("-", 1)
            try:
                a_i = int(a)
                b_i = int(b)
            except Exception:
                raise ValueError(f"Invalid range: {p}")
            if a_i <= b_i:
                runs.update(range(a_i, b_i + 1))
            else:
                runs.update(range(b_i, a_i + 1))
        elif ":" in p:
            a, b = p.split(":", 1)
            try:
                a_i = int(a)
                b_i = int(b)
            except Exception:
                raise ValueError(f"Invalid range: {p}")
            if a_i <= b_i:
                runs.update(range(a_i, b_i + 1))
            else:
                runs.update(range(b_i, a_i + 1))
        else:
            try:
                runs.add(int(p))
            except Exception:
                raise ValueError(f"Invalid run value: {p}")
    return sorted(runs)


# Inline path overrides are shown on the main tab. The previous pop-up dialog
# was removed to avoid blocking UI and to simplify flow.


class TemplateBatchTab(QWidget):
    """Tab for batch reduction using a template (reduce_from_template)

    Fields:
    - run list (supports comma and ranges with - or :)
    - experiment id
    - template path (text + browse)
    - plot checkbox
    - advanced override paths dialog
    - process button
    """

    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle("Batch template reduction")
        layout = QGridLayout()
        layout.setColumnStretch(1, 1)
        self.setLayout(layout)

        self.settings = QtCore.QSettings()
        # Experiment id and update defaults
        layout.addWidget(QLabel("Experiment id (IPTS-...):"), 0, 0)
        self.experiment_edit = QLineEdit()
        layout.addWidget(self.experiment_edit, 0, 1)
        self.update_defaults_btn = QPushButton("Update defaults")
        layout.addWidget(self.update_defaults_btn, 0, 2)

        # Inline paths (datapath, DBpath, Spath)
        layout.addWidget(QLabel("Datapath (NEXUS):"), 1, 0)
        self.datapath_edit = QLineEdit()
        layout.addWidget(self.datapath_edit, 1, 1)
        self.datapath_btn = QPushButton("Browse")
        layout.addWidget(self.datapath_btn, 1, 2)

        layout.addWidget(QLabel("DBpath:"), 2, 0)
        self.dbpath_edit = QLineEdit()
        layout.addWidget(self.dbpath_edit, 2, 1)
        self.dbpath_btn = QPushButton("Browse")
        layout.addWidget(self.dbpath_btn, 2, 2)

        layout.addWidget(QLabel("Spath (output):"), 3, 0)
        self.spath_edit = QLineEdit()
        layout.addWidget(self.spath_edit, 3, 1)
        self.spath_btn = QPushButton("Browse")
        layout.addWidget(self.spath_btn, 3, 2)

        # Template path
        layout.addWidget(QLabel("Template file:"), 4, 0)
        self.template_edit = QLineEdit()
        layout.addWidget(self.template_edit, 4, 1)
        self.template_btn = QPushButton("Browse")
        layout.addWidget(self.template_btn, 4, 2)

        # Runs input
        layout.addWidget(QLabel("Runs (e.g. 1234-1236,1240,1250:1252):"), 5, 0)
        self.runs_edit = QLineEdit()
        layout.addWidget(self.runs_edit, 5, 1)

        # Plot checkbox
        self.plot_checkbox = QCheckBox("Plot results")
        self.plot_checkbox.setChecked(False)
        layout.addWidget(self.plot_checkbox, 6, 1)

        # Process button
        self.process_btn = QPushButton("Process")
        self.process_btn.setStyleSheet("background-color : green")
        layout.addWidget(self.process_btn, 7, 1)

        # Log area
        self.log_edit = QtWidgets.QTextEdit()
        self.log_edit.setReadOnly(True)
        layout.addWidget(self.log_edit, 8, 0, 1, 3)

        # Connections
        self.template_btn.clicked.connect(self._browse_template)
        self.datapath_btn.clicked.connect(self._browse_datapath)
        self.dbpath_btn.clicked.connect(self._browse_dbpath)
        self.spath_btn.clicked.connect(self._browse_spath)
        self.update_defaults_btn.clicked.connect(self.update_defaults_from_experiment)
        self.process_btn.clicked.connect(self._on_process_clicked)

        # Load settings
        self.read_settings()

    def _browse_template(self):
        opts = QFileDialog.Options()
        opts |= QFileDialog.DontUseNativeDialog
        _file, _ = QFileDialog.getOpenFileName(self, "Select template file:", self.template_edit.text(), "", options=opts)
        if _file:
            self.template_edit.setText(_file)

    def _browse_datapath(self):
        opts = QFileDialog.Options()
        opts |= QFileDialog.DontUseNativeDialog
        _dir = QFileDialog.getExistingDirectory(self, "Select datapath (NEXUS):", self.datapath_edit.text(), options=opts)
        if os.path.isdir(_dir):
            self.datapath_edit.setText(_dir)

    def _browse_dbpath(self):
        opts = QFileDialog.Options()
        opts |= QFileDialog.DontUseNativeDialog
        _dir = QFileDialog.getExistingDirectory(self, "Select DBpath:", self.dbpath_edit.text(), options=opts)
        if os.path.isdir(_dir):
            self.dbpath_edit.setText(_dir)

    def _browse_spath(self):
        opts = QFileDialog.Options()
        opts |= QFileDialog.DontUseNativeDialog
        _dir = QFileDialog.getExistingDirectory(self, "Select Spath (output):", self.spath_edit.text(), options=opts)
        if os.path.isdir(_dir):
            self.spath_edit.setText(_dir)

    def update_defaults_from_experiment(self):
        expt = self.experiment_edit.text().strip()
        if not expt:
            QMessageBox.warning(self, "No experiment id", "Please enter an experiment id (e.g. IPTS-36119) before updating defaults")
            return
        try:
            from lr_reduction.nr_reduction_config import NRReductionConfig

            cfg = NRReductionConfig()
            cfg.experiment_id = expt
            # populate inline fields
            try:
                self.datapath_edit.setText(str(cfg.NEXUSpathRB))
            except Exception:
                pass
            try:
                self.dbpath_edit.setText(str(cfg.DBpath))
            except Exception:
                pass
            try:
                self.spath_edit.setText(str(cfg.Spath))
            except Exception:
                pass
        except Exception as e:
            QMessageBox.warning(self, "Defaults failed", f"Could not load defaults from NRReductionConfig: {e}")
        # normalize experiment id to start with IPTS-
        if not expt.upper().startswith("IPTS-"):
            expt = f"IPTS-{expt}"
            self.experiment_edit.setText(expt)

    def open_override_dialog(self):
        # kept for backwards compat but just opens a simple info box
        QMessageBox.information(self, "Note", "Path overrides are now inline on the main tab.")

    def read_settings(self):
        _runs = self.settings.value("template_runs", "")
        self.runs_edit.setText(_runs)

        _expt = self.settings.value("template_experiment_id", "")
        self.experiment_edit.setText(_expt)

        _template = self.settings.value("template_template_path", "")
        self.template_edit.setText(_template)

        _plot = self.settings.value("template_plot", False)
        self.plot_checkbox.setChecked(bool(_plot))
        # Inline path fields
        _datapath = self.settings.value("template_datapath", "")
        _DBpath = self.settings.value("template_DBpath", "")
        _Spath = self.settings.value("template_Spath", "")
        if _datapath:
            self.datapath_edit.setText(_datapath)
        if _DBpath:
            self.dbpath_edit.setText(_DBpath)
        if _Spath:
            self.spath_edit.setText(_Spath)

    def save_settings(self):
        self.settings.setValue("template_runs", self.runs_edit.text())
        self.settings.setValue("template_experiment_id", self.experiment_edit.text())
        self.settings.setValue("template_template_path", self.template_edit.text())
        self.settings.setValue("template_plot", bool(self.plot_checkbox.isChecked()))
        # persist inline overrides
        self.settings.setValue("template_datapath", self.datapath_edit.text())
        self.settings.setValue("template_DBpath", self.dbpath_edit.text())
        self.settings.setValue("template_Spath", self.spath_edit.text())

    def _on_process_clicked(self):
        # persist settings on main thread (avoids invokeMethod/slot lookup issues)
        try:
            self.save_settings()
        except Exception:
            pass

        plot_flag = bool(self.plot_checkbox.isChecked())
        # If plotting is requested, run synchronously on main thread so plotting
        # backends (which require the main interpreter thread) work correctly.
        if plot_flag:
            self.process_btn.setEnabled(False)
            try:
                self._process_sync()
            except Exception as e:
                self._show_error(str(e))
            finally:
                self.process_btn.setEnabled(True)
            return

        # otherwise run in worker thread
        self.process_btn.setEnabled(False)
        thread = threading.Thread(target=self._process)
        thread.daemon = True
        thread.start()

    def _process(self):
        try:
            self._process_sync()
        except Exception as e:
            QtCore.QMetaObject.invokeMethod(self, "_show_error", QtCore.Qt.QueuedConnection, QtCore.Q_ARG(str, str(e)))
        finally:
            QtCore.QMetaObject.invokeMethod(self.process_btn, "setEnabled", QtCore.Qt.QueuedConnection, QtCore.Q_ARG(bool, True))

    @QtCore.Slot(str)
    def _show_error(self, msg):
        QMessageBox.critical(None, "Error", msg)

    @QtCore.Slot(str)
    def _append_log(self, msg: str):
        # append text to the log widget (runs on main thread)
        try:
            self.log_edit.append(msg)
        except Exception:
            # best-effort
            pass

    def _process_sync(self):
        # validate inputs
        runs_text = self.runs_edit.text().strip()
        if not runs_text:
            raise ValueError("Please specify runs to process")
        try:
            runs = parse_run_list(runs_text)
        except Exception as e:
            raise

        expt = self.experiment_edit.text().strip()
        if not expt:
            raise ValueError("Please specify an experiment id (IPTS-...)')")

        template_text = self.template_edit.text().strip()
        if not template_text:
            raise ValueError("Please specify the template file path")
        template_path_obj = Path(template_text)
        # Require the template to be a file; if a directory has been provided, treat as error
        if template_path_obj.is_file():
            template_file = template_path_obj
            template_dir = template_file.parent
        elif template_path_obj.is_dir():
            raise ValueError("The template field must point to a template XML file, not a directory")
        else:
            raise ValueError(f"Template file not found: {template_text}")

        # deduce datapath and overrides from inline fields, falling back to NRReductionConfig
        datapath = self.datapath_edit.text().strip() or None
        DBpath = self.dbpath_edit.text().strip() or None
        Spath = self.spath_edit.text().strip() or None

        # as a final fallback, use NRReductionConfig defaults (set with the experiment id)
        try:
            from lr_reduction.nr_reduction_config import NRReductionConfig

            cfg = NRReductionConfig()
            cfg.experiment_id = expt
            if datapath is None:
                datapath = str(cfg.NEXUSpathRB)
            if DBpath is None:
                DBpath = str(cfg.DBpath)
            if Spath is None:
                Spath = str(cfg.Spath)
        except Exception:
            # leave None if not available
            pass

        plot_flag = bool(self.plot_checkbox.isChecked())

        # Save settings was already called on the main thread before the worker started.

        # verify new_template is importable
        if new_template is None:
            raise ImportError("Could not import new_reduction_from_template module")

        # iterate runs
        for run in runs:
            try:
                override_params = {}
                if DBpath:
                    override_params["DBpath"] = Path(DBpath)
                if Spath:
                    override_params["Spath"] = Path(Spath)
                # If datapath provided, also pass as NEXUSpathRB in override_params for completeness
                if datapath:
                    override_params["NEXUSpathRB"] = Path(datapath)

                # Log start
                QtCore.QMetaObject.invokeMethod(self, "_append_log", QtCore.Qt.QueuedConnection, QtCore.Q_ARG(str, f"Starting run {run}..."))

                # call reduce_from_template similar to examples
                new_template.reduce_from_template(
                    run,
                    Path(template_file),
                    expt,
                    datapath=Path(datapath) if datapath else None,
                    template_path=Path(template_dir),
                    override_params=override_params,
                    plot=plot_flag,
                )

                QtCore.QMetaObject.invokeMethod(self, "_append_log", QtCore.Qt.QueuedConnection, QtCore.Q_ARG(str, f"Run {run} completed"))
            except Exception as e:
                # continue processing other runs but surface the error
                tb = str(e)
                try:
                    import traceback

                    tb = traceback.format_exc()
                except Exception:
                    pass
                QtCore.QMetaObject.invokeMethod(self, "_append_log", QtCore.Qt.QueuedConnection, QtCore.Q_ARG(str, f"Run {run} failed: {e}\n{tb}"))
                QtCore.QMetaObject.invokeMethod(self, "_show_error", QtCore.Qt.QueuedConnection, QtCore.Q_ARG(str, f"Run {run} failed: {e}"))

        # finished
        QtCore.QMetaObject.invokeMethod(self, "_show_info", QtCore.Qt.QueuedConnection, QtCore.Q_ARG(str, "Processing complete"))

    @QtCore.Slot(str)
    def _show_info(self, msg):
        QMessageBox.information(None, "Info", msg)


__all__ = ["TemplateBatchTab"]
