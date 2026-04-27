#!/usr/bin/python3
import os
import threading
from pathlib import Path

from qtpy import QtCore, QtWidgets
from qtpy.QtWidgets import (
    QCheckBox,
    QFileDialog,
    QGridLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QWidget,
)

try:
    # Use QtAgg canvas for embedding matplotlib figures
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
except Exception:
    FigureCanvas = None
    NavigationToolbar = None
    try:
        import matplotlib.pyplot as plt
    except Exception:
        plt = None

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
        # Experiment id and update defaults — keep left-justified with other inputs
        ipts_label = QLabel("Experiment id (IPTS-...):")
        self.experiment_edit = QLineEdit()
        # limit width so the row stays compact and the button sits nearer the left
        try:
            self.experiment_edit.setMaximumWidth(220)
        except Exception:
            pass
        self.update_defaults_btn = QPushButton("Update defaults")
        try:
            self.update_defaults_btn.setMaximumWidth(120)
        except Exception:
            pass
        # Place into the same grid columns as other rows (0,1,2) so it's left-justified
        layout.addWidget(ipts_label, 0, 0)
        layout.addWidget(self.experiment_edit, 0, 1)
        layout.addWidget(self.update_defaults_btn, 0, 2)

        # Inline paths (datapath, DBpath, Spath)
        layout.addWidget(QLabel("Datapath (NEXUS):"), 1, 0)
        self.datapath_edit = QLineEdit()
        layout.addWidget(self.datapath_edit, 1, 1)
        self.datapath_btn = QPushButton("Browse")
        layout.addWidget(self.datapath_btn, 1, 2)

        layout.addWidget(QLabel("Direct Beam path:"), 2, 0)
        self.dbpath_edit = QLineEdit()
        layout.addWidget(self.dbpath_edit, 2, 1)
        self.dbpath_btn = QPushButton("Browse")
        layout.addWidget(self.dbpath_btn, 2, 2)

        layout.addWidget(QLabel("Save path (output):"), 3, 0)
        self.spath_edit = QLineEdit()
        layout.addWidget(self.spath_edit, 3, 1)
        self.spath_btn = QPushButton("Browse")
        layout.addWidget(self.spath_btn, 3, 2)

        # Template directory and file
        layout.addWidget(QLabel("Template directory:"), 4, 0)
        self.template_dir_edit = QLineEdit()
        layout.addWidget(self.template_dir_edit, 4, 1)
        self.template_dir_btn = QPushButton("Browse")
        layout.addWidget(self.template_dir_btn, 4, 2)

        layout.addWidget(QLabel("Template file name (xml):"), 5, 0)
        self.template_file_edit = QLineEdit()
        layout.addWidget(self.template_file_edit, 5, 1)
        # file browse is optional and can be disabled if causes hangs
        self.template_file_btn = QPushButton("Browse")
        layout.addWidget(self.template_file_btn, 5, 2)

        # Runs input
        layout.addWidget(QLabel("Runs (e.g. 1234-1236,1240,1250:1252):"), 6, 0)
        self.runs_edit = QLineEdit()
        layout.addWidget(self.runs_edit, 6, 1)

        layout.addWidget(QLabel("File Name Suffix (optional):"), 7, 0)
        self.subname_edit = QLineEdit()
        layout.addWidget(self.subname_edit, 7, 1)

        # Plot checkbox (separate row to avoid overlap)
        self.plot_checkbox = QCheckBox("Plot results")
        self.plot_checkbox.setChecked(False)
        layout.addWidget(self.plot_checkbox, 8, 1)

        # Enable browse buttons toggle (default disabled to avoid hangs)
        self.enable_browse_chk = QCheckBox("Enable Browse buttons")
        self.enable_browse_chk.setChecked(False)
        layout.addWidget(self.enable_browse_chk, 8, 2)

        # Process button
        self.process_btn = QPushButton("Process")
        self.process_btn.setStyleSheet("background-color : green")
        layout.addWidget(self.process_btn, 9, 1)

        # Log area (left column)
        self.log_edit = QtWidgets.QTextEdit()
        self.log_edit.setReadOnly(True)
        layout.addWidget(self.log_edit, 10, 0, 3, 1)

        # Embedded plot viewer (single canvas with prev/next navigation)
        # Keep a capped list of Figures to avoid creating too many tabs/windows.
        self._MAX_STORED_FIGURES = 200
        self.figures = []  # list of tuples: (Figure, meta_dict)
        self.current_fig_index = -1

        # Plot controls: Prev / Next / Clear + Max plots spinbox
        plot_ctrl_widget = QtWidgets.QWidget()
        plot_ctrl_layout = QtWidgets.QHBoxLayout()
        plot_ctrl_widget.setLayout(plot_ctrl_layout)
        self.prev_btn = QPushButton("Prev")
        self.next_btn = QPushButton("Next")
        self.clear_plots_btn = QPushButton("Clear plots")
        self.plot_index_label = QLabel("No plots")
        # user-configurable max plots (warning threshold)
        self.max_plots_label = QLabel("Max stored plots:")
        self.max_plots_spin = QSpinBox()
        self.max_plots_spin.setMinimum(1)
        self.max_plots_spin.setMaximum(1000)
        self.max_plots_spin.setValue(80)
        plot_ctrl_layout.addWidget(self.prev_btn)
        plot_ctrl_layout.addWidget(self.next_btn)
        plot_ctrl_layout.addWidget(self.plot_index_label)
        plot_ctrl_layout.addStretch()
        plot_ctrl_layout.addWidget(self.max_plots_label)
        plot_ctrl_layout.addWidget(self.max_plots_spin)
        plot_ctrl_layout.addWidget(self.clear_plots_btn)
        # place plot controls in right column (col 1..2)
        layout.addWidget(plot_ctrl_widget, 10, 1, 1, 2)

        # Canvas placeholder (right column)
        self.canvas_container = QtWidgets.QWidget()
        self.canvas_layout = QtWidgets.QVBoxLayout()
        self.canvas_container.setLayout(self.canvas_layout)
        layout.addWidget(self.canvas_container, 11, 1, 2, 2)

        # Connect plot controls
        self.prev_btn.clicked.connect(self._show_prev_figure)
        self.next_btn.clicked.connect(self._show_next_figure)
        self.clear_plots_btn.clicked.connect(self._clear_plots)

        # Connections
        self.template_dir_btn.clicked.connect(self._browse_template_dir)
        self.template_file_btn.clicked.connect(self._browse_template_file)
        self.datapath_btn.clicked.connect(self._browse_datapath)
        self.dbpath_btn.clicked.connect(self._browse_dbpath)
        self.spath_btn.clicked.connect(self._browse_spath)
        self.update_defaults_btn.clicked.connect(self.update_defaults_from_experiment)
        self.process_btn.clicked.connect(self._on_process_clicked)
        self.enable_browse_chk.toggled.connect(self._toggle_browse_buttons)
        # if matplotlib is available, ensure Prev/Next are enabled state matches availability
        self._update_plot_controls()

        # Load settings (this will also set the max plots spinbox)
        self.read_settings()

        # track user-configured max plots convenience property
        try:
            self.user_max_plots = int(self.max_plots_spin.value())
        except Exception:
            self.user_max_plots = 80
        self.max_plots_spin.valueChanged.connect(self._on_max_plots_changed)

    def _browse_template(self):
        # kept for backward compatibility; prefer using the separate template dir + file fields
        opts = QFileDialog.Options()
        opts |= QFileDialog.DontUseNativeDialog
        _file, _ = QFileDialog.getOpenFileName(self, "Select template file:", self.template_file_edit.text(), "", options=opts)
        if _file:
            p = Path(_file)
            # set directory and filename immediately, validate directory asynchronously
            self._set_path_async(self.template_dir_edit, str(p.parent), check_isdir=True)
            try:
                self.template_file_edit.setText(str(p.name))
            except Exception:
                pass

    def _browse_template_dir(self):
        opts = QFileDialog.Options()
        opts |= QFileDialog.DontUseNativeDialog
        _dir = QFileDialog.getExistingDirectory(self, "Select template directory:", self.template_dir_edit.text(), options=opts)
        if _dir:
            # set immediately and validate asynchronously
            self._set_path_async(self.template_dir_edit, str(_dir), check_isdir=True)

    def _browse_template_file(self):
        # Opens file dialog rooted at template_dir if available
        start = self.template_dir_edit.text() or os.path.expanduser("~")
        opts = QFileDialog.Options()
        opts |= QFileDialog.DontUseNativeDialog
        _file, _ = QFileDialog.getOpenFileName(self, "Select template file:", start, "", options=opts)
        if _file:
            p = Path(_file)
            self._set_path_async(self.template_dir_edit, str(p.parent), check_isdir=True)
            try:
                self.template_file_edit.setText(str(p.name))
            except Exception:
                pass

    def _browse_datapath(self):
        opts = QFileDialog.Options()
        opts |= QFileDialog.DontUseNativeDialog
        _dir = QFileDialog.getExistingDirectory(self, "Select datapath (NEXUS):", self.datapath_edit.text(), options=opts)
        if _dir:
            self._set_path_async(self.datapath_edit, str(_dir), check_isdir=True)

    def _browse_dbpath(self):
        opts = QFileDialog.Options()
        opts |= QFileDialog.DontUseNativeDialog
        _dir = QFileDialog.getExistingDirectory(self, "Select DBpath:", self.dbpath_edit.text(), options=opts)
        if _dir:
            self._set_path_async(self.dbpath_edit, str(_dir), check_isdir=True)

    def _browse_spath(self):
        opts = QFileDialog.Options()
        opts |= QFileDialog.DontUseNativeDialog
        _dir = QFileDialog.getExistingDirectory(self, "Select Spath (output):", self.spath_edit.text(), options=opts)
        if _dir:
            self._set_path_async(self.spath_edit, str(_dir), check_isdir=True)

    def update_defaults_from_experiment(self):
        expt = self.experiment_edit.text().strip()
        if not expt:
            QMessageBox.warning(self, "No experiment id", "Please enter an experiment id (e.g. IPTS-36119) before updating defaults")
            return
        # normalize experiment id to start with IPTS-
        if not expt.upper().startswith("IPTS-"):
            expt = f"IPTS-{expt}"
            self.experiment_edit.setText(expt)
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
            # template directory default
            try:
                self.template_dir_edit.setText(str(Path("/SNS/REF_L") / expt / "shared"))
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

        _subname = self.settings.value("template_subname", None)
        self.subname_edit.setText(_subname)

        _expt = self.settings.value("template_experiment_id", "")
        self.experiment_edit.setText(_expt)

        # template dir/file
        _template_dir = self.settings.value("template_dir", "")
        _template_file = self.settings.value("template_file", "")
        if _template_dir:
            self.template_dir_edit.setText(_template_dir)
        if _template_file:
            self.template_file_edit.setText(_template_file)

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
        # browse enable
        _browse = self.settings.value("template_enable_browse", True)
        self.enable_browse_chk.setChecked(bool(_browse))
        self._toggle_browse_buttons(bool(_browse))
        # max plots
        try:
            _max = int(self.settings.value("template_max_plots", 80))
            self.max_plots_spin.setValue(_max)
        except Exception:
            pass

    def save_settings(self):
        self.settings.setValue("template_runs", self.runs_edit.text())
        self.settings.setValue("template_experiment_id", self.experiment_edit.text())
        subname_text = self.subname_edit.text() if self.subname_edit else None
        self.settings.setValue("template_subname", subname_text or None)
        # persist template dir and file separately
        self.settings.setValue("template_dir", self.template_dir_edit.text())
        self.settings.setValue("template_file", self.template_file_edit.text())
        self.settings.setValue("template_plot", bool(self.plot_checkbox.isChecked()))
        # persist inline overrides
        self.settings.setValue("template_datapath", self.datapath_edit.text())
        self.settings.setValue("template_DBpath", self.dbpath_edit.text())
        self.settings.setValue("template_Spath", self.spath_edit.text())
        self.settings.setValue("template_enable_browse", bool(self.enable_browse_chk.isChecked()))
        try:
            self.settings.setValue("template_max_plots", int(self.max_plots_spin.value()))
        except Exception:
            pass

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

    def _set_path_async(self, line_edit: QLineEdit, path: str, check_isdir: bool = True):
        """Set the QLineEdit immediately, then validate the path in a background
        thread so slow network filesystems (or blocking os.stat) do not hang the
        GUI. If the path appears inaccessible, log a warning on the UI thread.
        """
        try:
            # show chosen path immediately so UI doesn't feel blocked
            line_edit.setText(str(path))
        except Exception:
            try:
                line_edit.setText(path)
            except Exception:
                pass

        if not check_isdir:
            return

        def _worker(p):
            ok = False
            try:
                ok = os.path.isdir(p)
            except Exception:
                ok = False
            if not ok:
                try:
                    # Post a message to the log on the GUI thread
                    QtCore.QMetaObject.invokeMethod(
                        self,
                        "_append_log",
                        QtCore.Qt.QueuedConnection,
                        QtCore.Q_ARG(str, f"Warning: chosen path appears inaccessible: {p}"),
                    )
                except Exception:
                    pass

        th = threading.Thread(target=_worker, args=(str(path),), daemon=True)
        th.start()

    def _toggle_browse_buttons(self, enabled: bool):
        # show/hide browse buttons depending on user preference
        try:
            self.template_dir_btn.setVisible(enabled)
            self.template_file_btn.setVisible(enabled)
            self.datapath_btn.setVisible(enabled)
            self.dbpath_btn.setVisible(enabled)
            self.spath_btn.setVisible(enabled)
        except Exception:
            pass

    def _on_max_plots_changed(self, val):
        try:
            self.user_max_plots = int(val)
        except Exception:
            self.user_max_plots = None
        # update controls (enable/disable prev/next depending on stored figs)
        self._update_plot_controls()

    def _update_plot_controls(self):
        enabled = FigureCanvas is not None and plt is not None
        self.prev_btn.setEnabled(enabled and len(self.figures) > 0)
        self.next_btn.setEnabled(enabled and len(self.figures) > 0)
        self.clear_plots_btn.setEnabled(enabled and len(self.figures) > 0)

    def _show_figure_at_index(self, idx: int):
        # Remove existing canvas widgets
        for i in reversed(range(self.canvas_layout.count())):
            w = self.canvas_layout.takeAt(i).widget()
            if w is not None:
                try:
                    w.setParent(None)
                    w.deleteLater()
                except Exception:
                    pass

        if idx < 0 or idx >= len(self.figures):
            self.current_fig_index = -1
            self.plot_index_label.setText("No plots")
            self._update_plot_controls()
            return

        fig, meta = self.figures[idx]
        # Create a canvas for the existing figure and insert into layout
        try:
            if FigureCanvas is None:
                # cannot embed
                self.plot_index_label.setText("Matplotlib not available")
                return
            canvas = FigureCanvas(fig)
            self.canvas_layout.addWidget(canvas)
            try:
                # optional navigation toolbar
                if NavigationToolbar is not None:
                    toolbar = NavigationToolbar(canvas, self.canvas_container)
                    self.canvas_layout.addWidget(toolbar)
            except Exception:
                pass
            canvas.draw()
            self.current_fig_index = idx
            self.plot_index_label.setText(f"Plot {idx+1}/{len(self.figures)} -- {meta.get('label','')}")
        except Exception as e:
            self.plot_index_label.setText(f"Could not show figure: {e}")
        self._update_plot_controls()

    def _show_prev_figure(self):
        if not self.figures:
            return
        idx = max(0, (self.current_fig_index - 1) % len(self.figures)) if self.current_fig_index >= 0 else 0
        self._show_figure_at_index(idx)

    def _show_next_figure(self):
        if not self.figures:
            return
        idx = (self.current_fig_index + 1) % len(self.figures) if self.current_fig_index >= 0 else 0
        self._show_figure_at_index(idx)

    def _clear_plots(self):
        # Close matplotlib figures and clear the store
        try:
            if plt is not None:
                for fig, _meta in self.figures:
                    try:
                        plt.close(fig)
                    except Exception:
                        pass
        except Exception:
            pass
        self.figures = []
        self.current_fig_index = -1
        self._show_figure_at_index(-1)


    def _process_sync(self):
        # validate inputs
        runs_text = self.runs_edit.text().strip()
        if not runs_text:
            raise ValueError("Please specify runs to process")
        try:
            runs = parse_run_list(runs_text)
        except Exception:
            raise

        expt = self.experiment_edit.text().strip()
        if not expt:
            raise ValueError("Please specify an experiment id (IPTS-...)')")

        # Build template path from directory + file fields
        template_dir_text = self.template_dir_edit.text().strip()
        template_file_text = self.template_file_edit.text().strip()
        if not template_file_text:
            raise ValueError("Please specify the template file name (xml)")
        if template_dir_text:
            template_file = Path(template_dir_text) / template_file_text
            template_dir = Path(template_dir_text)
        else:
            # If no directory provided, treat template_file_text as an absolute or relative path
            template_file = Path(template_file_text)
            template_dir = template_file.parent
        if not template_file.exists() or not template_file.is_file():
            raise ValueError(f"Template file not found: {str(template_file)}")

        template_subname = self.subname_edit.text() if self.subname_edit else None

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
                if template_subname:
                    override_params["subname"] = template_subname

                # Log start
                QtCore.QMetaObject.invokeMethod(self, "_append_log", QtCore.Qt.QueuedConnection, QtCore.Q_ARG(str, f"Starting run {run}..."))

                # call reduce_from_template similar to examples
                try:
                    # If plotting requested and matplotlib is present, capture any figures
                    if plot_flag and plt is not None:
                        try:
                            old_nums = set(plt.get_fignums())
                        except Exception:
                            old_nums = set()
                        # prevent the reduction from calling a blocking plt.show()
                        orig_show = getattr(plt, 'show', None)
                        try:
                            plt.show = lambda *args, **kwargs: None
                        except Exception:
                            orig_show = None

                        # run the reduction synchronously (figures will be created but not shown)
                        new_template.reduce_from_template(
                            run,
                            Path(template_file),
                            expt,
                            datapath=Path(datapath) if datapath else None,
                            template_path=Path(template_dir),
                            override_params=override_params,
                            plot=plot_flag,
                        )

                        # collect newly created figures
                        try:
                            new_nums = [n for n in plt.get_fignums() if n not in old_nums]
                        except Exception:
                            new_nums = []

                        for num in new_nums:
                            try:
                                fig = plt.figure(num)
                                # Ensure the figure clearly indicates which run it
                                # belongs to: prefer a suptitle, otherwise prefix
                                # existing axis titles.
                                try:
                                    fig.suptitle(f"Run {run}", fontsize=10)
                                except Exception:
                                    try:
                                        for ax in getattr(fig, 'axes', []):
                                            try:
                                                t = ax.get_title()
                                                if t:
                                                    ax.set_title(f"Run {run} -- {t}")
                                                else:
                                                    ax.set_title(f"Run {run}")
                                            except Exception:
                                                pass
                                    except Exception:
                                        pass
                                meta = {"run": run, "label": f"run {run}"}
                                self.figures.append((fig, meta))
                                # warn if the user-configured max is exceeded
                                try:
                                    user_max = int(self.max_plots_spin.value())
                                except Exception:
                                    user_max = None
                                if user_max is not None and len(self.figures) > user_max:
                                    QMessageBox.warning(self, "Plot history large", f"Stored plot count ({len(self.figures)}) exceeded your configured maximum ({user_max}).\nConsider clearing the plot history.")
                            except Exception:
                                pass

                        # enforce cap
                        while len(self.figures) > self._MAX_STORED_FIGURES:
                            fig_to_close, _m = self.figures.pop(0)
                            try:
                                plt.close(fig_to_close)
                            except Exception:
                                pass

                        # restore original show
                        try:
                            if orig_show is not None:
                                plt.show = orig_show
                        except Exception:
                            pass

                        # show the most recently added figure
                        if self.figures:
                            self._show_figure_at_index(len(self.figures) - 1)
                    else:
                        # no plotting requested or matplotlib missing: just run normally
                        new_template.reduce_from_template(
                            run,
                            Path(template_file),
                            expt,
                            datapath=Path(datapath) if datapath else None,
                            template_path=Path(template_dir),
                            override_params=override_params,
                            plot=plot_flag,
                        )
                finally:
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
