#!/usr/bin/python3
import os
from qtpy import QtCore
from qtpy.QtWidgets import (
    QFileDialog,
    QGridLayout,
    QLabel,
    QPushButton,
    QWidget,
    QLineEdit,
    QMessageBox,
    QHBoxLayout,
    QVBoxLayout,
    QCheckBox,
    QGroupBox,
    QFormLayout,
    QSpinBox,
    QDoubleSpinBox,
)

import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

# Try matplotlib Qt backends
FigureCanvas = None
NavigationToolbar = None
try:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
except Exception:
    try:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
    except Exception:
        FigureCanvas = None
        NavigationToolbar = None

# Attempt to import the Direct_Beam class
try:
    from lr_reduction.direct_beam_maker import Direct_Beam
except Exception:
    # Try to add project's src/ to sys.path so relative imports work when running from launcher/
    try:
        import sys
        this_dir = os.path.dirname(os.path.abspath(__file__))
        # launcher/apps -> go up two levels to new_workflow, then to src
        src_dir = os.path.abspath(os.path.join(this_dir, '..', '..', 'src/lr_reduction'))
        if os.path.isdir(src_dir) and src_dir not in sys.path:
            sys.path.insert(0, src_dir)
        from lr_reduction.direct_beam_maker import Direct_Beam
    except Exception:
        Direct_Beam = None


class CdSettingsDialog(QMessageBox):
    def __init__(self, parent=None, defaults=None, initial_defaults=None):
        super().__init__(parent)
        # Using simple dialog via QMessageBox with custom widget is cumbersome; build a QDialog-like widget
        from qtpy.QtWidgets import QDialog, QFormLayout, QDialogButtonBox
        self.dlg = QDialog(parent)
        self.dlg.setWindowTitle("Cd settings")
        layout = QVBoxLayout()
        form = QFormLayout()

        # Widgets
        self.mu_file_edit = QLineEdit(self.dlg)
        self.cd_foils_edit = QLineEdit(self.dlg)
        self.cd_edit = QLineEdit(self.dlg)
        self.flip_check = QCheckBox('Flip attenuator mapping', self.dlg)

        # store working defaults (what the dialog initially shows) and the canonical
        # initial defaults used for the Reset action
        self._defaults = defaults or {}
        self._initial_defaults = initial_defaults or defaults or {}
        if defaults:
            self.mu_file_edit.setText(defaults.get('mu_file', ''))
            self.cd_foils_edit.setText(','.join([str(x) for x in defaults.get('Cd_foils', [])]))
            self.cd_edit.setText(','.join([str(x) for x in defaults.get('Cd', [])]))
            # flip_atten default if provided
            try:
                self.flip_check.setChecked(bool(defaults.get('flip_atten', False)))
            except Exception:
                pass

        form.addRow('mu_file:', self.mu_file_edit)
        form.addRow('Cd_foils (comma separated):', self.cd_foils_edit)
        form.addRow('Cd (comma separated):', self.cd_edit)
        form.addRow(self.flip_check)

        layout.addLayout(form)

        # Reset defaults button
        self.reset_btn = QPushButton('Reset defaults')
        self.reset_btn.clicked.connect(self._reset_defaults)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.dlg.accept)
        buttons.rejected.connect(self.dlg.reject)

        # buttons layout: reset button then ok/cancel
        hb = QHBoxLayout()
        hb.addWidget(self.reset_btn)
        hb.addWidget(buttons)
        layout.addLayout(hb)
        self.dlg.setLayout(layout)

    def _reset_defaults(self):
        # restore widgets to the canonical initial defaults (not the current working defaults)
        d = self._initial_defaults or {}
        self.mu_file_edit.setText(d.get('mu_file', ''))
        # MUpath removed; keep mu_file reset only
        self.cd_foils_edit.setText(','.join([str(x) for x in d.get('Cd_foils', [])]))
        self.cd_edit.setText(','.join([str(x) for x in d.get('Cd', [])]))
        try:
            self.flip_check.setChecked(bool(d.get('flip_atten', False)))
        except Exception:
            pass

    def exec_(self):
        return self.dlg.exec_()

    def get_values(self):
        def parse_list(s):
            try:
                parts = [float(x.strip()) for x in s.split(',') if x.strip()]
                return parts
            except Exception:
                return []
        return {
            'mu_file': self.mu_file_edit.text().strip(),
            'Cd_foils': parse_list(self.cd_foils_edit.text()),
            'Cd': parse_list(self.cd_edit.text()),
            'flip_atten': bool(self.flip_check.isChecked()),
        }


class ModeratorDialog(QMessageBox):
    def __init__(self, parent=None, defaults=None, initial_defaults=None):
        from qtpy.QtWidgets import QDialog, QFormLayout, QDialogButtonBox
        super().__init__(parent)
        self.dlg = QDialog(parent)
        self.dlg.setWindowTitle('Moderator settings')
        layout = QVBoxLayout()
        form = QFormLayout()

        self.chop2_edit = QLineEdit(self.dlg)
        self.dMod_edit = QDoubleSpinBox(self.dlg)
        self.dMod_edit.setRange(0.0, 1e6)
        self.t0_edit = QLineEdit(self.dlg)

        # store working defaults and canonical initial defaults for reset
        self._defaults = defaults or {}
        self._initial_defaults = initial_defaults or defaults or {}
        if defaults:
            self.chop2_edit.setText(','.join([str(x) for x in defaults.get('Chop2_cut_fn', [])]))
            self.dMod_edit.setValue(defaults.get('dMod', 15500))
            self.t0_edit.setText(','.join([str(x) for x in defaults.get('t0', [])]))

        form.addRow('Chop2_cut_fn (a,b):', self.chop2_edit)
        form.addRow('dMod:', self.dMod_edit)
        form.addRow('t0 (t0_0,t0_1):', self.t0_edit)
        layout.addLayout(form)

        # Reset defaults button
        self.reset_btn = QPushButton('Reset defaults')
        self.reset_btn.clicked.connect(self._reset_defaults)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.dlg.accept)
        buttons.rejected.connect(self.dlg.reject)
        hb = QHBoxLayout()
        hb.addWidget(self.reset_btn)
        hb.addWidget(buttons)
        layout.addLayout(hb)
        self.dlg.setLayout(layout)

    def _reset_defaults(self):
        d = self._initial_defaults or {}
        try:
            self.chop2_edit.setText(','.join([str(x) for x in d.get('Chop2_cut_fn', [])]))
        except Exception:
            self.chop2_edit.setText('')
        try:
            self.dMod_edit.setValue(d.get('dMod', 15500))
        except Exception:
            pass
        try:
            self.t0_edit.setText(','.join([str(x) for x in d.get('t0', [])]))
        except Exception:
            self.t0_edit.setText('')

    def exec_(self):
        return self.dlg.exec_()

    def get_values(self):
        def parse_pair(s):
            try:
                parts = [float(x.strip()) for x in s.split(',') if x.strip()]
                return parts
            except Exception:
                return []
        return {
            'Chop2_cut_fn': parse_pair(self.chop2_edit.text()),
            'dMod': float(self.dMod_edit.value()),
            't0': parse_pair(self.t0_edit.text()),
        }


class DirectBeamTab(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle('Direct beam')
        self.settings = QtCore.QSettings()

        # UI layout: left controls, right canvas
        main_layout = QHBoxLayout()
        self.setLayout(main_layout)

        controls = QVBoxLayout()

        formbox = QGroupBox('Parameters')
        form = QFormLayout()
        formbox.setLayout(form)

        # Run list
        self.run_list_edit = QLineEdit(self)
        self.run_list_edit.setPlaceholderText('e.g. 1001,1002,1005')
        form.addRow('Run list:', self.run_list_edit)

        # IPTS
        self.ipts_edit = QLineEdit(self)
        form.addRow('IPTS:', self.ipts_edit)

        # IPTS path toggle
        self.ipts_toggle = QCheckBox('Use IPTS path structure')
        form.addRow(self.ipts_toggle)

        # Nexus path and save path
        self.nexus_edit = QLineEdit(self)
        form.addRow('NEXUS path:', self.nexus_edit)
        self.savepath_edit = QLineEdit(self)
        form.addRow('Save path:', self.savepath_edit)

        # Save name
        self.savename_edit = QLineEdit(self)
        self.savename_edit.setText('direct_beam.txt')
        form.addRow('Save name:', self.savename_edit)

        # Numeric parameters
        self.DTCcut_spin = QDoubleSpinBox(self)
        self.DTCcut_spin.setRange(0.0, 100.0)
        self.DTCcut_spin.setValue(1.25)
        form.addRow('DTCcut:', self.DTCcut_spin)

        self.DTCcut1_spin = QDoubleSpinBox(self)
        self.DTCcut1_spin.setRange(0.0, 100.0)
        self.DTCcut1_spin.setValue(1.5)
        form.addRow('DTCcut_config1:', self.DTCcut1_spin)

        self.Icut_spin = QDoubleSpinBox(self)
        self.Icut_spin.setRange(0.0, 1e6)
        self.Icut_spin.setDecimals(12)
        self.Icut_spin.setValue(1e-10)
        form.addRow('Icut:', self.Icut_spin)

        self.CutOffset_spin = QDoubleSpinBox(self)
        self.CutOffset_spin.setRange(-1e6, 1e6)
        self.CutOffset_spin.setValue(1)
        form.addRow('CutOffset:', self.CutOffset_spin)

        # y_ROI and low_res as simple text inputs (two ints comma separated)
        self.yroi_edit = QLineEdit(self)
        self.yroi_edit.setText('130,170')
        form.addRow('y_ROI (min,max):', self.yroi_edit)

        self.lowres_edit = QLineEdit(self)
        self.lowres_edit.setText('75,190')
        form.addRow('low_res (min,max):', self.lowres_edit)

        # Plot toggle
        self.plot_cb = QCheckBox('Plot')
        self.plot_cb.setChecked(True)
        form.addRow(self.plot_cb)

        # Cd settings and moderator settings popouts
        pop_row = QHBoxLayout()
        self.cd_btn = QPushButton('Cd settings')
        self.mod_btn = QPushButton('Moderator settings')
        pop_row.addWidget(self.cd_btn)
        pop_row.addWidget(self.mod_btn)
        form.addRow(pop_row)

        # Execute button
        self.run_btn = QPushButton('Create DB')
        form.addRow(self.run_btn)

        controls.addWidget(formbox)
        controls.addStretch()
        main_layout.addLayout(controls, 1)

        # Canvas
        canvas_layout = QVBoxLayout()
        if FigureCanvas is not None:
            self.figure = Figure(figsize=(5, 4))
            self.canvas = FigureCanvas(self.figure)
            self.toolbar = NavigationToolbar(self.canvas, self)
            canvas_layout.addWidget(self.toolbar)
            canvas_layout.addWidget(self.canvas)
        else:
            self.figure = None
            self.canvas = None
            self.toolbar = None
            canvas_layout.addWidget(QLabel('Embedded plotting unavailable'))

        main_layout.addLayout(canvas_layout, 3)

        # Defaults for dialogs taken from Direct_Beam if available
        default_vals = {}
        if Direct_Beam is not None:
            try:
                db = Direct_Beam()
                default_vals = {
                    # default mu_file left blank so Direct_Beam can discover via settings
                    'mu_file': '',
                    'Cd_foils': getattr(db, 'Cd_foils', []),
                    'Cd': getattr(db, 'Cd', []),
                    'flip_atten': getattr(db, 'flip_atten', False),
                    'Chop2_cut_fn': getattr(db, 'Chop2_cut_fn', []),
                    'dMod': getattr(db, 'dMod', 15500),
                    't0': getattr(db, 't0', []),
                }
            except Exception:
                default_vals = {}

        self.cd_defaults = default_vals
        self.mod_defaults = default_vals

        # Connections
        self.cd_btn.clicked.connect(self._open_cd_settings)
        self.mod_btn.clicked.connect(self._open_mod_settings)
        self.run_btn.clicked.connect(self._run_create_db)
        self.ipts_toggle.toggled.connect(self._ipts_toggled)

    def _open_cd_settings(self):
        # If the user has previously set cd values, use them as defaults so the dialog
        # retains the flip_atten checkbox and other values between openings.
        defaults = getattr(self, 'cd_vals', None) or self.cd_defaults
        dlg = CdSettingsDialog(self, defaults=defaults, initial_defaults=self.cd_defaults)
        if dlg.exec_() == 1:
            vals = dlg.get_values()
            self.cd_vals = vals

    def _open_mod_settings(self):
        # Use previous moderator values if available so the dialog retains edits between opens
        defaults = getattr(self, 'mod_vals', None) or self.mod_defaults
        dlg = ModeratorDialog(self, defaults=defaults, initial_defaults=self.mod_defaults)
        if dlg.exec_() == 1:
            vals = dlg.get_values()
            self.mod_vals = vals

    def _ipts_toggled(self, state):
        if state:
            ipts = self.ipts_edit.text().strip()
            if ipts:
                self.nexus_edit.setText(f"/SNS/REF_L/IPTS-{ipts}/nexus/")
                self.savepath_edit.setText(f"/SNS/REF_L/IPTS-{ipts}/shared/transmission/")
            # disable manual edits
            self.nexus_edit.setEnabled(False)
            self.savepath_edit.setEnabled(False)
        else:
            self.nexus_edit.setEnabled(True)
            self.savepath_edit.setEnabled(True)

    def _parse_run_list(self, text):
        import re
        tokens = re.split(r"[,\s]+", text.strip())
        runs = []
        for t in tokens:
            if not t:
                continue
            # allow ranges using '-' or ':' (inclusive)
            if '-' in t or ':' in t:
                try:
                    if '-' in t:
                        a, b = t.split('-', 1)
                    else:
                        a, b = t.split(':', 1)
                    a = int(a); b = int(b)
                    if b >= a:
                        runs.extend(list(range(a, b+1)))
                except Exception:
                    continue
            else:
                try:
                    runs.append(int(t))
                except Exception:
                    continue
        return runs

    def _run_create_db(self):
        # Validate Direct_Beam availability
        if Direct_Beam is None:
            QMessageBox.critical(self, 'Missing module', 'lr_reduction.direct_beam_maker.Direct_Beam not available')
            return
        runs = self._parse_run_list(self.run_list_edit.text())
        if not runs:
            QMessageBox.warning(self, 'No runs', 'Please provide a run list')
            return

        # Instantiate and set parameters
        try:
            db = Direct_Beam()
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Failed to create Direct_Beam: {e}')
            return

        # IPTS handling: prefer passing IPTS per-call to create_db so the
        # same Direct_Beam instance can be reused for multiple IPTS values.
        experiment_id_for_call = None
        if self.ipts_toggle.isChecked():
            ipts = self.ipts_edit.text().strip()
            if not ipts:
                QMessageBox.warning(self, 'IPTS missing', 'Please provide IPTS when using IPTS path structure')
                return
            # do not mutate the Direct_Beam instance; pass IPTS at call time
            experiment_id_for_call = ipts
        else:
            nx = self.nexus_edit.text().strip()
            sv = self.savepath_edit.text().strip()
            if nx:
                db.NEXUSpath = nx
            if sv:
                db.savepath = sv

        # other parameters
        db.DTCcut = float(self.DTCcut_spin.value())
        db.DTCcut_config1 = float(self.DTCcut1_spin.value())
        db.Icut = float(self.Icut_spin.value())
        db.CutOffset = float(self.CutOffset_spin.value())
        # parse y_roi and low_res
        def parse_pair(s, default):
            try:
                a, b = [int(x.strip()) for x in s.split(',')]
                return [a, b]
            except Exception:
                return default
        db.y_ROI = parse_pair(self.yroi_edit.text(), db.y_ROI if hasattr(db, 'y_ROI') else [130,170])
        db.low_res = parse_pair(self.lowres_edit.text(), db.low_res if hasattr(db, 'low_res') else [75,190])

        # cd settings
        mu_file_arg = None
        flip_atten = False
        try:
            cd = getattr(self, 'cd_vals', None)
            if cd:
                mu = cd.get('mu_file', '').strip()
                # prefer an explicit path if it exists; otherwise let Direct_Beam consult settings
                if mu:
                    # if mu already looks like an absolute or relative path, check it
                    if os.path.isabs(mu) or os.path.sep in mu:
                        if os.path.isfile(mu):
                            mu_file_arg = mu
                        else:
                            # provided path doesn't exist; don't pass it
                            mu_file_arg = None
                    else:
                        # just a filename was provided; leave None so the library will
                        # search known locations (settings/IPTS/BIN paths)
                        mu_file_arg = None
                # apply other cd overrides if provided
                if cd.get('Cd_foils'):
                    db.Cd_foils = cd.get('Cd_foils')
                if cd.get('Cd'):
                    db.Cd = cd.get('Cd')
                flip_atten = bool(cd.get('flip_atten', False))
        except Exception:
            pass

        # moderator settings
        try:
            mod = getattr(self, 'mod_vals', None)
            if mod:
                if mod.get('Chop2_cut_fn'):
                    db.Chop2_cut_fn = mod.get('Chop2_cut_fn')
                if mod.get('dMod'):
                    db.dMod = mod.get('dMod')
                if mod.get('t0'):
                    db.t0 = mod.get('t0')
        except Exception:
            pass

        # call create_db; request intermediate traces when plotting so we can show them
        return_traces = bool(self.plot_cb.isChecked())
        try:
            if return_traces:
                lam, inten, err, traces = db.create_db(
                    runs,
                    self.savename_edit.text().strip(),
                    plot=False,
                    mu_file=mu_file_arg,
                    flip_atten=flip_atten,
                    return_traces=True,
                    experiment_id=experiment_id_for_call,
                )
            else:
                lam, inten, err = db.create_db(
                    runs,
                    self.savename_edit.text().strip(),
                    plot=False,
                    mu_file=mu_file_arg,
                    flip_atten=flip_atten,
                    experiment_id=experiment_id_for_call,
                )
                traces = []
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'create_db failed: {e}')
            return

        # plot if requested
        if self.plot_cb.isChecked():
            if self.canvas is not None and self.figure is not None:
                try:
                    self.figure.clear()
                    ax = self.figure.add_subplot(111)
                    # plot intermediate per-run traces if present
                    if traces:
                        cmap = plt.get_cmap('tab10')
                        for idx, (L, I_trans, run) in enumerate(traces):
                            try:
                                color = cmap(idx % cmap.N)
                                ax.plot(L, I_trans, '-', color=color, alpha=0.6, linewidth=0.9, label=str(run))
                            except Exception:
                                pass
                    # plot aggregated result on top
                    ax.errorbar(lam, inten, yerr=err, fmt='o', markersize=2, color='k', label='aggregate')
                    # add legend for runs
                    try:
                        ax.legend(title='run', fontsize='small', ncol=2)
                    except Exception:
                        pass
                    ax.set_xscale('linear')
                    ax.set_yscale('log')
                    ax.set_xlabel('Lambda [A]')
                    ax.set_ylabel('I')
                    self.figure.tight_layout()
                    self.canvas.draw()
                except Exception as e:
                    QMessageBox.warning(self, 'Plot error', f'Failed to plot: {e}')
            else:
                # fallback to pyplot
                plt.figure()
                if traces:
                    cmap = plt.get_cmap('tab10')
                    for idx, (L, I_trans, run) in enumerate(traces):
                        try:
                            color = cmap(idx % cmap.N)
                            plt.plot(L, I_trans, '-', color=color, alpha=0.6, linewidth=0.9, label=str(run))
                        except Exception:
                            pass
                plt.errorbar(lam, inten, yerr=err, fmt='o', markersize=2, color='k', label='aggregate')
                try:
                    plt.legend(title='run', fontsize='small', ncol=2)
                except Exception:
                    pass
                plt.xscale('linear')
                plt.yscale('log')
                plt.xlabel('Lambda [A]')
                plt.ylabel('I')
                plt.tight_layout()
                plt.show()
 