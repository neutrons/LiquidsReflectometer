#!/usr/bin/python3
import json
import os
import subprocess

from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import QFileDialog, QGridLayout, QLabel, QMessageBox, QPushButton, QWidget

REFERENCE_DIRECTIVE = "Click to choose a 60Hz reference R(Q) file"
TEMPLATE_DIRECTIVE = "Click to choose a 30Hz template"
OUTPUT_DIR_DIRECTIVE = "Click to choose an output directory"


class Dynamic30Hz(QWidget):

    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle('Time-resolved reduction')
        layout = QGridLayout()
        self.setLayout(layout)

        self.settings = QtCore.QSettings()

        # 30Hz template file
        self.choose_template = QPushButton('30Hz template')
        layout.addWidget(self.choose_template, 1, 1)

        self.template_path = QLabel(self)
        layout.addWidget(self.template_path, 1, 2)

        # 60Hz reference file
        self.choose_ref = QPushButton('60Hz R(Q) reference')
        layout.addWidget(self.choose_ref, 2, 1)

        self.ref_path = QLabel(self)
        layout.addWidget(self.ref_path, 2, 2)

        # 30Hz reference run number
        self.ref_run_number_ledit = QtWidgets.QLineEdit()
        self.ref_run_number_ledit.setValidator(QtGui.QIntValidator())
        layout.addWidget(self.ref_run_number_ledit, 3, 1)
        self.ref_run_number_label = QLabel(self)
        self.ref_run_number_label.setText("Enter a 30Hz reference run number")
        layout.addWidget(self.ref_run_number_label, 3, 2)

        # 30Hz data to process
        self.data_run_number_ledit = QtWidgets.QLineEdit()
        #self.data_run_number_ledit.setValidator(QtGui.QIntValidator())
        layout.addWidget(self.data_run_number_ledit, 4, 1)
        self.data_run_number_label = QLabel(self)
        self.data_run_number_label.setText("Enter a 30Hz run number to reduce")
        layout.addWidget(self.data_run_number_label, 4, 2)

        # Time slice
        self.time_slice_ledit = QtWidgets.QLineEdit()
        self.time_slice_ledit.setValidator(QtGui.QDoubleValidator())
        layout.addWidget(self.time_slice_ledit, 5, 1)
        self.time_slice_label = QLabel(self)
        self.time_slice_label.setText("Enter time interval in seconds")
        layout.addWidget(self.time_slice_label, 5, 2)

        # Output directory
        self.choose_output_dir = QPushButton('Output directory')
        layout.addWidget(self.choose_output_dir, 6, 1)

        self.output_dir_label = QLabel(self)
        layout.addWidget(self.output_dir_label, 6, 2)

        # Process button
        self.perform_reduction = QPushButton('Reduce')
        layout.addWidget(self.perform_reduction, 7, 1)
        #self.perform_reduction_q = QPushButton('Reduce [Const-Q binning]')
        #layout.addWidget(self.perform_reduction_q, 8, 2)
        self.load_settings = QPushButton('Load settings')
        layout.addWidget(self.load_settings, 8, 1)

        # connections
        self.choose_template.clicked.connect(self.template_selection)
        self.choose_ref.clicked.connect(self.ref_selection)
        self.choose_output_dir.clicked.connect(self.output_dir_selection)
        self.perform_reduction.clicked.connect(self.reduce_new)
        #self.perform_reduction_q.clicked.connect(self.reduce_q)
        self.load_settings.clicked.connect(self.load_settings_from_file)

        # Populate from previous session
        self.read_settings()

    def load_settings_from_file(self):
        """
            Load the reduction options from a json file produced as the output
            of an earlier reduction. This file is saved in the same directory
            as the time-resolved reflectivity curves.
        """
        _settings_file, _ = QFileDialog.getOpenFileName(self, 'Open file',
                                                        self.output_dir_label.text(),
                                                        'Settings file (*.json)')

        if os.path.isfile(_settings_file):
            with open(_settings_file, 'r') as fd:
                options = json.load(fd)

                if 'template_30Hz' in options:
                    self.template_path.setText(options['template_30Hz'])

                if 'ref_data_60Hz' in options:
                    self.ref_path.setText(options['ref_data_60Hz'])

                if 'output_dir' in options:
                    self.output_dir_label.setText(options['output_dir'])

                if 'ref_run_30Hz' in options:
                    self.ref_run_number_ledit.setText(str(options['ref_run_30Hz']))

                if 'meas_run_30Hz' in options:
                    self.data_run_number_ledit.setText(str(options['meas_run_30Hz']))

                if 'time_interval' in options:
                    self.time_slice_ledit.setText(str(options['time_interval']))

    def template_selection(self):
        _template_file, _ = QFileDialog.getOpenFileName(self, 'Open file',
                                                        self.template_path.text(),
                                                        'Template file (*.xml)')
        if os.path.isfile(_template_file):
            self.template_path.setText(_template_file)

    def ref_selection(self):
        _ref_file, _ = QFileDialog.getOpenFileName(self, 'Open file',
                                                   self.ref_path.text(),
                                                   '60Hz reference file (*.txt)')
        if os.path.isfile(_ref_file):
            self.ref_path.setText(_ref_file)

    def output_dir_selection(self):
        _dir = QFileDialog.getExistingDirectory(None, 'Select a folder:',
                                                self.output_dir_label.text(),
                                                QFileDialog.ShowDirsOnly)
        if os.path.isdir(_dir):
            self.output_dir_label.setText(_dir)

    def read_settings(self):
        _template_file = self.settings.value("30Hz_template", TEMPLATE_DIRECTIVE)
        if len(_template_file.strip()) == 0:
            _template_file = TEMPLATE_DIRECTIVE
        self.template_path.setText(_template_file)

        _ref_file = self.settings.value("30Hz_reference", REFERENCE_DIRECTIVE)
        if len(_ref_file.strip()) == 0:
            _ref_file = REFERENCE_DIRECTIVE
        self.ref_path.setText(_ref_file)

        _out_dir = self.settings.value("30Hz_output_dir", OUTPUT_DIR_DIRECTIVE)
        if len(_out_dir.strip()) == 0:
            _out_dir = OUTPUT_DIR_DIRECTIVE
        self.output_dir_label.setText(_out_dir)

        _ref_run = self.settings.value("30Hz_ref_run_number", '')
        self.ref_run_number_ledit.setText(_ref_run)

        _data_run = self.settings.value("30Hz_data_run_number", '')
        self.data_run_number_ledit.setText(_data_run)

        _interval = self.settings.value("30Hz_time_slice", '')
        self.time_slice_ledit.setText(_interval)

    def save_settings(self):
        self.settings.setValue('30Hz_template', self.template_path.text())
        self.settings.setValue('30Hz_reference', self.ref_path.text())
        self.settings.setValue('30Hz_ref_run_number', self.ref_run_number_ledit.text())
        self.settings.setValue('30Hz_data_run_number', self.data_run_number_ledit.text())
        self.settings.setValue('30Hz_time_slice', self.time_slice_ledit.text())
        self.settings.setValue('30Hz_output_dir', self.output_dir_label.text())

    def check_inputs(self):
        error = None
        # Check files and output directory
        if not os.path.isfile(self.template_path.text()):
            error = "The chosen template file could not be found"
        elif not os.path.isfile(self.ref_path.text()):
            error = "The chosen reference reflectivity file could not be found"
        elif not os.path.isdir(self.output_dir_label.text()):
            error = "The chosen output directory could not be found"

        # Pop up a dialog if there were invalid inputs
        if error:
            self.show_dialog(error)
            return False
        return True

    def show_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Critical)
        msgBox.setText(text)
        msgBox.setWindowTitle("Invalid inputs")
        msgBox.setStandardButtons(QMessageBox.Ok)

    def parse_run_list(self, text):
        """
            Parse the run list string and expand it.
        """
        run_list = []
        for _r in text.split(','):
            try:
                run_list.append(int(_r))
            except:
                sub_toks = _r.split('-')
                if len(sub_toks) == 2:
                    run_list.extend(range(int(sub_toks[0]), int(sub_toks[1])+1))

        return run_list

    def reduce_new(self):
        return self.reduce(reduction_script='scripts/time_resolved_reduction.py')

    def reduce_q(self):
        return self.reduce(reduction_script='scripts/time_resolved_reduction.py', q_summing=True)

    def reduce(self, reduction_script='scripts/time_resolved_reduction.py', q_summing=False):
        if not self.check_inputs():
            print("Invalid inputs found")
            return

        self.save_settings()

        print("Reduce!")

        run_list = self.parse_run_list(self.data_run_number_ledit.text())
        for run in run_list:
            # python3 time_resolved_reduction.py dynamic30Hz
            # <meas_run_30Hz> <ref_run_30Hz> <ref_data_60Hz> <template_30Hz> <time_interval> <output_dir>
            args = ['python3', reduction_script, 'dynamic30Hz',
                    str(run),
                    self.ref_run_number_ledit.text(),
                    self.ref_path.text(),
                    self.template_path.text(),
                    self.time_slice_ledit.text(),
                    self.output_dir_label.text()]
            if not run == run_list[-1]:
                args.append('--no-plot')
            if q_summing:
                args.append('--qsumming')
            subprocess.run(args)
