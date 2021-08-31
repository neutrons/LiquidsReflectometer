#!/usr/bin/python3
import sys
import os
import subprocess

from qtpy import QtWidgets, QtGui, QtCore

from qtpy.QtWidgets import (QWidget, QGridLayout,
                            QFileDialog, QLabel,
                            QPushButton, QMessageBox)


TEMPLATE_DIRECTIVE = "Click to choose a 60Hz template"
OUTPUT_DIR_DIRECTIVE = "Click to choose an output directory"


class Dynamic60Hz(QWidget):

    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle('Time-resolved reduction')
        layout = QGridLayout()
        self.setLayout(layout)

        self.settings = QtCore.QSettings()

        # 30Hz template file
        self.choose_template = QPushButton('template')
        layout.addWidget(self.choose_template, 1, 1)

        self.template_path = QLabel(self)
        layout.addWidget(self.template_path, 1, 2)

        # 30Hz data to process
        self.data_run_number_ledit = QtWidgets.QLineEdit()
        self.data_run_number_ledit.setValidator(QtGui.QIntValidator())
        layout.addWidget(self.data_run_number_ledit, 4, 1)
        self.data_run_number_label = QLabel(self)
        self.data_run_number_label.setText("Enter a run number to reduce")
        layout.addWidget(self.data_run_number_label, 4, 2)

        # Time slice
        self.time_slice_ledit = QtWidgets.QLineEdit()
        self.time_slice_ledit.setValidator(QtGui.QDoubleValidator())
        layout.addWidget(self.time_slice_ledit, 5, 1)
        self.time_slice_label = QLabel(self)
        self.time_slice_label.setText("Enter time interval in seconds")
        layout.addWidget(self.time_slice_label, 5, 2)

        # Scan index
        self.idx_ledit = QtWidgets.QLineEdit()
        self.idx_ledit.setValidator(QtGui.QIntValidator())
        layout.addWidget(self.idx_ledit, 6, 1)
        self.idx_label = QLabel(self)
        self.idx_label.setText("Enter time the scan index to use within the template")
        layout.addWidget(self.idx_label, 6, 2)

        # Output directory
        self.choose_output_dir = QPushButton('Output directory')
        layout.addWidget(self.choose_output_dir, 7, 1)

        self.output_dir_label = QLabel(self)
        layout.addWidget(self.output_dir_label, 7, 2)

        # Process button
        self.perform_reduction = QPushButton('Reduce')
        layout.addWidget(self.perform_reduction, 8, 1)

        # connections
        self.choose_template.clicked.connect(self.template_selection)
        self.choose_output_dir.clicked.connect(self.output_dir_selection)
        self.perform_reduction.clicked.connect(self.reduce)

        # Populate from previous session
        self.read_settings()

    def template_selection(self):
        _template_file, _ = QFileDialog.getOpenFileName(self, 'Open file',
                                                        '', 'Template file (*.xml)')
        if os.path.isfile(_template_file):
            self.template_path.setText(_template_file)

    def output_dir_selection(self):
        _dir = QFileDialog.getExistingDirectory(None, 'Select a folder:', '', QFileDialog.ShowDirsOnly)
        if os.path.isdir(_dir):
            self.output_dir_label.setText(_dir)

    def read_settings(self):
        _template_file = self.settings.value("60Hz_template", TEMPLATE_DIRECTIVE)
        if len(_template_file.strip()) == 0:
            _template_file = TEMPLATE_DIRECTIVE
        self.template_path.setText(_template_file)

        _out_dir = self.settings.value("output_dir", OUTPUT_DIR_DIRECTIVE)
        if len(_out_dir.strip()) == 0:
            _out_dir = OUTPUT_DIR_DIRECTIVE
        self.output_dir_label.setText(_out_dir)

        _data_run = self.settings.value("60Hz_data_run_number", '')
        self.data_run_number_ledit.setText(_data_run)

        _interval = self.settings.value("60Hz_time_slice", '')
        self.time_slice_ledit.setText(_interval)

        _interval = self.settings.value("60Hz_scan_index", '')
        self.idx_ledit.setText(_interval)

    def save_settings(self):
        self.settings.setValue('60Hz_template', self.template_path.text())
        self.settings.setValue('60Hz_data_run_number', self.data_run_number_ledit.text())
        self.settings.setValue('60Hz_time_slice', self.time_slice_ledit.text())
        self.settings.setValue('60Hz_scan_index', self.idx_ledit.text())
        self.settings.setValue('output_dir', self.output_dir_label.text())

    def check_inputs(self):
        error = None
        # Check files and output directory
        if not os.path.isfile(self.template_path.text()):
            error = "The chosen template file could not be found"
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

    def reduce(self):
        if not self.check_inputs():
            print("Invalid inputs found")
            return

        self.save_settings()

        print("Reduce!")

        # python3 template_reduction.py dynamic60Hz <meas_run_60Hz> <template_60Hz> <time_interval> <output_dir>
        subprocess.run(['python3', 'scripts/template_reduction.py', 'dynamic60Hz',
                        self.data_run_number_ledit.text(),
                        self.template_path.text(),
                        self.time_slice_ledit.text(), self.output_dir_label.text(),
                        '--scan_index', self.idx_ledit.text()])

