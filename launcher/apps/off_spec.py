#!/usr/bin/python3
import os
import subprocess

from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import QFileDialog, QGridLayout, QLabel, QMessageBox, QPushButton, QSpacerItem, QWidget

DATA_FILE_DIRECTIVE = "Click to choose a file to process"
OUTPUT_DIR_DIRECTIVE = "Click to choose an output directory"


class OffSpec(QWidget):

    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle('Export off-spec data')
        layout = QGridLayout()
        layout.setColumnStretch(1, 0)
        layout.setColumnStretch(2, 1)
        self.setLayout(layout)

        self.settings = QtCore.QSettings()

        # Run number
        self.run_number_ledit = QtWidgets.QLineEdit()
        self.run_number_ledit.setValidator(QtGui.QIntValidator())
        layout.addWidget(self.run_number_ledit, 1, 1)
        self.run_number_label = QLabel(self)
        self.run_number_label.setText("Run number")
        layout.addWidget(self.run_number_label, 1, 2)

        # Time slice
        self.wl_step_ledit = QtWidgets.QLineEdit()
        self.wl_step_ledit.setValidator(QtGui.QDoubleValidator())
        layout.addWidget(self.wl_step_ledit, 2, 1)
        self.wl_step_label = QLabel(self)
        self.wl_step_label.setText("Enter wavelength step [Angstrom]")
        layout.addWidget(self.wl_step_label, 2, 2)

        # Output directory
        self.choose_output_dir = QPushButton('Output directory')
        layout.addWidget(self.choose_output_dir, 3, 1)

        self.output_dir_label = QLabel(self)
        layout.addWidget(self.output_dir_label, 3, 2)

        spacer = QSpacerItem(20, 20, QtWidgets.QSizePolicy.Minimum,
                             QtWidgets.QSizePolicy.Minimum)
        layout.addItem(spacer, 4, 1)

        # Process button
        self.perform_reduction = QPushButton('Process')
        self.perform_reduction.setStyleSheet("background-color : green")
        layout.addWidget(self.perform_reduction, 5, 1)

        spacer = QSpacerItem(10, 10, QtWidgets.QSizePolicy.Minimum,
                             QtWidgets.QSizePolicy.Expanding)
        layout.addItem(spacer, 6, 1)

        # connections
        self.choose_output_dir.clicked.connect(self.output_dir_selection)
        self.perform_reduction.clicked.connect(self.reduce)

        # Populate from previous session
        self.read_settings()

    def output_dir_selection(self):
        _dir = QFileDialog.getExistingDirectory(None, 'Select a folder:',
                                                self.output_dir_label.text(),
                                                QFileDialog.ShowDirsOnly)
        if os.path.isdir(_dir):
            self.output_dir_label.setText(_dir)

    def read_settings(self):
        _run_number = self.settings.value("offspec_run_number", '')
        self.run_number_ledit.setText(_run_number)

        _wl_step = self.settings.value("offspec_wl_step", '')
        self.wl_step_ledit.setText(_wl_step)

        _out_dir = self.settings.value("offspec_output_dir", OUTPUT_DIR_DIRECTIVE)
        if len(_out_dir.strip()) == 0:
            _out_dir = OUTPUT_DIR_DIRECTIVE
        self.output_dir_label.setText(_out_dir)

    def save_settings(self):
        self.settings.setValue('offspec_run_number', self.run_number_ledit.text())
        self.settings.setValue('offspec_wl_step', self.wl_step_ledit.text())
        self.settings.setValue('offspec_output_dir', self.output_dir_label.text())

    def check_inputs(self):
        error = None
        # Check files and output directory
        if not os.path.isdir(self.output_dir_label.text()):
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

        print("Processing!")

        #subprocess.run(['python3', 'scripts/off_spec.py',
        #                self.run_number_ledit.text(), self.wl_step_ledit.text(), self.output_dir_label.text()])

        reduction_script='scripts/off_spec.py'
        args = ['nsd-conda-wrap.sh', 'mantid',
                '--classic',
                reduction_script,
                self.run_number_ledit.text(),
                self.wl_step_ledit.text(),
                self.output_dir_label.text()
                ]

        subprocess.run(args)
