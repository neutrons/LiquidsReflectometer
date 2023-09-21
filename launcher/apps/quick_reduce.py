#!/usr/bin/python3
import os
import subprocess

from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import QFileDialog, QGridLayout, QLabel, QMessageBox, QPushButton, QSpacerItem, QWidget

DATA_FILE_DIRECTIVE = "Click to choose a file to process"
OUTPUT_DIR_DIRECTIVE = os.path.expanduser("~")


class QuickReduce(QWidget):

    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle('Quick reduce')
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

        # Estimated pixel
        self.peak_pixel_ledit = QtWidgets.QLineEdit()
        layout.addWidget(self.peak_pixel_ledit, 2, 1)
        self.peak_pixel_label = QLabel(self)
        self.peak_pixel_label.setText("Peak pixel estimate")
        layout.addWidget(self.peak_pixel_label, 2, 2)

        # DB run number
        self.db_run_number_ledit = QtWidgets.QLineEdit()
        self.db_run_number_ledit.setValidator(QtGui.QIntValidator())
        layout.addWidget(self.db_run_number_ledit, 3, 1)
        self.db_run_number_label = QLabel(self)
        self.db_run_number_label.setText("DB run number")
        layout.addWidget(self.db_run_number_label, 3, 2)

        # Estimated pixel
        self.db_peak_pixel_ledit = QtWidgets.QLineEdit()
        layout.addWidget(self.db_peak_pixel_ledit, 4, 1)
        self.db_peak_pixel_label = QLabel(self)
        self.db_peak_pixel_label.setText("DB peak pixel estimate")
        layout.addWidget(self.db_peak_pixel_label, 4, 2)

        # Output directory
        self.choose_output_dir = QPushButton('Output directory')
        layout.addWidget(self.choose_output_dir, 5, 1)

        self.output_dir_label = QLabel(self)
        layout.addWidget(self.output_dir_label, 5, 2)

        spacer = QSpacerItem(20, 20, QtWidgets.QSizePolicy.Minimum,
                             QtWidgets.QSizePolicy.Minimum)
        layout.addItem(spacer, 6, 1)

        # Process button
        self.perform_reduction = QPushButton('Process')
        self.perform_reduction.setStyleSheet("background-color : green")
        layout.addWidget(self.perform_reduction, 7, 1)

        spacer = QSpacerItem(10, 10, QtWidgets.QSizePolicy.Minimum,
                             QtWidgets.QSizePolicy.Expanding)
        layout.addItem(spacer, 8, 1)

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
        _run_number = self.settings.value("quick_run_number", '')
        self.run_number_ledit.setText(_run_number)

        _db_run_number = self.settings.value("quick_db_run_number", '')
        self.db_run_number_ledit.setText(_db_run_number)

        _pixel = self.settings.value("quick_pixel", 145)
        self.peak_pixel_ledit.setText(str(_pixel))

        _db_pixel = self.settings.value("quick_db_pixel", 145)
        self.db_peak_pixel_ledit.setText(str(_db_pixel))

        _out_dir = self.settings.value("db_output_dir", OUTPUT_DIR_DIRECTIVE)
        if len(_out_dir.strip()) == 0:
            _out_dir = OUTPUT_DIR_DIRECTIVE
        self.output_dir_label.setText(_out_dir)

    def save_settings(self):
        self.settings.setValue('quick_run_number', self.run_number_ledit.text())
        self.settings.setValue('quick_db_run_number', self.db_run_number_ledit.text())
        self.settings.setValue('quick_pixel', self.peak_pixel_ledit.text())
        self.settings.setValue('quick_db_pixel', self.db_peak_pixel_ledit.text())

        self.settings.setValue('quick_output_dir', self.output_dir_label.text())

    def check_inputs(self):
        error = None
        # Check files and output directory
        if not os.path.isdir(self.output_dir_label.text()):
            error = "The chosen output directory could not be found"

        try:
            int(self.run_number_ledit.text())
            int(self.db_run_number_ledit.text())
        except:
            error = "Check your run numbers"
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

        subprocess.run(['nsd-conda-wrap.sh', 'mantid', 'scripts/quick_reduce.py',
                        self.run_number_ledit.text(), self.db_run_number_ledit.text(),
                        self.peak_pixel_ledit.text(), self.db_peak_pixel_ledit.text(), self.output_dir_label.text()])
