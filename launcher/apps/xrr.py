#!/usr/bin/python3
import os
import subprocess

from qtpy import QtCore
from qtpy.QtWidgets import QFileDialog, QGridLayout, QLabel, QMessageBox, QPushButton, QWidget

DATA_FILE_DIRECTIVE = "Click to choose a file to process"
OUTPUT_DIR_DIRECTIVE = "Click to choose an output directory"


class XRR(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle("Rigaku XRR reduction")
        layout = QGridLayout()
        self.setLayout(layout)

        self.settings = QtCore.QSettings()

        # Data file
        self.choose_data = QPushButton("Data file")
        layout.addWidget(self.choose_data, 1, 1)

        self.data_path = QLabel(self)
        layout.addWidget(self.data_path, 1, 2)

        # Output directory
        self.choose_output_dir = QPushButton("Output directory")
        layout.addWidget(self.choose_output_dir, 2, 1)

        self.output_dir_label = QLabel(self)
        layout.addWidget(self.output_dir_label, 2, 2)

        # Process button
        self.perform_reduction = QPushButton("Process")
        layout.addWidget(self.perform_reduction, 3, 1)

        # connections
        self.choose_data.clicked.connect(self.data_selection)
        self.choose_output_dir.clicked.connect(self.output_dir_selection)
        self.perform_reduction.clicked.connect(self.reduce)

        # Populate from previous session
        self.read_settings()

    def data_selection(self):
        _data_file, _ = QFileDialog.getOpenFileName(
            self, "Open file", self.data_path.text(), "Rigaku data file (*.ras)"
        )
        if os.path.isfile(_data_file):
            self.data_path.setText(_data_file)

    def output_dir_selection(self):
        _dir = QFileDialog.getExistingDirectory(
            None, "Select a folder:", self.output_dir_label.text(), QFileDialog.ShowDirsOnly
        )
        if os.path.isdir(_dir):
            self.output_dir_label.setText(_dir)

    def read_settings(self):
        _data_file = self.settings.value("xrr_data_file", DATA_FILE_DIRECTIVE)
        if len(_data_file.strip()) == 0:
            _data_file = DATA_FILE_DIRECTIVE
        self.data_path.setText(_data_file)

        _out_dir = self.settings.value("output_dir", OUTPUT_DIR_DIRECTIVE)
        if len(_out_dir.strip()) == 0:
            _out_dir = OUTPUT_DIR_DIRECTIVE
        self.output_dir_label.setText(_out_dir)

    def save_settings(self):
        self.settings.setValue("xrr_data_file", self.data_path.text())
        self.settings.setValue("output_dir", self.output_dir_label.text())

    def check_inputs(self):
        error = None
        # Check files and output directory
        if not os.path.isfile(self.data_path.text()):
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

        subprocess.run(["python3", "scripts/xrr_reduction.py", self.data_path.text(), self.output_dir_label.text()])
