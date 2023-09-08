#!/usr/bin/python3
import os
import subprocess
import sys

from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import (
    QFileDialog,
    QGridLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QWidget,
)

OUTPUT_DIR_DIRECTIVE = "Click to choose an output directory"
DEFAULT_MATERIAL = "Si"


class Refracted(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle("Refracted beam")
        layout = QGridLayout()
        self.setLayout(layout)

        self.settings = QtCore.QSettings()

        # Run to process
        self.run_number_ledit = QtWidgets.QLineEdit()
        self.run_number_ledit.setValidator(QtGui.QIntValidator())
        layout.addWidget(self.run_number_ledit, 1, 1)
        self.run_number_label = QLabel(self)
        self.run_number_label.setText("Run to process")
        layout.addWidget(self.run_number_label, 1, 2)

        # Output directory
        self.choose_output_dir = QPushButton("Output directory")
        layout.addWidget(self.choose_output_dir, 2, 1)

        self.output_dir_label = QLabel(self)
        layout.addWidget(self.output_dir_label, 2, 2)

        # Run to process
        self.material_ledit = QtWidgets.QLineEdit()
        layout.addWidget(self.material_ledit, 3, 1)
        self.material_label = QLabel(self)
        self.material_label.setText("Material [Si or Quartz]")
        layout.addWidget(self.material_label, 3, 2)

        # Process button
        self.perform_reduction = QPushButton("Process")
        layout.addWidget(self.perform_reduction, 4, 1)

        # connections
        self.choose_output_dir.clicked.connect(self.output_dir_selection)
        self.perform_reduction.clicked.connect(self.reduce)

        # Populate from previous session
        self.read_settings()

    def output_dir_selection(self):
        _dir = QFileDialog.getExistingDirectory(
            None,
            "Select a folder:",
            self.output_dir_label.text(),
            QFileDialog.ShowDirsOnly,
        )
        if os.path.isdir(_dir):
            self.output_dir_label.setText(_dir)

    def read_settings(self):
        _out_dir = self.settings.value("refracted_output_dir", OUTPUT_DIR_DIRECTIVE)
        if len(_out_dir.strip()) == 0:
            _out_dir = OUTPUT_DIR_DIRECTIVE
        self.output_dir_label.setText(_out_dir)

        _run = self.settings.value("refracted_run_number", "")
        self.run_number_ledit.setText(_run)

        _material = self.settings.value("refracted_material", "Si")
        if len(_material.strip()) == 0:
            _material = "Si"
        self.material_ledit.setText(_material)

    def save_settings(self):
        self.settings.setValue("refracted_output_dir", self.output_dir_label.text())
        self.settings.setValue("refracted_run_number", self.run_number_ledit.text())
        self.settings.setValue("refracted_material", self.material_ledit.text())

    def check_inputs(self):
        error = None
        # Check output directory
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

        print("Process!")

        # python3 template_reduction.py dynamic60Hz <meas_run_60Hz> <template_60Hz> <time_interval> <output_dir>
        subprocess.run(
            [
                "python3",
                "scripts/refracted_beam.py",
                self.run_number_ledit.text(),
                self.output_dir_label.text(),
                "--material",
                self.material_ledit.text(),
            ]
        )
