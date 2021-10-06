#!/usr/bin/python3
import sys
import os
import subprocess

from qtpy import QtWidgets, QtGui, QtCore

from qtpy.QtWidgets import (QWidget, QGridLayout,
                            QFileDialog, QLabel,
                            QPushButton, QMessageBox)


TEMPLATE_DIRECTIVE = "Click to choose a 60Hz template"


class Reduction(QWidget):

    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle('Batch reduction')
        layout = QGridLayout()
        self.setLayout(layout)

        self.settings = QtCore.QSettings()

        # Standard template file
        self.choose_template = QPushButton('Template')
        layout.addWidget(self.choose_template, 1, 1)

        self.template_path = QLabel(self)
        layout.addWidget(self.template_path, 1, 2)

        # First run to process
        self.first_run_number_ledit = QtWidgets.QLineEdit()
        self.first_run_number_ledit.setValidator(QtGui.QIntValidator())
        layout.addWidget(self.first_run_number_ledit, 2, 1)
        self.first_run_number_label = QLabel(self)
        self.first_run_number_label.setText("First run to process")
        layout.addWidget(self.first_run_number_label, 2, 2)

        # Last run to process
        self.last_run_number_ledit = QtWidgets.QLineEdit()
        self.last_run_number_ledit.setValidator(QtGui.QIntValidator())
        layout.addWidget(self.last_run_number_ledit, 3, 1)
        self.last_run_number_label = QLabel(self)
        self.last_run_number_label.setText("Last run to process")
        layout.addWidget(self.last_run_number_label, 3, 2)

        # Process button
        self.perform_reduction = QPushButton('Reduce')
        layout.addWidget(self.perform_reduction, 5, 1)

        # connections
        self.choose_template.clicked.connect(self.template_selection)
        self.perform_reduction.clicked.connect(self.reduce)

        # Populate from previous session
        self.read_settings()

    def template_selection(self):
        _template_file, _ = QFileDialog.getOpenFileName(self, 'Open file',
                                                        self.template_path.text(),
                                                        'Template file (*.xml)')
        if os.path.isfile(_template_file):
            self.template_path.setText(_template_file)

    def read_settings(self):
        _template_file = self.settings.value("reduction_template", TEMPLATE_DIRECTIVE)
        if len(_template_file.strip()) == 0:
            _template_file = TEMPLATE_DIRECTIVE
        self.template_path.setText(_template_file)

        _first_run = self.settings.value("reduction_first_run_number", '')
        self.first_run_number_ledit.setText(_first_run)

        _interval = self.settings.value("reduction_last_run_number", '')
        self.last_run_number_ledit.setText(_interval)

    def save_settings(self):
        self.settings.setValue('reduction_template', self.template_path.text())
        self.settings.setValue('reduction_first_run_number', self.first_run_number_ledit.text())
        self.settings.setValue('reduction_last_run_number', self.last_run_number_ledit.text())

    def check_inputs(self):
        error = None
        # Check files and output directory
        if not os.path.isfile(self.template_path.text()):
            error = "The chosen template file could not be found"

        # Pop up a dialog if there were invalid inputs
        if error:
            self.show_dialog(error)
            return False
        return True

    def show_dialog(self, text, title="Invalid inputs"):
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Critical)
        msgBox.setText(text)
        msgBox.setWindowTitle(title)
        msgBox.setStandardButtons(QMessageBox.Ok)

    def reduce(self):
        if not self.check_inputs():
            print("Invalid inputs found")
            return

        # Get the IPTS from the template file
        toks = self.template_path.text().split('/')
        if toks[3].startswith('IPTS'):
            ipts = toks[3]
        else:
            self.show_dialog("The chosen template is not in an IPTS folder, please select another one.")
            return

        self.save_settings()

        print("Reduce!")

        # python3 batch_reduce.py <IPTS> <first run> <last run>
        subprocess.run(['python3', '/SNS/REF_L/shared/batch_reduce.py',
                        ipts,
                        self.first_run_number_ledit.text(),
                        self.last_run_number_ledit.text()])

        self.show_dialog("Task completed: please verify!", "Task completed")
