#!/usr/bin/python3
import os
import subprocess

from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import QFileDialog, QGridLayout, QLabel, QMessageBox, QPushButton, QWidget

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
        layout.addWidget(self.template_path, 1, 3)

        # First run to process
        self.first_run_number_ledit = QtWidgets.QLineEdit()
        self.first_run_number_ledit.setValidator(QtGui.QIntValidator())
        layout.addWidget(self.first_run_number_ledit, 3, 1)
        self.first_run_number_label = QLabel(self)
        self.first_run_number_label.setText("First run to process")
        layout.addWidget(self.first_run_number_label, 3, 3)

        # Last run to process
        self.last_run_number_ledit = QtWidgets.QLineEdit()
        self.last_run_number_ledit.setValidator(QtGui.QIntValidator())
        layout.addWidget(self.last_run_number_ledit, 4, 1)
        self.last_run_number_label = QLabel(self)
        self.last_run_number_label.setText("Last run to process")
        layout.addWidget(self.last_run_number_label, 4, 3)

        # Select const-q binning
        self.const_q_label = QLabel(self)
        self.const_q_label.setText("Constant-Q binning")
        self.const_q_label.setAlignment(QtCore.Qt.AlignRight)
        layout.addWidget(self.const_q_label, 6, 1)
        self.const_q_check = QtWidgets.QCheckBox()
        self.const_q_check.setChecked(False)
        layout.addWidget(self.const_q_check, 6, 2)

        # Select how to treat overlap
        self.average_overlapp_label = QLabel(self)
        self.average_overlapp_label.setText("Average overlapping points")
        self.average_overlapp_label.setAlignment(QtCore.Qt.AlignRight)
        layout.addWidget(self.average_overlapp_label, 7, 1)
        self.average_overlapp_check = QtWidgets.QCheckBox()
        self.average_overlapp_check.setChecked(True)
        layout.addWidget(self.average_overlapp_check, 7, 2)

        # Theta offset
        self.fix_offset_label = QLabel(self)
        self.fix_offset_label.setText("Theta offset")
        layout.addWidget(self.fix_offset_label, 9, 3)
        self.fix_offset_check = QtWidgets.QCheckBox()
        self.fix_offset_check.setChecked(False)
        layout.addWidget(self.fix_offset_check, 9, 2)

        self.fix_offset_ledit = QtWidgets.QLineEdit()
        self.fix_offset_ledit.setValidator(QtGui.QDoubleValidator())
        layout.addWidget(self.fix_offset_ledit, 9, 1)

        # Process button
        self.perform_reduction = QPushButton('Reduce')
        self.perform_reduction.setStyleSheet("background-color : green")
        layout.addWidget(self.perform_reduction, 10, 1)

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

        _avg = self.settings.value("reduction_avg_overlap", "true")
        self.average_overlapp_check.setChecked(_avg=='true')
        _const_q = self.settings.value("reduction_const_q", "false")
        self.const_q_check.setChecked(_const_q=='true')

        _fix_offset = self.settings.value("reduction_use_fix_offset", "false")
        self.fix_offset_check.setChecked(_fix_offset=='true')
        _fix_offset = self.settings.value("reduction_fix_offset", "0")
        self.fix_offset_ledit.setText(_fix_offset)


    def save_settings(self):
        self.settings.setValue('reduction_template', self.template_path.text())
        self.settings.setValue('reduction_first_run_number', self.first_run_number_ledit.text())
        self.settings.setValue('reduction_last_run_number', self.last_run_number_ledit.text())
        self.settings.setValue('reduction_use_old', False)
        self.settings.setValue('reduction_avg_overlap', self.average_overlapp_check.isChecked())
        self.settings.setValue('reduction_const_q', self.const_q_check.isChecked())

        self.settings.setValue('fit_first_peak', False)
        self.settings.setValue('reduction_use_fix_offset', self.fix_offset_check.isChecked())
        self.settings.setValue('reduction_fix_offset', self.fix_offset_ledit.text())

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

        options = ['python3', '/SNS/REF_L/shared/batch_reduce.py',
                   ipts,
                   self.first_run_number_ledit.text(),
                   self.last_run_number_ledit.text()]

        options.append('new')
        options.append(self.template_path.text())
        options.append(str(self.average_overlapp_check.isChecked()))
        options.append(str(self.const_q_check.isChecked()))
        options.append(False) # Was fit first peak option
        if self.fix_offset_check.isChecked():
            options.append(str(self.fix_offset_ledit.text()))

        # python3 batch_reduce.py <IPTS> <first run> <last run>
        print(' '.join(options))
        subprocess.run(options)

        self.show_dialog("Task completed: please verify!", "Task completed")
