#!/usr/bin/python3
import sys
import periodictable.nsf as nsf
import periodictable.xsf as xsf

from qtpy import QtWidgets, QtGui, QtCore

from qtpy.QtWidgets import (QWidget, QGridLayout,
                            QLabel, QPushButton, QMessageBox)


class SLD(QWidget):

    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle('Rigaku XRR reduction')
        layout = QGridLayout()
        self.setLayout(layout)

        self.settings = QtCore.QSettings()

        self.composition_ledit = QtWidgets.QLineEdit()
        layout.addWidget(self.composition_ledit, 1, 1)
        self.composition_label = QLabel(self)
        self.composition_label.setText("Enter composition")
        layout.addWidget(self.composition_label, 1, 2)

        # Density
        self.density_ledit = QtWidgets.QLineEdit()
        self.density_ledit.setValidator(QtGui.QDoubleValidator())
        layout.addWidget(self.density_ledit, 2, 1)
        self.density_label = QLabel(self)
        self.density_label.setText("Enter density [g/cm3]")
        layout.addWidget(self.density_label, 2, 2)

        # Wavelength [A]
        self.wl_ledit = QtWidgets.QLineEdit()
        self.wl_ledit.setValidator(QtGui.QDoubleValidator())
        layout.addWidget(self.wl_ledit, 3, 1)
        self.wl_label = QLabel(self)
        self.wl_label.setText("Enter wavelength [A]")
        layout.addWidget(self.wl_label, 3, 2)

        # Wavelength [A]
        self.output = QtWidgets.QTextEdit()
        layout.addWidget(self.output, 4, 1)
        self.output.setReadOnly(True)

        # Process button
        self.calculate = QPushButton('Calculate')
        layout.addWidget(self.calculate, 5, 1)

        # connections
        self.calculate.clicked.connect(self.compute_sld)

        # Populate from previous session
        self.read_settings()

    def read_settings(self):
        _composition = self.settings.value("sld_composition", 'Si')
        self.composition_ledit.setText(_composition)
        _wl = self.settings.value("sld_wavelength", '1.54')
        self.wl_ledit.setText(_wl)

    def save_settings(self):
        self.settings.setValue('sld_composition', self.composition_ledit.text())
        self.settings.setValue('sld_wavelength', self.wl_ledit.text())

    def show_dialog(self, text):
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Critical)
        msgBox.setText(text)
        msgBox.setWindowTitle("Invalid inputs")
        msgBox.setStandardButtons(QMessageBox.Ok)

    def compute_sld(self):
        self.save_settings()
        density = self.density_ledit.text()
        try:
            density = float(density)
        except:
            density = None
        wavelength = float(self.wl_ledit.text())
        composition = self.composition_ledit.text()

        try:
            sld, im_sld, incoh = nsf.neutron_sld(compound=composition, wavelength=wavelength, density=density)
            x_sld, x_im_sld = xsf.xray_sld(compound=composition, wavelength=wavelength, density=density)

            output_text =      "%-15s %6.6f\n" % ("Neutron SLD:", sld)
            output_text +=     "%-15s %6.6f\n" % ("   Imag SLD:", im_sld)
            output_text += "%-15s %6.6f\n\n" % ("   Incoh SLD:", incoh)

            output_text += "%-15s %6.6f\n" % ("X-ray SLD:", x_sld)
            output_text += "%-15s %6.6f\n\n" % (" Imag SLD:", x_im_sld)
            output_text += "All units in 10^-6 A^-2"
        except:
            output_text = str(sys.exc_info()[1])

        self.output.setText(output_text)
