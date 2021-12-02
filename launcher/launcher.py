#!/usr/bin/python3
import sys
import os
import subprocess

from qtpy import QtWidgets, QtGui, QtCore

from qtpy.QtWidgets import (QApplication, QWidget, QGridLayout,
                            QTabWidget, QFileDialog, QLabel,
                            QPushButton, QMessageBox)

from apps.dynamic_60Hz import Dynamic60Hz
from apps.dynamic_30Hz import Dynamic30Hz
from apps.xrr import XRR
from apps.reduction import Reduction
from apps.refracted import Refracted

REFERENCE_DIRECTIVE = "Click to choose a 60Hz reference R(Q) file"
TEMPLATE_DIRECTIVE = "Click to choose a 30Hz template"
OUTPUT_DIR_DIRECTIVE = "Click to choose an output directory"


class ReductionInterface(QTabWidget):

    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle('Reflectometry Launcher')
        layout = QGridLayout()
        self.setLayout(layout)

        self.settings = QtCore.QSettings()

        # Batch reduction
        tab_id = 0
        self.time_60Hz_tab = Reduction()
        self.addTab(self.time_60Hz_tab, "Batch reduction")
        self.setTabText(tab_id, "Batch reduction")

        # 60Hz time-resolved
        tab_id += 1
        self.time_60Hz_tab = Dynamic60Hz()
        self.addTab(self.time_60Hz_tab, "Time-resolved 60Hz")
        self.setTabText(tab_id, "Time-resolved 60Hz")

        # 30 Hz time-resolved
        tab_id += 1
        self.time_30Hz_tab = Dynamic30Hz()
        self.addTab(self.time_30Hz_tab, "Time-resolved 30Hz")
        self.setTabText(tab_id, "Time-resolved 30Hz")

        # XRR reduction
        tab_id += 1
        self.xrr_tab = XRR()
        self.addTab(self.xrr_tab, "XRR processing")
        self.setTabText(tab_id, "XRR processing")

        # Refracted beam analysis
        tab_id += 1
        self.refracted_tab = Refracted()
        self.addTab(self.refracted_tab, "Refraction analysis")
        self.setTabText(tab_id, "Refraction analysis")

if __name__ == '__main__':
    app = QApplication([])
    window = ReductionInterface()
    window.show()
    sys.exit(app.exec_())
