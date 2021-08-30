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

        # 60Hz time-resolved
        self.time_60Hz_tab = Dynamic60Hz()
        self.addTab(self.time_60Hz_tab,"Time-resolved 60Hz")
        self.setTabText(0,"Time-resolved 60Hz")

        self.time_30Hz_tab = Dynamic30Hz()
        self.addTab(self.time_30Hz_tab,"Time-resolved 30Hz")
        self.setTabText(1,"Time-resolved 30Hz")


if __name__ == '__main__':
    app = QApplication([])
    window = ReductionInterface()
    window.show()
    sys.exit(app.exec_())
