#!/usr/bin/env python
import sys

from apps.direct_beam import DirectBeamTab
from apps.overplot import Overplot
from apps.roi_selector import ROISelector
from apps.sld_calculator import SLD
from apps.template_batch import TemplateBatchTab
from apps.file_batch import FileBatchTab
from qtpy import QtCore
from qtpy.QtWidgets import QApplication, QGridLayout, QTabWidget, QWidget

#REFERENCE_DIRECTIVE = "Click to choose a 60Hz reference R(Q) file"
#TEMPLATE_DIRECTIVE = "Click to choose a 30Hz template"
#OUTPUT_DIR_DIRECTIVE = "Click to choose an output directory"

class ReductionInterface(QTabWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle("New Reflectometry Launcher")
        layout = QGridLayout()
        self.setLayout(layout)

        self.settings = QtCore.QSettings()

        # Overplot tab
        tab_id = 0
        self.overplot_tab = Overplot()
        self.addTab(self.overplot_tab, "Overplot")
        self.setTabText(tab_id, "Overplot")

        # Direct beam processing tab
        tab_id += 1
        self.direct_beam_tab = DirectBeamTab()
        self.addTab(self.direct_beam_tab, "Direct beam")
        self.setTabText(tab_id, "Direct beam")

        # Batch file-driven reduction tab (DAT/JSON)
        tab_id += 1
        self.file_batch_tab = FileBatchTab()
        self.addTab(self.file_batch_tab, "Batch file")
        self.setTabText(tab_id, "Batch file")
        
        
        # ROI selector
        tab_id += 1
        self.roi_tab = ROISelector()
        self.addTab(self.roi_tab, "ROI selector")
        self.setTabText(tab_id, "ROI selector")

        # SLD calculator
        tab_id += 1
        self.sld_tab = SLD()
        self.addTab(self.sld_tab, "SLD calculator")
        self.setTabText(tab_id, "SLD calculator")

        # Batch template reduction tab
        tab_id += 1
        self.template_batch_tab = TemplateBatchTab()
        self.addTab(self.template_batch_tab, "Batch template")
        self.setTabText(tab_id, "Batch template")

if __name__ == "__main__":
    app = QApplication([])
    window = ReductionInterface()
    window.show()
    sys.exit(app.exec_())
