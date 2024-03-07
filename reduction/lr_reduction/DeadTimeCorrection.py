"""
    Dead time correction algorithm for single-readout detectors.
"""
import time
import math
import os
from mantid.api import *
from mantid.simpleapi import *
from mantid.kernel import *
import numpy as np
import scipy

def call(InputWorkspace, DeadTime=4.2, TOFStep=100, Paralyzable=False, TOFRange=[0, 0], OutputWorkspace='correction'):
    """
        Function to make the algorithm call similar to a normal Mantid call
    """
    algo = SingleReadoutDeadTimeCorrection()
    algo.PyInit()
    algo.setProperty("InputWorkspace", InputWorkspace)
    algo.setProperty("DeadTime", DeadTime)
    algo.setProperty("TOFStep", TOFStep)
    algo.setProperty("Paralyzable", Paralyzable)
    algo.setProperty("TOFRange", TOFRange)
    algo.setProperty("OutputWorkspace", OutputWorkspace)
    algo.PyExec()
    return algo.getProperty('OutputWorkspace').value


class SingleReadoutDeadTimeCorrection(PythonAlgorithm):

    def category(self):
        return "Reflectometry\\SNS"

    def name(self):
        return "SingleReadoutDeadTimeCorrection"

    def version(self):
        return 1

    def summary(self):
        return "Single read-out dead time correction calculation"

    def PyInit(self):
        self.declareProperty(WorkspaceProperty("InputWorkspace", "", Direction.Input),
                             "Input workspace use to compute dead time correction")
        self.declareProperty("DeadTime", 4.2, doc="Dead time in microseconds")
        self.declareProperty("TOFStep", 100,
                             doc="TOF bins to compute deadtime correction for, in microseconds")
        self.declareProperty("Paralyzable", False,
                             doc="If true, paralyzable correction will be applied, non-paralyzing otherwise")
        self.declareProperty(FloatArrayProperty("TOFRange", [0., 0.],
                                                FloatArrayLengthValidator(2), direction=Direction.Input),
                             "TOF range to use")
        self.declareProperty(MatrixWorkspaceProperty("OutputWorkspace", "", Direction.Output), "Output workspace")

    def PyExec(self):
        # Event data must include error events (all triggers on the detector)
        ws_event_data = self.getProperty("InputWorkspace").value
        dead_time = self.getProperty("DeadTime").value
        tof_step = self.getProperty("TOFStep").value
        paralyzing = self.getProperty("Paralyzable").value

        # Rebin the data according to the tof_step we want to compute the correction with
        tof_min, tof_max = self.getProperty("TOFRange").value
        if tof_min == 0 and tof_max == 0:
            tof_min = ws_event_data.getTofMin()
            tof_max = ws_event_data.getTofMax()
        logger.notice("TOF range: %f %f" % (tof_min, tof_max))
        _ws_sc = Rebin(InputWorkspace=ws_event_data, Params="%s,%s,%s" % (tof_min, tof_step, tof_max), PreserveEvents=False)

        # Get the total number of counts on the detector for each TOF bin per pulse
        counts_ws = SumSpectra(_ws_sc)
        t_series = np.asarray(_ws_sc.getRun()['proton_charge'].value)
        non_zero = t_series > 0
        n_pulses = np.count_nonzero(non_zero)
        rate = counts_ws.readY(0) / n_pulses
        tof_bins = counts_ws.readX(0)

         # Compute the dead time correction for each TOF bin
        if paralyzing:
            true_rate = -scipy.special.lambertw(-rate * dead_time / tof_step).real / dead_time
            corr = true_rate / (rate / tof_step)
            # If we have no events, set the correction to 1 otherwise we will get a nan
            # from the equation above.
            corr[rate==0] = 1
        else:
            corr = 1/(1-rate * dead_time / tof_step)

        if np.min(corr) < 0:
            error = ( "Corrupted dead time correction:\n"
                      +"  Reflected: %s\n" % corr )
            logger.error(error)

        counts_ws.setY(0, corr)
        self.setProperty('OutputWorkspace', counts_ws)


AlgorithmFactory.subscribe(SingleReadoutDeadTimeCorrection)