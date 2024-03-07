import os
import numpy as np
from lr_reduction.DeadTimeCorrection import SingleReadoutDeadTimeCorrection

import mantid
import mantid.simpleapi as mtd_api
mtd_api.config["default.facility"] = "SNS"
mtd_api.config["default.instrument"] = "REF_L"


def test_deadtime():
    """
        Test the time-resolved reduction that uses a measured reference.
        It is generally used at 30 Hz but it also works at 60 Hz.
    """
    ws = mtd_api.Load("REF_L_198409")

    algo = SingleReadoutDeadTimeCorrection()
    algo.PyInit()
    algo.setProperty("InputWorkspace", ws)
    algo.setProperty("OutputWorkspace", "dead_time_corr")
    algo.PyExec()
    corr_ws = algo.getProperty('OutputWorkspace').value
    corr = corr_ws.readY(0)
    for c in corr:
        assert(c>0)
        assert(c<1.001)

def test_deadtime_paralyzable():
    """
        Test the time-resolved reduction that uses a measured reference.
        It is generally used at 30 Hz but it also works at 60 Hz.
    """
    ws = mtd_api.Load("REF_L_198409")

    algo = SingleReadoutDeadTimeCorrection()
    algo.PyInit()
    algo.setProperty("InputWorkspace", ws)
    algo.setProperty("Paralyzable", True)
    algo.setProperty("OutputWorkspace", "dead_time_corr")
    algo.PyExec()
    corr_ws = algo.getProperty('OutputWorkspace').value
    corr = corr_ws.readY(0)
    with open("dc.txt", 'w') as fd:
        fd.write(str(corr))
    for c in corr:
        assert(c>0)
        assert(c<1.001)
