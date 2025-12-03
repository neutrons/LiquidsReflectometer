import os

import mantid.simpleapi as mtd_api

from lr_reduction import template
from lr_reduction.dead_time_correction import SingleReadoutDeadTimeCorrection
from lr_reduction.utils import amend_config

mtd_api.config["default.facility"] = "SNS"
mtd_api.config["default.instrument"] = "REF_L"


def test_deadtime(nexus_dir):
    """
    Test the time-resolved reduction that uses a measured reference.
    It is generally used at 30 Hz but it also works at 60 Hz.
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_198409")

    algo = SingleReadoutDeadTimeCorrection()
    algo.PyInit()
    algo.setProperty("InputWorkspace", ws)
    algo.setProperty("OutputWorkspace", "dead_time_corr")
    algo.PyExec()
    corr_ws = algo.getProperty("OutputWorkspace").value
    corr = corr_ws.readY(0)
    for c in corr:
        assert c > 0
        assert c < 1.001


def test_deadtime_paralyzable(nexus_dir):
    """
    Test the time-resolved reduction that uses a measured reference.
    It is generally used at 30 Hz but it also works at 60 Hz.
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_198409")

    algo = SingleReadoutDeadTimeCorrection()
    algo.PyInit()
    algo.setProperty("InputWorkspace", ws)
    algo.setProperty("Paralyzable", True)
    algo.setProperty("OutputWorkspace", "dead_time_corr")
    algo.PyExec()
    corr_ws = algo.getProperty("OutputWorkspace").value
    corr = corr_ws.readY(0)
    for c in corr:
        assert c > 0
        assert c < 1.001

def test_deadtime_threshold(nexus_dir):
    """
    Test using the threshold. Here the threshold is set to 0,
    so all corrections will be 0
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_198409")

    algo = SingleReadoutDeadTimeCorrection()
    algo.PyInit()
    algo.setProperty("InputWorkspace", ws)
    algo.setProperty("OutputWorkspace", "dead_time_corr")

    # set the deadtime threshold
    # the mean corr is 1.0003, so use that as test value
    algo.setProperty("UseDeadTimeThreshold", True)
    algo.setProperty("DeadTimeThreshold", 1.0003)

    algo.PyExec()
    corr_ws = algo.getProperty("OutputWorkspace").value
    corr = corr_ws.readY(0)

    # manual inspection showed 88 corrections were above the threshold
    assert len(corr[corr == 0]) == 88

    for c in corr:
        assert c <= 1.0003

def test_full_reduction(nexus_dir, template_dir):
    """
    Test dead time from the reduction workflow
    """
    template_path = os.path.join(template_dir, "template.xml")
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_198409")

    sequence_number = ws.getRun().getProperty("sequence_number").value[0]
    template_data = template.read_template(template_path, sequence_number)
    template_data.dead_time = True

    _, r1, _ = template.process_from_template_ws(ws, template_data)

    template_data.dead_time = False
    _, r2, _ = template.process_from_template_ws(ws, template_data)

    corr = r1 / r2
    for c in corr:
        assert c > 0
        assert c < 1.001
