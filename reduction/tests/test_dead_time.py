import mantid.simpleapi as mtd_api

from lr_reduction import template
from lr_reduction.DeadTimeCorrection import SingleReadoutDeadTimeCorrection
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

def test_deadtime_threshold_ratio(nexus_dir):
    """
    Test using the threshold ratio. Here the threshold is set to 0,
    so all corrections will be 0
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_198409")

    algo = SingleReadoutDeadTimeCorrection()
    algo.PyInit()
    algo.setProperty("InputWorkspace", ws)
    algo.setProperty("OutputWorkspace", "dead_time_corr")
    algo.setProperty("ThresholdRatio", 0.0)
    algo.PyExec()
    corr_ws = algo.getProperty("OutputWorkspace").value
    corr = corr_ws.readY(0)
    for c in corr:
        assert c == 0

def test_full_reduction(nexus_dir):
    """
    Test dead time from the reduction workflow
    """
    template_path = "data/template.xml"
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


def test_full_reduction_with_threshold(nexus_dir):
    """
    Test dead time from the reduction workflow using a
    threshold ratio. Here the threshold is zero,
    so all corrections will be set to zero
    """
    template_path = "data/template.xml"
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_198409")

    sequence_number = ws.getRun().getProperty("sequence_number").value[0]
    template_data = template.read_template(template_path, sequence_number)
    template_data.dead_time = True
    template_data.dead_time_threshold_ratio = 0.5

    _, r1, _ = template.process_from_template_ws(ws, template_data)

    template_data.dead_time = False
    _, r2, _ = template.process_from_template_ws(ws, template_data)

    corr = r1 / r2
    for c in corr:
        assert c == 0.5
