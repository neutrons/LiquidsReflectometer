"""
    Scaling factors calculation workflow
"""
import os

from . import LRDirectBeamSort


def process_scaling_factors(ws, output_dir, tof_step=200., order_by_runs=True,
                            incident_medium='air', slit_tolerance=0.06, wait=False,
                            postfix='_auto', use_deadtime=True, paralyzable=True,
                            deadtime=4.2, deadtime_tof_step=200):
    """
        Compute scaling factors given a DB run, assumed to be the last
        one of a set.
        :param workspace ws: Mantid workspace for one of the direct beams to use.
    """
    # Read in the sequence information
    meta_data_run = ws.getRun()
    run_number = ws.getRunNumber()
    sequence_number = meta_data_run.getProperty("sequence_number").value[0]
    first_run_of_set = meta_data_run.getProperty("sequence_id").value[0]
    sequence_total = meta_data_run.getProperty("sequence_total").value[0]

    print(f"Run {run_number} - Sequence {first_run_of_set} [{first_run_of_set}/{sequence_total}]")

    # For an automated reduction situation, we have to wait until all the direct
    # beams are acquired.
    if wait and sequence_number < sequence_total:
        print(f"Waiting for at least {sequence_total} runs to compute scaling factors")
        return False

    # The medium for these direct beam runs may not be what was set in the template,
    # so either use the medium in the data file or a default name

    if meta_data_run.hasProperty("incident_medium"):
        incident_medium = meta_data_run.getProperty("incident_medium").value[0]

    file_id = incident_medium.replace("medium", "")

    output_cfg = os.path.join(output_dir, "sf_%s_%s%s.cfg" % (first_run_of_set, file_id, postfix))

    algo = LRDirectBeamSort.LRDirectBeamSort()
    algo.PyInit()
    algo.setProperty("RunList", list(range(first_run_of_set, first_run_of_set + sequence_total)))
    algo.setProperty("UseLowResCut", True)
    algo.setProperty("ComputeScalingFactors", True)
    algo.setProperty("TOFSteps", tof_step)
    algo.setProperty("IncidentMedium", incident_medium)
    algo.setProperty("SlitTolerance", slit_tolerance)
    algo.setProperty("OrderDirectBeamsByRunNumber", order_by_runs)
    algo.setProperty("UseDeadTimeCorrection", use_deadtime)
    algo.setProperty("ParalyzableDeadTime", paralyzable)
    algo.setProperty("DeadTime", deadtime)
    algo.setProperty("DeadTimeTOFStep", deadtime_tof_step)
    algo.setProperty("ScalingFactorFile", output_cfg)
    algo.PyExec()

    return True
