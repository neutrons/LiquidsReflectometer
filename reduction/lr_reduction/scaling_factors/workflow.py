"""
    Scaling factors calculation workflow
"""
import os

from lr_reduction.scaling_factors import LRDirectBeamSort
from lr_reduction.utils import mantid_algorithm_exec


def process_scaling_factors(ws, output_dir, tof_step=200., order_by_runs=True,
                            incident_medium='air', slit_tolerance=0.06, wait=False,
                            postfix='_auto', use_deadtime=True, paralyzable=True,
                            deadtime=4.2, deadtime_tof_step=200):
    """
        Compute scaling factors given a DB run, assumed to be the last
        one of a set.
        :param workspace ws: Mantid workspace for one of the direct beams to use.
        :param output_dir: path the the output directory
        :param tof_step: TOF binning for the scaling factor calculation
        :param order_by_runs: if True, the runs will be ordered by run number instead of meta data
        :param incident_medium: name of the incident medium
        :param slit_tolerance: tolerance to use when matching slits between runs
        :param wait: if True, scaling factors will only be processed if the workspace
                     given corresponds to the last run of the complete set
        :param postfix: string to add at the end of the output file
        :param use_deadtime: if True, a dead time correction will be applied
        :param paralyzable: if True, a paralyzable dead time correction will be applied
        :param deadtime: value of the dead time
        :param deadtime_tof_step: TOF binning to use when computing the dead time
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

    mantid_algorithm_exec(LRDirectBeamSort.LRDirectBeamSort,
                          RunList=list(range(first_run_of_set, first_run_of_set + sequence_total)),
                          UseLowResCut=True,
                          ComputeScalingFactors=True,
                          TOFSteps=tof_step,
                          IncidentMedium=incident_medium,
                          SlitTolerance=slit_tolerance,
                          OrderDirectBeamsByRunNumber=order_by_runs,
                          UseDeadTimeCorrection=use_deadtime,
                          ParalyzableDeadTime=paralyzable,
                          DeadTime=deadtime,
                          DeadTimeTOFStep=deadtime_tof_step,
                          ScalingFactorFile=output_cfg)

    return True
