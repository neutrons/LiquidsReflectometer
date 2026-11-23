"""
Example usage of time-sliced reduction (new_reduction_time_resolved)

Example usage only, not a test (plan time-slicing-reconcile §10 A5): the settings file, run and output folder are
one experiment's on the facility file system, and running this script writes the slices there. Point them at your
own data first. Nothing runs on import.
"""

from pathlib import Path

import lr_reduction.new_reduction_time_resolved as reduction


def slices_from_json(slices=2):
    run_list = 227164
    setting_file = "/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/REFL_227158_settings.json"
    experiment_id = "IPTS-36776"
    Spath = Path('/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/')

    output = reduction.reduce_time_slices(run_list, setting_file, experiment_id, slices,
                                          savepath=Spath, plot_time=True)

    return output

def times_from_json(starts=[0,80,120], stops=[80,120,136]):
    run_list = 227164
    setting_file = "/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/REFL_227158_settings.json"
    experiment_id = "IPTS-36776"
    Spath = Path('/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/')

    output = reduction.reduce_time_list(run_list, setting_file, experiment_id, starts, stops,
                                        savepath=Spath, plot_time=True, subname_input=None)

    return output

if __name__ == '__main__':
    # Run examples
    slices_from_json(slices=2)
    slices_from_json(slices=5)
    times_from_json()
