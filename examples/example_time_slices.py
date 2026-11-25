"""
Example usage of time-sliced reduction (new_reduction_time_resolved)

Inputs are read below the example data folder and the direct beams are written to the example output folder
(see example_paths). Nothing runs on import.
"""

import os
import numpy as np

import lr_reduction.new_reduction_time_resolved as reduction

from example_paths import data_folder, output_folder, run_examples

def settings_json():
    """The settings file most examples start from: the one example_from_dat writes."""
    return data_folder() / "shared" / "lr_reduction" / "new_workflow_test_outputs" / "REFL_227158_settings.json"

def _slices_from_json(slices=2):
    run_list = [227164]
    setting_file = settings_json()
    experiment_id = "IPTS-36776"

    output = reduction.reduce_time_slices(run_list, setting_file, experiment_id, slices,
                                          savepath=output_folder(), plot_time=True)

    return output

def example_slices_from_json_2():
    return _slices_from_json(slices=2)

def example_slices_from_json_5():
    return _slices_from_json(slices=5)

def example_times_from_json(starts=[0,80,120], stops=[80,120,136]):
    run_list = [227164]
    setting_file = settings_json()
    experiment_id = "IPTS-36776"

    output = reduction.reduce_time_list(run_list, setting_file, experiment_id, starts, stops,
                                        savepath=output_folder(), plot_time=True, subname_input=None)

    return output

if __name__ == '__main__':
    run_examples([example_slices_from_json_2, example_slices_from_json_5, example_times_from_json])
