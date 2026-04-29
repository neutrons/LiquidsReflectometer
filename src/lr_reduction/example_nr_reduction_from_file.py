# -*- coding: utf-8 -*-
"""
Example usage of the unified NR reduction class

Demonstrates how to configure and run reductions using both constantQ, constantTOF and MeanTheta methods
"""

import numpy as np
from pathlib import Path
import lr_reduction.new_reduction_from_file as reduction

def example_from_dat():
    run_list = [227158,227159,227160]
    run_list = [227158,227159,227160, 227161, 227162]
    setting_file = "/SNS/REF_L/IPTS-36776/shared/EBW_reduced/NR_meanTheta_8test_combined_8col.dat"
    experiment_id = "IPTS-36776"
    Spath = Path('/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/')
    override_params = {'Spath': Spath,}
    output = reduction.reduce_from_file(run_list, setting_file, experiment_id, override_params=override_params, plot=True, save_json=True)

    return output


def example_from_json():
    run_list = [227164,227165,227166]
    setting_file = "/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/REFL_227158_settings.json"
    experiment_id = "IPTS-36776"
    Spath = Path('/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/')
    override_params = {'Spath': Spath,}
    output = reduction.reduce_from_file(run_list, setting_file, experiment_id, override_params=override_params, plot=True, save_json=False)

    return output


if __name__ == '__main__':
    # Run examples
    example_from_dat()
    example_from_json()



