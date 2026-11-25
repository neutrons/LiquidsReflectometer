# -*- coding: utf-8 -*-
"""
Example usage of the reduction from a settings file: a .json of settings, or a reduced .dat whose header holds them

Inputs are read below the example data folder and results are written to the example output folder
(see example_paths). The settings files name the folders of the experiment they were written for, so each example
gives the folders to use (the raw data, the direct beams, the output) and they replace the file's.
Nothing runs on import.
"""

import lr_reduction.new_reduction_from_file as reduction

from example_paths import data_folder, output_folder, run_examples

EXPERIMENT = "IPTS-36776"


def raw_data(experiment_id=EXPERIMENT):
    return data_folder() / experiment_id / "nexus"


def direct_beams(experiment_id=EXPERIMENT):
    return data_folder() / experiment_id / "shared" / "transmission"


def settings_json():
    """The settings file most examples start from: the one example_from_dat writes."""
    return data_folder() / "shared" / "lr_reduction" / "new_workflow_test_outputs" / "REFL_227158_settings.json"


def example_from_dat():
    run_list = [227158,227159,227160, 227161, 227162]
    setting_file = data_folder() / EXPERIMENT / "shared" / "EBW_reduced" / "NR_meanTheta_8test_combined_8col.dat"
    override_params = {'Spath': output_folder(), 'DBpath': direct_beams()}
    output = reduction.reduce_from_file(run_list, setting_file, EXPERIMENT, datapath=raw_data(), override_params=override_params, plot=True, save_json=True)

    return output


def example_from_json():
    run_list = [227164,227165,227166]
    override_params = {'Spath': output_folder(), 'DBpath': direct_beams()}
    output = reduction.reduce_from_file(run_list, settings_json(), EXPERIMENT, datapath=raw_data(), override_params=override_params, plot=True, save_json=False)

    return output

def example_wrong_order():
    run_list = [227165,227164,227166]
    Sname = "REFL_227164_wrong_order_test"
    override_params = {'Spath': output_folder(), 'DBpath': direct_beams(), 'Sname': Sname}
    output = reduction.reduce_from_file(run_list, settings_json(), EXPERIMENT, datapath=raw_data(), override_params=override_params, plot=True, save_json=False)

    return output

def example_one_missing():
    run_list = [227165,227166]
    subname = "missing_test"
    override_params = {'Spath': output_folder(), 'DBpath': direct_beams(), "subname": subname}
    output = reduction.reduce_from_file(run_list, settings_json(), EXPERIMENT, datapath=raw_data(), override_params=override_params, plot=True, save_json=False)

    return output

def example_with_search():
    # check_for_prior looks in the output folder for the other runs of the set, reduced earlier: reduce them first.
    example_one_missing()
    run_list = [227164]
    subname = "missing_test"
    override_params = {'Spath': output_folder(), 'DBpath': direct_beams(), "subname": subname}
    output = reduction.reduce_from_file(run_list, settings_json(), EXPERIMENT, datapath=raw_data(), override_params=override_params, plot=True, save_json=False, check_for_prior=True)

    return output

def example_with_savepdf():
    run_list = [227164,227165,227166]
    Sname = "REFL_227164_save_test"
    override_params = {'Spath': output_folder(), 'DBpath': direct_beams(), "Sname": Sname}
    output = reduction.reduce_from_file(run_list, settings_json(), EXPERIMENT, datapath=raw_data(), override_params=override_params, plot=True, save_json=False, save_pdf_summary=True)

    return output

def example_with_savepdf_noshow():
    run_list = [227164,227165,227166]
    subname = "save_test_v2"
    override_params = {'Spath': output_folder(), 'DBpath': direct_beams(), "subname": subname}
    output = reduction.reduce_from_file(run_list, settings_json(), EXPERIMENT, datapath=raw_data(), override_params=override_params, plot=False, save_json=False, save_pdf_summary=True)

    return output

def example_with_search_savepdf():
    # check_for_prior looks in the output folder for the other runs of the set, reduced earlier: reduce them first.
    example_one_missing()
    run_list = [227164]
    subname = "missing_test"
    override_params = {'Spath': output_folder(), 'DBpath': direct_beams(), "subname": subname}
    output = reduction.reduce_from_file(run_list, settings_json(), EXPERIMENT, datapath=raw_data(), override_params=override_params, plot=True, save_json=False, check_for_prior=True, save_pdf_summary=True)

    return output

def example_with_8col():
    run_list = [227164,227165,227166]
    subname = "eight_test"
    override_params = {'Spath': output_folder(), 'DBpath': direct_beams(), "subname": subname, "save8col": True}
    output = reduction.reduce_from_file(run_list, settings_json(), EXPERIMENT, datapath=raw_data(), override_params=override_params, plot=True, save_json=False, check_for_prior=False, save_pdf_summary=True)

    return output

def example_with_8col_search():
    # check_for_prior looks in the output folder for the other runs of the set, reduced earlier: reduce them first.
    example_with_8col()
    run_list = [227165]
    subname = "eight_test"
    override_params = {'Spath': output_folder(), 'DBpath': direct_beams(), "subname": subname, "save8col": True}
    output = reduction.reduce_from_file(run_list, settings_json(), EXPERIMENT, datapath=raw_data(), override_params=override_params, plot=True, save_json=False, check_for_prior=True, save_pdf_summary=True)

    return output

def example_for_autoreduction():
    experiment_id = "IPTS-35331"
    run_list = [230395]
    setting_file = data_folder() / experiment_id / "shared" / "autoreduce" / "reduce_settings_down.json"
    subname = "autoreduction"
    override_params = {'Spath': output_folder(), 'DBpath': direct_beams(experiment_id), "subname": subname}
    output = reduction.reduce_from_file(run_list, setting_file, experiment_id, datapath=raw_data(experiment_id), override_params=override_params, plot=True, save_json=False, check_for_prior=True, save_pdf_summary=False)

    return output

def example_incidenttheta():
    experiment_id = "IPTS-36119"
    run_list = [231035]
    setting_file = data_folder() / experiment_id / "shared" / "reduced" / "Aug2026" / "REFL_231035_settings.json"
    DBpath = data_folder() / experiment_id / "shared" / "transmission" / "Aug2026"
    override_params = {'Spath': output_folder(), 'DBpath': DBpath}
    output = reduction.reduce_from_file(run_list, setting_file, experiment_id, datapath=raw_data(experiment_id), override_params=override_params, plot=True, save_json=False)

    return output

if __name__ == '__main__':
    run_examples([
        example_from_dat,
        example_from_json,
        example_wrong_order,
        example_one_missing,
        example_with_search,
        example_with_savepdf,
        example_with_savepdf_noshow,
        example_with_search_savepdf,
        example_with_8col,
        example_with_8col_search,
        example_for_autoreduction,
        example_incidenttheta,
    ])
