"""
Example script for running new reduction with config from template (e.g. for autoreduction workflow)

Inputs are read below the example data folder and results are written to the example output folder
(see example_paths). The reduction also writes the template it updates, REFL_<sequence>_template_new.xml, to
template_path: the examples give the output folder for it. Nothing runs on import.
"""

import lr_reduction.new_reduction_from_template as new_template

from example_paths import data_folder, output_folder, run_examples

EXPERIMENT = "IPTS-36119"


def raw_data():
    # The runs of these examples were measured in another experiment than the one named to the reduction.
    return data_folder() / "IPTS-30101" / "nexus"


def direct_beams():
    # The direct beams measured with these runs (runs 210927 to 210941, beam at pixel 161). Files of the same names
    # are kept in other folders for later detector positions; with those the reduction looks for the reflected beam
    # in the wrong place and stops.
    return data_folder() / "shared" / "EBW" / "Cd_DB_processing" / "DBs"


def templates():
    return data_folder() / "shared" / "lr_reduction" / "new_workflow_test_outputs"


def example_template_reduction():
    # Example for just a single angle.

    print("\n" + "="*60)
    print("TEMPLATE REDUCTION EXAMPLE")
    print("="*60)

    # This needs to include anything you want to set that isn't in the template file or you want to deviate from the template file.
    override_params = {
        'DBname': ['A1_air_div1_Cd_DB.dat'],
        'Spath': output_folder(),
        'NEXUSpathRB': raw_data(),
        'DBpath': direct_beams(),
        'useCalcTheta': True,
        'AutoScale': False
    }

    results, group_num = new_template.reduce_from_template(211029, templates() / "test_template.xml", EXPERIMENT, datapath=raw_data(), template_path=output_folder(), override_params=override_params, plot=True)

    print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")

    return results

def example_template_reduction_3ang():
    # Example with 3 angles in a set

    print("\n" + "="*60)
    print("TEMPLATE REDUCTION EXAMPLE")
    print("="*60)

    DB_list = ['A1_air_div1_Cd_DB.dat','A2_air_div1_Cd_DB.dat','A3_air_div1_Cd_DB.dat']

    runs = [211029, 211030, 211031]
    use_BS = [0,1,1]
    for idx, run in enumerate(runs):
        # This needs to include anything you want to set that isn't in the template file or you want to deviate from the template file.
        print(use_BS[idx])
        override_params = {
            'DBname': [DB_list[idx]],
            'Spath': output_folder(),
            'NEXUSpathRB': raw_data(),
            'DBpath': direct_beams(),
            'useCalcTheta': True,
            'AutoScale': True,
            'useBS': [use_BS[idx]]
        }

        results, group_num = new_template.reduce_from_template(run, templates() / "test_template_3ang.xml", EXPERIMENT, datapath=raw_data(), template_path=output_folder(), override_params=override_params, plot=True)

        print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")


    return results

def example_template_reduction_new():
    # Example to redo the prior reduction but reading from the new template.

    print("\n" + "="*60)
    print("TEMPLATE REDUCTION EXAMPLE")
    print("="*60)

    runs = [211029, 211030, 211031]
    for run in runs:
        # This needs to include anything you want to set that isn't in the template file or you want to deviate from the template file.
        override_params = {
            'Spath': output_folder(),
            'NEXUSpathRB': raw_data(),
            'DBpath': direct_beams()
        }

        results, group_num = new_template.reduce_from_template(run, templates() / "REFL_211029_template_new.xml", EXPERIMENT,
                                                    datapath=raw_data(), template_path=output_folder(), override_params=override_params, plot=True)

        print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")


    return results

if __name__ == '__main__':
    run_examples([example_template_reduction, example_template_reduction_3ang, example_template_reduction_new])
