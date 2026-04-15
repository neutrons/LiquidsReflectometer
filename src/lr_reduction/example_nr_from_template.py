"""
Example script for running new reduction with config from template (e.g. for autoreduction workflow)
"""

from pathlib import Path

import new_reduction_from_template as new_template


def example_template_reduction():
    # Example for just a single angle.

    print("\n" + "="*60)
    print("TEMPLATE REDUCTION EXAMPLE")
    print("="*60)

    # These shouldn't be needed in the longer term.
    datapath = Path('/SNS/REF_L/IPTS-30101/nexus')
    DBpath = Path('/SNS/REF_L/shared/Cd_DB_processing/DBs/')
    template_path = Path('/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/')
    Spath = Path('/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/')

    # This needs to include anything you want to set that isn't in the template file or you want to deviate from the template file.
    override_params = {
        'DBname': ['A1_air_div1_Cd_DB.dat'],
        'Spath': Spath,
        'NEXUSpathRB': datapath,
        'DBpath': DBpath,
        'useCalcTheta': True,
        'AutoScale': False
    }

    results = new_template.reduce_from_template(211029, template_path / "test_template.xml", "IPTS-36119", datapath=datapath, template_path=template_path, override_params=override_params, plot=True)

    #print(f"\nReduced {len(config.RBnum)} runs")
    print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")

    return results

def example_template_reduction_3ang():
    # Example with 3 angles in a set

    print("\n" + "="*60)
    print("TEMPLATE REDUCTION EXAMPLE")
    print("="*60)

    datapath = Path('/SNS/REF_L/IPTS-30101/nexus')
    DBpath = Path('/SNS/REF_L/shared/Cd_DB_processing/DBs/')
    template_path = Path('/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/')
    Spath = Path('/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/')

    DB_list = ['A1_air_div1_Cd_DB.dat','A2_air_div1_Cd_DB.dat','A3_air_div1_Cd_DB.dat']

    runs = [211029, 211030, 211031]
    use_BS = [0,1,1]
    for idx, run in enumerate(runs):
        # This needs to include anything you want to set that isn't in the template file or you want to deviate from the template file.
        print(use_BS[idx])
        override_params = {
            'DBname': [DB_list[idx]],
            'Spath': Spath,
            'NEXUSpathRB': datapath,
            'DBpath': DBpath,
            'useCalcTheta': True,
            'AutoScale': True,
            'useBS': [use_BS[idx]]
        }

        results = new_template.reduce_from_template(run, template_path / "test_template_3ang.xml", "IPTS-36119", datapath=datapath,template_path=template_path,override_params=override_params, plot=True)

        #print(f"\nReduced {len(config.RBnum)} runs")
        print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")


    return results

def example_template_reduction_new():
    # Example to redo the prior reduction but reading from the new template.

    print("\n" + "="*60)
    print("TEMPLATE REDUCTION EXAMPLE")
    print("="*60)

    datapath = Path('/SNS/REF_L/IPTS-30101/nexus')
    Spath = Path('/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/')
    DBpath = Path('/SNS/REF_L/shared/Cd_DB_processing/DBs/')
    template_path = Path('/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/')

    #DB_list = ['A1_air_div1_Cd_DB.dat','A1_air_div1_Cd_DB.dat','A1_air_div1_Cd_DB.dat']

    runs = [211029, 211030, 211031]
    for idx, run in enumerate(runs):
    #    # This needs to include anything you want to set that isn't in the template file or you want to deviate from the template file.
        override_params = {
            'Spath': Spath,
            'NEXUSpathRB': datapath,
            'DBpath': DBpath
        }

        results = new_template.reduce_from_template(run, template_path / "REFL_211029_template_new.xml", "IPTS-36119",
                                                    datapath=datapath, template_path=template_path, override_params=override_params, plot=True)

        #print(f"\nReduced {len(config.RBnum)} runs")
        print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")


    return results

if __name__ == '__main__':
    # Run examples
    #example_template_reduction()
    example_template_reduction_3ang()
    #example_template_reduction_new()
