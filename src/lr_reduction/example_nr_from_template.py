"""
Example script for running new reduction with config from template (e.g. for autoreduction workflow)
"""

import numpy as np
import new_reduction_from_template as new_template
from pathlib import Path

def example_template_reduction():
    
    print("\n" + "="*60)
    print("TEMPLATE REDUCTION EXAMPLE")
    print("="*60)
    
    datapath = Path('/Users/r2i/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Documents/4B/Reduction/GitHub/Workflowtesting/')
    DBpath = Path('/Users/r2i/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Documents/4B/Reduction/GitHub/Cd_DB_processing_013026/DBs/')

    # This needs to include anything you want to set that isn't in the template file or you want to deviate from the template file.
    override_params = {
        'DBname': ['A1_air_div1_Cd_DB.dat'],
        'Spath': datapath,
        'NEXUSpathRB': datapath,
        'DBpath': DBpath,
        'useCalcTheta': True,
        'AutoScale': False
    }

    results = new_template.reduce_from_template(211029, datapath / "test_template.xml", "IPTS-36119", datapath=datapath, override_params=override_params, plot=True)
   
    #print(f"\nReduced {len(config.RBnum)} runs")
    print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")
    
    return results

def example_template_reduction_3ang():
    
    print("\n" + "="*60)
    print("TEMPLATE REDUCTION EXAMPLE")
    print("="*60)
    
    datapath = Path('/Users/r2i/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Documents/4B/Reduction/GitHub/Workflowtesting/')
    DBpath = Path('/Users/r2i/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Documents/4B/Reduction/GitHub/Cd_DB_processing_013026/DBs/')

    DB_list = ['A1_air_div1_Cd_DB.dat','A1_air_div1_Cd_DB.dat','A1_air_div1_Cd_DB.dat']

    runs = [211029, 211030, 211031]
    for idx, run in enumerate(runs):
        # This needs to include anything you want to set that isn't in the template file or you want to deviate from the template file.
        override_params = {
            'DBname': [DB_list[idx]],
            'Spath': datapath,
            'NEXUSpathRB': datapath,
            'DBpath': DBpath,
            'useCalcTheta': True,
            'AutoScale': True
        }

        results = new_template.reduce_from_template(run, datapath / "test_template_3ang.xml", "IPTS-36119", datapath=datapath, override_params=override_params, plot=True)
   
        #print(f"\nReduced {len(config.RBnum)} runs")
        print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")
    

    return results

# THIS ONE DOESN'T WORK YET AS IT DOESN'T SAVE THE TEMPLATE CORRECTLY
def example_template_reduction_new():
    
    print("\n" + "="*60)
    print("TEMPLATE REDUCTION EXAMPLE")
    print("="*60)
    
    datapath = Path('/Users/r2i/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Documents/4B/Reduction/GitHub/Workflowtesting/')
    DBpath = Path('/Users/r2i/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/Documents/4B/Reduction/GitHub/Cd_DB_processing_013026/DBs/')

    #DB_list = ['A1_air_div1_Cd_DB.dat','A1_air_div1_Cd_DB.dat','A1_air_div1_Cd_DB.dat']

    runs = [211029, 211030, 211031]
    for idx, run in enumerate(runs):
    #    # This needs to include anything you want to set that isn't in the template file or you want to deviate from the template file.
    #    override_params = {
    #        'DBname': [DB_list[idx]],
    #        'Spath': datapath,
    #        'NEXUSpathRB': datapath,
    #        'DBpath': DBpath,
    #        'useCalcTheta': True,
    #        'AutoScale': True
    #    }

        results = new_template.reduce_from_template(run, datapath / "REFL_211029_template_new.xml", "IPTS-36119", datapath=datapath, plot=True)
   
        #print(f"\nReduced {len(config.RBnum)} runs")
        print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")
    

    return results

if __name__ == '__main__':
    # Run examples
    #example_template_reduction()
    example_template_reduction_3ang()
    example_template_reduction_new()