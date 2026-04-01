# -*- coding: utf-8 -*-
"""
Example usage of the unified NR reduction class

Demonstrates how to configure and run reductions using both constantQ, constantTOF and MeanTheta methods
"""

from nr_reduction_calc import NR_Reduction
from nr_reduction_config import NRReductionConfig
import numpy as np


def example_constant_q_reduction():
    """Example: Constant Q-line reduction"""
    print("\n" + "="*60)
    print("CONSTANT Q-LINE REDUCTION EXAMPLE")
    print("="*60)
    
    # Create configuration for constantQ method
    config = NRReductionConfig()
    config.method_per_run = ['constantQ']

    # Path needs to be setup for these tests but should be able to turn off later.    
    config.Spath = '/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/'
    #config.NEXUSpathRB = # Append path if required.
    config.DBpath = '/SNS/REF_L/shared/Cd_DB_processing/DBs/'

    config.experiment_id = "IPTS-30101"

    # Set data
    config.DBname = ['A1_air_div1_Cd_DB.dat', 'A2_air_div1_Cd_DB.dat', 'A3_air_div1_Cd_DB.dat']
    config.RBnum=[211029,211030,211031]
    config.RB_Ymin=[157,157,155]
    config.RB_Ymax=[165,166,169]

    config.ThetaShift = [0, 0, 0]
    config.useBS = [1, 1, 1]
    config.ScaleFactor = [1, 1, 1]

    RB_Ymin = np.array(config.RB_Ymin)
    RB_Ymax = np.array(config.RB_Ymax)
    config.BkgROI = np.column_stack((
        RB_Ymin - 8,
        RB_Ymin - 3,
        RB_Ymax + 3,
        RB_Ymax + 8
    ))

    # Output name
    config.Sname = "NR_constantQ"
    
    # Processing options
    config.Normalize = False
    config.AutoScale = False
    config.useCalcTheta = True
    config.plotON = True  # Set to False for batch processing
    
    # Method-specific parameters
    config.peak_pad = 5
    config.Qline_threshold = 0.66
    
    # Run reduction
    reducer = NR_Reduction(config)
    results = reducer.reduce()
    
    print(f"\nReduced {len(config.RBnum)} runs")
    print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")
    
    return results


def example_mean_theta_reduction():
    """Example: Mean theta reduction"""
    print("\n" + "="*60)
    print("MEAN THETA REDUCTION EXAMPLE")
    print("="*60)
    
    # Create configuration for MeanTheta method
    config = NRReductionConfig()
    config.method_per_run = ['meanTheta']
    
    # Path needs to be setup for these tests but should be able to turn off later.    
    config.Spath = '/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/'
    #config.NEXUSpathRB = # Append path if required.
    config.DBpath = '/SNS/REF_L/shared/Cd_DB_processing/DBs/'

    config.experiment_id = "IPTS-30101"
    # Set data
    config.DBname = ['A1_air_div10_Cd_DB.dat', 'A2_air_div10_Cd_DB.dat', 'A3_air_div10_Cd_DB.dat']
    config.RBnum = np.array([210975, 210976, 210977])
    config.RB_Ymin = [157, 155, 149]
    config.RB_Ymax = [165, 167, 174]
    config.ThetaShift = [0, 0, 0]
    config.useBS = [1, 1, 1]
    config.ScaleFactor = [1, 1, 1]

    RB_Ymin = np.array(config.RB_Ymin)
    RB_Ymax = np.array(config.RB_Ymax)
    config.BkgROI = np.column_stack((
        RB_Ymin - 8,
        RB_Ymin - 3,
        RB_Ymax + 3,
        RB_Ymax + 8
    ))

    # Output name
    config.Sname = "NR_meanTheta"
    
    # Processing options
    config.Normalize = False
    config.AutoScale = True 
    config.useCalcTheta = False
    config.plotON = True  # Set to False for batch processing
    
    # Method-specific parameters
    config.peak_pad = 1
    config.Qline_threshold = 1.0
    config.DetSigma = 0.8
    config.DetResFn = 'rectangular'
    
    # Run reduction
    reducer = NR_Reduction(config)
    results = reducer.reduce()
    
    print(f"\nReduced {len(config.RBnum)} runs with repeat averaging")
    print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")
    
    return results


def example_constant_tof_reduction():
    """Example: Constant TOF reduction"""
    print("\n" + "="*60)
    print("CONSTANT TOF REDUCTION EXAMPLE")
    print("="*60)
    
    # Create configuration for constantTOF method
    config = NRReductionConfig()
    config.method_per_run = ['constantTOF']
    
    # Path needs to be setup for these tests but should be able to turn off later.    
    config.Spath = '/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/'
    #config.NEXUSpathRB = # Append path if required.
    config.DBpath = '/SNS/REF_L/shared/Cd_DB_processing/DBs/'

    config.experiment_id = "IPTS-30101"
    # Set data
    config.DBname = ['A1_air_div1_Cd_DB.dat', 'A2_air_div1_Cd_DB.dat', 'A3_air_div1_Cd_DB.dat']
    config.RBnum=[211029,211030,211031]
    config.RB_Ymin=[157,157,155]
    config.RB_Ymax=[165,166,169]

    config.ThetaShift = [0, 0, 0]
    config.useBS = [1, 1, 1]
    config.ScaleFactor = [1, 1, 1]
    RB_Ymin = np.array(config.RB_Ymin)
    RB_Ymax = np.array(config.RB_Ymax)
    config.BkgROI = np.column_stack((
        RB_Ymin - 8,
        RB_Ymin - 3,
        RB_Ymax + 3,
        RB_Ymax + 8
    ))
    
    # Output name
    config.Sname = "NR_constantTOF"
    
    # Processing options
    config.Normalize = False
    config.AutoScale = True
    config.useCalcTheta = True
    config.plotON = True
    
    # Method-specific parameters
    config.peak_pad = 1
    config.plotQ4 = False  # ConstantTOF typically doesn't plot Q4
    
    # Run reduction
    reducer = NR_Reduction(config)
    results = reducer.reduce()
    
    print(f"\nReduced {len(config.RBnum)} runs (TOF-binned)")
    print(f"Q range: {results['Q'].min():.4f} - {results['Q'].max():.4f} Å⁻¹")
    
    return results


if __name__ == '__main__':
    # Run examples
    example_constant_q_reduction()
    example_mean_theta_reduction()
    example_constant_tof_reduction()


    
