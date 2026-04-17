# -*- coding: utf-8 -*-
"""
Example usage of the unified NR reduction class

Demonstrates how to configure and run reductions using both constantQ, constantTOF and MeanTheta methods
"""

import numpy as np
from nr_reduction_calc import NR_Reduction
from nr_reduction_config import NRReductionConfig



def example_mean_theta_reduction():
    """Example: Mean theta reduction"""
    print("\n" + "="*60)
    print("MEAN THETA REDUCTION EXAMPLE")
    print("="*60)

    # Create configuration for MeanTheta method
    config = NRReductionConfig()
    config.method_per_run = ['meanTheta']

    # Path needs to be setup for these tests but should be able to turn off later.
    config.Spath = '/SNS/REF_L/IPTS-36776/shared/EBW_reduced/'
    #config.NEXUSpathRB = # Append path if required.
    config.DBpath = '/SNS/REF_L/IPTS-36776/shared/transmission/'

    config.experiment_id = "IPTS-36776"
    # Set data
    config.DBname = ['DB_A1_Cd.txt', 'DB_A2_Cd.txt', 'DB_A3_Cd.txt']
    config.RBnum=[227147,227148,227149]
    config.RB_Ymin=[143,142,138]
    config.RB_Ymax=[157,157,159]
    
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
    config.Spath = '/SNS/REF_L/IPTS-36776/shared/EBW_reduced/'
    #config.NEXUSpathRB = # Append path if required.
    config.DBpath = '/SNS/REF_L/IPTS-36776/shared/transmission/'

    config.experiment_id = "IPTS-36776"
    # Set data
    config.DBname = ['DB_A1_Cd.txt', 'DB_A2_Cd.txt', 'DB_A3_Cd.txt']
    config.RBnum=[227147,227148,227149]
    config.RB_Ymin=[143,142,138]
    config.RB_Ymax=[157,157,159]

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
    config.useCalcTheta = False
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
    #example_mean_theta_reduction()
    example_constant_tof_reduction()



