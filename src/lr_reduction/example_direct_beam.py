"""
Example usage of the Direct_Beam class

Example usage only, not a test (plan time-slicing-reconcile §10 A5): the paths are experiments' folders on the
facility file system, and running this script writes the direct beams there. Point them at your own data first.
Nothing runs on import.
"""

import numpy as np

from lr_reduction.direct_beam_maker import Direct_Beam


def main():
    db = Direct_Beam()

    # Example usage
    db.savepath = '/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/'  # Set the directory where you want to save the output
    db.NEXUSpath = '/SNS/REF_L/IPTS-36119/nexus/'  # Set the path to your Nexus file
    db.Chop2_cut_fn = [2.077, -16818.0]
    db.DTCcut = 1.25
    db.DTCcut_config1 = 1.8
    DB_out = db.create_db(run_list=[226473,226474,226475,226476], save_name='direct_beam_test.dat',
                          plot=True)

    run_list = np.arange(225318,225321+1)
    DB_out = db.create_db(run_list=run_list, save_name='direct_beam_test.dat',
                          plot=True, flip_atten=True)


    # Example usage
    db = Direct_Beam()
    db.savepath = '/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/'  # Set the directory where you want to save the output
    db.NEXUSpath = '/SNS/REF_L/IPTS-36119/nexus/'  # Set the path to your Nexus file

    db.DTCcut_config1 = 1.8

    run_list = [231187]
    DB_out = db.create_db(run_list=run_list, save_name='DB_231187_Cd_test.dat',
                          plot=True, cd_list=[[1, 0, 0, 0]])


    run_list = [231189]
    DB_out = db.create_db(run_list=run_list, save_name='DB_231189_Cd_test.dat',
                          plot=True, cd_list=[[0, 1, 0, 0]])

    run_list = [231191]
    DB_out = db.create_db(run_list=run_list, save_name='DB_231191_Cd_test.dat',
                          plot=True, cd_list=[[1, 0, 0, 0]])
    return DB_out


if __name__ == '__main__':
    main()
