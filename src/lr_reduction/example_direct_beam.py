"""
Example usage of the Direct_Beam class
"""

from lr_reduction.direct_beam_maker import Direct_Beam
import numpy as np

db = Direct_Beam()  

# Example usage
db.savepath = '/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/'  # Set the directory where you want to save the output
db.NEXUSpath = '/SNS/REF_L/IPTS-36119/nexus/'  # Set the path to your Nexus file
db.Chop2_cut_fn = [2.077, -16818.0]
db.DTCcut = 1.25
db.DTCcut_config1 = 1.8
#DB_out = db.create_db(run_list=[226473,226474,226475,226476], save_name='direct_beam_test.dat', 
#                      plot=True)

run_list = np.arange(225318,225321+1)
DB_out = db.create_db(run_list=run_list, save_name='direct_beam_test.dat', 
                      plot=True, flip_atten=True)
