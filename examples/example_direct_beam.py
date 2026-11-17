"""
Example usage of the Direct_Beam class

Inputs are read below the example data folder and the direct beams are written to the example output folder
(see example_paths). Nothing runs on import.
"""

import numpy as np

from lr_reduction.direct_beam_maker import Direct_Beam

from example_paths import data_folder, output_folder, run_examples


def cd_attenuation():
    # The attenuation of Cd by wavelength. Left out (mu_file=None), the reduction reads the file its settings name,
    # which is this one at the facility.
    return data_folder() / "shared" / "Cd_transmission" / "Cd_mu_2025.dat"


def example_direct_beam():
    """A direct beam from four runs with different Cd attenuation"""
    db = Direct_Beam()
    db.savepath = output_folder()  # Set the directory where you want to save the output
    db.NEXUSpath = data_folder() / "IPTS-36119" / "nexus"  # Set the path to your Nexus file
    db.Chop2_cut_fn = [2.077, -16818.0]
    db.DTCcut = 1.25
    db.DTCcut_config1 = 1.8

    run_list = np.arange(225318,225321+1)
    DB_out = db.create_db(run_list=run_list, save_name='direct_beam_test.dat',
                          plot=True, mu_file=cd_attenuation(), flip_atten=True)
    return DB_out


def example_direct_beam_pixel_offset():
    """A direct beam whose position on the detector is offset: the y pixel range is set"""
    db = Direct_Beam()
    db.savepath = output_folder()  # Set the directory where you want to save the output
    db.NEXUSpath = data_folder() / "IPTS-36119" / "nexus"  # Set the path to your Nexus file
    db.Chop2_cut_fn = [2.077, -16818.0]
    db.DTCcut = 1.25
    db.DTCcut_config1 = 1.8
    db.y_ROI = [200,250]

    run_list = [235850]
    DB_out = db.create_db(run_list=run_list, save_name='direct_beam_test_pixel_offset.dat',
                          plot=True, mu_file=cd_attenuation(), flip_atten=False)
    return DB_out


if __name__ == '__main__':
    run_examples([example_direct_beam, example_direct_beam_pixel_offset])
