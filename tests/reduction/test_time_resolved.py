import unittest
import pytest
import os
import numpy as np
from reduction.lr_reduction import time_resolved
from mantid import config


class TimeResolvedTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.getcwd().endswith("LiquidsReflectometer"):
            os.chdir("tests")

        cwd = os.getcwd()
        dirs = [26010, 26776, 28662, 29196, 31279]
        for dir_num in dirs:
            config.appendDataSearchDir(str(os.path.join(cwd, f"data/liquidsreflectometer-data/SNS/REF_L/IPTS-{dir_num}/nexus")))

    @pytest.mark.datarepo()
    def test_reduce_workflow(self):
        """
        Test the time-resolved reduction that uses a measured reference.
        It is generally used at 30 Hz but it also works at 60 Hz.
        """
        template_path = "data/template.xml"
        output_dir = "/tmp"
        reduced_path = "data/reference_rq_avg_overlap.txt"
        ref_data = np.loadtxt(reduced_path).T

        reduced = time_resolved.reduce_30Hz_slices(
            198413,
            198413,
            ref_data_60Hz=reduced_path,
            template_30Hz=template_path,
            time_interval=300,
            output_dir=output_dir,
            scan_index=5,
            create_plot=False,
        )

        q_long = len(ref_data[0])
        q_short = len(reduced[0][0])
        n_match = 0
        n_pts = 0
        for i in range(q_long):
            if ref_data[0][i] > 0.03 and ref_data[0][i] < 0.047:
                n_pts += 1
                for k in range(q_short):
                    if np.fabs(reduced[0][0][k] - ref_data[0][i]) < 0.0001:
                        assert np.fabs(reduced[0][1][k] - ref_data[1][i]) < 1e-10
                        n_match += 1
        assert n_pts == n_match
