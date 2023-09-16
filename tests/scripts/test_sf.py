import pytest
import unittest
import os

from mantid import config
from scripts.autoreduce.sf_calculator import ScalingFactor


@pytest.mark.scripts()
class TestScalingFactor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.getcwd().endswith("LiquidsReflectometer"):
            os.chdir("tests")

        cwd = os.getcwd()
        dirs = [26010, 26776, 28662, 29196, 31279]
        for dir_num in dirs:
            config.appendDataSearchDir(str(os.path.join(cwd, f"data/liquidsreflectometer-data/SNS/REF_L/IPTS-{dir_num}/nexus")))

    @classmethod
    def test_scaling_factor(self):
        sf = ScalingFactor(run_list=range(184975, 184990), sf_file="/tmp/sf_184975_air.cfg", medium="air")

        sf.execute()
