import pytest
import unittest
import os

from mantid import config
from scripts.autoreduce.sf_calculator import ScalingFactor


@pytest.mark.scripts()
class TestScalingFactor(unittest.TestCase):
    def setUpClass(cls):
        if os.getcwd().endswith("LiquidsReflectometer"):
            os.chdir("tests")

        cwd = os.getcwd()
        config.appendDataSearchDir(str(os.path.join(cwd, "data/liquidsreflectometer-data/nexus")))
        print(config.getDataSearchDirs())

    def test_scaling_factor(self):
        sf = ScalingFactor(run_list=range(184975, 184990), sf_file="results/sf_184975_air.cfg", medium="air")

        sf.execute()
