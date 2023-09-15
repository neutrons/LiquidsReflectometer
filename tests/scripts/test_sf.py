import pytest
import unittest

from scripts.autoreduce.sf_calculator import ScalingFactor


class TestScalingFactor(unittest.TestCase):
    @pytest.mark.scripts()
    def test_scaling_factor(self):
        sf = ScalingFactor(run_list=range(184975, 184990), sf_file="/tmp/sf_184975_air.cfg", medium="air")

        sf.execute()
