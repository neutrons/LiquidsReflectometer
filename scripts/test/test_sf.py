import sys
import os
sys.path.append(os.path.expanduser('~/git/LiquidsReflectometer/reduction'))

from lr_reduction.sf_calculator import ScalingFactor

sf = ScalingFactor(run_list=range(197912, 197932),
                   sf_file="/tmp/sf_197912_si.cfg",
                   medium='Si')

sf.execute()

