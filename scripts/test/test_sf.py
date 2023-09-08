import sys

from autoreduce.sf_calculator import ScalingFactor

sf = ScalingFactor(
    run_list=range(184975, 184990), sf_file="/tmp/sf_184975_air.cfg", medium="air"
)

sf.execute()
