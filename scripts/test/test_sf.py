import sys
sys.path.append('../autoreduce')

from sf_calculator import ScalingFactor

sf = ScalingFactor(run_list=range(180459, 180483),
                   sf_file="/tmp/sf_180459_Si_md.cfg")
sf.execute()

