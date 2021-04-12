import sys
sys.path.append('../autoreduce')

from sf_calculator import ScalingFactor

# Testing subset
#    183112 wl= 4.25 thi=-8.216382133519986e-05 att=1    56,   60    -4,  259
#    183113 wl= 4.25 thi=-8.216382133519986e-05 att=2    60,   64    -4,  259
#    183114 wl= 4.25 thi=-8.216382133519986e-05 att=2    96,  141    -4,  259
#    183115 wl= 4.25 thi=-8.216382133519986e-05 att=2   137,  146   122,  132
# 

sf = ScalingFactor(run_list=range(183102, 183122),
                   sf_file="/tmp/sf_180459_Si_md.cfg",
                   medium='Air')
sf.execute()

