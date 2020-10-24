#!/usr/bin/python
"""
Example:

process_sf.py Si 178195 178216 sf_178195_Si2InDiam.cfg

"""
import sys
sys.path.append('/SNS/REF_L/shared/autoreduce')

import os
if len(sys.argv) < 5:
    print("\nUsage: python process_sf.py <incident medium> <first run> <last run> <cfg file name>")
    print("\nExample:\n   python process_sf.py Si 178195 178216 /SNS/REF_L/shared/autoreduce/sf_178195_Si2InDiam.cfg")
    sys.exit(0)

print("Incident medium: %s" % sys.argv[1])

try:
    first_run = int(sys.argv[2])
    last_run = int(sys.argv[3])
except:
    print("Your run range looks wrong: %s %s" % (sys.argv[2], sys.argv[3]))
    sys.exit(0)
print("Run range: %g - %g" % (first_run, last_run))


_fpath = os.path.abspath(sys.argv[4])
_dir, _file = os.path.split(_fpath)

if len(_dir)>0 and not os.path.isdir(_dir):
    print("The file you asked for is not in an existing directory: %s" % _dir)
    sys.exit(0)
print("Output file: %s" % _fpath)


from sf_calculator import ScalingFactor

sf = ScalingFactor(run_list=range(first_run,last_run+1),
                   sort_by_runs=True,
                   sf_file=_fpath,
                   tof_step=200,
                   medium=sys.argv[1],
                   slit_tolerance=0.06)
sf.execute()

