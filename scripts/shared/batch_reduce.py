#!/usr/bin/python3
import os
import subprocess
import sys
import time

PYTHON_CMD = "nsd-conda-wrap.sh mantid"

if len(sys.argv) < 4:
    print("\nUsage: python3 batch_reduce.py <IPTS> <first run> <last run>")
    print("\nExample:\n   python3 batch_reduce.py IPTS-20406 178195 178216")
    sys.exit(0)

try:
    first_run = int(sys.argv[2])
    last_run = int(sys.argv[3])
except:
    print("Your run range looks wrong: %s %s" % (sys.argv[2], sys.argv[3]))
    sys.exit(0)
print("Run range: %g - %g" % (first_run, last_run))

ipts = sys.argv[1]

# Look for additional options for the new reduction process
new_version = False
if len(sys.argv) > 4 and sys.argv[4] == "new":
    new_version = True

template_file = None
if len(sys.argv) > 5:
    template_file = sys.argv[5]

avg_overlap = True
if len(sys.argv) > 6:
    avg_overlap = sys.argv[6]

const_q = False
if len(sys.argv) > 7:
    const_q = sys.argv[7]

print("Using new version: %s" % new_version)
print("Using template: %s" % template_file)
print("  Average overlap: %s" % avg_overlap)
print("  Constant-Q binning: %s" % const_q)

t_0 = time.time()
for r in range(first_run, last_run + 1):
    _data_file_path = os.path.join("/SNS", "REF_L", ipts, "nexus", "REF_L_%d.nxs.h5" % r)
    _output_dir = os.path.join("/SNS", "REF_L", ipts, "shared", "autoreduce")
    if not os.path.isfile(_data_file_path):
        print("File does not exist: %s" % _data_file_path)
    else:
        print("Processing %s" % _data_file_path)
        if new_version:
            cmd = "%s /SNS/REF_L/shared/autoreduce/reduce_REF_L.py %s %s new %s %s %s" % (
                PYTHON_CMD,
                _data_file_path,
                _output_dir,
                template_file,
                avg_overlap,
                const_q,
            )
        else:
            if template_file is not None:
                cmd = "%s /SNS/REF_L/shared/autoreduce/reduce_REF_L.py %s %s old %s" % (
                    PYTHON_CMD,
                    _data_file_path,
                    _output_dir,
                    template_file,
                )
            else:
                cmd = "%s /SNS/REF_L/shared/autoreduce/reduce_REF_L.py %s %s" % (PYTHON_CMD, _data_file_path, _output_dir)

        out_log = os.path.join("/SNS", "REF_L", ipts, "shared", "autoreduce", "reduction_log", "REF_L_%d.nxs.h5.log" % r)
        out_err = os.path.join("/SNS", "REF_L", ipts, "shared", "autoreduce", "reduction_log", "REF_L_%d.nxs.h5.err" % r)
        logFile = open(out_log, "w")
        errFile = open(out_err, "w")
        proc = subprocess.Popen(
            cmd, shell=True, stdin=subprocess.PIPE, stdout=logFile, stderr=errFile, universal_newlines=True, cwd=_output_dir
        )
        proc.communicate()
        logFile.close()
        errFile.close()

        success = not os.path.isfile(out_err) or os.stat(out_err).st_size == 0
        if success:
            if os.path.isfile(out_err):
                os.remove(out_err)
        else:
            print("   Errors were found. Check %s" % out_log)

t_1 = time.time()
print("\nElapsed time: %g s [%g s/run]" % ((t_1 - t_0), ((t_1 - t_0) / (last_run - first_run + 1))))
