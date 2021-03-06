#!/usr/bin/python3
import sys
import os
import time
import subprocess

def local_submission(cmd, out_log, out_err):
    cmd = "%s %s %s %s/" % (configuration.python_executable, script, input_file, output_dir)
    logFile=open(out_log, "w")
    errFile=open(out_err, "w")
    if configuration.comm_only is False:
        proc = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                                stdout=logFile, stderr=errFile, universal_newlines = True,
                                cwd=output_dir)
        proc.communicate()
    logFile.close()
    errFile.close()

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

t_0 = time.time()
for r in range(first_run, last_run+1):
    _data_file_path = os.path.join('/SNS', 'REF_L', ipts, 'nexus', 'REF_L_%d.nxs.h5' % r)
    _output_dir = os.path.join('/SNS', 'REF_L', ipts, 'shared', 'autoreduce')
    if not os.path.isfile(_data_file_path):
        print("File does not exist: %s" % _data_file_path)
    else:
        print("Processing %s" % _data_file_path)
        cmd = "python3 /SNS/REF_L/shared/autoreduce/reduce_REF_L.py %s %s" % (_data_file_path, _output_dir)
        out_log = os.path.join('/SNS', 'REF_L', ipts, 'shared', 'autoreduce', 'reduction_log', 'REF_L_%d.nxs.h5.log' % r)
        out_err = os.path.join('/SNS', 'REF_L', ipts, 'shared', 'autoreduce', 'reduction_log', 'REF_L_%d.nxs.h5.err' % r)
        logFile=open(out_log, "w")
        errFile=open(out_err, "w")
        proc = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                                stdout=logFile, stderr=errFile, universal_newlines = True,
                                cwd=_output_dir)
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
print("\nElapsed time: %g s [%g s/run]" % ((t_1-t_0), ((t_1-t_0)/(last_run-first_run+1))))
