"""
    Acquire direct beams by decomposing each DB into NxN components by scanning the centers of Si and S1

    Adapted from Composite_DB_scan_LOOPS_60Hz_std.py by ESW, 2024-07-13
"""
import sys
import time
import json
import argparse


sys.path.append('/home/controls/var/tmp/scripts')

import instrument
import db_collector

# Instrument configurations to acquire DBs for.
# This should be read from a scan.csv file.
SCAN_60Hz = [[100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 15,     's1:Y:Gap': 0.39, 'si:Y:Gap': 0.25, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (1,1), 'C-DB 60Hz 15A'],
             [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 12.386, 's1:Y:Gap': 0.39, 'si:Y:Gap': 0.25, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (1,1), 'C-DB 60Hz 12.39A'],
             [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 9.74,   's1:Y:Gap': 0.39, 'si:Y:Gap': 0.25, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (1,1), 'C-DB 60Hz 9.74A'],
             [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 7.043,  's1:Y:Gap': 0.39, 'si:Y:Gap': 0.25, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (1,1), 'C-DB 60Hz 7.04A'],
             [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.25,  's1:Y:Gap': 0.39, 'si:Y:Gap': 0.25, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (2,2), 'C-DB 60Hz 4.25A'],
             [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.25,  's1:Y:Gap': 0.769, 'si:Y:Gap': 0.493, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (4,4), 'C-DB 60Hz 4.25A'],
             [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.25,  's1:Y:Gap': 1.523, 'si:Y:Gap': 0.976, 's3:Y:Gap': 20, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (5,5), 'C-DB 60Hz 4.25A'],
             [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.25,  's1:Y:Gap': 3.015, 'si:Y:Gap': 1.932, 's3:Y:Gap': 20, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (10,10), 'C-DB 60Hz 4.25A'],
             ]

SCAN_60Hz_New = [[1000, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 9.2,     's1:Y:Gap': 0.259, 'si:Y:Gap': 0.166, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (1,1), 'C-DB 60Hz-new 9.2A'],
                 [300, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 6.7, 's1:Y:Gap': 0.259, 'si:Y:Gap': 0.166, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (1,1), 'C-DB 60Hz-new 6.7A'],
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 6.7, 's1:Y:Gap': 0.259, 'si:Y:Gap': 0.166, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (2,2), 'C-DB 60Hz-new 6.7A'],
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,   's1:Y:Gap': 0.259, 'si:Y:Gap': 0.166, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (2,2), 'C-DB 60Hz-new 4.2A'],
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,   's1:Y:Gap': 0.259, 'si:Y:Gap': 0.166, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (4,4), 'C-DB 60Hz-new 4.2A'],
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 0.517, 'si:Y:Gap': 0.332, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (2,2), 'C-DB 60Hz-new 4.2A'],
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 0.517, 'si:Y:Gap': 0.332, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (4,4), 'C-DB 60Hz-new 4.2A'],
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 0.517, 'si:Y:Gap': 0.332, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (8,8), 'C-DB 60Hz-new 4.2A'],
                 
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 1.035, 'si:Y:Gap': 0.663, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (8,8), 'C-DB 60Hz-new 4.2A'],
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 1.035, 'si:Y:Gap': 0.663, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (10,10), 'C-DB 60Hz-new 4.2A'],
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 1.035, 'si:Y:Gap': 0.663, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (15,15), 'C-DB 60Hz-new 4.2A'],
                 
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 2.069, 'si:Y:Gap': 1.326, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (8,8), 'C-DB 60Hz-new 4.2A'],
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 2.069, 'si:Y:Gap': 1.326, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (10,10), 'C-DB 60Hz-new 4.2A'],
                 [100, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 2.069, 'si:Y:Gap': 1.326, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (15,15), 'C-DB 60Hz-new 4.2A'],
                 ]    

SCAN_20Hz = [   
             [100, {'BL4B:Chop:Gbl:SpeedReq': 20, 'BL4B:Chop:Gbl:WavelengthReq': 6., 'thm':-0.272, 's1:Y:Gap': 0.205, 'si:Y:Gap': 0.205, 's3:Y:Gap': 10, 'ths': 0.0, 'tthd': 0.0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (1,1), 'C-DB 20Hz 6A'],  # A1 DIV 0
             [100, {'BL4B:Chop:Gbl:SpeedReq': 20, 'BL4B:Chop:Gbl:WavelengthReq': 6., 'thm':-0.272, 's1:Y:Gap': 0.205, 'si:Y:Gap': 0.205, 's3:Y:Gap': 10, 'ths': 0.0, 'tthd': 0.0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (2,2), 'C-DB 20Hz 6A'],  # A1 DIV 0
             [100, {'BL4B:Chop:Gbl:SpeedReq': 20, 'BL4B:Chop:Gbl:WavelengthReq': 6., 'thm':-0.272, 's1:Y:Gap': 0.205, 'si:Y:Gap': 0.205, 's3:Y:Gap': 10, 'ths': 0.0, 'tthd': 0.0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (3,3), 'C-DB 20Hz 6A'],  # A1 DIV 0
             
             [100, {'BL4B:Chop:Gbl:SpeedReq': 20, 'BL4B:Chop:Gbl:WavelengthReq': 6., 'thm':-0.272, 's1:Y:Gap': 0.912, 'si:Y:Gap': 0.912, 's3:Y:Gap': 20, 'ths': 0.0, 'tthd': 0.0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (3,3), 'C-DB 20Hz 6A'],   # A2 DIV 0
             [100, {'BL4B:Chop:Gbl:SpeedReq': 20, 'BL4B:Chop:Gbl:WavelengthReq': 6., 'thm':-0.272, 's1:Y:Gap': 0.912, 'si:Y:Gap': 0.912, 's3:Y:Gap': 20, 'ths': 0.0, 'tthd': 0.0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (8,8), 'C-DB 20Hz 6A'],   # A2 DIV 0
             [100, {'BL4B:Chop:Gbl:SpeedReq': 20, 'BL4B:Chop:Gbl:WavelengthReq': 6., 'thm':-0.272, 's1:Y:Gap': 0.912, 'si:Y:Gap': 0.912, 's3:Y:Gap': 20, 'ths': 0.0, 'tthd': 0.0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (10,10), 'C-DB 20Hz 6A'], # A2 DIV 0
             
             [100, {'BL4B:Chop:Gbl:SpeedReq': 20, 'BL4B:Chop:Gbl:WavelengthReq': 6., 'thm':-0.272, 's1:Y:Gap': 3.624, 'si:Y:Gap': 0.725, 's3:Y:Gap': 10, 'ths': 0.0, 'tthd': 0.0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (10,10), 'C-DB 20Hz 6A'],  # A2 DIV 5
             [100, {'BL4B:Chop:Gbl:SpeedReq': 20, 'BL4B:Chop:Gbl:WavelengthReq': 6., 'thm':-0.272, 's1:Y:Gap': 3.624, 'si:Y:Gap': 0.725, 's3:Y:Gap': 10, 'ths': 0.0, 'tthd': 0.0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (20,20), 'C-DB 20Hz 6A'],  # A2 DIV 5
             [100, {'BL4B:Chop:Gbl:SpeedReq': 20, 'BL4B:Chop:Gbl:WavelengthReq': 6., 'thm':-0.272, 's1:Y:Gap': 5.769, 'si:Y:Gap': 0.577, 's3:Y:Gap': 20, 'ths': 0.0, 'tthd': 0.0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (10,10), 'C-DB 20Hz 6A'],  # A2 DIV 10
             [100, {'BL4B:Chop:Gbl:SpeedReq': 20, 'BL4B:Chop:Gbl:WavelengthReq': 6., 'thm':-0.272, 's1:Y:Gap': 5.769, 'si:Y:Gap': 0.577, 's3:Y:Gap': 20, 'ths': 0.0, 'tthd': 0.0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (20,20), 'C-DB 20Hz 6A'],  # A2 DIV 10
             ]



SCAN_TO_RUN = [
    [100, {'BL4B:Chop:Gbl:SpeedReq': 30, 'BL4B:Chop:Gbl:WavelengthReq': 6, 's1:Y:Gap': 1.202, 'si:Y:Gap': 0.129, 's3:Y:Gap': 30, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20, 'BL4B:Actuator:50M':1, 'BL4B:Actuator:100M':1, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':1}, (2,2), 'C-DB 30Hz angle 1'],
    [100, {'BL4B:Chop:Gbl:SpeedReq': 30, 'BL4B:Chop:Gbl:WavelengthReq': 6, 's1:Y:Gap': 3.557, 'si:Y:Gap': 0.381, 's3:Y:Gap': 30, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20, 'BL4B:Actuator:50M':0, 'BL4B:Actuator:100M':1, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':1}, (4,4), 'C-DB 30Hz angle 2'],
    [100, {'BL4B:Chop:Gbl:SpeedReq': 30, 'BL4B:Chop:Gbl:WavelengthReq': 6, 's1:Y:Gap': 10.528, 'si:Y:Gap': 1.128, 's3:Y:Gap': 30, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20, 'BL4B:Actuator:50M':1, 'BL4B:Actuator:100M':0, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':1}, (4,4), 'C-DB 30Hz angle 3'],
    [1, {'BL4B:Chop:Gbl:SpeedReq': 30, 'BL4B:Chop:Gbl:WavelengthReq': 6, 's1:Y:Gap': 10.528, 'si:Y:Gap': 1.128, 's3:Y:Gap': 30, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20, 'BL4B:Actuator:50M':1, 'BL4B:Actuator:100M':1, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':1}, (4,4), 'Take Cd out'],
                        
               ]


# Example usage
if __name__ == "__main__":
    # TODO Read the scan configuration from a file
    parser = argparse.ArgumentParser(description="Acquire direct beams for the reflectometer.")
    parser.add_argument("--scan", help="The scan configuration file to use.")
    parser.add_argument("--title", help="The title of the scan [like the medium].")
    parser.add_argument("--charge", help="Amount of charge to collect [replaces default in file].", default=100, type=float)
    args = parser.parse_args()

    if args.scan:
        with open(args.scan, 'r') as f:
            SCAN_TO_RUN = json.load(f)

    collector = db_collector.DBCollector(SCAN_TO_RUN, charge=args.charge)
    collector.collect()
