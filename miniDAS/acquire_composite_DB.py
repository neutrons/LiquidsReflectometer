"""
    Acquire direct beams by decomposing each DB into NxN components by scanning the centers of Si and S1

    Adapted from Composite_DB_scan_LOOPS_60Hz_std.py by ESW, 2024-07-13
"""
import sys
import time
import argparse


sys.path.append('/home/controls/var/tmp/scripts')

import instrument

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
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 9.2,  's1:Y:Gap': 0.259, 'si:Y:Gap': 0.166, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (1,1), 'C-DB A1 60Hz-new 9.2A'],
                 
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 6.7,  's1:Y:Gap': 0.259, 'si:Y:Gap': 0.166, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (1,1), 'C-DB A2 60Hz-new 6.7A'],
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 6.7,  's1:Y:Gap': 0.259, 'si:Y:Gap': 0.166, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (2,2), 'C-DB A2 60Hz-new 6.7A'],
                 
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 0.259, 'si:Y:Gap': 0.166, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (2,2), 'C-DB A3 60Hz-new 4.2A'],
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 0.259, 'si:Y:Gap': 0.166, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (4,4), 'C-DB A3 60Hz-new 4.2A'],
                 
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 0.517, 'si:Y:Gap': 0.332, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (2,2), 'C-DB A4 60Hz-new 4.2A'],
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 0.517, 'si:Y:Gap': 0.332, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (4,4), 'C-DB A4 60Hz-new 4.2A'],
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 0.517, 'si:Y:Gap': 0.332, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (8,8), 'C-DB A4 60Hz-new 4.2A'],
                 
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 1.035, 'si:Y:Gap': 0.663, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (8,8), 'C-DB A5 60Hz-new 4.2A'],
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 1.035, 'si:Y:Gap': 0.663, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (10,10), 'C-DB A5 60Hz-new 4.2A'],
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 1.035, 'si:Y:Gap': 0.663, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (15,15), 'C-DB A5 60Hz-new 4.2A'],
                 
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 2.069, 'si:Y:Gap': 1.326, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (8,8), 'C-DB A6 60Hz-new 4.2A'],
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 2.069, 'si:Y:Gap': 1.326, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (10,10), 'C-DB A6 60Hz-new 4.2A'],
                 [200, {'BL4B:Chop:Gbl:SpeedReq': 60, 'BL4B:Chop:Gbl:WavelengthReq': 4.2,  's1:Y:Gap': 2.069, 'si:Y:Gap': 1.326, 's3:Y:Gap': 10, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20}, (15,15), 'C-DB A6 60Hz-new 4.2A'],
                 
                 ]   


class DBCollector:
    def __init__(self, db_list: list):
        """
        Initializes the collector with a list of direct beam configurations.
        
        :param db_list: List of direct beam configurations.
        """
        self.db_list = db_list

        # Initialize the reflectometer
        self.lr = instrument.LiquidsReflectometer()

    def collect(self):
        """
        Collects the direct beams.
        """
        self.lr.initialize_series(length=len(self.db_list))

        for i, db in enumerate(self.db_list):
            print("Direct beam configuration: ", i)
            si_width = self.db_list[i][1]['si:X:Gap']
            s1_width = self.db_list[i][1]['s1:X:Gap']
            div = s1_width/si_width
            run_title = "%s: N=%sx%s Div=%g " % (self.db_list[i][3], self.db_list[i][2][0], self.db_list[i][2][1], div)
            
            self.lr.increment_sequence(title=run_title)
            try:
                scanner = CompositeDBScanner(self.db_list[i], grid_size=self.db_list[i][2])
                scanner.scan_centers()
            except:
                print("Error occurred while scanning direct beams", i)
                raise

        self.lr.stop()

class CompositeDBScanner:
    def __init__(self, positions: list, grid_size=None,
                 instr: instrument.LiquidsReflectometer = None):
        """
        Initializes the scanner with a frequency and grid size.
        
        :param positions: Instrument configuration to acquire a direct beam for.
        :param grid_size: The size of the grid for NxN decomposition (Si x S1)
        """
        self.positions = positions

        # Default grid size to 10x10 if not provided
        self.grid_size = grid_size if grid_size is not None else (10, 10)

        # Verify that the position contain the relevant information
        print(self.positions[1].keys())
        if 's1:X:Gap' not in self.positions[1] or 'si:X:Gap' not in self.positions[1]:
            raise ValueError("Invalid position configuration. Missing Si or S1 gap information.")

        if instr is not None:
            self.lr = instr
        else:
            self.lr = instrument.LiquidsReflectometer()

    def scan_centers(self):
        """
        Scans the centers of Si and S1.
        """
        # Move to the nominal position for all motors
        self.lr.move({'si:X:Center': 0, 's1:X:Center': 0})
        self.lr.move(self.positions[1])

        # Determine positions to raster over
        si_width = self.positions[1]['si:X:Gap']
        s1_width = self.positions[1]['s1:X:Gap']
        self.lr.move({'si:X:Gap': si_width/self.grid_size[0], 's1:X:Gap': s1_width/self.grid_size[1]})

        # The starting center should half a step from the left-most position
        si_start = si_width * (-1 + 1 / self.grid_size[0]) / 2
        s1_start = si_width * (-1 + 1 / self.grid_size[1]) / 2

        si_positions = [si_start + i * si_width / self.grid_size[0] for i in range(self.grid_size[0])]
        s1_positions = [s1_start + i * s1_width / self.grid_size[1] for i in range(self.grid_size[1])]

        
        charge_to_acquire_per_point = self.positions[0] / (self.grid_size[0] * self.grid_size[1])
        print("Charge to acquire per configuration:", charge_to_acquire_per_point)

        # Iterate over Si
        for si in si_positions:
            # Iterate over S1
            for s1 in s1_positions:
                print(f"Si X center: {si}\tS1 X center: {s1}")
                # Move motors to the specified positions
                self.lr.move({'si:X:Center': si, 's1:X:Center': s1})
                time.sleep(0.5)
                # Acquire neutrons
                self.lr.start_or_resume(charge=charge_to_acquire_per_point)

                # Pause to allow the next move
                self.lr.pause()

                rate = self.lr.get_rate()
                print(f"    Rate: {rate}")
                time.sleep(0.5)
        self.lr.stop()
        time.sleep(2)

# Example usage
if __name__ == "__main__":
    # TODO Read the scan configuration from a file
    parser = argparse.ArgumentParser(description="Acquire direct beams for the reflectometer.")
    parser.add_argument("--scan", help="The scan configuration file to use.")
    parser.add_argument("--title", help="The title of the scan [like the medium].")
    args = parser.parse_args()

    collector = DBCollector(SCAN_TO_RUN)
    collector.collect()
