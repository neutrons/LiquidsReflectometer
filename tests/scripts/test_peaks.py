# import mantid algorithms, numpy and matplotlib
import unittest
import pytest
import os

from scipy import ndimage
import numpy as np

import mantid.simpleapi as mtd_api
from mantid import config
from scripts.autoreduce.peak_finding import find_peaks, peak_prominences, peak_widths


@pytest.mark.scripts()
class ScanPeaksTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.getcwd().endswith("LiquidsReflectometer"):
            os.chdir("tests")

        cwd = os.getcwd()
        dirs = [26010, 26776, 28662, 29196, 31279]
        for dir_num in dirs:
            config.appendDataSearchDir(str(os.path.join(cwd, f"data/liquidsreflectometer-data/SNS/REF_L/IPTS-{dir_num}/nexus")))

    @classmethod
    def test_peak_finding(self):
        def scan_peaks(x):
            f1 = ndimage.gaussian_filter(x, 3)
            peaks, _ = find_peaks(f1)
            prom, _, _ = peak_prominences(f1, peaks)
            peaks_w, _, _, _ = peak_widths(f1, peaks)

            # The quality factor is the size of the peak (height*width).
            quality = -peaks_w * prom

            zipped = zip(peaks, peaks_w, quality, prom)
            ordered = sorted(zipped, key=lambda a: a[2])
            found_peaks = [[p[0], p[1]] for p in ordered]

            return found_peaks

        ws = mtd_api.Load("REF_L_179932")

        n_x = int(ws.getInstrument().getNumberParameter("number-of-x-pixels")[0])
        n_y = int(ws.getInstrument().getNumberParameter("number-of-y-pixels")[0])

        _integrated = mtd_api.Integration(InputWorkspace=ws)
        signal = _integrated.extractY()
        z = np.reshape(signal, (n_x, n_y))
        x = z.sum(axis=0)

        peaks = scan_peaks(x)

        print(peaks[0][0], peaks[0][1])
