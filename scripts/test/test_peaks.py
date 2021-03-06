# import mantid algorithms, numpy and matplotlib
import sys
sys.path.append('../autoreduce')
from mantid.simpleapi import *
import numpy as np
from scipy import ndimage
from peak_finding import find_peaks, peak_prominences, peak_widths


def scan_peaks(x):
    f1 = ndimage.gaussian_filter(x, 3)
    peaks, _ = find_peaks(f1)
    prom, _, _ = peak_prominences(f1, peaks)
    peaks_w, _, _, _ = peak_widths(f1, peaks)

    # The quality factor is the size of the peak (height*width).
    quality = -peaks_w * prom

    zipped = zip(peaks, peaks_w, quality, prom)
    ordered = sorted(zipped, key=lambda a:a[2])
    found_peaks = [[p[0],p[1]] for p in ordered]

    return found_peaks

ws = Load("REF_L_179932")

n_x = int(ws.getInstrument().getNumberParameter("number-of-x-pixels")[0])
n_y = int(ws.getInstrument().getNumberParameter("number-of-y-pixels")[0])

_integrated = Integration(InputWorkspace=ws)
signal = _integrated.extractY()
z=np.reshape(signal, (n_x, n_y))
x = z.sum(axis=0)

peaks = scan_peaks(x)

print(peaks[0][0], peaks[0][1])