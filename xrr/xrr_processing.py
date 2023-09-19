import os
import sys
import warnings

import numpy as np
from matplotlib import pyplot as plt

warnings.filterwarnings("ignore", module="numpy")
warnings.filterwarnings("ignore")


WAVELENGTH_META = "HW_XG_WAVE_LENGTH_ALPHA1"


def process_xrr(data_file, output_dir=None):
    """
    Process Rigaku .ras files to produce R(Q).

    data_file: full file path of the data file to process
    output_dir: optional output directory
    """
    data = np.loadtxt(data_file, comments=["#", "*"]).T

    # If no output directory was provided, use the location of the data file
    if output_dir is None:
        output_dir, _ = os.path.split(data_file)

    # Read meta data
    meta_data = dict()
    with open(data_file) as fd:
        for line in fd:
            if line.startswith("*"):
                toks = line.split()
                if len(toks) < 2:
                    # Single keywords are used to define meta data sections, skip them
                    pass
                else:
                    value = toks[1].replace('"', "")
                    try:
                        value = float(value)
                    except:
                        # keep value as a string
                        pass
                    meta_data[toks[0][1:]] = value

    if WAVELENGTH_META in meta_data:
        wl = float(meta_data[WAVELENGTH_META])
        print("X-ray wavelength found: %g A" % wl)
    else:
        wl = 1.5406
        print("Could not find wavelength, using %g A" % wl)

    # Points below q_min and above q_max will be cut
    q_min = 0.03
    q_max = 0.5

    ttheta = data[0]
    counts = data[1]

    q = 4 * np.pi / wl * np.sin(ttheta / 2 * np.pi / 180)

    # Select only points in a useful range
    _q_idx = (q > q_min) & (q < q_max)
    q = q[_q_idx]
    r = counts[_q_idx]

    # R(q) will be normalized to the average between q_min and norm_q_max
    norm_q_max = 0.04
    _q_idx = (q > q_min) & (q < norm_q_max)
    _norm = np.sum(r[_q_idx]) / len(q[_q_idx])
    print("Norm %g" % _norm)
    r /= _norm
    err = r * 0.1
    dq = np.ones(len(q)) * 0.002
    _rq_data = np.asarray([q, r, err, dq]).T

    # Choose a base name for the output files
    _name, _ext = os.path.splitext(data_file)

    np.savetxt(os.path.join(output_dir, "%s-Rq.txt" % _name), _rq_data)

    plt.figure(figsize=(10, 6))
    plt.plot(q, r)
    plt.xlabel(r"q [$1/\AA$]")
    plt.ylabel("R(q)")
    plt.yscale("log")
    plt.xscale("linear")
    plt.savefig(os.path.join(output_dir, "%s-Rq.png" % _name))


if len(sys.argv) < 2:
    print("\nUsage: python3 xrr_processing.py <data file in RAS format>")
    print("\nExample:\n   python3 xrr_processing.py xrr_data_file.ras")
    sys.exit(0)

# Data file to reduce. Must be in .ras format
data_file = sys.argv[1]

process_xrr(data_file)
