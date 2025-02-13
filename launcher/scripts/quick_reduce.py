import argparse
import os
import sys

import mantid
import mantid.simpleapi as api
import numpy as np

mantid.kernel.config.setLogLevel(3)

from matplotlib import pyplot as plt

sys.path.append("/SNS/REF_L/shared/reduction")

from lr_reduction import workflow


def process(run_number, db_run_number, peak_pixel, db_peak_pixel, output_dir=None):
    def _tof(_ws):
        tof_min = _ws.getTofMin()
        tof_max = _ws.getTofMax()
        _ws = api.Rebin(InputWorkspace=_ws, Params="%s,100,%s" % (tof_min, tof_max))
        y = _ws.extractY()
        x = _ws.extractX()
        charge = _ws.getRun().getProtonCharge()
        tof = (x[0][1:] + x[0][:-1]) / 2
        counts = np.sum(y, axis=0) / charge
        return tof, counts

    ws = api.Load("REF_L_%s" % run_number)
    ws_db = api.Load("REF_L_%s" % db_run_number)

    fig, ax = plt.subplots(2, 1, figsize=(10, 10))
    ax = plt.subplot(2, 1, 1)

    tof, counts = _tof(ws)
    plt.plot(tof, counts, label="SC %s" % run_number)

    tof, counts = _tof(ws_db)
    plt.plot(tof, counts, label="DB %s" % db_run_number)

    plt.xlabel("TOF")
    plt.ylabel("Counts")
    ax.set_yscale("linear")
    ax.set_xscale("linear")
    plt.legend()

    qz, r, dr = workflow.reduce_explorer(ws, ws_db, center_pixel=peak_pixel, db_center_pixel=db_peak_pixel)

    # Save reduced data
    np.savetxt(os.path.join(output_dir, "r%s_quick_reduce.txt" % run_number), np.asarray([qz, r, dr]).T)

    ax = plt.subplot(2, 1, 2)
    plt.errorbar(qz, r, yerr=dr)
    plt.xlabel("q [$1/\AA$]")
    plt.ylabel("R(q)")
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument("run_number", type=str, help="Run number to process")
    parser.add_argument("db_run_number", type=str, help="Direct beam run number")
    parser.add_argument("pixel", type=float, help="Peak pixel")
    parser.add_argument("db_pixel", type=float, help="Direct beam peak pixel")
    parser.add_argument("output_dir", type=str, help="Output directory")

    # Parse arguments
    args = parser.parse_args()
    process(
        args.run_number, args.db_run_number, np.rint(args.pixel).astype(int), np.rint(args.db_pixel).astype(int), output_dir=args.output_dir
    )
