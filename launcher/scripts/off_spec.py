import sys
import os
import argparse
import numpy as np

from matplotlib import pyplot as plt

import mantid
import mantid.simpleapi as api
mantid.kernel.config.setLogLevel(3)


DEFAULT_4B_SAMPLE_DET_DISTANCE = 1.83
DEFAULT_4B_SOURCE_DET_DISTANCE = 15.75


class PixelData():
    source_detector_distance = DEFAULT_4B_SOURCE_DET_DISTANCE
    det_distance = DEFAULT_4B_SAMPLE_DET_DISTANCE

    def __init__(self, run_number):
        self.ws = api.LoadEventNexus("REF_L_%s" % run_number, OutputWorkspace='r%s' % run_number)
        self.run_number = run_number
        self.get_parameters()

    def __repr__(self):
        _repr = "# Instrument parameters. All distance in meters\n"
        _repr += "# For more info about off-specular data, see for example:\n"
        _repr += "#   DOI: https://doi.org/10.1103/PhysRevB.38.2297\n"
        _repr += "# Pixel 0 is at the bottom of the detector\n"
        _repr += "# sample-detector-distance: %s\n" % self.det_distance
        _repr += "# source-detector-distance: %s\n" % self.source_detector_distance
        _repr += "# detector-size: %sx%s\n" % (self.n_x, self.n_y)
        _repr += "# pixel-width: %s\n" % self.pixel_width
        _repr += "# two-theta: %s\n" % self.tthd
        return _repr

    def get_parameters(self):
        if self.ws.getInstrument().hasParameter("sample-det-distance"):
            self.det_distance = self.ws.getInstrument().getNumberParameter("sample-det-distance")[0]

        if self.ws.getInstrument().hasParameter("source-det-distance"):
            self.source_detector_distance = self.ws.getInstrument().getNumberParameter("source-det-distance")[0]

        # Set up basic data
        self.n_x = int(self.ws.getInstrument().getNumberParameter("number-of-x-pixels")[0])
        self.n_y = int(self.ws.getInstrument().getNumberParameter("number-of-y-pixels")[0])

        self.pixel_width = float(self.ws.getInstrument().getNumberParameter("pixel-width")[0]) / 1000.0

        self.tthd = self.ws.getRun()['tthd'].value[0]

        h = 6.626e-34  # m^2 kg s^-1
        m = 1.675e-27  # kg
        self.constant = 1e-4 * m * self.source_detector_distance / h

    def process(self, output_dir, wavelength_step=1):
        print("Processing: %s" % self.run_number)

        tof_step = wavelength_step * self.constant

        tof_min = self.ws.getTofMin()
        tof_max = self.ws.getTofMax()
        print("%s %s %s" % (tof_min, tof_step, tof_max))
        _ws = api.Rebin(InputWorkspace=self.ws, Params="%s,%s,%s" % (tof_min, tof_step, tof_max))
        y = _ws.extractY()
        x = _ws.extractX()

        y = np.reshape(y, (256, 304, y.shape[1]))
        p_vs_t = np.sum(y, axis=0)

        err = np.sqrt(p_vs_t)

        # Use center of wl bins
        wl = x[0] / self.constant
        wl = (wl[1:] + wl[:-1])/2.0

        with open(os.path.join(output_dir, 'r%s-wl.txt' % self.run_number), 'w') as fd:
            fd.write(str(self))
            fd.write("# pixel\t wavelength\t signal\t error\n")
            for i in np.arange(p_vs_t.shape[0]):
                for i_wl in range(len(wl)):
                    fd.write("%8g %8g %8g %8g\n" % (i, wl[i_wl], p_vs_t[i][i_wl], err[i][i_wl]))

        plot_data(p_vs_t, wl, 'run %s' % self.run_number,
                  os.path.join(output_dir, 'r%s-counts.png' % self.run_number))


def plot_data(counts, wl, title, file_path, show=True):
    counts_vs_wl = np.sum(counts, axis=0)
    counts_vs_pixel = np.sum(counts, axis=1)

    fig, ax = plt.subplots(2, 1, figsize=(6,10))

    plt.subplot(2, 1, 1)
    plt.plot(np.arange(counts_vs_pixel.shape[0]), counts_vs_pixel)

    plt.title('Total counts per pixel - %s' % title)
    plt.xlabel('pixel number')
    plt.ylabel('Counts')
    ax[0].set_yscale('log')

    plt.subplot(2, 1, 2)
    plt.plot(wl, counts_vs_wl)

    plt.title('Total counts vs wavelength - %s' % title)
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Counts')
    ax[1].set_yscale('log')

    if show:
        plt.show()
    plt.savefig(file_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument('run_number', type=str,
                        help='Run number to process')
    parser.add_argument('wl_step', type=float,
                        help='Wavelength bin width', default=0.1)
    parser.add_argument('output_dir', type=str,
                        help='Output directory')

    # Parse arguments
    args = parser.parse_args()

    p = PixelData(args.run_number)
    p.get_parameters()
    p.process(args.output_dir, args.wl_step)
