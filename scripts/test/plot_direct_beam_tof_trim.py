"""
Manual check of TrimDirectBeamTof on one direct-beam run.

Runs the steps before the trim in the order of the direct-beam composition, prints the three
bounds and the kept TOF range, and plots the ROI intensity before and after the dead time
correction and after the trim.

Usage:
    pixi run python scripts/test/plot_direct_beam_tof_trim.py [NEXUS_FILE] [options]

Run with --help for the options. Without NEXUS_FILE the run REF_L_234266 of the test data
repository is used.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from mantid.kernel import config

REPOSITORY = Path(__file__).resolve().parents[2]

# mantid.simpleapi only wraps the algorithms it finds in the plugin directories when first imported,
# so the repository's algorithms are added for this session (not saved to the user's properties).
_PLUGINS = "user.python.plugins.directories"
config[_PLUGINS] = ";".join(filter(None, [config[_PLUGINS], str(REPOSITORY / "src/lr_reduction/operations/algorithms")]))

from mantid.simpleapi import (  # noqa: E402

    ApplyDeadTimeCorrection,
    ExtractDetectorROI,
    LoadErrorEventsNexus,
    LoadEventNexus,
    NormaliseByCurrent,
    Rebin,
    SingleReadoutDeadTimeCorrection,
    SumSpectra,
    TrimDirectBeamTof,
)

DEFAULT_RUN = REPOSITORY / "tests/data/liquidsreflectometer-data/nexus/REF_L_234266.nxs.h5"
PLOT_BINS = "8000,100,56000"  # tof_min, step, tof_max in µs


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("nexus_file", nargs="?", default=str(DEFAULT_RUN), help="direct-beam event NeXus file")
    parser.add_argument("--x-range", type=int, nargs=2, default=[70, 190], help="inclusive ROI x pixels")
    parser.add_argument("--y-range", type=int, nargs=2, default=[147, 155], help="inclusive ROI y pixels")
    parser.add_argument("--dead-time-threshold", type=float, default=1.25, help="largest reliable DTC factor")
    parser.add_argument(
        "--intensity-threshold", type=float, default=0.36, help="intensity floor per µA·h per 50 µs bin"
    )
    return parser.parse_args()


def summed_histogram(workspace, name):
    summed = SumSpectra(InputWorkspace=workspace, OutputWorkspace=name)
    return Rebin(InputWorkspace=summed, Params=PLOT_BINS, PreserveEvents=False, OutputWorkspace=name)


def main():
    args = parse_arguments()

    # The curve is computed on the full detector, before the ROI extraction and the normalization.
    LoadEventNexus(Filename=args.nexus_file, OutputWorkspace="tof")
    LoadErrorEventsNexus(Filename=args.nexus_file, OutputWorkspace="tof_error_events")
    SingleReadoutDeadTimeCorrection(
        InputWorkspace="tof",
        InputErrorEventsWorkspace="tof_error_events",
        UseDeadTimeThreshold=False,
        OutputWorkspace="dead_time_curve",
    )

    ExtractDetectorROI(InputWorkspace="tof", XRange=args.x_range, YRange=args.y_range, OutputWorkspace="roi")
    NormaliseByCurrent(InputWorkspace="roi", OutputWorkspace="roi")
    ApplyDeadTimeCorrection(
        InputWorkspace="roi", DeadTimeCorrectionWorkspace="dead_time_curve", OutputWorkspace="roi_dead_time"
    )
    _, tof_min, tof_max, chopper_cut, dead_time_cut, intensity_cut = TrimDirectBeamTof(
        InputWorkspace="roi_dead_time",
        CorrectionWorkspace="dead_time_curve",
        DeadTimeThreshold=args.dead_time_threshold,
        IntensityThreshold=args.intensity_threshold,
        OutputWorkspace="trimmed",
    )

    print(f"chopper cut:    {chopper_cut:.1f} µs")
    print(f"dead time cut:  {dead_time_cut:.1f} µs (0 when no bin exceeds the threshold)")
    print(f"intensity cut:  {intensity_cut:.1f} µs")
    print(f"kept TOF range: [{tof_min:.1f}, {tof_max:.1f}) µs")

    fig, ax = plt.subplots(subplot_kw={"projection": "mantid"})
    for workspace, label in (
        ("roi", "ROI, normalized"),
        ("roi_dead_time", "+ dead time correction"),
        ("trimmed", "+ TOF trim"),
    ):
        ax.plot(summed_histogram(workspace, f"{workspace}_sum"), wkspIndex=0, label=label)
    ax.axvspan(tof_min, tof_max, alpha=0.1, color="gray", label="kept TOF range")
    ax.set_title(Path(args.nexus_file).name)
    ax.legend()
    plt.show()


if __name__ == "__main__":
    main()
