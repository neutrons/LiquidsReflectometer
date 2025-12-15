"""
Autoreduction script for the Liquids Reflectometer (BL-4B, REF_L) instrument

CLI Arguments
-------------
    ``events_file``
        Path to the Nexus events file.
    ``output_dir``
        Path to the output directory.
    ``old_version_flag``
        (Optional) Not used, kept for compatibility with existing scripts.
    ``template_file``
        (Optional) Path to the template XML file.
    ``avg_overlap``
        (Optional) If 'true', use average overlap method.
    ``const_q``
        (Optional) If 'true', use constant Q summing.
    ``fit_first_peak``
        (Optional) Not used, kept for compatibility with existing scripts.
    ``theta_offset``
        (Optional) Theta offset value.
    ``--no_publish``
        (Optional) If provided, do not upload HTML report to monitor.sns.gov.

Usage
-----

.. code-block:: console

    python reduce_REF_L.py <events_file> <output_dir> [old_version_flag] \
[template_file] [avg_overlap] [const_q] [fit_first_peak] [theta_offset] \
[--no_publish]

"""

import argparse
import os
import subprocess
from pathlib import Path

import numpy as np
from mantid import logger
from mantid.simpleapi import LoadEventNexus
from plot_publisher import plot1d

from lr_reduction import workflow
from lr_reduction.mantid_utils import SampleLogValues
from lr_reduction.template import get_default_template_file
from lr_reduction.typing import MantidWorkspace

# Name of the conda environment to use - required by autoreduction
CONDA_ENV = 'lr_reduction'


def parse_command_arguments():
    """
    Parse command line arguments

    Returns
    -------
    argparse.Namespace
        Object holding the arguments as attributes

    Notes
    -------
    The command line arguments are used during batch reduction in ``nr_launcher``.
    """
    parser = argparse.ArgumentParser(description='Autoreduction script for REF_L')
    # Mandatory positional arguments
    parser.add_argument('events_file', type=str, help='Path to the Nexus events file.')
    parser.add_argument('output_dir', type=str, help='Output directory path.')
    # Existing behavior: optional 3rd positional arg
    parser.add_argument("old_version_flag", nargs="?", default=None)
    # Existing behavior: optional positional 4-8 parameters
    parser.add_argument("template_file", nargs="?", default=None, type=str, help='Path to the template XML file.')
    parser.add_argument("avg_overlap", nargs="?", default="false")
    parser.add_argument("const_q", nargs="?", default="false")
    parser.add_argument("fit_first_peak", nargs="?", default="false")
    parser.add_argument("theta_offset", nargs="?", default="0")
    # Optional arguments
    parser.add_argument('--no_publish', action='store_true', help='Do not upload HTML report to the livedata server.')

    return parser.parse_args()


def autoreduce(events_file: str, output_dir: str, template_file: str = None, avg_overlap: bool = False, const_q: bool = False, theta_offset: float = 0.0, publish: bool = True) -> None:
    """
    Autoreduce a single events file and upload the HTML report to the livedata server

    Parameters
    ----------
    events_file
        Path to the Nexus events file.
    output_dir
        Output directory path.
    template_file
        Template XML file.
    avg_overlap
        If True, use average overlap method when stitching together runs.
    const_q
        If True, use constant Q summing.
    theta_offset
        Theta offset.
    publish
        If True, upload the report to the livedata server.
    """
    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load workspace and sample logs
    ws = LoadEventNexus(Filename=events_file)
    sample_logs = SampleLogValues(ws)

    # Determine which template to use
    if template_file is None:
        template_file =  get_default_template_file(output_dir, sample_logs["tthd"])
    logger.notice(f"Using template: {template_file}")

    # Run the reduction
    workflow.reduce(ws, template_file, output_dir,
                    average_overlap=avg_overlap, theta_offset=theta_offset,
                    q_summing=const_q, bck_in_q=False)

    # Plot and publish results
    upload_report(output_dir, sample_logs, ws, publish)

    # Confirm data availability
    confirm_data_availability(sample_logs)


def upload_report(output_dir: str, sample_logs: SampleLogValues, ws: MantidWorkspace, publish: bool = True) -> None:
    """
    Creates and uploads HTML report to the livedata server.

    Parameters
    ----------
    output_dir
        Output directory path.
    sample_logs
        SampleLogValues object.
    ws
        Workspace object.
    publish
        If True, upload the report to the livedata server.
    """
    sequence_id = int(sample_logs["sequence_id"])
    sequence_number = int(sample_logs["sequence_number"])

    default_file_name = 'REFL_%s_combined_data_auto.txt' % sequence_id
    default_file_path = os.path.join(output_dir, default_file_name)
    if not os.path.isfile(default_file_path):
        raise ValueError("Combined data output file not found")

    logger.notice("Loading %s" % os.path.join(output_dir, default_file_name))
    x, y, dy, dx = np.loadtxt(os.path.join(output_dir, default_file_name)).T

    run_number = int(ws.getRunNumber())
    offset = sequence_id - run_number + sequence_number - 1

    multiplot = []
    run_position = int(run_number) - sequence_id
    logger.notice(f'run position: {run_position}')
    if run_position < 10:
        _run = int(run_number)
        for i in range(0, run_position + 1):
            _id = i + offset
            _run = sequence_id + i
            reduced_file_name = 'REFL_%s_%s_%s_partial.txt' % (sequence_id, _id + 1, _run)
            reduced_file_path = os.path.join(output_dir, reduced_file_name)
            if not os.path.isfile(reduced_file_path):
                logger.notice(f"File {reduced_file_name} not found, skipping run in sequence")
                continue
            try:
                xi, yi, dyi, dxi = np.loadtxt(reduced_file_path).T
            except Exception as e:  # noqa: BLE001
                logger.error(f"Failed to load {reduced_file_name}: {e}")
                continue
            multiplot.append([xi, yi, dyi, dxi])

        plot1d(_run, multiplot, instrument='REF_L',
               x_title=u"Q (1/A)", x_log=True,
               y_title="Reflectivity", y_log=True, show_dx=False, publish=publish)
    else:
        plot1d(run_number, [[x, y, dy, dx]], instrument='REF_L',
               x_title=u"Q (1/A)", x_log=True,
               y_title="Reflectivity", y_log=True, show_dx=False, publish=publish)


def confirm_data_availability(sample_logs: SampleLogValues) -> None:
    """Notify the external data confirmation utility that data are available.

    Raises subprocess exceptions on failure so the caller can handle/log them.
    """
    try:
        ipts = sample_logs["experiment_identifier"]
        ipts_number = ipts.split("-")[1]
        cmd = ["/SNS/software/nses/bin/confirm-data", "-s", "Yes", "BL-4B", ipts_number, "1",
               "Auto"]
        subprocess.run(cmd, check=True, timeout=30)
    except Exception:  # noqa: BLE001  # deliberately broad
        logger.notice("Could not set data availability")


def str_to_bool(s):
    return s.lower() == "true"


if __name__ == "__main__":
    args = parse_command_arguments()

    events_file_arg = args.events_file
    output_dir_arg = args.output_dir

    template_file_arg = args.template_file
    avg_overlap_arg = str_to_bool(args.avg_overlap)
    const_q_arg = str_to_bool(args.const_q)
    theta_offset_arg = float(args.theta_offset)

    publish_arg = True
    if args.no_publish:
        publish_arg = False

    autoreduce(events_file_arg, output_dir_arg, template_file_arg, avg_overlap_arg, const_q_arg, theta_offset_arg, publish_arg)
