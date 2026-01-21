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
import re
import subprocess
from pathlib import Path

from mantid import logger
from mantid.simpleapi import LoadEventNexus

from lr_reduction import workflow
from lr_reduction.data_info import DataType
from lr_reduction.mantid_utils import SampleLogValues
from lr_reduction.template import get_default_template_file
from lr_reduction.web_report import assemble_report, generate_report_sections, save_report, upload_report

# Name of the conda environment to use - required by autoreduction
CONDA_ENV = "lr_reduction"


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
    parser = argparse.ArgumentParser(description="Autoreduction script for REF_L")
    # Mandatory positional arguments
    parser.add_argument("events_file", type=str, help="Path to the Nexus events file.")
    parser.add_argument("output_dir", type=str, help="Output directory path.")
    # Existing behavior: optional 3rd positional arg
    parser.add_argument("old_version_flag", nargs="?", default=None)
    # Existing behavior: optional positional 4-8 parameters
    parser.add_argument("template_file", nargs="?", default=None, type=str, help="Path to the template XML file.")
    parser.add_argument("avg_overlap", nargs="?", default="false")
    parser.add_argument("const_q", nargs="?", default="false")
    parser.add_argument("fit_first_peak", nargs="?", default="false")
    parser.add_argument("theta_offset", nargs="?", default="0")
    # Optional arguments
    parser.add_argument("--no_publish", action="store_true", help="Do not upload HTML report to the livedata server.")

    return parser.parse_args()


def autoreduce(
    events_file: str,
    output_dir: str,
    template_file: str | None = None,
    avg_overlap: bool = False,
    const_q: bool = False,
    theta_offset: float = 0.0,
    publish: bool = True,
) -> None:
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
    data_type = DataType.from_workspace(ws)

    # Determine which template to use
    if template_file is None:
        template_file = get_default_template_file(output_dir, sample_logs["tthd"])
    logger.notice(f"Using template: {template_file}")

    if data_type == DataType.REFLECTED_BEAM:
        # Run the reduction
        _, report = workflow.reduce(ws, template_file, output_dir,
                        average_overlap=avg_overlap, theta_offset=theta_offset,
                        q_summing=const_q, bck_in_q=False, return_report=True)
    elif data_type == DataType.DIRECT_BEAM:
        # Generate simple report
        report_sections = generate_report_sections(ws, template_file)
        report = assemble_report(None, report_sections)
    elif data_type == DataType.UNKNOWN:
        logger.notice(f"Data type {data_type} not supported for autoreduction.")
        return
    else:
        raise ValueError(f"Unhandled data type: {data_type}")

    # Save to disk and (optionally) upload the HTML report
    match = re.search(r'REF_L_(\d+)', events_file)
    if not match:
        raise ValueError(f"Could not extract run number from events file: {events_file}")
    run_number = match.group(1)
    save_report(report, os.path.join(output_dir, f'REF_L_{run_number}.html'))
    if publish:
        upload_report(report, run_number=run_number)

    # Confirm data availability
    confirm_data_availability(sample_logs)


def confirm_data_availability(sample_logs: SampleLogValues) -> None:
    """Notify the external data confirmation utility that data are available.

    Raises subprocess exceptions on failure so the caller can handle/log them.
    """
    try:
        ipts = sample_logs["experiment_identifier"]
        ipts_number = ipts.split("-")[1]
        cmd = ["/SNS/software/nses/bin/confirm-data", "-s", "Yes", "BL-4B", ipts_number, "1", "Auto"]
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

    autoreduce(
        events_file_arg, output_dir_arg, template_file_arg, avg_overlap_arg, const_q_arg, theta_offset_arg, publish_arg
    )
