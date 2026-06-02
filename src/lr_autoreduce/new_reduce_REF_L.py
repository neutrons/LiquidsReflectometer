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

#from lr_reduction import workflow
from lr_reduction.data_info import DataType
from lr_reduction.mantid_utils import SampleLogValues
from lr_reduction.template import get_default_template_file
from lr_reduction.web_report import assemble_report, generate_report_sections, save_report, upload_report
import lr_reduction.new_reduction_from_file as nrff
from plot_publisher import plot1d

# Name of the conda environment to use - required by autoreduction
#CONDA_ENV = "lr_reduction"


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
    parser.add_argument("settings_file",nargs="?", default=None, type=str, help="Path to the settings JSON file.")

    return parser.parse_args()

def get_default_setting_file(output_dir: str, tthd: float) -> str:
    """
    Determine the default setting file based on the tthd value.
    ADAPTED FROM TEMPLATE EQUIVALENT.

    @param output_dir: Directory with setting files
    @param tthd: Two-theta detector value
    @return: Path to the setting file
    """
    output_dir = Path(output_dir)
    setting_file = None

    # Read tthd to determine the geometry
    if tthd > 0:
        setting_path = output_dir / "reduce_settings_up.json"
    else:
        setting_path = output_dir / "reduce_settings_down.json"

    default_setting_path = output_dir / "reduce_settings.json"

    if setting_path.exists():
        setting_file = str(setting_path)
    elif default_setting_path.exists():
        setting_file = str(default_setting_path)

    if setting_file is None:
        raise ValueError("No settings file found: place a settings file in shared/autoreduce")

    return setting_file

def autoreduce_new(
	events_file,
	output_dir,
	setting_file = None,
    template_file = None,
	publish = True):

    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # KEEP MANTID LOGIC HERE FOR NOW
    # Load workspace and sample logs
    ws = LoadEventNexus(Filename=events_file)
    sample_logs = SampleLogValues(ws)
    data_type = DataType.from_workspace(ws)
    
    experiment_id = sample_logs["experiment_identifier"]

    match = re.search(r'REF_L_(\d+)', events_file)
    if not match:
        raise ValueError(f"Could not extract run number from events file: {events_file}")
    run_number = match.group(1)
    #run_number = sample_logs["entry_identifier"]

    # Determine which setting file to use
    if setting_file is None:
        setting_file = get_default_setting_file(output_dir, sample_logs["tthd"])
    logger.notice(f"Using settings from: {setting_file}")

    # TODO: separate from template if can?
    # Determine which template to use
    if template_file is None:
        template_file = get_default_template_file(output_dir, sample_logs["tthd"])
    logger.notice(f"Using template: {template_file}")

    if data_type == DataType.REFLECTED_BEAM:
        # Run the reduction
        subname = "autoreduction"
        savepath = Path(output_dir) / "new_reduction"
        # Ensure output directory exists
        Path(savepath).mkdir(parents=True, exist_ok=True)
        override_params = {"Spath": savepath, "subname": subname}
        output, plots, run_number_list, config_out = nrff.reduce_from_file([run_number], setting_file, experiment_id, plot=False, 
                override_params = override_params, check_for_prior=True, save_pdf_summary=False)
        print('config out', vars(config_out))
        # Now have the config so can feed that in later...
        plot_out = generate_ref_plot(output[0], run_number_list)
        # TODO: Add the meta_data part
        # TODO: here for the config input...
        report_sections = generate_report_sections(ws,template_file, meta_data = None, config_in=config_out)
        report = assemble_report(plot_out, report_sections)
    elif data_type == DataType.DIRECT_BEAM:
        # TODO: Fix the DB part...

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
    save_report(report, os.path.join(output_dir, f'REF_L_{run_number}_new.html'))
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


def generate_ref_plot(output, run_nums):
        """
        Plot the combined reflectivity curve for a collection of runs

        Returns
        -------
        str
            HTML div containing the combined reflectivity curve plot
        """
        refl_curves = []
        run_names = []

        for i in range(len(output)):
            refl_curves.append([output[i]["Q"], output[i]["R"], output[i]["dR"], output[i]["dQ"]])
            # TODO: get the SF pulled through.
            #run_names.append(f"Run: {run_num}  SF: {output[i]['config'].ScaleFactor[i]:.3f}")
            run_names.append(f"Run: {run_nums[i]}")

        # run_number parameter is only used when publish=True
        return plot1d(run_number="dummy_run", data_list=refl_curves, data_names=run_names,
                      instrument='REF_L',
                      x_title=u"Q (1/A)", x_log=True,
                      y_title="Reflectivity", y_log=True, show_dx=False, publish=False)

if __name__ == "__main__":
    args = parse_command_arguments()

    events_file_arg = args.events_file
    output_dir_arg = args.output_dir

    settings_file_arg = args.settings_file
    template_file_arg = args.template_file

    publish_arg = True
    if args.no_publish:
        publish_arg = False

    print("Running new reduction")
    autoreduce_new(
        events_file_arg, output_dir_arg, settings_file_arg,template_file_arg, publish_arg 
    )


