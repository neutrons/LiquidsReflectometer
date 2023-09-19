"""
    Auto-reduction script for the Liquids Reflectometer
    For reference:
        Type 0: Normal sample data
        Type 1: Direct beams for scaling factors
        Type 2: Zero-attenuator direct beams
        Type 3: Data that we don't need to treat
"""
import os
import sys
import warnings

import numpy as np

warnings.simplefilter("ignore")

# New reduction code
sys.path.append("/SNS/REF_L/shared/reduction")

CONDA_ENV = "mantid"

from mantid.simpleapi import *

event_file_path = sys.argv[1]
output_dir = sys.argv[2]

template_file = None
if len(sys.argv) > 4:
    template_file = sys.argv[4]

avg_overlap = False
if len(sys.argv) > 5:
    avg_overlap = sys.argv[5].lower() == "true"

const_q = False
if len(sys.argv) > 6:
    const_q = sys.argv[6].lower() == "true"

event_file = os.path.split(event_file_path)[-1]
# The legacy format is REF_L_xyz_event.nxs
# The new format is REF_L_xyz.nxs.h5
run_number = event_file.split("_")[2]
run_number = int(run_number.replace(".nxs.h5", ""))

# The new reduction will be used by default starting in June 2023
old_version = False
if len(sys.argv) > 3 and sys.argv[3] == "old":
    old_version = True

# Load data for auto-reduction
ws = LoadEventNexus(Filename=event_file_path)

# Locate the template file
if template_file is None:
    default_template_path = os.path.join(output_dir, "template.xml")

    # Read tthd to determine the geometry
    tthd = ws.getRun().getProperty("tthd").value[0]
    if tthd > 0:
        template_path = os.path.join(output_dir, "template_up.xml")
    else:
        template_path = os.path.join(output_dir, "template_down.xml")

    if os.path.isfile(template_path):
        template_file = template_path
    elif os.path.isfile(default_template_path):
        template_file = default_template_path

print("Using template: %s" % template_file)

# Check the measurement geometry. This will be useful later
# if ws.getRun().getProperty('BL4B:CS:ExpPl:OperatingMode').value[0] == 'Free Liquid':

# Determine whether we have to go through the legacy reduction to
# process direct beams
if not old_version and ws.getRun().hasProperty("data_type"):
    data_type = ws.getRun().getProperty("data_type").value[0]
    old_version = not data_type == 0

if old_version:
    SetInstrumentParameter(ws, ParameterName="dq-constant", Value="0.0", ParameterType="Number")
    output = LRAutoReduction(  # Filename=event_file_path,
        InputWorkspace=ws,
        ScaleToUnity=False,
        ScalingWavelengthCutoff=10,
        OutputDirectory=output_dir,
        SlitTolerance=0.07,
        ReadSequenceFromFile=True,
        OrderDirectBeamsByRunNumber=True,
        TemplateFile=template_file,
    )
    first_run_of_set = int(output[1])
else:
    print("Average overlap: %s" % avg_overlap)
    print("Constant-Q binning: %s" % const_q)
    from lr_reduction import workflow

    if True:
        first_run_of_set = workflow.reduce(ws, template_file, output_dir, average_overlap=avg_overlap, q_summing=const_q, bck_in_q=False)
    else:
        first_run_of_set = workflow.reduce_fixed_two_theta(
            ws, template_file, output_dir, average_overlap=avg_overlap, q_summing=const_q, bck_in_q=False
        )

# -------------------------------------------------------------------------
# Produce plot for the web monitor
# Wait 30 seconds in order to avoid race condition with live reduction
# time.sleep(30)

sequence_id = int(ws.getRun().getProperty("sequence_id").value[0])
sequence_number = int(ws.getRun().getProperty("sequence_number").value[0])

default_file_name = "REFL_%s_combined_data_auto.txt" % sequence_id
default_file_path = os.path.join(output_dir, default_file_name)
if os.path.isfile(default_file_path):
    # Set flag to announce that the data is available
    try:
        ipts = ws.getRun().getProperty("experiment_identifier").value
        ipts_number = ipts.split("-")[1]
        os.system("/SNS/software/nses/bin/confirm-data -s Yes BL-4B %s 1 Auto" % ipts_number)
    except:
        logger.notice("Could not set data availability")
    print("Loading %s" % os.path.join(output_dir, default_file_name))

    x, y, dy, dx = np.loadtxt(os.path.join(output_dir, default_file_name)).T

    plotting_ready = True
    try:
        from postprocessing.publish_plot import plot1d
    except ImportError:
        from finddata.publish_plot import _determine_config_file, plot1d, publish_plot

        if _determine_config_file(None) is None:
            plotting_ready = False

    offset = sequence_id - run_number + sequence_number - 1

    if int(run_number) - sequence_id < 10:
        for i in range(0, 10):
            _id = i + offset
            _run = sequence_id + i
            reduced_file_name = "REFL_%s_%s_%s_auto.nxs" % (sequence_id, _id + 1, _run)
            reduced_file_path = os.path.join(output_dir, reduced_file_name)
            reduced_file_name2 = "REFL_%s_%s_%s_partial.txt" % (sequence_id, _id + 1, _run)
            reduced_file_path2 = os.path.join(output_dir, reduced_file_name2)
            if os.path.isfile(reduced_file_path) or os.path.isfile(reduced_file_path2):
                # Look to see whether submitting the plot is enabled
                if plotting_ready:
                    plot1d(
                        _run,
                        [[x, y, dy, dx]],
                        instrument="REF_L",
                        x_title="Q (1/A)",
                        x_log=True,
                        y_title="Reflectivity",
                        y_log=True,
                        show_dx=False,
                    )
                else:
                    plot_div = plot1d(
                        _run,
                        [[x, y, dy, dx]],
                        instrument="REF_L",
                        x_title="q (1/A)",
                        x_log=True,
                        y_title="Reflectivity",
                        y_log=True,
                        show_dx=False,
                        publish=False,
                    )
                    publish_plot("REF_L", _run, files={"file": plot_div}, config="/SNS/REF_L/shared/.livedata.conf")
    else:
        plot1d(
            run_number,
            [[x, y, dy, dx]],
            instrument="REF_L",
            x_title="Q (1/A)",
            x_log=True,
            y_title="Reflectivity",
            y_log=True,
            show_dx=False,
        )
