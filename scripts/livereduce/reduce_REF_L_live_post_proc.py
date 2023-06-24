import sys
import os
import json

import mantid.simpleapi as mtd_api

import numpy as np
import time

sys.path.append("/SNS/REF_L/shared/reduction")
from lr_reduction import workflow


DEBUG = True
if DEBUG:
    logfile = open("/SNS/REF_L/shared/livereduce/LR_live_outer.log", 'a')
    logfile.write("\nStarting post-proc: %s\n" % time.ctime())

def logthis(msg):
    if DEBUG:
        logfile.write(msg)

plotting_ready = True
LIVE_DATA_WS = 'accumulation'

try:
    from finddata.publish_plot import plot1d, publish_plot
except ImportError:
    plotting_ready = False


def reduction():
    """
        Perform reduction on live data
    """
    ws = mtd_api.mtd[LIVE_DATA_WS]
    run_number = ws.getRunNumber()
    ws.getRun().integrateProtonCharge()

    # Find the template to use
    expt = ws.getRun().getProperty('experiment_identifier').value
    if len(expt)>0:
        output_dir = '/SNS/REF_L/%s/shared/autoreduce' % expt
        logthis("IPTS %s [%s]\n" % (expt, output_dir))
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
        else:
            logthis("No template found\n")
            return ''

        first_run_of_set = workflow.reduce(ws, template_file,
                                           output_dir, pre_cut=1, post_cut=1,
                                           average_overlap=False, q_summing=False,
                                           bck_in_q=False, is_live=True)

        reduced_data = os.path.join(output_dir, 'REFL_%s_live_estimate.txt' % first_run_of_set)
        r = _data = np.loadtxt(reduced_data).T
        if plotting_ready:
            plot_div = plot1d(run_number, [[r[0], r[1], r[2], r[3]]], instrument='REF_L', 
                                      x_title=u"Q (1/A)", x_log=True,
                                      y_title="Reflectivity", y_log=True, show_dx=False, publish=False)
            return plot_div
    else:
        logthis("No experiment ID\n")
    return ''


def time_resolved():
    logthis("\nStarting time-resolved processing\n")
    ws = mtd_api.mtd[LIVE_DATA_WS]
    run_number = ws.getRunNumber()
    ws.getRun().integrateProtonCharge()
    charge = ws.getRun().getProtonCharge()

    plot_data = []
    data_names = []
    ws = mtd_api.SumSpectra(LIVE_DATA_WS)
    tof = mtd_api.Rebin(ws, [ws.getTofMin(), 300, ws.getTofMax()], OutputWorkspace='tof_')
    x = tof.readX(0)
    x = (x[1:]+x[:-1])/2.0
    y = tof.readY(0)

    time_data = get_live_data(run_number)
    # If we changed run, we should not use the previous data
    if len(time_data) == 0:
        logthis("New run: clearing previous data\n")
        if 'previous_data' in mtd_api.mtd:
            mtd_api.DeleteWorkspace("previous_data")

    if "previous_data" in mtd_api.mtd and len(time_data)>0:
        _previous_data = mtd_api.mtd["previous_data"]
        _previous_charge = _previous_data.getRun().getProtonCharge()

        tof_previous_data = mtd_api.Rebin(_previous_data,
                                          [ws.getTofMin(), 300, ws.getTofMax()],
                                          OutputWorkspace='tof_previous_data_')

        x_prev = tof_previous_data.readX(0)
        x_prev = (x_prev[1:]+x_prev[:-1])/2.0
        y_prev = tof_previous_data.readY(0)
        _nevts = int(np.sum(y-y_prev))
        signal = (y-y_prev) / (charge - _previous_charge)
        plot_data.append([x_prev, signal])
        data_names.append('Last 30s [%s events]' % _nevts)

        # Append to time data
        time_data.append([time.time(), [list(x), list(signal)]])
        save_live_data(run_number, time_data)

        # A minute ago
        if len(time_data) > 1:
            plot_data.append(time_data[-2][1])
            data_names.append('Previous 30s')

        # Five minutes
        if len(time_data) > 11:
            plot_data.append(time_data[-11][1])
            d_time = int(time.time() - time_data[-11][0])
            data_names.append('5 minutes ago')

        # Very first 30s that is only for this run
        if len(time_data) > 1:
            plot_data.append(time_data[1][1])
            first_time = int(time.time() - time_data[1][0])
            data_names.append('First 30s [%gs ago]' % first_time)

    else:
        time_data.append([time.time(), [list(x), list(y)]])
        save_live_data(run_number, time_data)
        plot_data.append([x, y])
        data_names.append('Last 30s')

    previous_data = mtd_api.CloneWorkspace(ws)
    logthis("plotting...\n")
    plot_div = plot1d(run_number, plot_data, data_names=data_names, instrument='REF_L',
                      x_title="TOF", x_log=True, title=time.ctime(),
                      y_title="Counts / charge", y_log=True, show_dx=False, publish=False)
    return plot_div


def get_live_data(run_number):
    """
        Load stored live data or return an empty set
    """
    live_data_path = '/SNS/REF_L/shared/livereduce/live_data.json'
    if os.path.isfile(live_data_path):
        with open(live_data_path, 'r') as fd:
            live_data = json.load(fd)
            if 'run_number' in live_data and run_number == live_data['run_number']:
                logthis("Found stored live data: %g\n" % len(live_data['time_data']))
                return live_data['time_data']
    return []


def save_live_data(run_number, time_data):
    """
        Save time data
    """
    live_data = dict(run_number=run_number, time_data=time_data)
    live_data_path = '/SNS/REF_L/shared/livereduce/live_data.json'
    with open(live_data_path, 'w') as fp:
        json.dump(live_data, fp)


html_div = ''
if LIVE_DATA_WS in mtd_api.mtd:
    try:
        ws = mtd_api.mtd[LIVE_DATA_WS]
        run_number = ws.getRunNumber()
        n_events = ws.getNumberEvents()
        ws.getRun().integrateProtonCharge()
        charge = ws.getRun().getProtonCharge()

        logthis("\nRun %s    Events: %g [charge=%s]\n" % (run_number, n_events, charge))

        # Call the reduction
        reduction_div = reduction()

        # Time-resolved plot
        plot_div = time_resolved()

        if plotting_ready:
            html_div += reduction_div
            html_div += ''
            html_div += plot_div
            publish_plot('REF_L', run_number, files={'file': html_div},
                         config="/SNS/REF_L/shared/.livedata.conf")
    except:
        logthis("failure: %s" % sys.exc_info()[1])

    # Creating a 'result' workspace is necessary for the service not to crash
    mtd_api.SumSpectra(LIVE_DATA_WS, OutputWorkspace='result')
else:
    logthis("No live data available")
