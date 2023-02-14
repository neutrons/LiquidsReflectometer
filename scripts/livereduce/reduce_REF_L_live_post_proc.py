import sys
import os
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
    from finddata.publish_plot import plot1d, _determine_config_file, publish_plot
except ImportError:
    plotting_ready = False

logthis(str(mtd_api.mtd.getObjectNames()))

html_div = ''
if False and LIVE_DATA_WS in mtd_api.mtd:
    try:
        ws = mtd_api.mtd[LIVE_DATA_WS]
        n_events = ws.getNumberEvents()
        run_number = ws.getRunNumber()
        ws.getRun().integrateProtonCharge()
        charge = ws.getRun().getProtonCharge()
        #charge = np.sum(charge)
        logthis("\nEvents: %g [charge=%s]\n" % (n_events, charge))

        expt = ws.getRun().getProperty('experiment_identifier').value
        if len(expt)>0:
            output_dir = '/SNS/REF_L/%s/shared/autoreduce' % expt
            logthis("IPTS %s [%s]" % (expt, output_dir))
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

            first_run_of_set = workflow.reduce(ws, template_file,
                                           output_dir, pre_cut=1, post_cut=1,
                                           average_overlap=False,
                                           q_summing=False, bck_in_q=False, is_live=True)

            reduced_data = os.path.join(output_dir, 'REFL_%s_combined_data_auto.txt' % first_run_of_set)
            #reduced_data = os.path.join(output_dir, 'REFL_%s_live_estimate.txt' % first_run_of_set)
            r = _data = np.loadtxt(reduced_data).T
            if plotting_ready:
                plot_div = plot1d(run_number, [[r[0], r[1], r[2], r[3]]], instrument='REF_L', 
                                          x_title=u"Q (1/A)", x_log=True,
                                          y_title="Reflectivity", y_log=True, show_dx=False, publish=False)
                html_div += plot_div
                publish_plot('REF_L', run_number, files={'file': html_div},
                             config="/SNS/REF_L/shared/.livedata.conf")

        else:
            logthis("No experiment ID")

        ws = mtd_api.SumSpectra(LIVE_DATA_WS, OutputWorkspace='result')

        if False:
            tof = mtd_api.Rebin(ws, [ws.getTofMin(), 300, ws.getTofMax()], OutputWorkspace='tof_')
            x = tof.readX(0)
            logthis(str(x))
            x = (x[1:]+x[:-1])/2.0
            y = tof.readY(0)

            previous_data = mtd_api.CloneWorkspace(ws)
    
            if plotting_ready:
                plot_div = plot1d(run_number, [[x, y, 0*x, 0*y]], instrument='REF_L',
                                          x_title=u"TOF", x_log=True,
                                          y_title="Counts", y_log=True, show_dx=False, publish=False)
                html_div += plot_div
                publish_plot('REF_L', run_number, files={'file': html_div},
                             config="/SNS/REF_L/shared/.livedata.conf")

    except:
        logthis("failure: %s" % sys.exc_info()[1])

else:
    logthis("No live data available")
