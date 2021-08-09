import sys
import os
import numpy as np

from matplotlib import pyplot as plt

import mantid
from mantid.api import *
import mantid.simpleapi as api
from mantid.kernel import *
from reduction_gui.reduction.reflectometer.refl_data_series import DataSeries
from reduction_gui.reduction.reflectometer.refl_data_script import DataSets


TOLERANCE = 0.02

def read_template(template_file, sequence_number):
    """
        Read template from file.
        @param sequence_number: the ID of the data set within the sequence of runs
    """
    fd = open(template_file, "r")
    xml_str = fd.read()
    s = DataSeries()
    s.from_xml(xml_str)

    if len(s.data_sets) >= sequence_number:
        data_set = s.data_sets[sequence_number - 1]
    elif len(s.data_sets) > 0:
        data_set = s.data_sets[0]
    else:
        raise RuntimeError("Invalid reduction template")

    return data_set


def reduce_with_mantid(ws, data_set):
    """
        @param ws: Mantid workspace
        @param data_set: template object
    """
    kwargs = {
        "InputWorkspace": ws,
        "NormalizationRunNumber": str(data_set.norm_file),
        "SignalPeakPixelRange": data_set.DataPeakPixels,
        "SubtractSignalBackground": data_set.DataBackgroundFlag,
        "SignalBackgroundPixelRange": data_set.DataBackgroundRoi[:2],
        "NormFlag": False,
        "LowResDataAxisPixelRangeFlag": data_set.data_x_range_flag,
        "LowResDataAxisPixelRange": data_set.data_x_range,
        "TOFRange": data_set.DataTofRange,
        "TOFSteps": 40,
        "ApplyScalingFactor": False,
        "GeometryCorrectionFlag": False,
        "QMin": data_set.q_min,
        "QStep": data_set.q_step,
        "AngleOffset": data_set.angle_offset,
        "AngleOffsetError": data_set.angle_offset_error,
        "SlitsWidthFlag": data_set.slits_width_flag,
        "ApplyPrimaryFraction": False,
        "SlitTolerance": 0.02,
        "OutputWorkspace": 'reflectivity_%s' % str(ws)
    }
    quartz_ws = api.LiquidsReflectometryReduction(**kwargs)
    _q_mtd = quartz_ws.readX(0)
    _r_mtd = quartz_ws.readY(0)
    _dr_mtd = quartz_ws.readE(0)
    return np.asarray([_q_mtd, _r_mtd, _dr_mtd])


def reduce_30Hz(meas_run_30Hz, ref_run_30Hz, ref_data_60Hz, template_30Hz, scan_index=1):
    """
        Perform 30Hz reduction
        @param meas_run_30Hz: run number of the data we want to reduce
        @param ref_run_30Hz: run number of the reference data, take with the same config
        @param ref_data_60Hz: file path of the reduce data file at 60Hz
        @param template_30Hz: file path of the template file for 30Hz
        @param scan_index: scan index to use within the template.
    """
    # Load the template
    template_data = read_template(template_30Hz, scan_index)

    # Reduce the quartz at 30Hz
    ref_ws_30Hz = api.LoadEventNexus("REF_L_%s" % ref_run_30Hz)

    # Reduce the sample data at 30Hz
    meas_ws_30Hz = api.LoadEventNexus("REF_L_%s" % meas_run_30Hz)

    # Load the 60Hz reference data
    data_60Hz = np.loadtxt(ref_data_60Hz).T

    return reduce_30Hz_from_ws(meas_ws_30Hz, ref_ws_30Hz, data_60Hz, template_data, scan_index=scan_index)


def reduce_30Hz_from_ws(meas_ws_30Hz, ref_ws_30Hz, data_60Hz, template_data, scan_index=1):
    """
        Perform 30Hz reduction
        @param meas_ws_30Hz: Mantid workspace of the data we want to reduce
        @param ref_ws_30Hz: Mantid workspace of the reference data, take with the same config
        @param data_60Hz: reduced reference data at 60Hz
        @param template_data: template data object (for 30Hz)
        @param scan_index: scan index to use within the template.
    """
    # Reduce the quartz at 30Hz
    r_ref = reduce_with_mantid(ref_ws_30Hz, template_data)

    # Reduce the sample data at 30Hz
    r_meas = reduce_with_mantid(meas_ws_30Hz, template_data)

    # Identify the bins we need to overlap with the 30Hz measurement
    # The assumption is that the binning is the same
    _q_idx_ref = np.asarray(np.where((data_60Hz[0] > r_ref[0][0]-0.1*(r_ref[0][1]-r_ref[0][0])) \
                                     & (data_60Hz[0] < r_ref[0][-1]+0.1*(r_ref[0][-1]-r_ref[0][-2]))))[0]
    _q_idx_meas = np.asarray(np.where(r_ref[0] <= data_60Hz[0][-1]))[0]
    
    # Confirm identical binning
    _sum = np.sum(data_60Hz[0][_q_idx_ref]-r_ref[0][_q_idx_meas])
    if _sum > r_ref[0][0]/100:
        print("Binning not identical!")

    r_q_final = r_meas[1][_q_idx_meas]/r_ref[1][_q_idx_meas]*data_60Hz[1][_q_idx_ref]

    dr_q_final = np.sqrt((r_meas[2][_q_idx_meas]/r_ref[1][_q_idx_meas]*data_60Hz[1][_q_idx_ref])**2 \
                         +(r_meas[1][_q_idx_meas]/r_ref[1][_q_idx_meas]*data_60Hz[2][_q_idx_ref])**2 \
                         +(r_meas[1][_q_idx_meas]/r_ref[1][_q_idx_meas]**2*data_60Hz[1][_q_idx_ref]*r_ref[2][_q_idx_meas])**2)

    print("Q range: %s - %s" % (r_meas[0][0], r_meas[0][_q_idx_meas][-1]))
    return np.asarray([r_meas[0][_q_idx_meas], r_q_final, dr_q_final])


def reduce_30Hz_slices(meas_run_30Hz, ref_run_30Hz, ref_data_60Hz, template_30Hz, 
                       time_interval, output_dir, scan_index=1, create_plot=True):

    meas_ws_30Hz = api.LoadEventNexus("REF_L_%s" % meas_run_30Hz)

    return reduce_30Hz_slices_ws(meas_ws_30Hz, ref_run_30Hz, ref_data_60Hz, template_30Hz, 
                       time_interval, output_dir, scan_index=scan_index, create_plot=create_plot)

def reduce_30Hz_slices_ws(meas_ws_30Hz, ref_run_30Hz, ref_data_60Hz, template_30Hz, 
                       time_interval, output_dir, scan_index=1, create_plot=True):
    """
        Perform 30Hz reduction
        @param meas_run_30Hz: run number of the data we want to reduce
        @param ref_run_30Hz: run number of the reference data, take with the same config
        @param ref_data_60Hz: file path of the reduce data file at 60Hz
        @param template_30Hz: file path of the template file for 30Hz
        @param scan_index: scan index to use within the template.
    """
    # Load the template
    print("Reading template")
    template_data = read_template(template_30Hz, scan_index)

    # Reduce the quartz at 30Hz
    print("Reading reference data at 30Hz")
    ref_ws_30Hz = api.LoadEventNexus("REF_L_%s" % ref_run_30Hz)

    # Reduce the sample data at 30Hz
    print("Reading sample data at 30Hz")
    #meas_ws_30Hz = api.LoadEventNexus("REF_L_%s" % meas_run_30Hz)

    # Some meta-data are not filled in for the live data stream
    # Use dummy values for those
    try:
        duration = meas_ws_30Hz.getRun()['duration'].value
    except:
        duration = 0
    try:
        meas_run_30Hz = meas_ws_30Hz.getRun()['run_number'].value
    except:
        meas_run_30Hz = 0

    # Time slices
    print("Slicing data")
    splitws, infows = api.GenerateEventsFilter(InputWorkspace=meas_ws_30Hz, TimeInterval=time_interval)

    api.FilterEvents(InputWorkspace=meas_ws_30Hz,
        SplitterWorkspace=splitws,
        InformationWorkspace=infows,
        OutputWorkspaceBaseName='time_ws',
        GroupWorkspaces=True,
        FilterByPulseTime = True,
        OutputWorkspaceIndexedFrom1 = True,
        CorrectionToSample = "None",
        SpectrumWithoutDetector = "Skip",
        SplitSampleLogs = False,
        OutputTOFCorrectionWorkspace='mock')
    wsgroup = api.mtd["time_ws"]
    wsnames = wsgroup.getNames()

    # Load the 60Hz reference data
    print("Loading reference R(Q)")
    data_60Hz = np.loadtxt(ref_data_60Hz).T

    reduced = []
    total_time = 0
    for name in sorted(wsnames):
        tmpws = api.mtd[name]
        print("workspace %s has %d events" % (name, tmpws.getNumberEvents()))
        try:
            _reduced = reduce_30Hz_from_ws(tmpws, ref_ws_30Hz, data_60Hz, template_data, scan_index=scan_index)
            reduced.append(_reduced)
            np.savetxt(os.path.join(output_dir, 'r%s_t%g.txt' % (meas_run_30Hz, total_time)), _reduced.T)
        except:
            print(sys.exc_info()[0])
        total_time += time_interval

    if create_plot:
        plot_slices(reduced, title='Duration: %g seconds' % duration,
                    time_interval=time_interval,
                    file_path=os.path.join(output_dir, 'r%s.png' % meas_run_30Hz))
        
    return reduced

def plot_slices(reduced, title, time_interval, file_path):
    fig, ax = plt.subplots(figsize=(6,6))

    total_time = 0

    for _data in reduced:
        qz, refl, d_refl = _data

        plt.errorbar(qz, refl, yerr=d_refl, markersize=4, marker='o',
                     label='T=%g s' % total_time)

        total_time += time_interval

    plt.legend()
    plt.title(title)
    plt.xlabel('q [$1/\AA$]')
    plt.ylabel('R(q)')
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.show()
    plt.savefig(file_path)

if __name__ == "__main__":
    if len(sys.argv) < 7:
        print("\nUsage: python3 template_reduction.py <meas_run_30Hz> <ref_run_30Hz> <ref_data_60Hz> <template_30Hz> <time_interval> <output_dir>")
        print("%s" % str(sys.argv))
        sys.exit(0)

    meas_30Hz = sys.argv[1]
    ref_run_30Hz = sys.argv[2]
    ref_data_60Hz = sys.argv[3]
    template_30Hz = sys.argv[4]
    time_interval = sys.argv[5]
    output_dir = sys.argv[6]

    print("Measured run to reduce: %s" % meas_30Hz)
    reduced = reduce_30Hz_slices(meas_30Hz, ref_run_30Hz, ref_data_60Hz, template_30Hz,
                                 time_interval=float(time_interval), output_dir=output_dir, scan_index=1)

    print("DONE")