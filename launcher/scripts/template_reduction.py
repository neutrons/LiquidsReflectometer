import sys
import os
import numpy as np
import argparse

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


def reduce_with_mantid(ws, data_set, apply_db=False, apply_scaling_factor=False):
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
        "NormFlag": apply_db,
        "NormPeakPixelRange": data_set.NormPeakPixels,
        "SubtractNormBackground": data_set.NormBackgroundFlag,
        "NormBackgroundPixelRange": data_set.NormBackgroundRoi[:2],
        "LowResDataAxisPixelRangeFlag": data_set.data_x_range_flag,
        "LowResDataAxisPixelRange": data_set.data_x_range,
        "LowResNormAxisPixelRangeFlag": data_set.norm_x_range_flag,
        "LowResNormAxisPixelRange": data_set.norm_x_range,
        "TOFRange": data_set.DataTofRange,
        "TOFSteps": 40,
        "ApplyScalingFactor": apply_scaling_factor,
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

    if apply_scaling_factor:
        kwargs["ScalingFactorFile"] = data_set.scaling_factor_file
        kwargs["IncidentMediumSelected"] = data_set.incident_medium_list[data_set.incident_medium_index_selected]

    quartz_ws = api.LiquidsReflectometryReduction(**kwargs)
    _q_mtd = quartz_ws.readX(0)
    _r_mtd = quartz_ws.readY(0)
    _dr_mtd = quartz_ws.readE(0)

    # Q resolution
    dq = 0.0004 + 0.025*_q_mtd

    return np.asarray([_q_mtd, _r_mtd, _dr_mtd, dq])


def reduce_30Hz(meas_run_30Hz, ref_run_30Hz, ref_data_60Hz, template_30Hz, scan_index=1, template_reference=None):
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
    ref_ws_30Hz = api.LoadEventNexus("REF_L_%s.nxs.h5" % ref_run_30Hz)

    # Reduce the sample data at 30Hz
    meas_ws_30Hz = api.LoadEventNexus("REF_L_%s.nxs.h5" % meas_run_30Hz)

    # Load the 60Hz reference data
    data_60Hz = np.loadtxt(ref_data_60Hz).T

    return reduce_30Hz_from_ws(meas_ws_30Hz, ref_ws_30Hz, data_60Hz, template_data,
                               scan_index=scan_index, template_reference=template_reference)


def compute_resolution(ws, default_dq=0.03):
    """
        Compute the Q resolution from the meta data.
        @param default_dq: resolution to use if we can't compute it
    """

    # We can't compute the resolution if the value of xi is not in the logs.
    # Since it was not always logged, check for it here.
    if not ws.getRun().hasProperty("BL4B:Mot:xi.RBV"):
        print("Could not find BL4B:Mot:xi.RBV: using supplied dQ/Q")
        return default_dq

    # Xi reference would be the position of xi if the si slit were to be positioned
    # at the sample. The distance from the sample to si is then xi_reference - xi.
    xi_reference = 445
    if ws.getInstrument().hasParameter("xi-reference"):
        ws.getInstrument().getNumberParameter("xi-reference")[0]

    # Distance between the s1 and the sample
    s1_sample_distance = 1485
    if ws.getInstrument().hasParameter("s1-sample-distance"):
        ws.getInstrument().getNumberParameter("s1-sample-distance")[0]

    s1h = abs(ws.getRun().getProperty("S1VHeight").value[0])
    ths = abs(ws.getRun().getProperty("ths").value[0])
    xi = abs(ws.getRun().getProperty("BL4B:Mot:xi.RBV").value[0])
    sample_si_distance = xi_reference - xi
    slit_distance = s1_sample_distance - sample_si_distance
    dq_over_q = s1h / slit_distance * 180 / 3.1416 / ths
    return dq_over_q


def reduce_30Hz_from_ws(meas_ws_30Hz, ref_ws_30Hz, data_60Hz, template_data, scan_index=1, template_reference=None):
    """
        Perform 30Hz reduction
        @param meas_ws_30Hz: Mantid workspace of the data we want to reduce
        @param ref_ws_30Hz: Mantid workspace of the reference data, take with the same config
        @param data_60Hz: reduced reference data at 60Hz
        @param template_data: template data object (for 30Hz)
        @param scan_index: scan index to use within the template.
    """
    # Reduce the quartz at 30Hz
    if template_reference is None:
        r_ref = reduce_with_mantid(ref_ws_30Hz, template_data)
    else:
        r_ref = reduce_with_mantid(ref_ws_30Hz, template_reference)

    # Reduce the sample data at 30Hz
    r_meas = reduce_with_mantid(meas_ws_30Hz, template_data)

    # Identify the bins we need to overlap with the 30Hz measurement
    # The assumption is that the binning is the same
    _tolerance = 0.0001
    _max_q = min(r_ref[0].max(), r_meas[0].max())
    _min_q = max(r_ref[0].min(), r_meas[0].min())
    
    print("60Hz:      %g %g" % (data_60Hz[0].min(), data_60Hz[0].max()))
    print("Ref 30Hz:  %g %g" % (r_meas[0].min(), r_meas[0].max()))
    print("Meas 30Hz: %g %g" % (r_ref[0].min(), r_ref[0].max()))
    
    _q_idx_60 = np.asarray(np.where((data_60Hz[0] > _min_q-_tolerance) & (data_60Hz[0] < _max_q+_tolerance)))[0]
    _q_idx_meas30 = np.asarray(np.where((r_meas[0] > _min_q-_tolerance) & (r_meas[0] < _max_q+_tolerance)))[0]
    _q_idx_ref30 = np.asarray(np.where((r_ref[0] > _min_q-_tolerance) & (r_ref[0] < _max_q+_tolerance)))[0]

    if not data_60Hz[0][_q_idx_60].shape[0] == r_meas[0][_q_idx_ref30].shape[0]:
        print("60Hz reference may have been reduced with different binning!")

    # Confirm identical binning
    _sum = np.sum(data_60Hz[0][_q_idx_60]-r_ref[0][_q_idx_ref30])
    if _sum > r_ref[0][0]/100:
        print("Binning 60Hz and ref 30Hz not identical!")
        

    _sum = np.sum(data_60Hz[0][_q_idx_60]-r_meas[0][_q_idx_meas30])
    if _sum > r_ref[0][0]/100:
        print("Binning 60Hz and meas 30Hz not identical!")
    
    r_q_final = r_meas[1][_q_idx_meas30]/r_ref[1][_q_idx_ref30]*data_60Hz[1][_q_idx_60]

    dr_q_final = np.sqrt((r_meas[2][_q_idx_meas30]/r_ref[1][_q_idx_ref30]*data_60Hz[1][_q_idx_60])**2 \
                         +(r_meas[1][_q_idx_meas30]/r_ref[1][_q_idx_ref30]*data_60Hz[2][_q_idx_60])**2 \
                         +(r_meas[1][_q_idx_meas30]/r_ref[1][_q_idx_ref30]**2*data_60Hz[1][_q_idx_60]*r_ref[2][_q_idx_ref30])**2)

    print("Q range: %s - %s" % (r_meas[0][0], r_meas[0][_q_idx_meas30][-1]))
    q = r_meas[0][_q_idx_meas30]
    _idx = (r_q_final > 0) & (r_q_final < np.inf)

    # Q resolution
    #   Assume a constant term of 0 unless it is specified
    dq0 = 0
    if meas_ws_30Hz.getInstrument().hasParameter("dq-constant"):
        dq0 = meas_ws_30Hz.getInstrument().getNumberParameter("dq-constant")[0]
    dq_slope = compute_resolution(meas_ws_30Hz)
    print("Resolution: %g + %g Q" % (dq0, dq_slope))
    dq = dq0 + dq_slope * q[_idx]
    return np.asarray([q[_idx], r_q_final[_idx], dr_q_final[_idx], dq])


def reduce_30Hz_slices(meas_run_30Hz, ref_run_30Hz, ref_data_60Hz, template_30Hz, 
                       time_interval, output_dir, scan_index=1, create_plot=True, template_reference=None):

    meas_ws_30Hz = api.LoadEventNexus("REF_L_%s.nxs.h5" % meas_run_30Hz)

    return reduce_30Hz_slices_ws(meas_ws_30Hz, ref_run_30Hz, ref_data_60Hz, template_30Hz, 
                                 time_interval, output_dir, scan_index=scan_index, create_plot=create_plot,
                                template_reference=template_reference)

def reduce_60Hz_slices(meas_run, template_file,
                       time_interval, output_dir, scan_index=1, create_plot=True):

    meas_ws = api.LoadEventNexus("REF_L_%s.nxs.h5" % meas_run)

    return reduce_60Hz_slices_ws(meas_ws, template_file,
                                 time_interval, output_dir, scan_index=scan_index, create_plot=create_plot)

def reduce_30Hz_slices_ws(meas_ws_30Hz, ref_run_30Hz, ref_data_60Hz, template_30Hz, 
                          time_interval, output_dir, scan_index=1, create_plot=True,
                          template_reference=None):
    """
        Perform 30Hz reduction
        @param meas_ws_30Hz: workspace of the data we want to reduce
        @param ref_ws_30Hz: workspace of the reference data, take with the same config
        @param ref_data_60Hz: file path of the reduce data file at 60Hz
        @param template_30Hz: file path of the template file for 30Hz
        @param time_interval: time step in seconds
        @param scan_index: scan index to use within the template.
    """
    # Load the template
    print("Reading template")
    template_data = read_template(template_30Hz, scan_index)

    # Reduce the quartz at 30Hz
    print("Reading reference data at 30Hz")
    if os.path.isfile(ref_run_30Hz):
        ref_ws_30Hz = api.LoadEventNexus(ref_run_30Hz)
    else:
        ref_ws_30Hz = api.LoadEventNexus("REF_L_%s.nxs.h5" % ref_run_30Hz)

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
    for name in wsnames:
        tmpws = api.mtd[name]
        print("workspace %s has %d events" % (name, tmpws.getNumberEvents()))
        try:
            _reduced = reduce_30Hz_from_ws(tmpws, ref_ws_30Hz, data_60Hz, template_data,
                                           scan_index=scan_index, template_reference=template_reference)
            reduced.append(_reduced)
            _filename = 'r{0}_t{1:06d}.txt'.format(meas_run_30Hz, int(total_time))
            np.savetxt(os.path.join(output_dir, _filename), _reduced.T)
        except:
            print(sys.exc_info()[0])
        total_time += time_interval

    if create_plot:
        plot_slices(reduced, title='Duration: %g seconds' % duration,
                    time_interval=time_interval,
                    file_path=os.path.join(output_dir, 'r%s.png' % meas_run_30Hz))

    return reduced

def reduce_60Hz_slices_ws(meas_ws, template_file, 
                          time_interval, output_dir, scan_index=1, create_plot=True):
    """
        Perform 30Hz reduction
        @param meas_ws: workspace of the data we want to reduce
        @param template_file: autoreduction template file
        @param time_interval: time step in seconds
        @param scan_index: scan index to use within the template.
    """

    # Load the template
    print("Reading template")
    template_data = read_template(template_file, scan_index)

    # Some meta-data are not filled in for the live data stream
    # Use dummy values for those
    try:
        duration = meas_ws.getRun()['duration'].value
    except:
        duration = 0
    try:
        meas_run = meas_ws.getRun()['run_number'].value
    except:
        meas_run = 0

    # Time slices
    print("Slicing data")
    splitws, infows = api.GenerateEventsFilter(InputWorkspace=meas_ws, TimeInterval=time_interval)

    api.FilterEvents(InputWorkspace=meas_ws,
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

    reduced = []
    total_time = 0
    for name in wsnames:
        tmpws = api.mtd[name]
        print("workspace %s has %d events" % (name, tmpws.getNumberEvents()))
        try:
            _reduced = reduce_with_mantid(tmpws, template_data, apply_scaling_factor=True)
            reduced.append(_reduced)
            _filename = 'r{0}_t{1:06d}.txt'.format(meas_run, int(total_time))
            np.savetxt(os.path.join(output_dir, _filename), _reduced.T)
        except:
            print(sys.exc_info()[0])
        total_time += time_interval

    if create_plot:
        plot_slices(reduced, title='Duration: %g seconds' % duration,
                    time_interval=time_interval,
                    file_path=os.path.join(output_dir, 'r%s.png' % meas_run))

    return reduced

def plot_slices(reduced, title, time_interval, file_path, offset=10):
    fig, ax = plt.subplots(figsize=(6,6))

    total_time = 0
    _running_offset = 1.
    for _data in reduced:
        qz, refl, d_refl, _ = _data

        plt.errorbar(qz, refl*_running_offset, yerr=d_refl*_running_offset, markersize=4, marker='o',
                     label='T=%g s' % total_time)

        total_time += time_interval
        _running_offset *= offset

    plt.legend()
    plt.title(title)
    plt.xlabel('q [$1/\AA$]')
    plt.ylabel('R(q)')
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.show()
    plt.savefig(file_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Time-resolved at 30Hz
    dynanic30_parser = subparsers.add_parser('dynamic30Hz', help='Reduce time-resolved 30Hz [-h for help]')
    dynanic30_parser.add_argument('meas_run_30Hz', type=int,
                             help='Run number for the data to be processed')
    dynanic30_parser.add_argument('ref_run_30Hz', type=str,
                             help='Run number for the reference 30Hz data, measured at the same settings as the data to be processed')
    dynanic30_parser.add_argument('ref_data_60Hz', type=str,
                             help='Reference R(Q), measured at 60Hz')
    dynanic30_parser.add_argument('template_30Hz', type=str,
                             help='File path for the 30Hz reduction template')
    dynanic30_parser.add_argument('time_interval', type=float,
                             help='Time interval to use, in seconds')
    dynanic30_parser.add_argument('output_dir', type=str,
                             help='Output directory')
    dynanic30_parser.add_argument('--scan_index', type=int, dest='scan_index',
                                  help='Template scan index', required=False, default=1)
    dynanic30_parser.add_argument('--no-plot', dest='create_plot', action='store_false')
    dynanic30_parser.set_defaults(create_plot=True)

    # Time-resolved at 60Hz
    dynanic60_parser = subparsers.add_parser('dynamic60Hz', help='Reduce time-resolved 60Hz [-h for help]')
    dynanic60_parser.add_argument('meas_run_60Hz', type=int,
                             help='Run number for the data to be processed')
    dynanic60_parser.add_argument('template_60Hz', type=str,
                             help='File path for the 60Hz reduction template')
    dynanic60_parser.add_argument('time_interval', type=float,
                             help='Time interval to use, in seconds')
    dynanic60_parser.add_argument('output_dir', type=str,
                             help='Output directory')
    dynanic60_parser.add_argument('--scan_index', type=int, dest='scan_index',
                                  help='Template scan index', required=True, default=1)

    # Parse arguments
    args = parser.parse_args()

    if args.command=='dynamic30Hz':
        print("Time-resolved reduction at 30Hz: run %s" % args.meas_run_30Hz)
        reduced = reduce_30Hz_slices(args.meas_run_30Hz, args.ref_run_30Hz, args.ref_data_60Hz, args.template_30Hz,
                                     time_interval=args.time_interval, output_dir=args.output_dir,
                                     scan_index=args.scan_index, create_plot=args.create_plot)
    elif args.command=='dynamic60Hz':
        print("Time-resolved reduction at 60Hz: run %s" % args.meas_run_60Hz)
        reduced = reduce_60Hz_slices(args.meas_run_60Hz, args.template_60Hz,
                                     time_interval=args.time_interval, output_dir=args.output_dir, scan_index=args.scan_index)
