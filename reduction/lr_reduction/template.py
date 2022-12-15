import sys
import os
import numpy as np

from matplotlib import pyplot as plt

import mantid
from mantid.api import *
import mantid.simpleapi as api
from mantid.kernel import *

from functools import reduce 

from . import event_reduction
from . import reduction_template_reader

TOLERANCE = 0.02


def read_template(template_file, sequence_number):
    """
        Read template from file.
        @param sequence_number: the ID of the data set within the sequence of runs
    """
    fd = open(template_file, "r")
    xml_str = fd.read()
    data_sets = reduction_template_reader.from_xml(xml_str)
    if len(data_sets) >= sequence_number:
        data_set = data_sets[sequence_number - 1]
    elif len(data_sets) > 0:
        data_set = data_sets[0]
    else:
        raise RuntimeError("Invalid reduction template")
    return data_set


def scaling_factor(scaling_factor_file, workspace, match_slit_width=True):
    """
        Apply scaling factor from reference scaling data
        @param workspace: Mantid workspace
    """
    if not os.path.isfile(scaling_factor_file):
        print("Could not find scaling factor file: %s" % scaling_factor_file)
        return workspace

    # Get the wavelength
    lr = workspace.getRun().getProperty('LambdaRequest').value[0]
    lr_value = float("{0:.2f}".format(lr))

    s1h = abs(workspace.getRun().getProperty("S1VHeight").value[0])
    s1w = abs(workspace.getRun().getProperty("S1HWidth").value[0])
    s2h = abs(workspace.getRun().getProperty("SiVHeight").value[0])
    s2w = abs(workspace.getRun().getProperty("SiHWidth").value[0])

    def _reduce(accumulation, item):
        """
            Reduce function that accumulates values in a dictionary
        """
        toks_item = item.split('=')
        if len(toks_item)!=2:
            return accumulation
        if isinstance(accumulation, dict):
            accumulation[toks_item[0].strip()] = toks_item[1].strip()
        else:
            toks_accum = accumulation.split('=')
            accumulation = {toks_item[0].strip(): toks_item[1].strip(),
                            toks_accum[0].strip(): toks_accum[1].strip()}
        return accumulation

    def _value_check(key, data, reference):
        """
            Check an entry against a reference value
        """
        if key in data:
            return abs(abs(float(data[key])) - abs(float(reference))) <= TOLERANCE
        return False

    with open(scaling_factor_file, 'r') as fd:
        file_content = fd.read()

    data_found = None
    for line in file_content.split('\n'):
        if line.startswith('#'):
            continue

        # Parse the line of data and produce a dict
        toks = line.split()
        data_dict = reduce(_reduce, toks, {})

        # Get ordered list of keys
        keys = []
        for token in toks:
            key_value = token.split('=')
            if len(key_value)==2:
                keys.append(key_value[0].strip())

        # Skip empty lines
        if len(keys)==0:
            continue
        # Complain if the format is non-standard
        elif len(keys)<10:
            print("Bad scaling factor entry\n  %s" % line)
            continue

        # Sanity check
        if keys[0] != 'IncidentMedium' and keys[1] != 'LambdaRequested' \
                and keys[2] != 'S1H':
            print("The scaling factor file isn't standard: bad keywords")
        # The S2H key has been changing in the earlier version of REFL reduction.
        # Get the key from the data to make sure we are backward compatible.
        s2h_key = keys[3]
        s2w_key = keys[5]
        if 'IncidentMedium' in data_dict \
                and _value_check('LambdaRequested', data_dict, lr_value) \
                and _value_check('S1H', data_dict, s1h) \
                and _value_check(s2h_key, data_dict, s2h):

            if not match_slit_width or (_value_check('S1W', data_dict, s1w)
                                        and _value_check(s2w_key, data_dict, s2w)):
                data_found = data_dict
                break

    if data_found is not None:
        a = float(data_found['a'])
        b = float(data_found['b'])
        a_error = float(data_found['error_a'])
        b_error = float(data_found['error_b'])
    else:
        return 1, 0, 0, 0
    return a, b, a_error, b_error


def process_from_template(run_number, template_path, q_summing=False,
                          tof_weighted=False, bck_in_q=False, clean=False):
    # Load data
    ws_sc = api.Load("REF_L_%s" % run_number)
    return process_from_template_ws(ws_sc, template_path, q_summing=q_summing,
                                    tof_weighted=tof_weighted, bck_in_q=bck_in_q,
                                    clean=clean)

def process_from_template_ws(ws_sc, template_path, q_summing=False,
                             tof_weighted=False, bck_in_q=False, clean=False):
    # Get the sequence number
    sequence_number = 1
    if ws_sc.getRun().hasProperty("sequence_number"):
        sequence_number = ws_sc.getRun().getProperty("sequence_number").value[0]

    # Load the template
    template_data = read_template(template_path, sequence_number)

    # Load normalization run
    ws_db = api.LoadEventNexus("REF_L_%s" % template_data.norm_file)

    # Get the angle offset
    offset = template_data.angle_offset

    thi_value = ws_sc.getRun()['thi'].value[0]
    theta = np.fabs(ws_sc.getRun()['ths'].value[0]) + offset
    _wl = ws_sc.getRun()['LambdaRequest'].value[0]
    print('wl=%g; ths=%g; offset=%g' % (_wl, theta, offset))

    theta = theta * np.pi / 180.

    # Get the reduction parameters from the template
    peak = template_data.data_peak_range
    peak_bck = [template_data.background_roi[0], template_data.background_roi[1]]
    peak_center = (peak[0]+peak[1])/2.0
    low_res = template_data.data_x_range

    norm_peak = template_data.norm_peak_range
    norm_low_res = template_data.norm_x_range
    norm_bck = [template_data.norm_background_roi[0], template_data.norm_background_roi[1]]

    [tof_min, tof_max] = template_data.tof_range
    q_min = template_data.q_min
    q_step = -template_data.q_step
    
    # Perform the reduction
    event_refl = event_reduction.EventReflectivity(ws_sc, ws_db,
                                                   signal_peak=peak, signal_bck=peak_bck,
                                                   norm_peak=norm_peak, norm_bck=norm_bck,
                                                   specular_pixel=peak_center,
                                                   signal_low_res=low_res, norm_low_res=norm_low_res,
                                                   q_min=q_min, q_step=q_step, q_max=None,
                                                   tof_range=[tof_min, tof_max],
                                                   theta=np.abs(theta),
                                                   instrument=event_reduction.EventReflectivity.INSTRUMENT_4B)

    # Get the scaling factors
    a, b, err_a, err_b = scaling_factor(template_data.scaling_factor_file, ws_sc)

    # R(Q)
    qz, refl, d_refl = event_refl.specular(q_summing=q_summing, tof_weighted=tof_weighted,
                                           bck_in_q=bck_in_q, clean=clean)
    qz_mid = (qz[:-1] + qz[1:])/2.0
    _tof = 4*np.pi*np.sin(event_refl.theta)*event_refl.constant/qz
    _tof_mid = (_tof[1:] + _tof[:-1])/2.0

    a_q = _tof_mid*b + a
    d_a_q = np.sqrt(_tof_mid**2 * err_b**2 + err_a**2)

    d_refl = np.sqrt(d_refl**2/a_q**2 + refl**2*d_a_q**2/a_q**4)
    refl /= a_q

    return qz_mid, refl, d_refl


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
    r_ref = process_from_template_ws(ref_ws_30Hz, template_data)

    # Reduce the sample data at 30Hz
    r_meas = process_from_template_ws(meas_ws_30Hz, template_data)

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
    return np.asarray([q[_idx], r_q_final[_idx], dr_q_final[_idx]])


def reduce_30Hz_slices(meas_run_30Hz, ref_run_30Hz, ref_data_60Hz, template_30Hz, 
                       time_interval, output_dir, scan_index=1, slice_range=None):
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
    meas_ws_30Hz = api.LoadEventNexus("REF_L_%s" % meas_run_30Hz)
    duration = meas_ws_30Hz.getRun()['duration'].value

    # Time slices
    print("Slicing data: total duration = %s s" % duration)
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
    wsnames = wsgroup.getNames()[:-1]

    # Load the 60Hz reference data
    print("Loading reference R(Q)")
    data_60Hz = np.loadtxt(ref_data_60Hz).T

    reduced = []
    total_time = 0

    if slice_range is not None:
        wsnames = wsnames[:slice_range]

    for name in wsnames:
        tmpws = api.mtd[name]
        print("workspace %s has %d events" % (name, tmpws.getNumberEvents()))
        try:
            _reduced = reduce_30Hz_from_ws(tmpws, ref_ws_30Hz, data_60Hz, template_data, scan_index=scan_index)
            reduced.append(_reduced)
            _filename = 'r{0}_t{1:06d}.txt'.format(meas_run_30Hz, int(total_time))
            np.savetxt(os.path.join(output_dir, _filename), _reduced.T)
        except:
            print(sys.exc_info())
        total_time += time_interval

    return reduced


def plot_slices(reduced, title, time_interval, file_path, offset=1):
    fig, ax = plt.subplots(figsize=(6,6))

    total_time = 0

    _running_offset = 1.
    for _data in reduced:
        qz, refl, d_refl = _data

        plt.errorbar(qz, refl*_running_offset, yerr=d_refl*_running_offset, markersize=4, marker='o', 
                     #linestyle='--', 
                     label='T=%g s' % total_time)

        total_time += time_interval
        _running_offset *= offset

    plt.legend()
    plt.title(title)
    plt.xlabel('q [$1/A$]')
    plt.ylabel('R(q)')
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.show()
    plt.savefig(file_path)
