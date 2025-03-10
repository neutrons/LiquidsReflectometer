"""
Time-resolved data reduction
"""

import json
import os
import sys

import mantid
import mantid.simpleapi as api
import numpy as np
from mantid.api import *
from mantid.kernel import *
from matplotlib import pyplot as plt

mantid.kernel.config.setLogLevel(3)

from . import template
from .event_reduction import apply_dead_time_correction, compute_resolution


def reduce_30Hz_from_ws(
    meas_ws_30Hz,
    ref_ws_30Hz,
    data_60Hz,
    template_data,
    scan_index=1,  # noqa ARG001
    template_reference=None,
    q_summing=False,
):
    """
    Perform 30Hz reduction
    @param meas_ws_30Hz: Mantid workspace of the data we want to reduce
    @param ref_ws_30Hz: Mantid workspace of the reference data, take with the same config
    @param data_60Hz: reduced reference data at 60Hz
    @param template_data: template data object (for 30Hz)
    @param scan_index: scan index to use within the template.
    """
    # Reduce the reference at 30Hz
    if template_reference is None:
        r_ref = template.process_from_template_ws(ref_ws_30Hz, template_data, q_summing=q_summing, normalize=False)
    else:
        r_ref = template.process_from_template_ws(ref_ws_30Hz, template_reference, q_summing=q_summing, normalize=False)

    # Reduce the sample data at 30Hz
    r_meas = template.process_from_template_ws(meas_ws_30Hz, template_data, q_summing=q_summing, normalize=False)

    # Identify the bins we need to overlap with the 30Hz measurement
    # The assumption is that the binning is the same
    _max_q = min(r_ref[0].max(), r_meas[0].max())
    _min_q = max(r_ref[0].min(), r_meas[0].min())

    _binning_60hz = (data_60Hz[0][1] - data_60Hz[0][0]) / data_60Hz[0][0]
    _binning_ref = (r_ref[0][1] - r_ref[0][0]) / r_ref[0][0]
    _binning_meas = (r_meas[0][1] - r_meas[0][0]) / r_meas[0][0]

    if not np.fabs(_binning_60hz - _binning_ref) < 1e-5:
        print("ERROR: The binning of the 60 Hz reference is not the same as the dynamic template")
        print("   %s <> %s" % (_binning_60hz, _binning_ref))

    # The tolerance will be half a bin
    _tolerance = _binning_ref / 2.0

    print("60Hz:      %g %g    [%g]" % (data_60Hz[0].min(), data_60Hz[0].max(), _binning_60hz))
    print("Ref 30Hz:  %g %g    [%g]" % (r_meas[0].min(), r_meas[0].max(), _binning_ref))
    print("Meas 30Hz: %g %g    [%g]" % (r_ref[0].min(), r_ref[0].max(), _binning_meas))

    _q_idx_60 = np.asarray(np.where((data_60Hz[0] > _min_q * (1 - _tolerance)) & (data_60Hz[0] < _max_q * (1 + _tolerance))))[0]
    _q_idx_meas30 = np.asarray(np.where((r_meas[0] > _min_q * (1 - _tolerance)) & (r_meas[0] < _max_q * (1 + _tolerance))))[0]
    _q_idx_ref30 = np.asarray(np.where((r_ref[0] > _min_q * (1 - _tolerance)) & (r_ref[0] < _max_q * (1 + _tolerance))))[0]

    if not data_60Hz[0][_q_idx_60].shape[0] == r_ref[0][_q_idx_ref30].shape[0]:
        print("\n\n60Hz reference may have been reduced with different binning!")
        print("Array sizes: %g %g %g" % (len(_q_idx_60), len(_q_idx_meas30), len(_q_idx_ref30)))
        print("Remember to average overlapping points!\n\n")

    r_q_final = r_meas[1][_q_idx_meas30] / r_ref[1][_q_idx_ref30] * data_60Hz[1][_q_idx_60]

    dr_q_final = np.sqrt(
        (r_meas[2][_q_idx_meas30] / r_ref[1][_q_idx_ref30] * data_60Hz[1][_q_idx_60]) ** 2
        + (r_meas[1][_q_idx_meas30] / r_ref[1][_q_idx_ref30] * data_60Hz[2][_q_idx_60]) ** 2
        + (r_meas[1][_q_idx_meas30] / r_ref[1][_q_idx_ref30] ** 2 * data_60Hz[1][_q_idx_60] * r_ref[2][_q_idx_ref30]) ** 2
    )

    print("Q range: %s - %s" % (r_meas[0][0], r_meas[0][_q_idx_meas30][-1]))
    print("Constant-Q binning: %s" % str(q_summing))
    q = r_meas[0][_q_idx_meas30]
    # Skip infinite points if they exist, and skip the first two points which
    # often has binning artifact
    _idx = (q > r_meas[0][1]) & (r_q_final < np.inf) & (q < r_meas[0][-2]) & (r_q_final > 0)

    # Q resolution
    #   Assume a constant term of 0 unless it is specified
    dq0 = 0
    if meas_ws_30Hz.getInstrument().hasParameter("dq-constant"):
        dq0 = meas_ws_30Hz.getInstrument().getNumberParameter("dq-constant")[0]
    dq_slope = compute_resolution(meas_ws_30Hz)
    print("Resolution: %g + %g Q" % (dq0, dq_slope))
    dq = dq0 + dq_slope * q[_idx]
    return np.asarray([q[_idx], r_q_final[_idx], dr_q_final[_idx], dq])


def reduce_30Hz_slices(
    meas_run_30Hz,
    ref_run_30Hz,
    ref_data_60Hz,
    template_30Hz,
    time_interval,
    output_dir,
    scan_index=1,
    create_plot=True,
    template_reference=None,
    q_summing=False,
):  # noqa ARG001
    meas_ws_30Hz = api.LoadEventNexus("REF_L_%s" % meas_run_30Hz)

    return reduce_30Hz_slices_ws(
        meas_ws_30Hz,
        ref_run_30Hz,
        ref_data_60Hz,
        template_30Hz,
        time_interval,
        output_dir,
        scan_index=scan_index,
        create_plot=create_plot,
        template_reference=template_reference,
        q_summing=False,
    )


def reduce_slices(meas_run, template_file, time_interval, output_dir, scan_index=1, theta_offset=None, create_plot=True):
    meas_ws = api.LoadEventNexus("REF_L_%s" % meas_run)

    return reduce_slices_ws(
        meas_ws, template_file, time_interval, output_dir, scan_index=scan_index, theta_offset=theta_offset, create_plot=create_plot
    )


def reduce_30Hz_slices_ws(
    meas_ws_30Hz,
    ref_run_30Hz,
    ref_data_60Hz,
    template_30Hz,
    time_interval,
    output_dir,
    scan_index=1,
    create_plot=True,
    template_reference=None,
    q_summing=False,
):
    """
    Perform 30Hz reduction
    @param meas_ws_30Hz: workspace of the data we want to reduce
    @param ref_ws_30Hz: workspace of the reference data, take with the same config
    @param ref_data_60Hz: file path of the reduce data file at 60Hz
    @param template_30Hz: file path of the template file for 30Hz
    @param time_interval: time step in seconds
    @param scan_index: scan index to use within the template.
    """
    # Save options
    options = dict(
        meas_run_30Hz=meas_ws_30Hz.getRun()["run_number"].value,
        ref_run_30Hz=ref_run_30Hz,
        ref_data_60Hz=ref_data_60Hz,
        template_30Hz=template_30Hz,
        time_interval=time_interval,
        output_dir=output_dir,
        scan_index=scan_index,
        template_reference=template_reference,
        q_summing=q_summing,
    )
    with open(os.path.join(output_dir, "options.json"), "w") as fp:
        json.dump(options, fp)

    # Load the template
    print("Reading template: %s" % template_30Hz)
    template_data = template.read_template(template_30Hz, scan_index)

    # Reduce the quartz at 30Hz
    print("Reading reference data at 30Hz")
    if os.path.isfile(ref_run_30Hz):
        ref_ws_30Hz = api.LoadEventNexus(ref_run_30Hz)
    else:
        ref_ws_30Hz = api.LoadEventNexus("REF_L_%s" % ref_run_30Hz)

    # Reduce the sample data at 30Hz
    print("Reading sample data at 30Hz")
    # meas_ws_30Hz = api.LoadEventNexus("REF_L_%s" % meas_run_30Hz)

    # Some meta-data are not filled in for the live data stream
    # Use dummy values for those
    try:
        duration = meas_ws_30Hz.getRun()["duration"].value
    except:
        duration = 0
    try:
        meas_run_30Hz = meas_ws_30Hz.getRun()["run_number"].value
    except:
        meas_run_30Hz = 0

    # Time slices
    print("Slicing data")
    splitws, infows = api.GenerateEventsFilter(InputWorkspace=meas_ws_30Hz, TimeInterval=time_interval)

    api.FilterEvents(
        InputWorkspace=meas_ws_30Hz,
        SplitterWorkspace=splitws,
        InformationWorkspace=infows,
        OutputWorkspaceBaseName="time_ws",
        GroupWorkspaces=True,
        FilterByPulseTime=True,
        OutputWorkspaceIndexedFrom1=True,
        CorrectionToSample="None",
        SpectrumWithoutDetector="Skip",
        SplitSampleLogs=False,
        OutputTOFCorrectionWorkspace="mock",
    )
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
            _reduced = reduce_30Hz_from_ws(
                tmpws,
                ref_ws_30Hz,
                data_60Hz,
                template_data,
                scan_index=scan_index,
                template_reference=template_reference,
                q_summing=q_summing,
            )
            # Remove first point
            reduced.append(_reduced)
            _filename = "r{0}_t{1:06d}.txt".format(meas_run_30Hz, int(total_time))
            np.savetxt(os.path.join(output_dir, _filename), _reduced.T)
        except:
            print("reduce_30Hz_slices_ws: %s" % sys.exc_info()[0])
        total_time += time_interval

    # Save output
    output_file = os.path.join(output_dir, "r%s-time-resolved.json" % meas_run_30Hz)
    print("Saving t-NR to %s" % output_file)
    package_json_data(meas_run_30Hz, output_dir, out_array=output_file)

    if create_plot:
        plot_slices(
            reduced,
            title="Duration: %g seconds" % duration,
            time_interval=time_interval,
            file_path=os.path.join(output_dir, "r%s.png" % meas_run_30Hz),
        )

    return reduced


def reduce_slices_ws(meas_ws, template_file, time_interval, output_dir, scan_index=1, theta_offset=0, theta_value=None, create_plot=True):
    """
    Perform time-resolved reduction
    :param meas_ws: workspace of the data we want to reduce
    :param template_file: autoreduction template file
    :param time_interval: time step in seconds
    :param scan_index: scan index to use within the template.
    :param theta_value: force theta value
    :param theta_offset: add a theta offset, defaults to zero
    """
    # Save options
    options = dict(
        meas_run=meas_ws.getRun()["run_number"].value,
        template=template_file,
        time_interval=time_interval,
        output_dir=output_dir,
        scan_index=scan_index,
        theta_value=theta_value,
        theta_offset=theta_offset,
    )
    with open(os.path.join(output_dir, "options.json"), "w") as fp:
        json.dump(options, fp)

    # Load the template
    print("Reading template")
    template_data = template.read_template(template_file, scan_index)

    # Add theta offset
    if theta_offset is not None:
        print("Theta offset: %s" % theta_offset)
        template_data.angle_offset = theta_offset

    # Some meta-data are not filled in for the live data stream
    # Use dummy values for those
    try:
        duration = meas_ws.getRun()["duration"].value
    except:
        duration = 0
    try:
        meas_run = meas_ws.getRun()["run_number"].value
    except:
        meas_run = 0

    # Apply dead time correction up front since we are using the error events but
    # not filtering them.
    if template_data.dead_time:
        apply_dead_time_correction(meas_ws, template_data)

    # Time slices
    print("Slicing data")
    splitws, infows = api.GenerateEventsFilter(InputWorkspace=meas_ws, TimeInterval=time_interval)

    api.FilterEvents(
        InputWorkspace=meas_ws,
        SplitterWorkspace=splitws,
        InformationWorkspace=infows,
        OutputWorkspaceBaseName="time_ws",
        GroupWorkspaces=True,
        FilterByPulseTime=True,
        OutputWorkspaceIndexedFrom1=True,
        CorrectionToSample="None",
        SpectrumWithoutDetector="Skip",
        SplitSampleLogs=False,
        OutputTOFCorrectionWorkspace="mock",
    )
    wsgroup = api.mtd["time_ws"]
    wsnames = wsgroup.getNames()

    # Load DB so we do it only once
    ws_db = api.LoadEventNexus("REF_L_%s" % template_data.norm_file)

    # Apply dead time correction up-front
    if template_data.dead_time:
        apply_dead_time_correction(ws_db, template_data)

    # Turn off dead time since we already did it
    template_data.dead_time = False

    reduced = []
    total_time = 0
    for name in wsnames:
        tmpws = api.mtd[name]
        print("workspace %s has %d events" % (name, tmpws.getNumberEvents()))
        try:
            _reduced = template.process_from_template_ws(tmpws, template_data, theta_value=theta_value, ws_db=ws_db)

            dq0 = 0
            dq_slope = compute_resolution(tmpws)
            dq = dq0 + dq_slope * _reduced[0]
            _reduced = [_reduced[0], _reduced[1], _reduced[2], dq]
            _reduced = np.asarray(_reduced)

            reduced.append(_reduced)
            _filename = "r{0}_t{1:06d}.txt".format(meas_run, int(total_time))
            np.savetxt(os.path.join(output_dir, _filename), _reduced.T)
        except:
            print("reduce_slices_ws: %s" % sys.exc_info()[0])
        total_time += time_interval

    if create_plot:
        plot_slices(
            reduced,
            title="Duration: %g seconds" % duration,
            time_interval=time_interval,
            file_path=os.path.join(output_dir, "r%s.png" % meas_run),
        )

    return reduced


def plot_slices(reduced, title, time_interval, file_path, offset=10, show=True):
    fig, ax = plt.subplots(figsize=(6, 6))

    total_time = 0
    _running_offset = 1.0
    for _data in reduced:
        qz, refl, d_refl, _ = _data

        plt.errorbar(qz, refl * _running_offset, yerr=d_refl * _running_offset, markersize=4, marker="o", label="T=%g s" % total_time)

        total_time += time_interval
        _running_offset *= offset

    plt.legend()
    plt.title(title)
    plt.xlabel(r"q [$1/\AA$]")
    plt.ylabel("R(q)")
    ax.set_yscale("log")
    ax.set_xscale("log")
    if show:
        plt.show()
    plt.savefig(file_path)


def package_json_data(dynamic_run, dyn_data_dir, out_array=None):
    compiled_array = []
    compiled_times = []

    _file_list = sorted(os.listdir(dyn_data_dir))

    # Get only the files for the run we're interested in
    _good_files = [_f for _f in _file_list if _f.startswith("r%s_t" % dynamic_run)]

    for i, _file in enumerate(_good_files):
        if _file.startswith("r%s_t" % dynamic_run):
            _data = np.loadtxt(os.path.join(dyn_data_dir, _file)).T
            _data_name, _ = os.path.splitext(_file)
            _time = int(_data_name.replace("r%s_t" % dynamic_run, ""))
            compiled_array.append(_data.tolist())
            compiled_times.append(_time)

    if out_array:
        with open(out_array, "w") as fp:
            json.dump(dict(times=compiled_times, data=compiled_array), fp)

    return compiled_times, compiled_array
