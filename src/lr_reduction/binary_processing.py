import logging
from pathlib import Path

import h5py
import numpy as np
from scipy.special import lambertw

## NOTES: This moves away from events after the dead-time correction has been applied.
## Various aspects of this are easier through h5py than mantid wksp. Here will include a few functions to
## work through the steps which might not all be needed in the final version but shows the idea.

# TODO: link up the parts that are self. from copying across.

logger = logging.getLogger(__name__)


def convert_to_binary(fname, lowres, collapse_x = True, tofbin=50, tofmax=100000, tofmin=0, deadtime=4.2, tof_step=100, n_y=304, n_x=256, start_times=None, end_times=None):  # noqa: ARG001 -- collapse_x documented as planned/not implemented; callers pass it
    '''
    Main function for converting to load the file, apply the dead-time correction and obtain y vs tof data (non-event).

    :param fname: File to load
    :param lowres: pixel range for low-res direction (i.e. x-pixels for LR). expects of format [min, max]
    :param collapse_x: planned option to keep the x-pixel direction but not implemented. True sums over x-pixels in the lowres range.
    :param tofbin: default bin size for tof histogramming
    :param tofmax: default max tof for histogramming (mainly important for non-standard chopper configurations)
    :param start_times: default None, the whole run. The start of each time window, in seconds from the run's first
        pulse (see load_and_extract and time_windows). Must match the length of end_times.
    :param end_times: default None. The end of each time window. Must match the length of start_times.
    :return: tof_array, y_tof_corr, error_array_corr, log_values, DTC
    :raises ValueError: when no pulse read has proton charge (a window after the run or between two pulses, or a
        run without beam): there is nothing to normalise by.
    '''
    # Include option to collapse along the x-pixel direction between the min/max bounds
    # Assume wksp has had the DTC applied.

    # Example using the h5py to extract info as can be easier to manipulate that mantid wksp
    e_offset, event_id, error_event_offset, pcharge, cPc, log_values = load_and_extract(fname, start_times=start_times, end_times=end_times)

    if not np.any(cPc != 0):
        raise ValueError(f"No proton charge in {Path(fname).name}{describe_windows(start_times, end_times)}: "
                         "nothing to reduce")

    tof_array, DTC, error_counts = get_deadtime_correction(error_event_offset, e_offset, cPc, tofbin, tofmax, tofmin, deadtime=deadtime, tof_step=tof_step)
    tof_array, y_tof, error_array = get_y_tof(tof_array, event_id, e_offset, lowres, pcharge, n_y, n_x)


    #y_tof_collapse = np.sum(y_tof, axis=0)
    # Apply the dead-time correction
    y_tof_corr = y_tof * DTC
    error_array_corr = error_array * DTC
    y_tof_corr = np.nan_to_num(y_tof_corr, nan=0)
    error_array_corr = np.nan_to_num(error_array_corr, nan=0)

    tof_array = tof_array / 1000
    return tof_array, y_tof_corr, error_array_corr, log_values, DTC


def time_windows(start_times, end_times):
    '''
    The time windows to read, as (start, end) pairs in the order given, or None for the whole run.

    Times are seconds from the run's first pulse (event_time_zero). A window is half-open, start <= t < end, at every
    end (pulse_mask); the slicing functions close the run's final window by moving its end past the last pulse
    (close_final_window). Two numbers (or 0-d arrays) are one window.
    Windows may be disjoint or overlap: the selection concatenates them and deduplicates nothing.

    :raises ValueError: one of the lists missing, lists of unequal length, an empty list (None is the spelling of
        "the whole run"), or a start not before its end. Raised before anything is read.
    '''
    if start_times is None and end_times is None:
        return None
    if start_times is None or end_times is None:
        raise ValueError("start_times and end_times must either both be provided or both be None.")
    starts = [start_times] if np.ndim(start_times) == 0 else list(start_times)  # a number, or a 0-d array: one window
    ends = [end_times] if np.ndim(end_times) == 0 else list(end_times)
    if len(starts) != len(ends):
        raise ValueError("Start and end times for time slices must be the same length.")
    if not starts:
        raise ValueError("Empty start and end time lists hold no window; None reads the whole run.")
    for start, end in zip(starts, ends):
        if start >= end:
            raise ValueError(f"Start time ({start}) must be less than end time ({end}).")
    return list(zip(starts, ends))


def describe_windows(start_times, end_times):
    '''" for the time window(s) [a, b), ..." for messages, or "" for the whole run.'''
    windows = time_windows(start_times, end_times)
    if windows is None:
        return ""
    return " for the time window(s) " + ", ".join(f"[{start}, {end})" for start, end in windows)


def load_and_extract(fname, start_times = None, end_times = None):
    '''
    Load the nexus file and extract the relevant arrays using h5py.

    With start_times and end_times, only the pulses inside the time windows are kept (time_windows says what a window
    is; event_time_filter selects). The windows are checked before the file is opened.

    :param fname: File to load
    :param start_times: default None, the whole run; else the windows' starts, seconds from the run's first pulse
    :param end_times: default None; else the windows' ends
    :return: e_offset, event_id, error_event_offset, pcharge, cPC, log_values. For the whole run pcharge is the run's
        total as the file records it (entry/proton_charge), as before time slicing; for windows it is the sum of the
        selected pulses' charge.
    '''
    windows = time_windows(start_times, end_times)
    with h5py.File(fname, 'r') as f:
        e_offset = np.array(f['entry/bank1_events/event_time_offset'][:]) # TOF
        event_id = np.array(f['entry/bank1_events/event_id'][:]) # position on detector
        error_event_offset = np.array(f['entry/bank_error_events/event_time_offset'][:]) # error TOF
        cPC=np.array(f['entry/DASlogs/proton_charge/value'][:]) # the charge of each pulse
        if windows is None:
            # This is single value. TODO: streamline so don't need this and the previous log.
            pcharge=np.array(f['entry/proton_charge'][:])
        else:
            pulses = {
                "event_time": np.array(f['entry/bank1_events/event_time_zero'][:]), # pulses, s
                "event_index": np.array(f['entry/bank1_events/event_index'][:]), # each pulse's first event
                "error_event_time": np.array(f['entry/bank_error_events/event_time_zero'][:]), # error pulses, s
                "error_event_index": np.array(f['entry/bank_error_events/event_index'][:]), # each error pulse's first event
                "charge_time": np.array(f['entry/DASlogs/proton_charge/time'][:]), # the charge log's pulses, s
            }

    log_values = get_log_values(fname)

    if windows is None:
        return e_offset, event_id, error_event_offset, pcharge, cPC, log_values

    logger.info("Processing with time filter: %s%s", Path(fname).name, describe_windows(start_times, end_times))
    masked_e_offset, masked_event_id, masked_e_offset_error, masked_cPC = event_time_filter(
        start_times, end_times, pulses["event_time"], event_id, e_offset, pulses["event_index"],
        pulses["error_event_time"], error_event_offset, pulses["error_event_index"], cPC, pulses["charge_time"])
    masked_pcharge = np.array([np.sum(masked_cPC)])
    return masked_e_offset, masked_event_id, masked_e_offset_error, masked_pcharge, masked_cPC, log_values


def pulse_mask(pulse_times, start, end):
    '''
    The pulses of a window: a mask over pulse_times, true for each pulse with start <= t < end. A pulse belongs to a
    window by its own time, wherever it sits in the file: real runs can record pulse times that go back near the end
    (REF_L_198410), so the times are never assumed sorted, and an out-of-order pulse is taken by the window its time is
    in. The window is half-open at every end, the last pulse's time included: a window is closed at the run's end only
    by moving its end past the run's last pulse, which the slicing functions do for the final window
    (close_final_window).
    '''
    pulse_times = np.asarray(pulse_times)
    return (pulse_times >= start) & (pulse_times < end)


def event_mask(event_index, n_events, pulses):
    '''
    The events of the selected pulses: a mask over the n_events events, true for each event of a pulse whose mask
    entry is true. Pulse i's events are event_index[i] : event_index[i + 1] (the last pulse's run to n_events); a pulse
    without events selects none.
    '''
    event_index = np.asarray(event_index, dtype=np.int64)
    counts = np.diff(np.append(event_index, n_events))
    return np.repeat(np.asarray(pulses, dtype=bool), counts)


def read_run_end(fname):
    '''
    The run's end, as close_final_window needs it: (last_pulse, duration). last_pulse is the latest pulse time over the
    three banks the selection reads (the detector events, the error events and the charge log), the latest, not the
    last recorded; None for a run with no pulse. duration is entry/duration as recorded (float32 on the real files,
    which can round below the last pulse time), or None.
    '''
    with h5py.File(fname, 'r') as f:
        times = [np.asarray(f[name][:]) for name in ('entry/bank1_events/event_time_zero',
                                                     'entry/bank_error_events/event_time_zero',
                                                     'entry/DASlogs/proton_charge/time') if name in f]
        duration = float(np.asarray(f['entry/duration'][:]).ravel()[0]) if 'entry/duration' in f else None
    latest = [float(np.max(t)) for t in times if len(t)]
    return (max(latest) if latest else None), duration


def final_windows(windows, last_pulse, duration):
    '''
    The indices of the run's final window among windows ((start, end) pairs): those that reach the run's end and end
    last. A window reaches the run's end when its end is at or past the earlier of the run's last pulse time and its
    recorded duration (the float32 duration can round below the last pulse). Of the windows that reach it, the final
    one ends last; windows that overlap and end together at that end are all final, each taking its share of the
    run's last pulses (a pulse in two windows is taken twice). A window list that stops short of the run's end has
    none. Contiguous windows have one.
    '''
    marks = [m for m in (last_pulse, duration) if m is not None]
    if not marks:
        return []
    reach = min(marks)
    ends = [end for _start, end in windows if end >= reach]
    if not ends:
        return []
    last = max(ends)
    return [i for i, (_start, end) in enumerate(windows) if end == last]


def close_final_window(windows, last_pulse, duration):
    '''
    The windows, with the run's final window (final_windows) closed at the run's end: if its end does not pass the
    last pulse, it is moved just past it (numpy.nextafter), so its half-open selection takes every remaining pulse in
    every bank. Every other window is returned as given and stays half-open, even when its end equals the last pulse
    time, so contiguous windows partition the run: no pulse lost at the end, none taken twice.
    '''
    windows = [(start, end) for start, end in windows]
    if last_pulse is None:
        return windows
    for i in final_windows(windows, last_pulse, duration):
        if windows[i][1] <= last_pulse:
            windows[i] = (windows[i][0], float(np.nextafter(last_pulse, np.inf)))
    return windows


def event_time_filter(start_times, end_times, event_time, event_id, e_offset, event_index,
                       error_event_time, error_event_offset, error_event_index, cPC, charge_time):
    '''
    The events, error events and pulse charges of the time windows, concatenated in the order of the windows.

    Each bank is selected by the same predicate on its own pulse times (pulse_mask: start <= t < end, whatever order
    the times are in): the detector events by event_time_zero and event_index, the error events by their bank's, and
    the charge by the proton-charge log's own times. A pulse in two overlapping windows is taken twice.

    :return: masked_e_offset, masked_event_id, masked_e_offset_error, masked_cPC
    '''
    masked_event_id = []
    masked_e_offset = []
    masked_e_offset_error = []
    masked_cPC = []

    for start, end in time_windows(start_times, end_times):
        events = event_mask(event_index, len(event_id), pulse_mask(event_time, start, end))
        masked_event_id.append(event_id[events])
        masked_e_offset.append(e_offset[events])

        errors = event_mask(error_event_index, len(error_event_offset), pulse_mask(error_event_time, start, end))
        masked_e_offset_error.append(error_event_offset[errors])

        masked_cPC.append(cPC[pulse_mask(charge_time, start, end)])

    # Concatenate them back into one.
    masked_event_id = np.concatenate(masked_event_id)
    masked_e_offset = np.concatenate(masked_e_offset)
    masked_e_offset_error = np.concatenate(masked_e_offset_error)
    masked_cPC = np.concatenate(masked_cPC)

    return masked_e_offset, masked_event_id, masked_e_offset_error, masked_cPC


def get_deadtime_correction(error_event_offset, e_offset, pcharge, tofbin=50, tofmax=50000, tofmin=0, use_bad_counts=True, deadtime=4.2, tof_step=100):  # noqa: ARG001 -- tof_step kept for API compatibility (callers pass it)
    '''
    Gets and applies the dead-time correction.

    :param error_event_offset: Description
    :param e_offset: Description
    :param pcharge: Description
    :param tofbin: Description
    :param tofmax: Description
    :param use_bad_counts: Description
    :return: tof_array, DTC, error_counts
    '''
    # Probably want tofbin to dfault to the tof_step value below.

    pGood=len(pcharge[pcharge != 0])

    # Setup arrays for histogramming
    tof_array = np.arange(tofmin, tofmax, tofbin)
    #d_tof = np.diff(tof_array)[0] # Step size
    bin_edges = np.concatenate([[tof_array[0] - tofbin/2], tof_array + tofbin/2])
    # Fill with the event time offsets
    counts, _ = np.histogram(e_offset, bins=bin_edges)

    if use_bad_counts:
        # Include the bad counts for dead-time correction
        bad_counts, _ = np.histogram(error_event_offset, bins=bin_edges)
        counts += bad_counts

    error_counts = np.sqrt(counts) / pGood
    # Normalise by proton charge
    counts_norm = counts / pGood

    # Calculate the dead-time correction - this links to existing expression and properties.
    with np.errstate(divide='ignore', invalid='ignore'):
        b = -lambertw(-counts_norm * deadtime / tofbin) / (deadtime / tofbin)
        DTC = np.real(b / counts_norm)
        DTC = np.nan_to_num(DTC, nan=1.0, posinf=1.0, neginf=1.0)

    return tof_array, DTC, error_counts

def get_y_tof(tof_array, event_id, e_offset, lowres, pcharge, n_y = 304, n_x = 256):  # noqa: ARG001 -- n_x kept for API compatibility (callers pass it)
    '''
    Collapses the event data into y vs tof histogram.

    :param tof_array: array of tof
    :param event_id: array of event ids
    :param e_offset: array of event time offsets
    :param lowres: low-res (i.e. x-direction) pixel range [min, max]. Uses lr_reduction notation.
    :param pcharge: proton charge for normalisation
    :return: tof_array, y_tof, error_array
    '''
    # Get the y vs tof:
    y_tof = np.zeros((n_y, len(tof_array)))
    # convert the event id into x and y pixel values
    xvals = event_id // n_y
    yvals = event_id % n_y

    x_good = (xvals >= lowres[0]) & (xvals <= lowres[1])
    e_offset_good = e_offset[x_good]
    y_good = yvals[x_good]

    # Compute TOF bin edges
    d_tof = np.diff(tof_array)[0] if len(tof_array) > 1 else 1
    bin_edges = np.concatenate([[tof_array[0] - d_tof / 2], tof_array + d_tof / 2])

    # Use np.digitize to assign TOF bins
    bin_indices = np.digitize(e_offset_good, bins=bin_edges) - 1
    bin_indices = np.clip(bin_indices, 0, len(tof_array) - 1)

    np.add.at(y_tof, (y_good, bin_indices), 1)

    error_array = np.sqrt(y_tof)

    pcharge = pcharge[0] if len(pcharge) == 1 else pcharge

    y_tof /= pcharge
    error_array /= pcharge

    return tof_array, y_tof, error_array

def get_log_values(fname):
    # Read any log values needed in the reduction process on a per-run basis.
    f = h5py.File(fname, 'r')
    log_values = {}
    log_values["thi"] = f['entry/DASlogs/BL4B:Mot:thi.RBV/value'][-1]
    log_values["ths"] = f['entry/DASlogs/BL4B:Mot:ths.RBV/value'][-1]
    log_values["tthd"] = f['entry/DASlogs/BL4B:Mot:tthd.RBV/value'][-1]
    # Autoreduction indicators
    log_values["seq_num"] = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Num/value'][0]
    log_values["seq_id"] = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Id/value'][0]
    try:
        #log_values["seq_total"] = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Total/value'][0]
        #log_values["center_pixel"] = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:CenterPixel/value'][0]
        #log_values["data_type"] = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:DataType/value'][0]
        log_values["scale_multiplier"] = f['entry/DASlogs/BL4B:CS:Autoreduce:ScaleMultiplier/value'][0]
        #log_values["distance_sample_detector"] = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:DistanceSampleDetector/value'][0]
        #log_values["template_id"] = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:TemplateId/value'][0]
    except Exception as e:  # noqa: BLE001 -- tolerate any missing/odd DASlog entry; partial log_values is the intended behavior
        logger.warning("Cannot process %s because %s fails to extract", fname, e)
        pass
    # Get slit gap openings.
    log_values["siY"]=np.array(f['entry/DASlogs/BL4B:Mot:si:Y:Gap:Readback/average_value'][0])
    log_values["s1Y"]=np.array(f['entry/DASlogs/BL4B:Mot:s1:Y:Gap:Readback/average_value'][0])
    log_values["siX"]=np.array(f['entry/DASlogs/BL4B:Mot:si:X:Gap:Readback/average_value'][0])
    log_values["s1X"]=np.array(f['entry/DASlogs/BL4B:Mot:s1:X:Gap:Readback/average_value'][0])
    log_values["xi"]=np.array(f['entry/DASlogs/BL4B:Mot:xi.RBV/average_value'][0])

    log_values['start_time'] = f['entry/start_time'].asstr()[0]
    log_values['title'] = f['entry/title'].asstr()[0]

    log_values["op_mode"] = f['entry/DASlogs/BL4B:CS:ExpPl:OperatingMode/value'][0] # This is one that can say "Free Liquid"
    log_values["op_mode"] = log_values["op_mode"][0].decode("utf-8")
    try:
        log_values["coordinates"] = f['entry/DASlogs/BL4B:CS:Mode:Coordinates/value'][0] # This is one that shows earth vs beam center 0=earth; 1=beam
    except:
        logger.info("Older run doesn't include coordinates PV")
        pass

    try:
        log_values["incident_theta"] = f['entry/DASlogs/BL4B:CS:BeamToEarthCenterAngleOffset/value'][0]
    except:
        logger.info("Older run doesn't include Beam-Earth centered angle, set to default 4.0deg")
        log_values["incident_theta"] = 4.0

    # TODO: check if we need any of the other chopper parts.
    log_values["frequency"]=np.array(f['entry/DASlogs/BL4B:Det:TH:BL:Frequency/value'][0])
    log_values["lam_request"]=np.array(f['entry/DASlogs/BL4B:Det:TH:BL:Lambda/value'][0])

    try:
        log_values["chopper_mod"]=np.array(f['entry/DASlogs/BL4B:Chop:Skf2:ChopperModerator/value'][0])
    except:
        logger.warning("Run missing the chopper moderator log value")
        pass
    log_values['emission_mod_distance'] = np.array(f['entry/DASlogs/BL4B:Det:TH:DlyDet:BasePath/value'][0]) * 1000
    off =np.array(f['entry/DASlogs/BL4B:Chop:Skf2:ChopperOffset/value'][0]) # 114.0
    mult =np.array(f['entry/DASlogs/BL4B:Chop:Skf2:ChopperMultiplier/value'][0])  # 29.5
    log_values["emission_coefficients"]=np.array([off/1000, mult/1000])

    log_values["chop2_PD"] = np.array(f['entry/DASlogs/BL4B:Chop:Skf2:PhaseAccuracy/average_value'][0])

    # TODO: Test these!!
    try:
        atten_menu = np.array(f['entry/DASlogs/BL4B:Actuator:Menu/value'][0])
        atten_lookup = {0: [0,0,0,0],
                        1: [1,0,0,0],
                        2: [0,1,0,0],
                        3: [1,1,0,0],
                        4: [0,0,1,0],
                        5: [1,0,1,0],
                        6: [0,1,1,0],
                        7: [1,1,1,0],
                        8: [0,0,0,1],
                        9: [1,0,0,1],
                        10: [0,1,0,1],
                        11: [1,1,0,1],
                        12: [0,0,1,1],
                        13: [1,0,1,1],
                        14: [0,1,1,1],
                        15: [1,1,1,1]
                        }
        log_values['Atten'] = np.array(atten_lookup.get(float(atten_menu)))

    except:
        Att1 = np.array(f['entry/DASlogs/BL4B:Actuator:50MRb/average_value'][0])
        Att2 = np.array(f['entry/DASlogs/BL4B:Actuator:100MRb/average_value'][0])
        Att3 = np.array(f['entry/DASlogs/BL4B:Actuator:200MRb/average_value'][0])
        Att4 = np.array(f['entry/DASlogs/BL4B:Actuator:400MRb/average_value'][0])
        log_values['Atten'] = np.array([Att1,Att2,Att3,Att4])

    f.close()

    return log_values
