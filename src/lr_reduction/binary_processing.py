import h5py
import numpy as np
from scipy.special import lambertw

## NOTES: This moves away from events after the dead-time correction has been applied.
## Various aspects of this are easier through h5py than mantid wksp. Here will include a few functions to 
## work through the steps which might not all be needed in the final version but shows the idea.

# TODO: link up the parts that are self. from copying across.

def convert_to_binary(fname, lowres, collapse_x = True, tofbin=50, tofmax=100000, tofmin=0, deadtime=4.2, tof_step=100, n_y=304, n_x=256):
    '''
    Main function for converting to load the file, apply the dead-time correction and obtain y vs tof data (non-event).
    
    :param fname: File to load
    :param lowres: pixel range for low-res direction (i.e. x-pixels for LR). expects of format [min, max]
    :param collapse_x: planned option to keep the x-pixel direction but not implemented. True sums over x-pixels in the lowres range.
    :param tofbin: default bin size for tof histogramming
    :param tofmax: default max tof for histogramming (mainly important for non-standard chopper configurations)
    :return: tof_array, y_tof_corr, error_array_corr
    '''
    # Include option to collapse along the x-pixel direction between the min/max bounds
    # Assume wksp has had the DTC applied.

    # Example using the h5py to extract info as can be easier to manipulate that mantid wksp
    e_offset, event_id, error_event_offset, pcharge, cPc, log_values = load_and_extract(fname)
    tof_array, DTC, error_counts = get_deadtime_correction(error_event_offset, e_offset, cPc, tofbin, tofmax, tofmin, deadtime=deadtime, tof_step=tof_step)
    tof_array, y_tof, error_array = get_y_tof(tof_array, event_id, e_offset, lowres, pcharge, n_y, n_x)


    #y_tof_collapse = np.sum(y_tof, axis=0)
    # Apply the dead-time correction
    y_tof_corr = y_tof * DTC
    error_array_corr = error_array * DTC
    y_tof_corr = np.nan_to_num(y_tof_corr, nan=0)
    error_array_corr = np.nan_to_num(error_array_corr, nan=0)

    tof_array = tof_array / 1000
    return tof_array, y_tof_corr, error_array_corr, log_values

def load_and_extract(fname):
    '''
    Load the nexus file and extract the relevant arrays using h5py.

    :param fname: File to load
    :return: e_offset, event_id, error_event_offset, pcharge
    '''
    # Example using the h5py to extract info as can be easier to manipulate that mantid wksp
    f = h5py.File(fname, 'r')

    e_offset = np.array(f['entry/bank1_events/event_time_offset'][:])
    event_id = np.array(f['entry/bank1_events/event_id'][:])
    error_event_offset = np.array(f['entry/bank_error_events/event_time_offset'][:])
    pcharge=np.array(f['entry/proton_charge'][:])
    # This is single value. TODO: streamline so don't need this and the previous log.
    cPC=np.array(f['entry/DASlogs/proton_charge/value'][:])
    
    log_values = get_log_values(fname)

    return e_offset, event_id, error_event_offset, pcharge, cPC, log_values


def get_deadtime_correction(error_event_offset, e_offset, pcharge, tofbin=50, tofmax=50000, tofmin=0, use_bad_counts=True, deadtime=4.2, tof_step=100):
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

def get_y_tof(tof_array, event_id, e_offset, lowres, pcharge, n_y = 304, n_x = 256):
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
    log_values["seq_num"] = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Num/value'][0]
    log_values["seq_id"] = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Id/value'][0]
    # Get slit gap openings.
    log_values["siY"]=np.array(f['entry/DASlogs/BL4B:Mot:si:Y:Gap:Readback/average_value'][0])
    log_values["s1Y"]=np.array(f['entry/DASlogs/BL4B:Mot:s1:Y:Gap:Readback/average_value'][0])
    log_values["siX"]=np.array(f['entry/DASlogs/BL4B:Mot:si:X:Gap:Readback/average_value'][0])
    log_values["s1X"]=np.array(f['entry/DASlogs/BL4B:Mot:s1:X:Gap:Readback/average_value'][0])
    log_values["xi"]=np.array(f['entry/DASlogs/BL4B:Mot:xi.RBV/average_value'][0])

    log_values['start_time'] = f['entry/start_time'].asstr()[0]

    log_values["op_mode"] = f['entry/DASlogs/BL4B:CS:ExpPl:OperatingMode/value'][0] # This is one that can say "Free Liquid"
    try:
        log_values["coordinates"] = f['entry/DASlogs/BL4B:CS:Mode:Coordinates/value'][0] # This is one that shows earth vs beam center 0=earth; 1=beam
    except:
        print("Older run doesn't include coordinates PV")
        pass

    # TODO: check if we need any of the other chopper parts.
    log_values["frequency"]=np.array(f['entry/DASlogs/BL4B:Det:TH:BL:Frequency/value'][0])
    log_values["lam_request"]=np.array(f['entry/DASlogs/BL4B:Det:TH:BL:Lambda/value'][0])
    off =np.array(f['entry/DASlogs/BL4B:Chop:Skf2:ChopperOffset/value'][0]) # 114.0
    mult =np.array(f['entry/DASlogs/BL4B:Chop:Skf2:ChopperMultiplier/value'][0])  # 29.5
    log_values["emission_coefficients"]=np.array([off/1000, mult/1000])

        
    f.close()
    return log_values

