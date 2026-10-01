# New reduction time slicing

from pathlib import Path

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

import lr_reduction.new_reduction_from_file as reduction

def reduce_time_slices(run, settings_file, experiment_id, num_slices, savepath=None, plot_time = True, plot_ref=False, subname_input=None):
    '''
    Function to reduce the data, splitting into the number of time slices
    '''

    # Get the duration of the run

    # TODO: add a part to read the nexus location from the setting file
    nexus_path = Path("/SNS/REF_L") / experiment_id / "nexus"
    fname = f"REF_L_{run}.nxs.h5"

    f = h5py.File(nexus_path / fname, 'r')
    duration = np.array(f['entry/duration'][0])

    # determine the list of start/stop values
    time_int = duration/num_slices
    starts = [0]
    stops = [duration]
    for ii in range(1,num_slices):
        starts.append(time_int*ii)
        stops.insert(-1,time_int*ii)

    all_outputs = []
    for ii in range(num_slices):
        if not subname_input:
            subname = f"slice_{ii+1}of{num_slices}"
        else:
            subname = f"{subname_input}_slice_{ii+1}of{num_slices}"
        print('starting number', ii+1, starts[ii], stops[ii])
        slice_outputs, _ = reduce_time_list(run, settings_file, experiment_id,
                                        starts=[starts[ii]], ends=[stops[ii]], savepath=savepath,
                                        plot_ref=plot_ref, plot_time=False, subname_input=subname)
        print('finished', ii+1, starts[ii], stops[ii])
        all_outputs.extend(slice_outputs)
    if plot_time:
        plots = plot_kinetic(all_outputs)
    else:
        plots = None

    return all_outputs, plots


def reduce_time_list(run, settings_file, experiment_id, starts, ends, savepath=None, plot_time = True, plot_ref=False, subname_input=None):

    # starts and ends are lists for each separate file. Each of these can be a nested list.

    if len(starts) != len(ends):
        raise ValueError("Length of starts and ends lists must be the same")

    run_list = [run]
    if savepath:
        Spath = Path(savepath)
    else:
        Spath = Path("/SNS/REF_L") / experiment_id / "shared" / "reduced"

    store_outputs = []
    # run the looped reduction
    for slice_idx in range(len(starts)):
        print("starting:", slice_idx, starts[slice_idx], ends[slice_idx])
        if not subname_input:
            subname = f"slice_{int(starts[slice_idx])}_{int(ends[slice_idx])}"
        else:
            subname = f"{subname_input}_slice_{int(starts[slice_idx])}_{int(ends[slice_idx])}"

        override_params = {'Spath': Spath, "subname": subname}
        output = reduction.reduce_from_file(run_list, settings_file, experiment_id,
                                    override_params=override_params, plot=plot_ref, save_json=False,
                                    start_times=starts[slice_idx], end_times=ends[slice_idx])
        print("finish:", slice_idx, starts[slice_idx], ends[slice_idx])

        all_results, _, _, _ = output
        flat_results = flatten_reduced_results(all_results)
        if not flat_results:
            raise ValueError(f"No reduced result data returned for slice {slice_idx} from run {run}")

        store_outputs.append(flat_results)
        print(len(store_outputs))
    # create plot of set
    if plot_time:
        plots = plot_kinetic(store_outputs)
    else:
        plots = None

    return store_outputs, plots

# TODO: add this one.
'''
def reduce_time_log_filter(run, settings, log_id, log_min, log_max):

    # Get the log values for the run
    # determine the list of start/stop values
    # run the looped reduction
    # create plot of set

    return full_set
'''


def plot_kinetic(output_list):
    # Plot an offset graph and a colour map.
    # output_list is expected to be a list of per-slice result packs, where each pack is a
    # list of reduced data dicts, e.g. [ {"Q":..., "R":..., "dR":..., "dQ":...}, ... ]

    fig, ax = plt.subplots()
    store_q = []
    store_r = []
    store_dr = []
    store_dq = []

    spacing = 0.5

    for ii, slice_results in enumerate(output_list):
        for dataset in slice_results:
            if not isinstance(dataset, dict):
                continue
            if not {"Q", "R", "dR", "dQ"}.issubset(dataset.keys()):
                continue

            Q_vals = np.asarray(dataset["Q"], dtype=float)
            R_vals = np.asarray(dataset["R"], dtype=float)
            dR = np.asarray(dataset["dR"], dtype=float)
            dQ = np.asarray(dataset["dQ"], dtype=float)

            offset = 10**(ii * spacing)
            ax.errorbar(Q_vals, R_vals * offset, yerr=dR * offset, xerr=dQ, fmt='o', markersize=1)

            store_q.append(Q_vals)
            store_r.append(R_vals)
            store_dr.append(dR)
            store_dq.append(dQ)

    ax.set_ylabel('R')
    ax.set_xscale('log')
    ax.set_title('Time_slices', fontsize=16)
    ax.set_yscale('log')
    Angstrom = '\u212B'
    ax.set_xlabel('Q [1/' + Angstrom + ']', fontsize=14)

    # create the colour map. Store the data into a 2D array and plot with imshow
    fig2, ax2 = plt.subplots()
    if not store_dr:
        raise ValueError("No reduced data available for kinetic plot")

    Z = np.array([np.asarray(arr) for arr in store_dr])
    mask = np.isfinite(Z) & (Z > 0)
    if not np.any(mask):
        raise ValueError("No positive dR values available for kinetic plot")

    vmin, vmax = np.percentile(Z[mask], [2, 98])
    im = ax2.imshow(Z, aspect='auto', origin='lower', extent=[store_q[0][0], store_q[0][-1], 0, len(store_dr)-1],
                    norm=LogNorm(vmin=vmin, vmax=vmax), cmap='viridis')

    fig2.colorbar(im, ax=ax2, label='R')
    ax2.set_xlabel('Q [1/' + Angstrom + ']', fontsize=14)
    ax2.set_xscale('log')
    ax2.set_ylabel('Time slice')
    ax2.set_title('Time_slices', fontsize=16)

    plt.show()

    return fig, fig2


def flatten_reduced_results(result):
    """Normalize a reduce_from_file() output into a flat list of result dicts.

    The same reduction call can return either a single dict, a list of dicts, or a nested
    list structure depending on whether a prior-combine step was used. This helper ensures
    the plotting code always gets a flat list of reduced datasets.
    """
    flat = []

    if isinstance(result, dict):
        return [result]

    if isinstance(result, (list, tuple)):
        for item in result:
            if isinstance(item, dict):
                flat.append(item)
            elif isinstance(item, (list, tuple)):
                flat.extend(flatten_reduced_results(item))

    return flat