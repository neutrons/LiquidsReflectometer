# New reduction time slicing

import logging
from pathlib import Path

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure

import lr_reduction.binary_processing as BP
import lr_reduction.new_reduction_from_file as reduction

logger = logging.getLogger(__name__)


def run_folders(settings_file, experiment_id, datapath=None, savepath=None):
    '''
    The run's NeXus folder and the output folder, as the reduction itself resolves them: from the config that
    reduce_from_file builds out of the settings file. NEXUSpathRB is datapath, else the settings' own, else the
    experiment's nexus folder; Spath is savepath, else the settings' own, else the experiment's reduced folder.
    '''
    config = reduction.json_to_config(reduction.load_from_file(settings_file)["config"])
    config.experiment_id = experiment_id
    if datapath:
        config.NEXUSpathRB = datapath
    if savepath:
        config.Spath = savepath
    return Path(config.NEXUSpathRB), Path(config.Spath)


def window_span(start, end):
    '''The first start and last end of a slice's windows (one window, or a list of windows read as one slice).'''
    windows = BP.time_windows(start, end)
    return min(s for s, _ in windows), max(e for _, e in windows)


def window_text(start, end):
    '''A slice's windows for a message: "[a, b)", or "[a, b), [c, d)"; windows that are not valid, as given.'''
    try:
        return ", ".join(f"[{s}, {e})" for s, e in BP.time_windows(start, end))
    except ValueError:
        return f"[{start}, {end})"


def close_final_entries(starts, ends, last_pulse, duration):
    '''
    ends, with the run's final window closed (binary_processing.close_final_window), chosen over every window of every
    entry: reduce_time_list reduces each entry by itself, so only it knows which window is the run's final one. An entry
    that is a list of windows stays a list. An entry that is not a valid window is left out of the choice, and is
    reported when it is reduced.
    '''
    flat = []
    for k, (start, end) in enumerate(zip(starts, ends)):
        try:
            windows = BP.time_windows(start, end)
        except ValueError:
            continue
        flat.extend((k, j, window) for j, window in enumerate(windows))
    closed = BP.close_final_window([window for _, _, window in flat], last_pulse, duration)
    ends = [list(end) if np.ndim(end) else end for end in ends]
    for (k, j, _), (_, end) in zip(flat, closed):
        if isinstance(ends[k], list):
            ends[k][j] = end
        else:
            ends[k] = end
    return ends


def reduce_time_slices(run, settings_file, experiment_id, num_slices, savepath=None, plot_time = True, plot_ref=False, subname_input=None, show_plots=True, datapath=None, figure=None):
    '''
    Function to reduce the data, splitting into the number of time slices

    The run is split into num_slices equal windows over its duration (entry/duration), and each is reduced through
    reduce_time_list with slice_<i>of<n> (or <subname_input>_slice_<i>of<n>) as its subname_input, so that slice's files
    carry slice_<i>of<n>_slice_<start>_<end>. The last window is the run's final one and takes every remaining pulse:
    entry/duration is float32 and can round below the last pulse time, so its end is moved past the run's last pulse
    when it does not pass it already (binary_processing.close_final_window). The slices partition the run.
    The run's file and the output folder are resolved as the reduction resolves them (run_folders).

    :return: (outputs, plots): one flat list of reduced data per slice, and the kinetic plot (None without plot_time),
        drawn into figure when one is given (plot_kinetic)
    :raises ValueError: num_slices below 1; or, once every other slice has been reduced and written, the slices that
        could not be reduced, each with its window and its error.
    '''
    if num_slices < 1:
        raise ValueError(f"The number of time slices must be at least 1, not {num_slices}.")

    nexus_path, _ = run_folders(settings_file, experiment_id, datapath)
    fname = f"REF_L_{run}.nxs.h5"

    with h5py.File(nexus_path / fname, 'r') as f:
        duration = np.array(f['entry/duration'][0])

    # determine the list of start/stop values
    time_int = duration/num_slices
    starts = [0]
    stops = [duration]
    for ii in range(1,num_slices):
        starts.append(time_int*ii)
        stops.insert(-1,time_int*ii)

    windows = BP.close_final_window(list(zip(starts, stops)), *BP.read_run_end(nexus_path / fname))

    mid_points = []
    all_outputs = []
    failures = []
    for ii in range(num_slices):
        if not subname_input:
            subname = f"slice_{ii+1}of{num_slices}"
        else:
            subname = f"{subname_input}_slice_{ii+1}of{num_slices}"
        logger.info("Slice %d of %d: %s to %s s", ii + 1, num_slices, starts[ii], stops[ii])
        try:
            slice_outputs, _ = reduce_time_list(run, settings_file, experiment_id,
                                            starts=[windows[ii][0]], ends=[windows[ii][1]], savepath=savepath,
                                            plot_ref=plot_ref, plot_time=False, subname_input=subname,
                                            datapath=nexus_path, close_final=False)
        except ValueError as error:  # the slice's own report, naming its window: the other slices still run
            failures.append(str(error))
            continue
        logger.info("Slice %d of %d reduced", ii + 1, num_slices)
        all_outputs.extend(slice_outputs)
        mid_point = (stops[ii] - starts[ii]) / 2 + starts[ii]
        mid_points.append(mid_point)

    if failures:
        raise ValueError("; ".join(failures))

    if plot_time:
        plots = plot_kinetic(all_outputs, run, times=mid_points, show=show_plots, figure=figure)
    else:
        plots = None

    return all_outputs, plots


def reduce_time_list(run, settings_file, experiment_id, starts, ends, savepath=None, plot_time = True, plot_ref=False, subname_input=None, show_plots=True, datapath=None, close_final=True, figure=None):
    '''
    Reduce the run once per entry of starts/ends, in order.

    starts and ends are lists for each separate file. Each entry is one window (two numbers) or a list of windows
    read as one slice (see binary_processing.time_windows), named by its span, slice_<earliest start>_<latest end> in
    whole seconds (truncated), or <subname_input>_slice_<...>. The run's file and the output folder are resolved as the
    reduction resolves them (run_folders).

    Every window is half-open, start <= t < end on each pulse's own time, except the run's final one: with close_final
    (the default), the window that reaches the run's end and ends last takes every remaining pulse
    (close_final_entries). A window that ends on the last pulse but is not the final one leaves it to the next, so
    contiguous windows partition the run. reduce_time_slices, which has closed its own final window, passes False.

    :return: (outputs, plots): one flat list of reduced data per entry, in order, and the kinetic plot (None without
        plot_time), drawn into figure when one is given (plot_kinetic)
    :raises ValueError: starts and ends of different lengths; or, once every other entry has been reduced and written,
        the entries that could not be reduced (the reduction raised, or returned no data), each with its window.
    '''

    if len(starts) != len(ends):
        raise ValueError("Length of starts and ends lists must be the same")

    nexus_path, Spath = run_folders(settings_file, experiment_id, datapath, savepath)
    run_list = [run]
    if close_final:
        ends = close_final_entries(starts, ends, *BP.read_run_end(nexus_path / f"REF_L_{run}.nxs.h5"))

    store_outputs = []
    mid_points = []
    failures = []
    # run the looped reduction
    for slice_idx in range(len(starts)):
        logger.info("Window %d: %s", slice_idx, window_text(starts[slice_idx], ends[slice_idx]))
        try:  # a window that is not valid is reported with the others below, as one that cannot be reduced
            span_start, span_end = window_span(starts[slice_idx], ends[slice_idx])
            if not subname_input:
                subname = f"slice_{int(span_start)}_{int(span_end)}"
            else:
                subname = f"{subname_input}_slice_{int(span_start)}_{int(span_end)}"

            override_params = {'Spath': Spath, "subname": subname}
            output = reduction.reduce_from_file(run_list, settings_file, experiment_id, datapath=nexus_path,
                                        override_params=override_params, plot=plot_ref, save_json=False,
                                        start_times=starts[slice_idx], end_times=ends[slice_idx])
        except (ValueError, RuntimeError) as error:  # this window only; reported with the others below
            failures.append(f"window {window_text(starts[slice_idx], ends[slice_idx])}: {error}")
            continue
        logger.info("Window %d reduced", slice_idx)

        all_results, _, _, _ = output
        flat_results = flatten_reduced_results(all_results)
        if not flat_results:
            failures.append(f"window {window_text(starts[slice_idx], ends[slice_idx])}: no reduced result data returned")
            continue

        mid_point = (span_end - span_start) / 2 + span_start
        mid_points.append(mid_point)
        store_outputs.append(flat_results)
        logger.debug("%d windows reduced so far", len(store_outputs))

    if failures:
        raise ValueError(f"No reduced result data for run {run} in " + "; ".join(failures))

    # create plot of set
    if plot_time:
        plots = plot_kinetic(store_outputs, run, times=mid_points, show=show_plots, figure=figure)
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


def plot_kinetic(output_list, run, times, show=True, figure=None):
    '''
    Plot an offset graph and a colour map of the slices' R, and return the figure.

    output_list is expected to be a list of per-slice result packs, where each pack is a list of reduced data dicts,
    e.g. [ {"Q":..., "R":..., "dR":..., "dQ":...}, ... ]. With a figure, the plot replaces its contents (the launcher's
    tab draws into its own canvas this way); without one, a new figure is made: through pyplot when show is true (a
    script's window), else outside pyplot, so that no figure is left registered there.
    '''

    # TODO: update the color plot to take time rather than slice index...

    if figure is not None:
        figure.clear()
        fig = figure
    elif show:
        fig = plt.figure(figsize=(15, 6))
    else:
        fig = Figure(figsize=(15, 6))
    ax = fig.subplots(nrows = 1, ncols = 2)
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
            ax[0].errorbar(Q_vals, R_vals * offset, yerr=dR * offset, xerr=dQ, fmt='o', markersize=1, label=f'mid_point: {times[ii]:.1f}')

            store_q.append(Q_vals)
            store_r.append(R_vals)
            store_dr.append(dR)
            store_dq.append(dQ)

    ax[0].set_ylabel('R')
    ax[0].set_xscale('log')
    ax[0].set_title(f'Time_slices for run {run}', fontsize=16)
    ax[0].set_yscale('log')
    ax[0].legend()
    Angstrom = '\u212B'
    ax[0].set_xlabel('Q [1/' + Angstrom + ']')

    # create the colour map. Store the data into a 2D array and plot with imshow
    if not store_r:
        raise ValueError("No reduced data available for kinetic plot")

    # The colour map is R, slice by slice, as its colour bar says (the contribution mapped dR under the label "R")
    Z = np.array([np.asarray(arr) for arr in store_r], dtype=float)
    mask = np.isfinite(Z) & (Z > 0)
    if not np.any(mask):
        raise ValueError("No positive R values available for kinetic plot")

    q_vals = np.asarray(store_q[0], dtype=float)
    y_vals = np.asarray(times, dtype=float)
    if y_vals.size != Z.shape[0]:
        raise ValueError(
            f"Time list length ({y_vals.size}) does not match the number of slices ({Z.shape[0]}). "
            "Pass one time value per slice."
        )

    # Convert the supplied midpoint times into full bin edges. The plotted data are
    # assigned to the row centers (mid_points), but the colour map should extend from
    # the previous edge to the next edge so it spans the entire slice widths.
    if y_vals.size == 1:
        y_edges = np.array([y_vals[0] - 0.5, y_vals[0] + 0.5], dtype=float)
    else:
        y_edges = np.empty(y_vals.size + 1, dtype=float)
        y_edges[1:-1] = 0.5 * (y_vals[:-1] + y_vals[1:])
        y_edges[0] = y_vals[0] - 0.5 * (y_vals[1] - y_vals[0])
        y_edges[-1] = y_vals[-1] + 0.5 * (y_vals[-1] - y_vals[-2])

    vmin, vmax = np.percentile(Z[mask], [2, 98])
    im = ax[1].imshow(
        Z,
        aspect='auto',
        origin='lower',
        extent=[q_vals.min(), q_vals.max(), y_edges[0], y_edges[-1]],
        norm=LogNorm(vmin=vmin, vmax=vmax),
        cmap='viridis',
    )

    fig.colorbar(im, ax=ax[1], label='R')
    ax[1].set_xlabel('Q [1/' + Angstrom + ']')
    ax[1].set_xscale('log')
    ax[1].set_ylabel('Time (s)')
    ax[1].set_title(f'Time_slices for run {run}', fontsize=16)

    fig.tight_layout()

    if show:
        plt.show()

    return fig


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
