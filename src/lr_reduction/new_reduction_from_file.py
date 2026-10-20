import ast
import json
import os
import re
from pathlib import Path

import h5py
import numpy as np
from matplotlib import pyplot as plt

import lr_reduction.nr_tools as tools
import lr_reduction.save_reduced_data as save_fn

#import template
from lr_reduction.nr_reduction_calc import NR_Reduction  # TODO: Fix names of files!!
from lr_reduction.nr_reduction_config import NRReductionConfig

# Per-run logs recorded in the output headers ("Run Title" and "Angles")
LOG_KEYS = ("title", "ths", "thi", "ThCen")


def reduce_from_file(run_array, setting_file, experiment_id, datapath: Path = None, override_params: dict = None, plot=True,
                     save_json=False, check_for_prior=False, save_pdf_summary=False):
    """
    Wrapper function to reduce a single run with reading of parameters from the header of a file or a saved json file, instead of the xml.
    Then collect like results within the save folder and combine them together.
    This is intended for use with e.g. autoreduction.

    run_array: set of reflectivity run number, which will be grouped by seq_id
    setting_file: .dat file to read header or .json setting file, including Path
    experiment_id: str IPTS number which appends to datapath if this isn't provided (e.g. "IPTS-36119")
    datapath: Path optional override of location to look up NEXUS file
    override_params: dict  Dictionary of config settings to override the defaults in either the template reader or the NRReduction config defaults.
    plot: bool Toggle to plot outputs during reduction steps.
    save_json: optional to save out a separate json setting file
    check_for_prior: optional to look in the save folder for other runs with the same prefix to join together. Important for autoreduction.
    save_pdf_summary: optional save of plots as a pdf output

    returns
        combined_results: dict of Q, R, dR, dQ. This is assembled with anything else on same seq num. #TODO: check if need individual one returned too.
    """

    all_results = []
    all_figures = []
    group_output = group_runs(run_array, experiment_id, datapath)
    output_figures = []

    # validation on the settings file can be .json or .dat
    # TODO: add validation step on file type
    # Loop over the groups
    for idx in range(len(group_output["run_nums"])):
        print('--------------------------------------------')
        print(f"Beginning run set {group_output['run_nums'][idx]}")

        # reset config
        file_load = load_from_file(setting_file) # moving here to try fix scale factor issue
        config_new = json_to_config(file_load["config"])

        # TODO: Sort out the sequence num part to make sure it runs the right indices of the config
        # This assumes want to process in increasing seq order
        group_output_sorted = sort_runs(group_output, idx)

        # TODO: validate check that config is long enough?

        config_new.RBnum = group_output_sorted["run_nums"]
        config_new.experiment_id = experiment_id

        config_new.Sname = f"REFL_{group_output_sorted['seq_ids'][0]}"
        config_new.plotON = plot

        if datapath:
            config_new.NEXUSpathRB = datapath

        # override template with anything provided
        if override_params:
            for key, value in override_params.items():
                if hasattr(config_new, key):
                    setattr(config_new, key, value)
                else:
                    raise AttributeError(f"{key} is not a valid config parameter")

        eight_col = config_new.save8col

        # Run reduction
        reducer = NR_Reduction(config_new)
        results = reducer.reduce(eight_col=eight_col, plot=plot, save_pdf_summary=save_pdf_summary)

        config_final = results["config"]
        figures_out = results["figures"]
        output_figures.extend([results["figures"][-1]]) # Choose to output only the overlapped plot of settings

        if check_for_prior:
            # Look in folder for files of correct format
            (dict_output, combine_results, scaling_factors, matched_files, sorted_run_nums, sorted_seq_nums,
             position_run_nums, used_theta_vals) = find_combine_priors(config_final, group_output_sorted["run_nums"], results, group_output_sorted, eight_col)

            # check dictionaries and arrays aren't empty - they're empty if no priors found
            if dict_output:
                # update the config scaling factors
                config_final.ScaleFactor = scaling_factors
                # the headers describe the whole set, one entry per sequence position
                config_final.RBnum = position_run_nums

                # save files
                # non-concatenated
                # TODO: this is resaving them. Think this is the best option.
                for i in range(len(dict_output)):
                    sname = f"{config_final.Sname}_{sorted_seq_nums[i]}_{sorted_run_nums[i]}{config_final.subname}"
                    save_fn.save_results(dict_output[i], config_final, used_theta_vals, sname=sname)
                    if eight_col:
                        save_fn.save_results(dict_output[i], config_final, used_theta_vals, sname=sname, eight_column=True)
                # Always make the final plot, only show it if plot is True
                new_plot = plot_reflectivity(dict_output, RQ4=False, show_fig=plot)
                figures_out.append(new_plot)
                output_figures.append(new_plot)
                # concatenated
                try:
                    save_fn.save_results(combine_results, config_final, used_theta_vals, full=True, sname=f"{config_final.Sname}_combined{config_final.subname}")
                except KeyError as e:
                    print(f"Warning: combined results missing expected key {e}; skipping save_results for combined output")
                if eight_col:
                    try:
                        save_fn.save_results(combine_results, config_final, used_theta_vals, eight_column=True, full=True, sname=f"{config_final.Sname}_combined{config_final.subname}")
                    except KeyError as e:
                        print(f"Warning: combined results missing expected key {e}; skipping save_results (8col) for combined output")

                if save_pdf_summary and plot:
                    # Overwrite output with new plot if created
                    save_fn.save_plot_pdf_summary(config_final.Spath, f"{config_final.Sname}{config_final.subname}", figures_out)

                all_results.append(dict_output)
            else:
                # when no priors found
                # TODO: handling of other missing return parts from no priors found
                print("No prior runs found")
                for val in range(len(results["Q_per_run"])):
                    to_add = {"Q": results["Q_per_run"][val],
                              "R": results["R_per_run"][val],
                              "dR": results["dR_per_run"][val],
                              "dQ": results["dQ_per_run"][val]}

                    all_results.append(to_add)
                sorted_run_nums = group_output['run_nums'][idx]
                all_results = [all_results]
        else:
            # TODO: This should be made consistent with the searching for prior return too
            all_results.append(results)
            sorted_run_nums = group_output['run_nums'][idx]
        all_figures.extend(figures_out)

        if save_json:
            filepath_out = Path(config_final.Spath / f"{config_final.Sname}_settings.json")
            # Through the helper, not a second inline copy. The duplicate is how
            # the public saver stayed broken and unnoticed: the flow that
            # actually saves never exercised it.
            save_config_json(filepath_out, config_final)

    # Might need to come back to which figures are output.
    return all_results, output_figures, sorted_run_nums, config_final

def highest_run(run_nums):
    """The run whose data and header entries a sequence position takes when it was measured more than once.

    D-5 (``[human, 2026-10-02]``): the highest run number wins, whatever the order in which the runs
    were reduced or named. Compared as integers. File names give ``int`` from the regex and callers
    pass ``int``, but text order ranks ``"99999"`` above ``"100000"``.
    """
    return max(run_nums, key=int)


def sort_runs(group_output, idx):
    seq_to_use = group_output["seq_nums"][idx]
    sort_seq = np.argsort(seq_to_use)

    group_output_new={}
    group_output_new["seq_nums"] = [group_output["seq_nums"][idx][i] for i in sort_seq]
    group_output_new["run_nums"] = [group_output["run_nums"][idx][i] for i in sort_seq]
    group_output_new["seq_ids"] = [group_output["seq_ids"][idx][i] for i in sort_seq]

    # Put None where the seq_id are missing
    full_range = list(range(1, max(group_output_new["seq_nums"]) + 1))
    # Map values to their positions. A call that names two runs of one position reduces only the
    # highest (D-5). A plain dict(zip(...)) kept whichever came last, so the argument order decided
    # and the other run was dropped without a word.
    by_position = {}
    for seq, run, seq_id in zip(group_output_new["seq_nums"], group_output_new["run_nums"],
                                group_output_new["seq_ids"]):
        by_position.setdefault(seq, []).append((run, seq_id))
    map2, map3 = {}, {}
    for seq, entries in by_position.items():
        runs = [run for run, _ in entries]
        map2[seq] = highest_run(runs)
        map3[seq] = next(seq_id for run, seq_id in entries if run == map2[seq])
        if len(set(runs)) > 1:
            print(f"Warning: this call names runs {sorted(set(runs), key=int)} at sequence position {seq}; "
                  f"reducing run {map2[seq]} only")
    # Build aligned lists
    base_aligned  = [x if x in group_output_new["seq_nums"] else None for x in full_range]
    list2_aligned = [map2.get(x, None) for x in full_range]
    list3_aligned = [map3.get(x, group_output_new["seq_ids"][-1]) for x in full_range]

    group_output_sorted = {}
    group_output_sorted = {"run_nums": list2_aligned,
                            "seq_nums": base_aligned,
                            "seq_ids": list3_aligned
                            }

    return group_output_sorted

def group_runs(run_array, experiment_id, datapath):
    current_seq_id = None

    # loop over array of runs and join together where seq_id changes
    run_nums = []
    seq_nums = []
    seq_ids = []
    for idx, runno in enumerate(run_array):
        # Get sequence number from file
        if not datapath:
            datapath = Path("/SNS/REF_L") / experiment_id / "nexus"
        fname = datapath / f"REF_L_{runno}.nxs.h5"
        f = h5py.File(fname, 'r')
        seq_num = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Num/value'][0]
        seq_id = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Id/value'][0]
        f.close()

        if seq_id != current_seq_id:
            if idx != 0:
                run_nums.append(run_num_subset)
                seq_nums.append(seq_num_subset)
                seq_ids.append(seq_id_subset)
            # clear the current arrays
            run_num_subset = []
            seq_num_subset = []
            seq_id_subset = []
            # add the first value
            run_num_subset.append(runno)
            seq_num_subset.append(seq_num)
            seq_id_subset.append(seq_id)
            current_seq_id = seq_id
        else:
            run_num_subset.append(runno)
            seq_num_subset.append(seq_num)
            seq_id_subset.append(seq_id)

            current_seq_id = seq_id

    # add the last ones
    run_nums.append(run_num_subset)
    seq_nums.append(seq_num_subset)
    seq_ids.append(seq_id_subset)

    # TODO work out if need to repeat the id etc.
    return {"run_nums": run_nums, "seq_nums": seq_nums, "seq_ids": seq_ids}

def find_priors(updated_config, eight_col, run_nums):
    pattern = re.compile(
                rf"^{re.escape(updated_config.Sname)}_(\d+)_(\d+){re.escape(updated_config.subname)}{'_8col' if eight_col else ''}.dat$"
                    )
    #print(f"Looking for files with pattern: {pattern.pattern} in {updated_config.Spath}")
    file_list = sorted(os.listdir(updated_config.Spath))
    matched_files = []
    for item in file_list:
        m = pattern.match(item)
        if m:
            num1, num2 = map(int, m.groups())
            # ignore any already in the currently processed set
            if num2 not in run_nums:
                matched_files.append((item, num1, num2))
    print("New files found:", len(matched_files)) # only finds those not re-processed.

    return matched_files

def load_prior_data(results, matched_files, updated_config, initial_seq, initial_run_nums):
    existing_data = []
    loaded_seq_nums = []
    loaded_run_nums = []
    #TODO: This needs cleaning up!

    # Logs of the current reduction, by sequence position like initial_seq (None where no run)
    current_logs = results["used_log_vals"]
    loaded_logs = [{k: current_logs[k][pos] for k in LOG_KEYS}
                   for pos, seq in enumerate(initial_seq) if seq is not None]

    # Remove None entries from initial set
    initial_seq = [x for x in initial_seq if x is not None]
    initial_run_nums = [x for x in initial_run_nums if x is not None]

    for val in range(len(results["Q_per_run"])):

        # arrange to columns
        if updated_config.save8col:
            data_to_add = np.vstack((results["Q_per_run"][val],
                                            results["R_per_run"][val],
                                            results["dR_per_run"][val],
                                            results["dQ_per_run"][val],
                                            results["L_per_run"][val],
                                            results["dL_per_run"][val],
                                            results["T_per_run"][val],
                                            results["dT_per_run"][val]
                                            ))
        else:
            data_to_add = np.vstack((results["Q_per_run"][val],
                                results["R_per_run"][val],
                                results["dR_per_run"][val],
                                results["dQ_per_run"][val]
                                ))

        existing_data.append(data_to_add)
        loaded_seq_nums.append(initial_seq[val])
        loaded_run_nums.append(initial_run_nums[val])

    # The candidates for each sequence position: this call's run there (no file name) and every prior
    # file found for it. Which candidate a position takes is decided before any prior file is read.
    candidates = {}
    for seq_num, run_num in zip(loaded_seq_nums, loaded_run_nums):
        candidates.setdefault(seq_num, []).append((run_num, None))
    for filename, seq_num, run_num in matched_files:
        candidates.setdefault(seq_num, []).append((run_num, filename))

    chosen = {}  # sequence position -> (run, prior file name or None for this call's run)
    for seq_num, position in sorted(candidates.items()):
        runs = sorted({run_num for run_num, _ in position}, key=int)
        if len(runs) > 1:
            # A prior file must belong to the position it names: the pre-fix merge named files by Q
            # order, and with the highest run winning, such a copy would displace the true run. Checked
            # only where runs compete, so a position with one run reads no NeXus file. A file whose
            # NeXus file cannot be read cannot be disproved and stays.
            belongs = []
            for run_num, filename in position:
                if filename is None:
                    belongs.append((run_num, filename))
                    continue
                nexus = Path(updated_config.NEXUSpathRB) / f"REF_L_{run_num}.nxs.h5"
                try:
                    with h5py.File(nexus, "r") as f:
                        true_seq = int(f["entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Num/value"][0])
                except (OSError, KeyError, IndexError, ValueError, TypeError) as e:
                    print(f"Warning: {filename} could not be checked against the NeXus file of run {run_num} "
                          f"({type(e).__name__}: {e}); kept for sequence position {seq_num}")
                    belongs.append((run_num, filename))
                    continue
                if true_seq == seq_num:
                    belongs.append((run_num, filename))
                else:
                    print(f"Warning: {filename} is not used: the NeXus file of run {run_num} puts it at "
                          f"sequence position {true_seq}, not {seq_num}")
            position = belongs
        if position:
            # D-5: the highest run number wins, whether its data come from this call or from a file
            winner = highest_run([run_num for run_num, _ in position])
            chosen[seq_num] = next(c for c in position if c[0] == winner)
        if len(runs) > 1:
            used = f"using run {chosen[seq_num][0]}" if seq_num in chosen else "none of them belongs to it"
            print(f"Warning: sequence position {seq_num} has files for runs {runs}; {used}")

    # Join by sequence position: the winner's data and logs (a position with no run stays None)
    highest_seq_num = max(chosen, default=0)
    print(highest_seq_num)
    combined_data = [None] * highest_seq_num
    combined_seq_nums = [None] * highest_seq_num
    combined_run_nums = [None] * highest_seq_num
    combined_logs = [None] * highest_seq_num
    for seq_num, (run_num, filename) in chosen.items():
        if filename is None:
            idx = loaded_seq_nums.index(seq_num)
            data, logs = existing_data[idx], loaded_logs[idx]
        else:
            filepath = Path(updated_config.Spath) / filename
            logs = own_logs_from_header(read_prior_header(filepath), seq_num, run_num)
            if logs is None:
                print(f"{filename}: header does not identify the logs of run {run_num}; reading title and angles "
                      f"from its NeXus file, ThCen unknown until run {run_num} is reduced again")
                logs = read_logs_from_nexus(run_num, updated_config.NEXUSpathRB)
            data = np.loadtxt(filepath, unpack=True)
        combined_data[seq_num - 1] = data
        combined_seq_nums[seq_num - 1] = seq_num
        combined_run_nums[seq_num - 1] = run_num
        combined_logs[seq_num - 1] = logs

    # sort the positions that have data (a gap in the sequence stays None)
    present = [i for i in range(len(combined_data)) if combined_data[i] is not None]
    indices = sorted(present, key=lambda i: combined_data[i][0][0])
    sorted_data = [combined_data[i] for i in indices]
    sorted_seq_num = [combined_seq_nums[i] for i in indices]
    sorted_run_num = [combined_run_nums[i] for i in indices]

    # logs stay by sequence position, like the config lists they are written next to
    position_logs = {k: [None if logs is None else logs[k] for logs in combined_logs] for k in LOG_KEYS}

    return sorted_data, sorted_seq_num, sorted_run_num, combined_run_nums, position_logs


def read_prior_header(filepath):
    """
    Read the run list, format marker and per-run logs from the header of a saved partial file.
    Missing or unreadable entries are None.
    """
    header = {"format": None, "NR_runs": None, **{k: None for k in LOG_KEYS}}
    with open(filepath, "r") as f:
        for line in f:
            if not line.startswith("#"):
                break
            try:
                if line.startswith("# NR_runs = "):
                    header["NR_runs"] = ast.literal_eval(line[len("# NR_runs = "):].strip())
                elif line.startswith("# Run Title: "):
                    header["title"] = json.loads(line[len("# Run Title: "):])["title"]
                elif line.startswith("# Angles: "):
                    angles = json.loads(line[len("# Angles: "):])
                    header["ths"], header["thi"], header["ThCen"] = angles["THS"], angles["THI"], angles["ThCen"]
                elif line.startswith("# Header format: "):
                    header["format"] = int(line[len("# Header format: "):].split()[0])
            except (ValueError, SyntaxError, KeyError, TypeError):
                pass  # an unreadable entry stays None, so the header cannot vouch for its logs
    return header


def own_logs_from_header(header, seq_num, run_num):
    """
    The logs of run `run_num` at sequence position `seq_num`, from the header of its own saved file,
    or None when the header cannot vouch for them.

    Header format 2 lists are indexed by sequence position. Older headers ("legacy") hold lists
    appended in reduction order, and NR_runs holds the runs of the call that wrote them. There the
    file's entry is known only when that call reduced this run alone (NR_runs has one run, at the
    file's position): it appended its own entry last, to the longest list among the other files
    (or to nothing, when there were none).
    """
    runs = header["NR_runs"]
    lists = [header[k] for k in LOG_KEYS]
    if not isinstance(runs, list) or not all(isinstance(values, list) and values for values in lists):
        return None
    if len({len(values) for values in lists}) != 1 or len(runs) < seq_num or runs[seq_num - 1] != run_num:
        return None
    if header["format"] is None:
        if [r for r in runs if r is not None] != [run_num]:
            return None
        return {k: header[k][-1] for k in LOG_KEYS}
    if len(lists[0]) != len(runs):
        return None
    return {k: header[k][seq_num - 1] for k in LOG_KEYS}


def read_logs_from_nexus(run_num, datapath):
    """
    Title, ths and thi of a run, read from its NeXus file as NR_Reduction.reduce() records them.
    ThCen is None: it is a result of the reduction (theta shift or calculated theta), not a log.
    """
    fname = Path(datapath) / f"REF_L_{run_num}.nxs.h5"
    try:
        with h5py.File(fname, "r") as f:
            return {"title": f["entry/title"].asstr()[0],
                    "ths": float(np.round(f["entry/DASlogs/BL4B:Mot:ths.RBV/value"][-1], 3)),
                    "thi": float(np.round(f["entry/DASlogs/BL4B:Mot:thi.RBV/value"][-1], 3)),
                    "ThCen": None}
    except (OSError, KeyError) as e:
        print(f"Could not read the logs of run {run_num} from {fname}: {e}")
        return {k: None for k in LOG_KEYS}

def find_combine_priors(updated_config, run_nums, results, group_output_sorted, eight_col=False):
    # Find the files
    matched_files = find_priors(updated_config, eight_col, run_nums)
    print(matched_files)

    initial_scalefactors = updated_config.ScaleFactor
    initial_seq = group_output_sorted["seq_nums"]
    initial_run_nums = group_output_sorted["run_nums"]

    if len(matched_files) > 0:

        print(f"Found {len(matched_files)} to combine")
        sorted_data, sorted_seq_num, sorted_run_num, position_run_nums, position_logs = load_prior_data(results, matched_files, updated_config, initial_seq, initial_run_nums)
        # TODO: quite a bit is a duplicate of in the calc file. Can be smarter here.
        Q, R, dR, dQ = [], [], [], []
        T, L, dT, dL = [], [], [], []
        dict_output = []
        scaling_factors = []
        print(sorted_seq_num, sorted_run_num, len(sorted_data))
        for run, result in enumerate(sorted_data):
            if updated_config.AutoScale and run != 0:
                #mask1 = Q[run-1] >= min(result[0, :])
                #mask2 = result[0, :] <= max(Q[run-1])

                #y1 = R[run-1][mask1]
                #y2 = result[1, :][mask2]
                #e1 = dR[run-1][mask1]
                #e2 = result[2, :][mask2]

                #scale, sigma_scale = tools.weighted_mean(y1, y2, e1, e2)

                #TODO: the function can take a series of arrays rather than stepwise.
                y1=R[run-1]
                y2=result[1, :]
                e1=dR[run-1]
                e2=result[2, :]

                scaling_factors_out = tools.calc_scaling_factors(xarrays = [Q[run-1], result[0,:]], yarrays = [y1, y2], yerrarrays=[e1, e2], method="weighted_mean", sigma_mask=3, set_first_to_one=False)
                scale = scaling_factors_out[0]

                if not np.isfinite(scale):
                    print(f"Unable to find scaling factor for run {sorted_run_num[run]}")
                    scale = 1

                result[1, :] *= scale
                result[2, :] *= scale

                print('Scaling factor:', np.round(scale, 3))
                scaling_factors.append(scale)
                position = sorted_seq_num[run] - 1
                while position >= len(initial_scalefactors):
                    initial_scalefactors.append(1)
                initial_scalefactors[position] *= scale

            Q.append(result[0, :])
            R.append(result[1, :])
            dR.append(result[2, :])
            dQ.append(result[3, :])
            if eight_col:
                L.append(result[4, :])
                dL.append(result[5, :])
                T.append(result[6, :])
                dT.append(result[7, :])

                dict_output.append({'Q': result[0,:], 'R': result[1,:], 'dR': result[2,:], 'dQ': result[3,:],
                                    'L': result[4,:], 'dL': result[5,:], 'T': result[6,:], 'dT': result[7,:]})
            else:
                # This is a bit muddled with a few things in arrays and dict. TODO: clean-up so don't need both.
                dict_output.append({'Q': result[0,:], 'R': result[1,:], 'dR': result[2,:], 'dQ': result[3,:]})

        # Combine results for all settings
        Q_combined = np.concatenate(Q)
        R_combined = np.concatenate(R)
        dR_combined = np.concatenate(dR)
        dQ_combined = np.concatenate(dQ)
        if eight_col:
            T_combined = np.concatenate(T)
            dT_combined = np.concatenate(dT)
            L_combined = np.concatenate(L)
            dL_combined = np.concatenate(dL)

        # Sort by Q for combined data
        idx = np.argsort(Q_combined)
        if eight_col:
            combine_results = {'Q': Q_combined[idx], 'R': R_combined[idx], 'dR': dR_combined[idx], 'dQ': dQ_combined[idx],
                            'L': L_combined[idx], 'dL': dL_combined[idx], 'T': T_combined[idx], 'dT': dT_combined[idx]}
        else:
            combine_results = {'Q': Q_combined[idx], 'R': R_combined[idx], 'dR': dR_combined[idx], 'dQ': dQ_combined[idx]}

        #print(combine_results)
        return (dict_output, combine_results, initial_scalefactors, matched_files, sorted_run_num, sorted_seq_num,
                position_run_nums, position_logs)

    else:
        # TODO: fix this...
        return {}, {}, [], [], [], [], [], {}

def json_to_config(json_input):
    config_init = NRReductionConfig()

    # override provided json input
    for key, value in json_input.items():
        if hasattr(config_init, key):
            setattr(config_init, key, value)
        else:
            raise AttributeError(f"{key} is not a valid config parameter")

    return config_init

def load_from_file(filepath):
    filepath = Path(filepath)

    if filepath.suffix == ".json":
        # JSON: config only (no data unless you choose to store it there)
        with open(filepath, "r") as f:
            config = json.load(f)
        data = None

    elif filepath.suffix == ".dat":
        config = None

        with open(filepath, "r") as f:
            for line in f:
                if line.startswith("# Config:"):
                    config = json.loads(line[len("# Config: "):])
                    break

            data = np.loadtxt(f, unpack=True)

    else:
        raise ValueError(f"Unsupported file type: {filepath.suffix}")

    return {"data": data, "config": config}

def save_config_json(filepath, config):
    """Write ``config`` as a settings JSON that :func:`json_to_config` reads back.

    This called `json_to_config` — the LOADER — where it needed the saver, so
    there was no input for which it worked: a dict raised ``TypeError`` (the
    config the loader builds is not serializable) and a config raised
    ``AttributeError`` (no ``.items()``). Nothing called it, because the working
    saver was three lines inlined in the reduction flow; the codebase carried a
    correct saver that could not be reused and a reusable saver that was broken,
    and the broken one was the public name a newcomer reaches for.

    ``__dict__`` is the right source rather than an implementation-detail
    shortcut, though not for the reason first given here: ``Spath`` **does**
    have a setter, so the public form would survive the ``hasattr`` gate. The
    real reasons are that the private ``_*_override`` names **re-derive** their
    paths from ``experiment_id`` on load, whereas the public form freezes the
    resolved absolute path into the file; and that ``base_path`` is a property
    with **no** setter, so any "save every public attribute" scheme hits
    ``AttributeError: property 'base_path' has no setter`` on reload.
    """
    # Serialize BEFORE opening. `json.dump` raising part-way through left a
    # truncated file behind, so a failed save destroyed the previous settings.
    payload = json.dumps(save_fn.make_json_safe(config.__dict__), indent=2)
    with open(filepath, "w") as f:
        f.write(payload)

def plot_reflectivity(data_array, RQ4=False, log_x = True, show_fig=True):
    """
    Plot reflectivity assuming data_array is an array of dictionaries with keys of Q, R, dQ and dR.
    # TODO: Add error handling for not being an array or a proper dictionary.

    :param data_array: Description
    :param RQ4: Description
    :param log_x: Description
    """
    fig, ax = plt.subplots()

    for vals in data_array:
        Q = vals["Q"]
        R = vals["R"]
        dQ = vals["dQ"]
        dR = vals["dR"]
        if RQ4:
            ax.errorbar(Q,R*Q**4, yerr=dR*Q**4, xerr=dQ,  fmt='o', markersize=1)
            ax.set_ylabel(r'$R \cdot Q^4$', fontsize=14)
        else:
            ax.errorbar(Q,R, yerr=dR, xerr=dQ,  fmt='o', markersize=1)
            ax.set_ylabel('R', fontsize=14)
        if log_x:
            ax.set_xscale('log')
        ax.set_yscale('log')
    Angstrom = '\u212B'
    ax.set_xlabel('Q [1/'+Angstrom+']', fontsize=14)
    ax.set_title('NR data') # TODO: add better title handling
    if show_fig:
        plt.show()
    else:
        plt.close()

    return fig
