import copy
import os
from pathlib import Path
import json
import re

import h5py
import numpy as np
from matplotlib import pyplot as plt

#import template
from lr_reduction.nr_reduction_calc import NR_Reduction  # TODO: Fix names of files!!
from lr_reduction.nr_reduction_config import NRReductionConfig
import lr_reduction.save_reduced_data as save_fn
import lr_reduction.nr_tools as tools

def reduce_from_file(run_array, setting_file, experiment_id, datapath: Path = None, override_params: dict = None, plot=True, save_json=False, check_for_prior=False, save_pdf_summary=False):
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
        logs_out = results["used_log_vals"]
        output_figures.extend([results["figures"][-1]]) # Choose to output only the overlapped plot of settings

        if check_for_prior:
            # Look in folder for files of correct format
            dict_output, combine_results, scaling_factors, matched_files, sorted_run_nums, angle_logs = find_combine_priors(config_final, group_output_sorted["run_nums"], results, group_output_sorted, eight_col)

            # check dictionaries and arrays aren't empty and pass if they are
            if not dict_output:
                continue
            if not combine_results:
                continue
            if not scaling_factors:
                continue
            if not matched_files:
                continue

            # update the config scaling factors
            config_final.ScaleFactor = scaling_factors

            # TODO: Need to read in the used_theta_vals
            #used_theta_vals = {"thi":[], "ths":[], "ThCen":[], "title": []}
            used_theta_vals = {k: angle_logs.get(k, []) + logs_out.get(k, []) for k in angle_logs.keys() | logs_out.keys()}
            # save files
            # non-concatenated
            # TODO: this is resaving them. Think this is the best option.
            for i in range(len(dict_output)):
                save_fn.save_results(dict_output[i], config_final, used_theta_vals, sname=f"{config_final.Sname}_{i+1}_{sorted_run_nums[i]}{config_final.subname}")
                if eight_col:
                    save_fn.save_results(dict_output[i], config_final, used_theta_vals, sname=f"{config_final.Sname}_{i+1}_{sorted_run_nums[i]}{config_final.subname}", eight_column=True)
            # Always make the final plot, only show it if plot is True
            new_plot = plot_reflectivity(dict_output, RQ4=False, show_fig=plot)
            figures_out.extend(new_plot)
            output_figures.extend(new_plot)
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

        all_results.append(results)
        all_figures.extend(figures_out)

        if save_json:
            filepath_out = Path(config_final.Spath / f"{config_final.Sname}_settings.json")
            with open(filepath_out, "w") as f:
                json.dump(
                    save_fn.make_json_safe(config_final.__dict__),
                    f,
                    indent=2
                    )

    num_figures = len(all_figures)

    # Might need to come back to which figures are output.
    return all_results, output_figures

def sort_runs(group_output, idx):
    seq_to_use = group_output["seq_nums"][idx]
    sort_seq = np.argsort(seq_to_use)

    group_output_new={}
    group_output_new["seq_nums"] = [group_output["seq_nums"][idx][i] for i in sort_seq]
    group_output_new["run_nums"] = [group_output["run_nums"][idx][i] for i in sort_seq]
    group_output_new["seq_ids"] = [group_output["seq_ids"][idx][i] for i in sort_seq]

    # Put None where the seq_id are missing
    full_range = list(range(1, max(group_output_new["seq_nums"]) + 1))
    # Map values to their positions
    map2 = dict(zip(group_output_new["seq_nums"], group_output_new["run_nums"]))
    map3 = dict(zip(group_output_new["seq_nums"], group_output_new["seq_ids"]))
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
                rf"^{re.escape(updated_config.Sname)}_(\d+)_(\d+){'_8col' if eight_col else ''}{re.escape(updated_config.subname)}.dat$"
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

    # Load, sort data order
    angle_logs_thi = []
    angle_logs_ths = []
    angle_logs_thcen = []
    title_log = []

    for item in matched_files:
        filepath = Path(updated_config.Spath) / item[0]
        with open(filepath, "r") as f:
            for line in f:
                if line.startswith("# Angles: "):
                    angles_out = json.loads(line[len("# Angles: "):])
                elif line.startswith("# Run Title:"):
                    title_out = json.loads(line[len("# Run Title: "):])

        if angles_out:
            angle_logs_ths.append(angles_out["THS"])
            angle_logs_thi.append(angles_out["THI"])
            angle_logs_thcen.append(angles_out["ThCen"])
        if title_out:
            title_log.append(title_out["title"])
        data = np.loadtxt(filepath, unpack=True)
        existing_data.append(data)
        loaded_seq_nums.append(item[1])
        loaded_run_nums.append(item[2])
    
    # sort
    indices = sorted(range(len(existing_data)), key=lambda i: existing_data[i][0][0])
    sorted_data = [existing_data[i] for i in indices]
    sorted_seq_num = [loaded_seq_nums[i] for i in indices]
    sorted_run_num = [loaded_run_nums[i] for i in indices]
    angle_logs_thcen = max(angle_logs_thcen, key=len) # THis is weird and messy and needs a fix but because it adds on in prior cycles...!!
    angle_logs_ths = max(angle_logs_ths, key=len)
    angle_logs_thi = max(angle_logs_thi, key=len)
    title_log = max(title_log, key=len)

    angle_logs = {"ths": angle_logs_ths, "thi": angle_logs_thi, "ThCen": angle_logs_thcen, "title": title_log}

    return sorted_data, sorted_seq_num, sorted_run_num, angle_logs

def find_combine_priors(updated_config, run_nums, results, group_output_sorted, eight_col=False):
    # Find the files
    matched_files = find_priors(updated_config, eight_col, run_nums)
    print(matched_files)

    initial_scalefactors = updated_config.ScaleFactor
    initial_seq = group_output_sorted["seq_nums"]
    initial_run_nums = group_output_sorted["run_nums"]

    if len(matched_files) > 0:
        
        print(f"Found {len(matched_files)} to combine")
        sorted_data, sorted_seq_num, sorted_run_num, angle_logs = load_prior_data(results, matched_files, updated_config, initial_seq, initial_run_nums)
        # TODO: quite a bit is a duplicate of in the calc file. Can be smarter here.
        Q, R, dR, dQ = [], [], [], []
        T, L, dT, dL = [], [], [], []
        dict_output = []
        scaling_factors = []
        for run, result in enumerate(sorted_data):
            if updated_config.AutoScale and run != 0:
                mask1 = Q[run-1] >= min(result[0, :])
                mask2 = result[0, :] <= max(Q[run-1])

                y1 = R[run-1][mask1]
                y2 = result[1, :][mask2]
                e1 = dR[run-1][mask1]
                e2 = result[2, :][mask2]

                scale, sigma_scale = tools.weighted_mean(y1, y2, e1, e2)

                if not np.isfinite(scale):
                    print(f"Unable to find scaling factor for {run}")
                    scale = 1

                result[1, :] *= scale
                result[2, :] *= scale

                print('Scaling factor:', np.round(scale, 3))
                scaling_factors.append(scale)
                position = sorted_seq_num[run] - 1
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
        return dict_output, combine_results, initial_scalefactors, matched_files, sorted_run_num, angle_logs
    
    else:
        # TODO: fix this...
        return {}, {}, [], [], [], {}

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
    with open(filepath, "w") as f:
        json.dump(json_to_config(config), f, indent=2)

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
