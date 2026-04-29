import copy
import os
from pathlib import Path
import json

import h5py
import numpy as np
from matplotlib import pyplot as plt

#import template
from lr_reduction.nr_reduction_calc import NR_Reduction  # TODO: Fix names of files!!
from lr_reduction.nr_reduction_config import NRReductionConfig
import lr_reduction.save_reduced_data as save_fn

def reduce_from_file(run_array, setting_file, experiment_id, datapath: Path = None, override_params: dict = None, plot=True, eight_col=None, save_json=False):
    """
    Wrapper function to reduce a single run with reading of parameters from the header of a file or a saved json file, instead of the xml.
    Then collect like results within the save folder and combine them together.
    This is intended for use with e.g. autoreduction.

    It loads the template and extracts values for that run, creates the reduction config and fills with parameters from the template
    file. It then overrides any parameters set within the override parameter dictionary to give full flexibility. It saves out the
    individual run in an equivalent format to prior autoreduction processes. Then loads and sorts any which have the same seq_id in
    the title and combines them into an output file. Equivalent to prior autoreduction processes.

    runno: reflectivity run number
    template_file: template file, including Path
    experiment_id: str IPTS number which appends to datapath if this isn't provided (e.g. "IPTS-36119")
    datapath: Path optional override of location to look up NEXUS file
    template_path: Path optional override of template location. Otherwise uses IPTS shared folder
    override_params: dict  Dictionary of config settings to override the defaults in either the template reader or the NRReduction config defaults.
    plot: bool Toggle to plot outputs during reduction steps.

    returns
        combined_results: dict of Q, R, dR, dQ. This is assembled with anything else on same seq num. #TODO: check if need individual one returned too.
    """

    all_results = []
    group_output = group_runs(run_array, experiment_id, datapath)

    # validation on the settings file can be .json or .dat
    # TODO: add validation step on file type
    # Loop over the groups
    for idx in range(len(group_output["run_nums"])):
        print('--------------------------------------------')
        print(f"Beginning run set {group_output['run_nums'][idx]}")

        # reset config
        file_load = load_from_file(setting_file) # moving here to try fix scale factor issue
        config_new = json_to_config(file_load["config"])
        print(config_new.ScaleFactor)
        config_new.RBnum = group_output["run_nums"][idx]
        config_new.experiment_id = experiment_id

        # TODO: work out the saving, does it need to match prior behavior?
        config_new.Sname = f"REFL_{group_output['seq_ids'][idx][0]}"
        print(config_new.Sname)
        config_new.plotON = plot

        if datapath:
            config_new.NEXUSpathRB = datapath

        # override the saving of 8 column if provided into the function. TODO: decide if needed.
        if eight_col:
            config_new.save8col = eight_col

        # override template with anything provided
        if override_params:
            for key, value in override_params.items():
                if hasattr(config_new, key):
                    setattr(config_new, key, value)
                else:
                    raise AttributeError(f"{key} is not a valid config parameter")

        #print(vars(config)) # Print to check the changes.
        #print(vars(template_data))

        # Run reduction
        reducer = NR_Reduction(config_new)
        results = reducer.reduce(eight_col=True)

        all_results.append(results)

        #TODO: decide about saving json etc.
        if save_json:
            config_final = results["config"]
            filepath_out = Path(config_final.Spath / f"{config_final.Sname}_settings.json")
            with open(filepath_out, "w") as f:
                json.dump(
                    save_fn.make_json_safe(config_final.__dict__),
                    f,
                    indent=2
                    )

    return all_results

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

# Work out if need to restore any file types back
def restore_config_types(cfg):
    if "ScaleFactor" in cfg:
        cfg["ScaleFactor"] = np.array(cfg["ScaleFactor"])
    return cfg

def save_config_json(filepath, config):
    with open(filepath, "w") as f:
        json.dump(json_to_config(config), f, indent=2)