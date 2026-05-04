import copy
import os
from pathlib import Path

import h5py
import lr_reduction.new_reduction_template_reader as reduction_template_reader
import lr_reduction.nr_tools as tools
import numpy as np
from matplotlib import pyplot as plt
from lr_reduction.new_reduction_template_reader import ReductionParameters

#import template
from lr_reduction.nr_reduction_calc import NR_Reduction  # TODO: Fix names of files!!
from lr_reduction.nr_reduction_config import NRReductionConfig
import lr_reduction.save_reduced_data as save_fn
import lr_reduction.new_reduction_from_file as nrff

def reduce_from_template(runno, template_file, experiment_id, datapath: Path = None, template_path: Path = None, override_params: dict = None, plot=True, eight_col=None):
    """
    Wrapper function to reduce a single run with reading of parameters from an xml template of the lr_reduction format.
    Then collect like results within the save folder and combine them together.
    This is intended for use with e.g. autoreduction.

    It loads the template and extracts values for that run, creates the reduction config and fills with parameters from the template
    file. It then overrides any parameters set within the override parameter dictionary to give full flexibility. It saves out the
    individual run in an equivalent format to prior autoreduction processes. Then loads and sorts any which have the same seq_id in
    the title and combines them into an output file. Equivalent to prior autoreduction processes.

    At the moment this goes between template file xml -> template dict -> config -> reduction -> new_template dict -> xml.
    Can be simplified in future.

    NOTE: new files for reading the template were made to not cause issues with prior setup but can be changed in future.

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

    # Get sequence number from file
    if not datapath:
        datapath = Path("/SNS/REF_L") / experiment_id
    fname = datapath / f"REF_L_{runno}.nxs.h5"
    f = h5py.File(fname, 'r')
    seq_num = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Num/value'][0]
    seq_id = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Id/value'][0]
    f.close()

    # Read the template file
    template_data = read_template(template_file, seq_num)
    # Apply template to config
    config = config_from_template(template_data)
    config.RBnum = [runno]
    config.experiment_id = experiment_id
    config.Sname = f"REFL_{seq_id}_{seq_num}_{runno}"
    config.plotON = plot

    if datapath:
        config.NEXUSpathRB = datapath

    # override the normalisation flag if isn't the first in the sequence.
    if seq_num != 1:
        config.Normalize = False

    # override the saving of 8 column if provided into the function. TODO: decide if needed.
    if eight_col:
        config.save8col = eight_col

    # override template with anything provided
    if override_params:
        for key, value in override_params.items():
            if hasattr(config, key):
                setattr(config, key, value)
            else:
                raise AttributeError(f"{key} is not a valid config parameter")

    #print(vars(config)) # Print to check the changes.
    #print(vars(template_data))

    # Run reduction
    reduce_calc = NR_Reduction(config)
    results, config_out, log_values = reduce_calc._reduce_single_run(i=0, rb_num=config.RBnum[0], save=True)

    if eight_col:
        # Save the single run output. TODO: this might need cleaning up between the different functions.
        result = {'Q': results['q'], 'R': results['r'], 'dR': results['dr'], 'dQ': results['dq'],
                'T': result['t'], 'L': result['l'], 'dT': result['dt'], 'dL': result['dl']}
    else:
        result = {'Q': results['q'], 'R': results['r'], 'dR': results['dr'], 'dQ': results['dq']}       

    # Save data
    if config.subname:
        partial_name = f"{config.subname}_partial"
    else:
        partial_name = "partial"

    save_fn.save_results(result, config_out, log_values, sname=f"{config.Sname}{partial_name}")
    if eight_col: #TODO: decide whether this is instead of prior save
        save_fn.save_results(result, config_out, log_values, sname=f"{config.Sname}{partial_name}", eight_column=True)

    # Collect "like" runs together
    seq_list, run_list, combined_results, scaling_factors, config_array = assemble_results(seq_id, config.Spath, autoscale=config.AutoScale, plot=plot, RQ4=config.plotQ4, eight_col=eight_col)
    # Add scaling factor to output
    scale_list = np.array([np.float64(1)] + scaling_factors)
    config.ScaleFactor *= scale_list

    # Save combined data
    if config.subname:
        combined_name = f"{config.subname}_combined_data"
    else:
        combined_name = "combined_data"

    save_fn.save_results(combined_results, config_out, log_values, sname=f"REFL_{seq_id}{combined_name}", full=True)
    if eight_col: #TODO: decide whether this is instead of prior save
        reduce_calc.save_results(combined_results, config_out, log_values, sname=f"REFL_{seq_id}{combined_name}", full=False, eight_column=True)

    # plot
    if plot:
        plot_reflectivity([combined_results], RQ4 = config.plotQ4)

    # Save new template. Including new flags for other aspects. #TODO: check if it covers all we need.
    if template_path is None:
        template_path = Path("/SNS/REF_L") / experiment_id / "shared"

    template_updated = template_to_config(config, template_data)
    # For now is set to save with a "new" appended file name to be separate from refred but this can be changed moving forward.
    template_save_name = f"REFL_{seq_id}_template_new.xml"
    prior_template = template_file
    if Path(template_path / template_save_name).exists():
        file_to_change = template_path / template_save_name # If later in sequence want to update the new one not the old refred one.
    else:
        file_to_change = template_path / prior_template
    write_template(seq_list, run_list, file_to_change, template_updated, seq_num, template_path, save_name=template_save_name, prior_template=prior_template)

    run_group = len(run_list)

    return combined_results, run_group


def config_from_template(template_data):
    """
    Create an NRReductionConfig from a template file.

    Parameters
    ----------
    template_data :

    Returns
    -------
    NRReductionConfig
        Configuration object ready for NR_Reduction

    """
    # NOTE: Various of the config items expect arrays so ensure is set as single item array here.

    # Determine reduction method from template
    if template_data.q_method != False and template_data.q_method is not None: # Should this be None or False? Not sure where the False comes from...
        method = template_data.q_method
    else:
        method = 'MeanTheta' if template_data.const_q else 'constantTOF' #TODO: decide if this should be MeanTheta or constantQ

    # Initialize configuration
    config = NRReductionConfig()
    config.method_per_run = [method]

    # Update other parameters from the template
    config.RB_Ymin = [template_data.data_peak_range[0]]
    config.RB_Ymax = [template_data.data_peak_range[1]]

    if template_data.subtract_background == True:
        config.useBS = [1]
    else:
        config.useBS = [0]
    config.BkgROI = [template_data.background_roi]

    config.tof_max = [template_data.tof_range[1]]    # TODO: check which tof one to use...
    config.tof_min = [template_data.tof_range[0]]    # TODO: check which tof one to use...
    config.data_x_range = template_data.data_x_range

    if template_data.lam_range is not None:
        config.LambdaMin = [template_data.lam_range[0]]
        config.LambdaMax = [template_data.lam_range[1]]

    config.qmin = template_data.q_min
    config.dqbin = template_data.q_step

    config.ThetaShift = [template_data.angle_offset]

    config.dead_time = template_data.dead_time_value
    config.dead_time_tof_step = template_data.dead_time_tof_step

    config.use_emission_time = template_data.use_emission_time

    # TODO: does the gravity direction part need adding? does the emission time use need adding?
    #       do the flags on instrument settings need to be added?

    # Update with new flags in template if they exist
    config.Normalize = getattr(template_data, "norm_scale", config.Normalize)
    config.DBname = [getattr(template_data, "DB_file", config.DBname)]
    config.AutoScale = getattr(template_data, "autoscale", config.AutoScale)
    config.useCalcTheta = getattr(template_data, "use_calc_theta", config.useCalcTheta)
    config.Qline_threshold = getattr(template_data, "qline_threshold", config.Qline_threshold)
    config.ScaleFactor = [getattr(template_data, "scale_factor", 1.0)]
    config.save8col = getattr(template_data, "save_eight_col", config.save8col)

    return config

def template_to_config(config_data, template_data):
    """
    Reverse of the config settings back into the template format. Needed whilst keeping the xml format.
    """
    template = template_data
    template.q_method = config_data.method_per_run[0]
    template.data_peak_range = [config_data.RB_Ymin[0], config_data.RB_Ymax[0]]
    if config_data.useBS[0] == 1:
        template.subtract_background = True
    else:
        template.subtract_background = False
    template.background_roi = config_data.BkgROI[0] # TODO: need to check what happens with the order of settings here.
    template.tof_range = [config_data.tof_min[0], config_data.tof_max[0]]
    template.data_x_range = config_data.data_x_range
    template.q_min = config_data.qmin
    template.q_step = config_data.dqbin
    template.angle_offset = config_data.ThetaShift[0]
    template.dead_time_value = config_data.dead_time
    template.dead_time_tof_step = config_data.dead_time_tof_step
    template.norm_scale = config_data.Normalize
    template.DB_file = config_data.DBname[0]
    template.autoscale = config_data.AutoScale
    template.use_calc_theta = config_data.useCalcTheta
    template.qline_threshold = config_data.Qline_threshold
    template.scale_factor = config_data.ScaleFactor[0]
    template.use_emission_time = config_data.use_emission_time
    template.save8col = config_data.save8col

    #template.source_detector_distance = config_data.source_detector_distance

    if (config_data.LambdaMin is not None) & (config_data.LambdaMax is not None):
        template.lam_range = [config_data.LambdaMin[0], config_data.LambdaMax[0]]

    return template


def assemble_results(seq_id, output_dir, autoscale = True, plot=True, RQ4=False, eight_col = False):
    """
    Assemble the results for any files in the saved directory that have the same sequence number.
    Finds files saved witht the "partial.dat" ending, sorts the data , autoscales if the flag applied,
    combines and returns the seq list, run list and combined data as a dict.

    seq_id: Sequence ID of set to combined
    output_dir: Path to save out combined file
    autoscale: Bool to apply scaling in R between settings, scales to first in series
    plot: Bool to plot combined output
    RQ4: Bool to plot as RQ4 instead of R.

    returns
        seq_list, run_list, combined results

    """

    # Keep track of sequence IDs and run numbers so we can make a new template
    seq_list = []
    run_list = []
    full_names = []

    # Find the files
    file_list = sorted(os.listdir(output_dir))
    print("Files found:", len(file_list))
    for item in file_list:
        if eight_col:
            search_flag = item.endswith("partial_8col.dat") #TODO: make sure it searches on the subname too.
        else:
            search_flag = item.endswith("partial.dat")
        if item.startswith("REFL_%s" % seq_id) and search_flag:
            toks = item.split("_")
            if not len(toks) == 5 or int(toks[2]) == 0: # TODO: check was using 3 as a check
                print('here')
                continue
            seq_list.append(int(toks[2]))
            run_list.append(int(toks[3]))
            full_names.append(item)

    # Load the data
    data_array = []
    config_array = []
    for idx, file in enumerate(full_names):
        seq_id_store = seq_list[idx]
        data = np.loadtxt(Path(output_dir) / file, unpack=True) # TODO: sort as Path.
        #output = nrff.load_from_file(Path(output_dir) / file)
        #data_load = output["data"]
        #config_load = output["config"]
        data_array.append(data)
        #config_array.append(config_load)
    print("Data loaded:", len(data_array))

    # Do a sort based on lowest q?? TODO: work out this sorting part... Hasn't been fully tested.
    to_sort = []
    for run in range(len(data_array)):
        first_q = data_array[run][0, 0]
        to_sort.append(first_q)
    sort_list = np.argsort(to_sort)

    sorted_data = [data_array[i] for i in sort_list]

    # TODO: add better autoscaling options. Make scaling a function in nr_tools?
    Q, R, dR, dQ = [], [], [], []
    T, L, dT, dL = [], [], [], []
    dict_output = []
    scaling_factors = []
    for run, result in enumerate(sorted_data):
        if autoscale and run != 0:
            mask1 = Q[run-1] >= min(result[0, :])
            mask2 = result[0, :] <= max(Q[run-1])

            y1 = R[run-1][mask1]
            y2 = result[1, :][mask2]
            e1 = dR[run-1][mask1]
            e2 = result[2, :][mask2]

            scale, sigma_scale = tools.weighted_mean(y1, y2, e1, e2)

            result[1, :] *= scale
            result[2, :] *= scale

            print('Scaling factor:', np.round(scale, 3))
            scaling_factors.append(scale)

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

    if len(Q) == 0:
        raise ValueError(f"No valid runs found for sequence {seq_id}")

    if plot:
        plot_reflectivity(dict_output, RQ4)

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

    return seq_list, run_list, combine_results, scaling_factors, config_array

def write_template(seq_list, run_list, file_to_change, template_data_updated, seq_updated, output_dir, save_name=None, prior_template=None):
    """
    Read the appropriate entry in a template file and save an updated
    copy with the updated run number.

    Parameters
    ----------
    seq_list : list
        The sequence identifiers
    run_list : list
        The run numbers
    file_to_change : Path
        Path to template to open and fill with new entry
    template_data_updated: template object
        updated settings to be added into the template file associated with the seq to update
    seq_updated : float
        seq num to be updated with the updated template settings
    output_dir : str
        Directory where the output files are saved
    save_name : bool (optional)
        ability to override the default template naming convention
    prior_template : Path (optional)
        path to prior template with greater entries in case file to change is too small. (This might need changing)
    """
    print("Reading file", file_to_change)
    with open(file_to_change, "r") as fd:
        xml_str = fd.read()
        # Read the template
        data_sets = reduction_template_reader.from_xml(xml_str)
        # Do a check on the length of the data_sets
        # This part should allow sequence to be handled out of order. #TODO: Needs fully checking. Small check complete.
        # TODO: This function is messy and needs sorting out!
        max_seq_num = np.max(seq_list)
        if len(data_sets) >= max_seq_num:
            data_sets = data_sets
        elif prior_template:
            # if prior template given, use the other entries from the prior template file
            with open(prior_template, "r") as pt:
                xml_str_prior = pt.read()
                data_sets = reduction_template_reader.from_xml(xml_str_prior)
                print("Using prior template data")
        else:
            print("Template might not have enough entries, will fill the ones it can.")

        # For each requested sequence number, produce an entry to save.
        # Only the sequence matching `seq_updated` should receive the full `template_data_updated`.
        to_save = []
        for i in range(len(seq_list)):
            seq_i = seq_list[i]
            run_i = run_list[i]
            if len(data_sets) >= seq_i:
                if seq_i == seq_updated:
                    # Use the updated template for this sequence
                    new_entry = copy.deepcopy(template_data_updated)
                    new_entry.data_files = [run_i]
                else:
                    # Preserve the original entry for this sequence but set the run number
                    new_entry = copy.deepcopy(data_sets[seq_i - 1])
                    new_entry.data_files = [run_i]
                to_save.append(new_entry)
            elif len(data_sets) == seq_i - 1:
                # extend the list by duplicating the last existing entry, then behave as above
                data_sets.append(copy.deepcopy(data_sets[-1]))
                if seq_i == seq_updated:
                    new_entry = copy.deepcopy(template_data_updated)
                    new_entry.data_files = [run_i]
                else:
                    new_entry = copy.deepcopy(data_sets[seq_i - 1])
                    new_entry.data_files = [run_i]
                to_save.append(new_entry)
            else:
                print("Too few entries [%s] in template for sequence number %s" % (len(data_sets), seq_i))

    if not save_name:
        save_name = "REF_L_%s_auto_template.xml" % run_list[0]

    xml_str = reduction_template_reader.to_xml(to_save)
    with open(os.path.join(output_dir, save_name), "w") as fd:
        fd.write(xml_str)
    print("Saving to file", save_name)

# TODO: Fix to use the one inside template.py.
def read_template(template_file: str, sequence_number: int) -> ReductionParameters:
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

def plot_reflectivity(data_array, RQ4=False, log_x = True):
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
    plt.show()
