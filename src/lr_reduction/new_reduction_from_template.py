import h5py
#import template
from nr_reduction_unified import NR_Reduction  # TODO: Fix names of files!!
from nr_reduction_config import NRReductionConfig
from pathlib import Path
import os
import numpy as np
import nr_tools as tools
from new_reduction_template_reader import ReductionParameters
import new_reduction_template_reader as reduction_template_reader
from matplotlib import pyplot as plt


def reduce_from_template(runno, template_file, experiment_id, datapath: Path = None, override_params: dict = None, plot=True):
    """
    Wrapper function to reduce a single run with reading of parameters from an xml template of the lr_reduction format.
    Then collect like results within the save folder and combine them together.
    This is intended for use with e.g. autoreduction.

    It loads the template and extracts values for that run, creates the reduction config and fills with parameters from the template
    file. It then overrides any parameters set within the override parameter dictionary to give full flexibility. It saves out the
    individual run in an equivalent format to prior autoreduction processes. Then loads and sorts any which have the same seq_id in
    the title and combines them into an output file. Equivalent to prior autoreduction processes.

    NOTE: new files for reading the template were made to separate from Mantid due to a typing.py bug. Once rectified might be better
    to move back to existing reader.
    
    runno: reflectivity run number
    template_file: template file, including Path
    experiment_id: str IPTS number which appends to datapath if this isn't provided (e.g. "IPTS-36119")
    datapath: Path optional override of location to look up NEXUS file
    override_params: dict  Dictionary of config settings to override the defaults in either the template reader or the NRReduction config defaults.
    plot: bool Toggle to plot outputs during reduction steps.
    """
    # NOTE: DBname will need to be in override_params for now.

    # Get sequence number from file
    if not datapath:
        datapath = Path("/SNS/REF_L") / experiment_id
    fname = datapath / f"REF_L_{runno}.nxs.h5"
    f = h5py.File(fname, 'r')
    seq_num = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Num/value'][0]
    seq_id = f['entry/DASlogs/BL4B:CS:Autoreduce:Sequence:Id/value'][0]
    f.close()

    # Read the template file
    # TODO: decide if this is a good reader or needs some altering. Does a load of mantid stuff!
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

    # override template with anything provided
    if override_params:
        for key, value in override_params.items():
            if hasattr(config, key):
                setattr(config, key, value)
            else:
                raise AttributeError(f"{key} is not a valid config parameter")

    #print(vars(config)) # Print to check the changes.

    # Run reduction
    reduce_calc = NR_Reduction(config)
    results = reduce_calc._reduce_single_run(i=0, rb_num=config.RBnum[0], save=True)

    # Save the single run output. TODO: this might need cleaning up between the different functions.
    result = {'Q': results['q'], 'R': results['r'], 'dR': results['dr'], 'dQ': results['dq']}
    reduce_calc.save_results(result, sname=f"{config.Sname}_partial")

    # Collect "like" runs together 
    seq_list, run_list, combined_results = assemble_results(seq_id, config.Spath, autoscale=config.AutoScale, plot=plot, RQ4=config.plotQ4)
    # Save combined data
    reduce_calc.save_results(combined_results, sname=f"REFL_{seq_id}_combined_data")

    # plot
    if plot:
        plot_reflectivity([combined_results], RQ4 = config.plotQ4)

    # Save new template. Including new flags for other aspects?!
    # Normalize, Autoscale, useCalcTheta, Qline_threshold, ScaleFactor, DetResFn, DetSigma, Thetamethod, DBname
    # If add these into the template also need to add them to the config_from_template.
    template_save_name = f"REFL_{seq_id}_template_new.xml"
    write_template(seq_list, run_list, template_file, datapath, save_name=template_save_name)
    # TODO: FIX THE template save path...

    # TODO: Issue is it doesn't populate the new fields correctly...

    return combined_results


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
    # const_q=True -> MeanTheta; const_q=False -> constantTOF
    # TODO: Add better/different flag for method. Expect other paramters needed elsewhere too...
    method = 'MeanTheta' if template_data.const_q else 'constantTOF'
    # attempt 1
    if template_data.q_method is not None:
        method = template_data.q_method
    else:
        method = 'MeanTheta' if template_data.const_q else 'constantTOF'
    
    # Initialize configuration
    config = NRReductionConfig(method=method)
    
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

    config.qmin = template_data.q_min
    config.dqbin = template_data.q_step

    config.ThetaShift = [template_data.angle_offset]

    config.dead_time = template_data.dead_time_value
    config.dead_time_tof_step = template_data.dead_time_tof_step

    # TODO: does the gravity direction part need adding? does the emission time use need adding?
    #       do the flags on instrument settings need to be added?

    # TODO: Update with new flags in template if they work!
    
    config.Normalize = getattr(template_data, "norm_scale", config.Normalize)
    config.DBname = getattr(template_data, "DB_file", config.DBname)
    config.AutoScale = getattr(template_data, "autoscale", config.AutoScale)
    config.useCalcTheta = getattr(template_data, "use_cal_theta", config.useCalcTheta)
    config.Qline_threshold = getattr(template_data, "qline_threshold", config.Qline_threshold)
    config.ScaleFactor = [getattr(template_data, "scale_factor", 1.0)]
    
    return config

def assemble_results(seq_id, output_dir, autoscale = True, plot=True, RQ4=False):
 
    # Keep track of sequence IDs and run numbers so we can make a new template
    seq_list = []
    run_list = []
    full_names = []

    # Find the files
    file_list = sorted(os.listdir(output_dir))
    print("Files found:", len(full_names))
    for item in file_list:
        if item.startswith("REFL_%s" % seq_id) and item.endswith("partial.dat"):
            toks = item.split("_")
            if not len(toks) == 5 or int(toks[2]) == 0:
                continue
            seq_list.append(int(toks[2]))
            run_list.append(int(toks[3]))
            full_names.append(item)

    # Load the data
    data_array = []
    for idx, file in enumerate(full_names):
        seq_id_store = seq_list[idx]
        data = np.loadtxt(Path(output_dir) / file, unpack=True) # TODO: sort as Path.
        data_array.append(data)
    print("Data loaded:", len(data_array))

    # Do a sort based on lowest q?? TODO: work out this sorting part...
    to_sort = []
    for run in range(len(data_array)):
        first_q = data_array[run][0,0]
        to_sort.append(first_q)
    sort_list = np.argsort(to_sort)

    sorted_data = [data_array[i] for i in sort_list]

    # TODO: add better autoscaling options. Make scaling a function in nr_tools?
    Q, R, dR, dQ = [], [], [], []
    dict_output = []
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

        Q.append(result[0, :])
        R.append(result[1, :])
        dR.append(result[2, :])
        dQ.append(result[3, :])  

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
    
    # Sort by Q for combined data
    idx = np.argsort(Q_combined)
    combine_results = {'Q': Q_combined[idx], 'R': R_combined[idx], 'dR': dR_combined[idx], 'dQ': dQ_combined[idx]}

    return seq_list, run_list, combine_results

# TODO: Fix and update this part!!
def write_template(seq_list, run_list, template_file, output_dir, save_name=None):
    """
    Read the appropriate entry in a template file and save an updated
    copy with the updated run number.

    Parameters
    ----------
    seq_list : list
        The sequence identifiers
    run_list : list
        The run numbers
    template_file : str
        Path to the template file
    output_dir : str
        Directory where the output files are saved
    """
    with open(template_file, "r") as fd:
        xml_str = fd.read()
        data_sets = reduction_template_reader.from_xml(xml_str)

        new_data_sets = []
        for i in range(len(seq_list)):
            if len(data_sets) >= seq_list[i]:
                data_sets[seq_list[i] - 1].data_files = [run_list[i]]
                new_data_sets.append(data_sets[seq_list[i] - 1])
            else:
                print("Too few entries [%s] in template for sequence number %s" % (len(data_sets), seq_list[i]))

    if not save_name:
        save_name = "REF_L_%s_auto_template.xml" % run_list[0]

    # Save the template that was used
    xml_str = reduction_template_reader.to_xml(new_data_sets)
    with open(os.path.join(output_dir, save_name), "w") as fd:
        fd.write(xml_str)

# TODO: Fix to use the one inside template.py. This is needed at the moment from issue with mantid/typing.py
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

