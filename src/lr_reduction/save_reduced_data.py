from pathlib import Path
import json

import lr_reduction.nr_tools as tools
import numpy as np
from lr_reduction.nr_reduction_config import NRReductionConfig
from matplotlib.backends.backend_pdf import PdfPages

def save_results(results, config_header, log_values, sname = None, full=True, eight_column=False, sequence=None):
    """
    Save results as .dat file with header
    results: results to save
    config_header: config to be used for the meta data in header
    log_values: log values to include in meta data header (specifically theta values)
    sname: optional save name to overwrite the default
    full: flag to include more information in the header
    eight_column: option to save out 8-column data with L, dL, T, dT in addition to the standard 4 column
    sequence: option to specify index within the set of runs for selecting settings from config into the header #TODO: should this be in the method?

    Parameters
    ----------
    results : dict
        Results from reduce() method
    """

    if eight_column:
        array = np.column_stack((results['Q'], results['R'], results['dR'], results['dQ'],
                                        results['L'], results['dL'], results['T'], results['dT']))
    else:
        array = np.column_stack((results['Q'], results['R'], results['dR'], results['dQ']))
        
    head = _build_header(full=full, eight_column=eight_column, config_header=config_header, sequence=sequence, log_values=log_values)
    
    if not sname:
        output_file = config_header.Spath / f"{config_header.Sname}"
    else:
        output_file = config_header.Spath / f"{sname}"

    if eight_column:
        output_file = f"{output_file}_8col.dat"
    else:
        output_file = f"{output_file}.dat"
    np.savetxt(output_file,
                array, header=head, delimiter='\t')
    print(f"Saved result to {output_file}")

def _build_header(config_header, log_values, full=True, eight_column=False, sequence=None):
    """
    Wrapper to handle assembly logic for the output file header.
    """

    if eight_column:
        col_label = "columns = Q, R, dR, dQ (sigma), L, dL, T, dT"
    else:
        col_label = "columns = Q, R, dR, dQ (sigma)"

    if sequence:
        sorted_config = {k: tools.maybe_index(v, sequence) for k, v in config_header.items()}
    else:
        sorted_config = config_header

    try:
        config_json = json.dumps(make_json_safe(sorted_config))
    except:
        config_json = sorted_config.__dict__
        config_json = json.dumps(make_json_safe(config_json))

    # TODO: Lambda Use values need to be arrays and stored angles need to be arrays.
    if full:
        head = (
            f"NR_runs = {sorted_config.RBnum}\n"
            f"Run Title = {log_values['title']}\n"
            f"DB = {sorted_config.DBname}\n"
            f"Method = {sorted_config.method_per_run}\n"
            f"Normalize = {sorted_config.Normalize}\n"
            f"Autoscale = {sorted_config.AutoScale}\n"
            f"Scaling factors = {sorted_config.ScaleFactor}\n"
            f"Lambda Range = {sorted_config.LambdaMinUse}\u212B to {sorted_config.LambdaMaxUse}\u212B\n"
            f"THS = {log_values['ths']}, THI = {log_values['thi']}, ThCen = {log_values['ThCen']}\n"
            f"{'---' * 20}\n"
            f"Config: {config_json}\n"
            f"{'---' * 20}\n"
            f"{col_label}\n"
            f"{'---' * 20}"
        )

    else:
        head = (
            f"NR_runs = {sorted_config.RBnum}\n"
            f"Run Title = {log_values['title']}\n"
            f"DB = {sorted_config.DBname}\n"
            f"Method = {sorted_config.method_per_run}\n"
            f"Normalize = {sorted_config.Normalize}\n"
            f"Autoscale = {sorted_config.AutoScale}\n"
            f"Scaling factors = {sorted_config.ScaleFactor}\n"
            f"Lambda Range = {sorted_config.LambdaMinUse}\u212B to {sorted_config.LambdaMaxUse}\u212B\n"
            f"THS = {log_values['ths']}, THI = {log_values['thi']}, ThCen = {log_values['ThCen']}\n"
            f"{'---' * 20}\n"
            f"Config: {config_json}\n"
            f"{'---' * 20}\n"
            f"{col_label}\n"
            f"{'---' * 20}"
        )

    return head


# Function to convert dictionary to json safe inputs
def make_json_safe(obj):
    if isinstance(obj, dict):
        return {k: make_json_safe(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [make_json_safe(v) for v in obj]
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    elif isinstance(obj, Path):
        return str(obj)
    elif callable(obj):
        # Convert functions/methods to a readable string to avoid JSON serialization errors
        try:
            return str(obj)
        except Exception:
            return repr(obj)
    else:
        return obj
    
# TODO: link this up to store the figures and save them out.
def save_plot_pdf_summary(savepath, savename, fig_list):

    output_path = Path(savepath.text() or ".") / f"plot_summary_{savename}.pdf"
    with PdfPages(output_path) as pdf:
        for fig, meta in fig_list:
            pdf.savefig(fig)
    print(f"Saved plot summary {output_path}.")





