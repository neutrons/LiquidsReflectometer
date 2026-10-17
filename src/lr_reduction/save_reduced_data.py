import json
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

import lr_reduction.nr_tools as tools

# Header format 2: "Run Title", "Angles" and "NR_runs" are lists indexed by sequence position
# (seq_num - 1), with null where no run was reduced. Files without the line predate it.
# Header format 3 (header-scale-factors-per-position, F17): "Scaling factors" and "Lambda Range" are
# per-position records too -- the factor applied to each position's data and the wavelength range its
# reduction used -- and Config's LambdaMinUse/LambdaMaxUse hold them as lists. Written only when the
# caller supplies the records; a caller without them (the template path) still writes format 2.
HEADER_FORMAT = 3
LOGS_ONLY_HEADER_FORMAT = 2

#: The per-position records a format-3 header carries, beside the logs (NR_Reduction.reduce() builds them).
RECORD_KEYS = ("scale", "lambda_min", "lambda_max")


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

def _header_runs(runs):
    """The run numbers as NR_runs writes them: plain ints, None for a gap (R6).

    A numpy integer, digits, or an integral float are written as the int they are, so read_prior_header's
    ast.literal_eval reads the line back and the file can vouch for its own entries (``np.int64(221472)`` cannot
    be read back). Anything else is refused, never written.
    """
    if not isinstance(runs, (list, tuple)):
        return runs
    plain = []
    for run in runs:
        if run is None:
            plain.append(None)
        elif isinstance(run, (bool, np.bool_)):
            raise ValueError(f"run number {run!r} is not an integer; NR_runs cannot hold it")
        elif isinstance(run, (int, np.integer)):
            plain.append(int(run))
        elif isinstance(run, str) and run.strip().isdigit():
            plain.append(int(run))
        elif isinstance(run, (float, np.floating)) and float(run).is_integer():
            plain.append(int(run))
        else:
            raise ValueError(f"run number {run!r} is not an integer; NR_runs cannot hold it")
    return plain


def _build_header(config_header, log_values, full=True, eight_column=False, sequence=None):
    """
    Wrapper to handle assembly logic for the output file header.

    With the per-position records (RECORD_KEYS) in ``log_values``, the header is format 3: the Scaling factors
    and Lambda Range lines, and Config's LambdaMinUse/LambdaMaxUse, are lists by sequence position (the config
    object itself keeps its scalars). Without them, today's lines and format 2.
    """

    if eight_column:
        col_label = "columns = Q, R, dR, dQ (sigma), L, dL, T, dT"
    else:
        col_label = "columns = Q, R, dR, dQ (sigma)"

    if sequence:
        sorted_config = {k: tools.maybe_index(v, sequence) for k, v in config_header.items()}
    else:
        sorted_config = config_header

    records = all(key in log_values for key in RECORD_KEYS)
    try:
        config_values = make_json_safe(sorted_config)
        json.dumps(config_values)
    except:
        config_values = make_json_safe(sorted_config.__dict__)
    if records:
        # Config's runtime-owned range is per position in the file; the config object keeps its scalars (R5)
        config_values = {**config_values,
                         "LambdaMinUse": tools.clean_log_value(log_values["lambda_min"]),
                         "LambdaMaxUse": tools.clean_log_value(log_values["lambda_max"])}
    config_json = json.dumps(config_values)

    angle_header = json.dumps({"THS": tools.clean_log_value(log_values['ths']), "THI": tools.clean_log_value(log_values['thi']), "ThCen": tools.clean_log_value(log_values['ThCen'])})
    title_header = json.dumps({"title": log_values['title']})
    if records:
        scale_factor_header = json.dumps({"scale_factor": tools.clean_log_value(log_values["scale"])})
        lambda_header = json.dumps({"lambda_min": tools.clean_log_value(log_values["lambda_min"]),
                                    "lambda_max": tools.clean_log_value(log_values["lambda_max"])})
        marker = (f"Header format: {HEADER_FORMAT} (Run Title, Angles, NR_runs, Scaling factors and Lambda Range "
                  f"are indexed by sequence position)")
    else:
        scale_factor_header = json.dumps({"scale_factor": tools.clean_log_value(sorted_config.ScaleFactor)})
        lambda_header = f"{sorted_config.LambdaMinUse}\u212B to {sorted_config.LambdaMaxUse}\u212B"
        marker = (f"Header format: {LOGS_ONLY_HEADER_FORMAT} (Run Title, Angles and NR_runs are indexed by "
                  f"sequence position)")
    nr_runs = _header_runs(sorted_config.RBnum)
    if full:
        head = (
            f"NR_runs = {nr_runs}\n"
            f"Run Title: {title_header}\n"
            f"DB = {sorted_config.DBname}\n"
            f"Method = {sorted_config.method_per_run}\n"
            f"Normalize = {sorted_config.Normalize}\n"
            f"Autoscale = {sorted_config.AutoScale}\n"
            f"Scaling factors = {scale_factor_header}\n"
            f"Lambda Range = {lambda_header}\n"
            f"Angles: {angle_header}\n"
            f"{marker}\n"
            f"{'---' * 20}\n"
            f"Config: {config_json}\n"
            f"{'---' * 20}\n"
            f"{col_label}\n"
            f"{'---' * 20}"
        )

    else:
        head = (
            f"NR_runs = {nr_runs}\n"
            f"Run Title: {title_header}\n"
            f"DB = {sorted_config.DBname}\n"
            f"Method = {sorted_config.method_per_run}\n"
            f"Normalize = {sorted_config.Normalize}\n"
            f"Autoscale = {sorted_config.AutoScale}\n"
            f"Scaling factors = {scale_factor_header}\n"
            f"Lambda Range = {lambda_header}\n"
            f"Angles: {angle_header}\n"
            f"{marker}\n"
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
        except Exception:  # noqa: BLE001 -- serialization fallback must absorb any str() failure
            return repr(obj)
    else:
        return obj

# TODO: link this up to store the figures and save them out.
def save_plot_pdf_summary(savepath, savename, fig_list):

    output_path = Path(savepath) / f"plot_summary_{savename}.pdf"
    with PdfPages(output_path) as pdf:
        for fig in fig_list:
            pdf.savefig(fig)
    print(f"Saved plot summary {output_path}.")
