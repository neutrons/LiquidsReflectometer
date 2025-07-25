"""
Reduce a data run using a template generated by RefRed
"""

import os
from functools import reduce

import mantid.simpleapi as api
import numpy as np
from mantid.api import *
from mantid.kernel import *

from lr_reduction import event_reduction, peak_finding, reduction_template_reader
from lr_reduction.instrument_settings import InstrumentSettings
from lr_reduction.reduction_template_reader import ReductionParameters

TOLERANCE = 0.07
OUTPUT_NORM_DATA = False


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


def scaling_factor(scaling_factor_file, workspace, match_slit_width=True):
    """
    Apply scaling factor from reference scaling data
    @param workspace: Mantid workspace
    """
    if not os.path.isfile(scaling_factor_file):
        print("Could not find scaling factor file: %s" % scaling_factor_file)
        return workspace

    # Get the wavelength
    lr = workspace.getRun().getProperty("LambdaRequest").value[0]
    lr_value = float(f"{lr:.2f}")

    s1h = abs(workspace.getRun().getProperty("S1VHeight").value[0])
    s1w = abs(workspace.getRun().getProperty("S1HWidth").value[0])
    s2h = abs(workspace.getRun().getProperty("SiVHeight").value[0])
    s2w = abs(workspace.getRun().getProperty("SiHWidth").value[0])

    def _reduce(accumulation, item):
        """
        Reduce function that accumulates values in a dictionary
        """
        toks_item = item.split("=")
        if len(toks_item) != 2:
            return accumulation
        if isinstance(accumulation, dict):
            accumulation[toks_item[0].strip()] = toks_item[1].strip()
        else:
            toks_accum = accumulation.split("=")
            accumulation = {toks_item[0].strip(): toks_item[1].strip(), toks_accum[0].strip(): toks_accum[1].strip()}
        return accumulation

    def _value_check(key, data, reference):
        """
        Check an entry against a reference value
        """
        if key in data:
            return abs(abs(float(data[key])) - abs(float(reference))) <= TOLERANCE
        return False

    with open(scaling_factor_file, "r") as fd:
        file_content = fd.read()

    data_found = None
    for line in file_content.split("\n"):
        if line.startswith("#"):
            continue

        # Parse the line of data and produce a dict
        toks = line.split()
        data_dict = reduce(_reduce, toks, {})

        # Get ordered list of keys
        keys = []
        for token in toks:
            key_value = token.split("=")
            if len(key_value) == 2:
                keys.append(key_value[0].strip())

        # Skip empty lines
        if len(keys) == 0:
            continue
        # Complain if the format is non-standard
        elif len(keys) < 10:
            print("Bad scaling factor entry\n  %s" % line)
            continue

        # Sanity check
        if keys[0] != "IncidentMedium" and keys[1] != "LambdaRequested" and keys[2] != "S1H":
            print("The scaling factor file isn't standard: bad keywords")
        # The S2H key has been changing in the earlier version of REFL reduction.
        # Get the key from the data to make sure we are backward compatible.
        s2h_key = keys[3]
        s2w_key = keys[5]
        if (
            "IncidentMedium" in data_dict
            and _value_check("LambdaRequested", data_dict, lr_value)
            and _value_check("S1H", data_dict, s1h)
            and _value_check(s2h_key, data_dict, s2h)
        ):
            if not match_slit_width or (_value_check("S1W", data_dict, s1w) and _value_check(s2w_key, data_dict, s2w)):
                data_found = data_dict
                break

    if data_found is not None:
        a = float(data_found["a"])
        b = float(data_found["b"])
        a_error = float(data_found["error_a"])
        b_error = float(data_found["error_b"])
    else:
        return 1, 0, 0, 0
    return a, b, a_error, b_error


def process_from_template(
    run_number, template_path, q_summing=False, normalize=True, tof_weighted=False, bck_in_q=False, clean=False, info=False
):
    """
    The clean option removes leading zeros and the drop when doing q-summing
    """
    # For backward compatibility, consider the case of a list of run numbers to be added
    if "," in str(run_number):
        list_of_runs = str(run_number).split(",")
        run_number = "+".join(list_of_runs)
    # Load data
    ws_sc = api.Load("REF_L_%s" % run_number, OutputWorkspace="REF_L_%s" % run_number)
    return process_from_template_ws(
        ws_sc, template_path, q_summing=q_summing, tof_weighted=tof_weighted, bck_in_q=bck_in_q, clean=clean, info=info, normalize=normalize
    )


def process_from_template_ws(
    ws_sc,
    template_data,
    q_summing=False,
    tof_weighted=False,
    bck_in_q=False,
    clean=False,
    info=False,
    normalize=True,
    theta_value=None,
    ws_db=None,
):
    # Get the sequence number
    sequence_number = 1
    if ws_sc.getRun().hasProperty("sequence_number"):
        sequence_number = ws_sc.getRun().getProperty("sequence_number").value[0]

    # Load the template
    if isinstance(template_data, str):
        template_data = read_template(template_data, sequence_number)

    if template_data.dead_time:
        ws_sc = event_reduction.apply_dead_time_correction(ws_sc, template_data)

    # Load normalization run
    normalize = normalize and template_data.apply_normalization
    if ws_db is None and normalize:
        ws_db = api.LoadEventNexus("REF_L_%s" % template_data.norm_file)
        attenuator_thickness = event_reduction.get_attenuation_info(ws_db)
        if attenuator_thickness > 0:
            ws_db = event_reduction.process_attenuation(ws_db, thickness=attenuator_thickness)

    # Apply dead time correction
    if normalize and template_data.dead_time:
        ws_db = event_reduction.apply_dead_time_correction(ws_db, template_data)

    # Apply instrument settings
    if template_data.apply_instrument_settings:
        instrument_settings = InstrumentSettings(
            template_data.apply_instrument_settings,
            template_data.source_detector_distance,
            template_data.sample_detector_distance,
            template_data.num_x_pixels,
            template_data.num_y_pixels,
            template_data.pixel_width,
            template_data.xi_reference,
            template_data.s1_sample_distance,
        )
    else:
        instrument_settings = None

    # If we run in theta-theta geometry, we'll need thi
    thi_value = ws_sc.getRun()["thi"].value[0]
    ths_value = ws_sc.getRun()["ths"].value[0]

    # NOTE: An offset is no longer used be default. To use itm we can use
    # the EventReflectivity directly.

    _wl = ws_sc.getRun()["LambdaRequest"].value[0]

    # Determine scattering angle. If a value was given, don't add the offset
    if theta_value is not None:
        theta = theta_value * np.pi / 180.0
    else:
        if (
            "BL4B:CS:ExpPl:OperatingMode" in ws_sc.getRun()
            and ws_sc.getRun().getProperty("BL4B:CS:ExpPl:OperatingMode").value[0] == "Free Liquid"
        ):
            theta = thi_value * np.pi / 180.0
        else:
            theta = ths_value * np.pi / 180.0

        # Add offset
        theta += template_data.angle_offset * np.pi / 180.0

    theta_degrees = theta * 180 / np.pi
    print("wl=%g; ths=%g; thi=%g; offset=%g; theta used=%g" % (_wl, ths_value, thi_value, template_data.angle_offset, theta_degrees))

    # Get the reduction parameters from the template
    peak = template_data.data_peak_range
    if template_data.subtract_background:
        peak_bck = template_data.background_roi
        if template_data.two_backgrounds is False:
            peak_bck = peak_bck[0:2]  # retain only the first background
    else:
        peak_bck = None

    peak_center = (peak[0] + peak[1]) / 2.0

    # Fit the reflected beam position, which may not be in the middle and is
    # used in the q-summing calculation
    if q_summing:
        x_min = template_data.data_peak_range[0]
        x_max = template_data.data_peak_range[1]
        _, _x, _y = peak_finding.process_data(ws_sc, summed=True, tof_step=200)
        peak_center = np.argmax(_y[x_min:x_max]) + x_min
        peak_center, sc_width, _ = peak_finding.fit_signal_flat_bck(_x, _y, x_min=x_min, x_max=x_max, center=peak_center, sigma=3.0)
        print("Peak center: %g" % peak_center)

    if template_data.data_x_range_flag:
        low_res = template_data.data_x_range
    else:
        low_res = None

    norm_peak = template_data.norm_peak_range
    if template_data.norm_x_range_flag:
        norm_low_res = template_data.norm_x_range
    else:
        norm_low_res = None

    # We are not subtracting background for the direct beam
    if template_data.subtract_norm_background:
        norm_bck = template_data.norm_background_roi
    else:
        norm_bck = None

    [tof_min, tof_max] = template_data.tof_range
    q_min = template_data.q_min
    q_step = -template_data.q_step

    # Perform the reduction
    event_refl = event_reduction.EventReflectivity(
        ws_sc,
        ws_db,
        signal_peak=peak,
        signal_bck=peak_bck,
        norm_peak=norm_peak,
        norm_bck=norm_bck,
        specular_pixel=peak_center,
        signal_low_res=low_res,
        norm_low_res=norm_low_res,
        q_min=q_min,
        q_step=q_step,
        q_max=None,
        tof_range=[tof_min, tof_max],
        theta=np.abs(theta),
        instrument=event_reduction.EventReflectivity.INSTRUMENT_4B,
        functional_background=template_data.two_backgrounds,
        dead_time=template_data.dead_time,
        paralyzable=template_data.paralyzable,
        dead_time_value=template_data.dead_time_value,
        dead_time_tof_step=template_data.dead_time_tof_step,
        use_dead_time_threshold=template_data.use_dead_time_threshold,
        dead_time_threshold=template_data.dead_time_threshold,
        instrument_settings=instrument_settings,
        use_emission_time=template_data.use_emission_time,
    )
    print(f"{'*'*88}\nevent_refl:\n{event_refl}\n{'*'*88}")

    # R(Q)
    qz, refl, d_refl = event_refl.specular(
        q_summing=q_summing, tof_weighted=tof_weighted, bck_in_q=bck_in_q, clean=clean, normalize=normalize
    )
    qz_mid = (qz[:-1] + qz[1:]) / 2.0

    # When using composite direct beam, we don't need a scaling
    # factor file if the multiplier is in the logs
    scaling_factor_PV = "BL4B:CS:Autoreduce:ScaleMultiplier"
    if normalize and scaling_factor_PV in ws_db.getRun():
        # The standard changed during early implementation...
        if int(template_data.norm_file) > 212005:
            a = ws_db.getRun()[scaling_factor_PV].value[0] ** 2
        else:
            a = ws_db.getRun()[scaling_factor_PV].value[0]
        print("Composite scaling factor: %s" % a)
        d_refl /= a
        refl /= a
        b = err_a = err_b = 0
    elif normalize and template_data.scaling_factor_flag:
        print("Normalization options: %s %s" % (normalize, template_data.scaling_factor_flag))
        # Get the scaling factors
        a, b, err_a, err_b = scaling_factor(template_data.scaling_factor_file, ws_sc)

        _tof = 4 * np.pi * np.sin(event_refl.theta) * event_refl.constant / qz
        _tof_mid = (_tof[1:] + _tof[:-1]) / 2.0

        a_q = _tof_mid * b + a
        d_a_q = np.sqrt(_tof_mid**2 * err_b**2 + err_a**2)

        d_refl = np.sqrt(d_refl**2 / a_q**2 + refl**2 * d_a_q**2 / a_q**4)
        refl /= a_q
    else:
        print("Skipping scaling factor")
        a = b = 1
        err_a = err_b = 0

    # Trim ends as needed
    npts = len(qz_mid)
    qz_mid = qz_mid[template_data.pre_cut : npts - template_data.post_cut]
    refl = refl[template_data.pre_cut : npts - template_data.post_cut]
    d_refl = d_refl[template_data.pre_cut : npts - template_data.post_cut]

    if normalize and OUTPUT_NORM_DATA:
        lr = ws_sc.getRun().getProperty("LambdaRequest").value[0]
        s1h = abs(ws_sc.getRun().getProperty("S1VHeight").value[0])
        s1w = abs(ws_sc.getRun().getProperty("S1HWidth").value[0])
        s2h = abs(ws_sc.getRun().getProperty("SiVHeight").value[0])
        s2w = abs(ws_sc.getRun().getProperty("SiHWidth").value[0])

        # Apply scaling factor
        _norm = (a_q * event_refl.norm)[template_data.pre_cut : npts - template_data.post_cut]
        _d_norm = (a_q * event_refl.d_norm)[template_data.pre_cut : npts - template_data.post_cut]

        wl = 4.0 * np.pi / qz_mid * np.sin(np.abs(theta))

        with open("%s_%s_%2.3g_counts.txt" % (template_data.norm_file, lr, np.abs(ths_value)), "w") as fd:
            fd.write("# Theta: %g\n" % np.abs(ths_value))
            fd.write("# Slit 1 (h x w): %g x %g\n" % (s1h, s1w))
            fd.write("# Slit 2 (h x w): %g x %g\n" % (s2h, s2w))
            fd.write("# Q, wavelength, counts/mC, error\n")
            for i in range(len(_norm)):
                fd.write("%g %g %g %g\n" % (qz_mid[i], wl[i], _norm[i], _d_norm[i]))

    # We can optionally return details about the reduction process
    if info:
        meta_data = event_refl.to_dict()
        meta_data["scaling_factors"] = dict(a=a, err_a=err_a, b=b, err_b=err_b)
        meta_data["q_summing"] = q_summing
        meta_data["tof_weighted"] = tof_weighted
        meta_data["bck_in_q"] = bck_in_q
        meta_data["theta_offset"] = template_data.angle_offset
        return qz_mid, refl, d_refl, meta_data

    return qz_mid, refl, d_refl
