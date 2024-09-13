"""
    Autoreduction process for the Liquids Reflectometer
"""
import json
import os

import mantid.simpleapi as mtd_api
import numpy as np

from . import event_reduction, output, reduction_template_reader, template


def reduce(ws, template_file, output_dir, average_overlap=False,
           theta_offset: float | None = 0,
           q_summing=False, bck_in_q=False, is_live=False):
    """
        Function called by reduce_REFL.py, which lives in /SNS/REF_L/shared/autoreduce
        and is called by the automated reduction workflow.

        If average_overlap is used, overlapping points will be averaged, otherwise they
        will be left in the final data file.

        :param average_overlap: if True, the overlapping points will be averaged
        :param q_summing: if True, constant-Q binning will be used
        :param bck_in_q: if True, and constant-Q binning is used, the background will be estimated
                         along constant-Q lines rather than along TOF/pixel boundaries.
        :param theta_offset: Theta offset to apply. If None, the template value will be used.
    """
    # Get the sequence number
    sequence_number = 1
    if ws.getRun().hasProperty("sequence_number"):
        sequence_number = ws.getRun().getProperty("sequence_number").value[0]
    # Read template if it was passed as a file path
    # It can be passsed as a dict directly
    if isinstance(template_file, str):
        template_data = template.read_template(template_file, sequence_number)
    else:
        template_data = template_file

    # Set the theta offset if it was provided
    if theta_offset is not None:
        template_data.angle_offset = theta_offset

    # Call the reduction using the template
    qz_mid, refl, d_refl, meta_data = template.process_from_template_ws(ws, template_data,
                                                                        q_summing=q_summing,
                                                                        tof_weighted=q_summing,
                                                                        clean=q_summing,
                                                                        bck_in_q=bck_in_q,
                                                                        info=True)

    # Save partial results
    coll = output.RunCollection()
    coll.add(qz_mid, refl, d_refl, meta_data=meta_data)

    # If this is live data, put it in a separate file to avoid conflict with auto-reduction
    if is_live:
        reduced_file = os.path.join(output_dir, 'REFL_live_partial.txt')
    else:
        reduced_file = os.path.join(output_dir, 'REFL_%s_%s_%s_partial.txt' % (meta_data['sequence_id'],
                                                                               meta_data['sequence_number'],
                                                                               meta_data['run_number']))
    coll.save_ascii(reduced_file, meta_as_json=True)

    # Assemble partial results into a single R(q)
    seq_list, run_list = assemble_results(meta_data['sequence_id'], output_dir,
                                          average_overlap, is_live=is_live)

    # Save template. This will not happen if the template_file input was
    # template data, which the template processing allows.
    if isinstance(template_file, str):
        write_template(seq_list, run_list, template_file, output_dir)
    else:
        print("Template data was passed instead of a file path: template data not saved")

    # Return the sequence identifier
    return run_list[0]


def assemble_results(first_run, output_dir, average_overlap=False, is_live=False):
    """
        Find related runs and assemble them in one R(q) data set
    """
    # Keep track of sequence IDs and run numbers so we can make a new template
    seq_list = []
    run_list = []
    coll = output.RunCollection(average_overlap=average_overlap)

    file_list = sorted(os.listdir(output_dir))
    for item in file_list:
        if item.startswith("REFL_%s" % first_run) and item.endswith('partial.txt'):
            toks = item.split('_')
            if not len(toks) == 5 or int(toks[2]) == 0:
                continue
            seq_list.append(int(toks[2]))
            run_list.append(int(toks[3]))

            # Read the partial data and add to a collection
            _, _, _, _, _meta = output.read_file(os.path.join(output_dir, item))
            if is_live or not _meta['start_time'] == "live":
                coll.add_from_file(os.path.join(output_dir, item))
        #elif item == "REFL_live_partial.txt":
        #    coll.add_from_file(os.path.join(output_dir, item))

    output_file_name = 'REFL_%s_combined_data_auto.txt' % first_run
    if is_live:
        output_file_name = 'REFL_%s_live_estimate.txt' % first_run
    coll.save_ascii(os.path.join(output_dir, output_file_name))

    return seq_list, run_list


def write_template(seq_list, run_list, template_file, output_dir):
    """
        Read the appropriate entry in a template file and save an updated
        copy with the updated run number.
    """
    with open(template_file, "r") as fd:
        xml_str = fd.read()
        data_sets = reduction_template_reader.from_xml(xml_str)

        new_data_sets = []
        for i in range(len(seq_list)):
            if len(data_sets) >= seq_list[i]:
                data_sets[seq_list[i]-1].data_files = [run_list[i]]
                new_data_sets.append(data_sets[seq_list[i]-1])
            else:
                print("Too few entries [%s] in template for sequence number %s" % (len(data_sets), seq_list[i]))

    # Save the template that was used
    xml_str = reduction_template_reader.to_xml(new_data_sets)
    with open(os.path.join(output_dir, 'REF_L_%s_auto_template.xml' % run_list[0]), 'w') as fd:
        fd.write(xml_str)


def offset_from_first_run(
        ws,
        template_file :str,
        output_dir: str,
        ):
    """
    Find a theta offset from the first peak.
    Used when sample is misaligned.
    """
    from . import peak_finding

    print("\nProcessing: %s" % ws.getRunNumber())
    # Get the sequence number
    sequence_number = 1
    sequence_id = 1
    if ws.getRun().hasProperty("sequence_number"):
        sequence_number = ws.getRun().getProperty("sequence_number").value[0]
        sequence_id = ws.getRun().getProperty("sequence_id").value[0]

    # Theta value that we are aiming for
    ths_value = ws.getRun()['ths'].value[-1]

    # Read template so we can load the direct beam run
    template_data = template.read_template(template_file, sequence_number)

    # Load normalization run
    print("    DB: %s" % template_data.norm_file)
    ws_db = mtd_api.LoadEventNexus("REF_L_%s" % template_data.norm_file)

    # Look for parameters that might have been determined earlier for this measurement
    options_file = os.path.join(output_dir, 'REFL_%s_options.json' % sequence_id)

    if sequence_number > 1 and os.path.isfile(options_file):
        with open(options_file, 'r') as fd:
            options = json.load(fd)
            return options['theta_offset']
    else:
        # Fit direct beam position
        x_min=template_data.norm_peak_range[0]
        x_max=template_data.norm_peak_range[1]
        _, _x, _y = peak_finding.process_data(ws_db, summed=True, tof_step=200)
        peak_center = np.argmax(_y)
        db_center, db_width, _ = peak_finding.fit_signal_flat_bck(_x, _y, x_min=x_min, x_max=x_max,
                                                           center=peak_center,
                                                           sigma=1.)
        print("    DB center: %g\t Width: %g from [%g %g]" % (db_center, db_width,
                                                              template_data.norm_peak_range[0],
                                                              template_data.norm_peak_range[1]))

        # Fit the reflected beam position
        x_min=template_data.data_peak_range[0]
        x_max=template_data.data_peak_range[1]
        _, _x, _y = peak_finding.process_data(ws, summed=True, tof_step=200)
        peak_center = np.argmax(_y[x_min:x_max]) + x_min
        sc_center, sc_width, _ = peak_finding.fit_signal_flat_bck(_x, _y, x_min=x_min, x_max=x_max, center=peak_center, sigma=3.)
        pixel_offset = sc_center - peak_center
        print("    SC center: %g\t Width: %g" % (sc_center, sc_width))

        settings = event_reduction.read_settings(ws)
        sample_det_distance = settings["sample-det-distance"]
        pixel_width = settings["pixel-width"] / 1000.0

        theta = np.arctan((sc_center-db_center) * pixel_width / sample_det_distance) / 2.0 * 180 / np.pi
        theta_offset = theta - ths_value

        # If this is the first angle, keep the value for later
        options = dict(theta_offset = theta_offset,
                       pixel_offset = pixel_offset)
        with open(options_file, 'w') as fp:
            json.dump(options, fp)

    return theta_offset


def reduce_explorer(ws, ws_db, theta_pv=None, center_pixel=145, db_center_pixel=145, peak_width=10):
    """
    """
    from . import peak_finding

    if theta_pv is None:
        if 'BL4B:CS:ExpPl:OperatingMode' in ws.getRun() \
            and ws.getRun().getProperty('BL4B:CS:ExpPl:OperatingMode').value[0] == 'Free Liquid':
            theta_pv = 'thi'
        else:
            theta_pv = 'ths'
    print("\nProcessing: %s" % ws.getRunNumber())

    # Theta value that we are aiming for
    theta_value = np.fabs(ws.getRun()[theta_pv].value[0])

    # Load normalization run
    tthd_value = ws.getRun()['tthd'].value[0]

    # Fit direct beam position
    x_min = center_pixel - 25
    x_max = center_pixel + 25
    tof, _x, _y = peak_finding.process_data(ws_db, summed=True, tof_step=200)
    peak_center = np.argmax(_y)
    db_center, db_width, _ = peak_finding.fit_signal_flat_bck(_x, _y, x_min=x_min, x_max=x_max, center=peak_center)
    print("    DB center: %g\t Width: %g" % (db_center, db_width))

    # Fit the reflected beam position
    x_min=db_center_pixel-peak_width
    x_max=db_center_pixel+peak_width
    tof, _x, _y = peak_finding.process_data(ws, summed=True, tof_step=200)
    peak_center = np.argmax(_y)
    sc_center, sc_width, _ = peak_finding.fit_signal_flat_bck(_x, _y, x_min=x_min, x_max=x_max, center=peak_center)
    print("    SC center: %g\t Width: %g" % (sc_center, sc_width))

    pixel_width = float(ws.getInstrument().getNumberParameter("pixel-width")[0]) / 1000.0
    sample_det_distance = event_reduction.EventReflectivity.DEFAULT_4B_SAMPLE_DET_DISTANCE
    twotheta = np.arctan((db_center-sc_center)*pixel_width / sample_det_distance) / 2.0 * 180 / np.pi

    # Store the tthd of the direct beam and account for the fact that it may be
    # different from our reflected beam for this calibration data.
    # This will allow us to be compatible with both fixed and moving detector arm.
    tthd_db = ws_db.getRun()['tthd'].value[0]
    twotheta = twotheta + tthd_value - tthd_db

    print("    Theta = %g   Two-theta = %g" % (theta_value, twotheta))

    # Perform the reduction
    width_mult = 2.5
    peak = [np.rint(sc_center - width_mult*sc_width).astype(int), np.rint(sc_center + width_mult*sc_width).astype(int)]
    norm_peak = [np.rint(db_center - width_mult*db_width).astype(int), np.rint(db_center + width_mult*db_width).astype(int)]
    peak_bck = [peak[0]-3, peak[1]+3]
    norm_bck = [norm_peak[0]-3, norm_peak[1]+3]

    tof_min = ws.getTofMin()
    tof_max = ws.getTofMax()

    theta = theta_value * np.pi / 180.

    #TODO: dead time correction should be applied here
    event_refl = event_reduction.EventReflectivity(ws, ws_db,
                                                   signal_peak=peak, signal_bck=peak_bck,
                                                   norm_peak=norm_peak, norm_bck=norm_bck,
                                                   specular_pixel=sc_center.value,
                                                   signal_low_res=[65,180], norm_low_res=[65,180],
                                                   q_min=None, q_max=None,
                                                   tof_range = [tof_min, tof_max],
                                                   theta=theta)

    # R(Q)
    qz, refl, d_refl = event_refl.specular(q_summing=False, tof_weighted=False,
                                           bck_in_q=False, clean=False, normalize=True)
    qz_mid = (qz[:-1] + qz[1:])/2.0

    return qz_mid, refl, d_refl
