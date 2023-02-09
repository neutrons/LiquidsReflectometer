"""
    Autoreduction process for the Liquids Reflectometer
"""
import sys
import os
import json
import numpy as np

import mantid
import mantid.simpleapi as mtd_api

from . import template
from . import reduction_template_reader
from . import output
from . import event_reduction


def reduce(ws, template_file, output_dir, pre_cut=1, post_cut=1, average_overlap=False,
           q_summing=False, bck_in_q=False):
    """
        Function called by reduce_REFL.py, which lives in /SNS/REF_L/shared/autoreduce
        and is called by the automated reduction workflow.

        If average_overlap is used, overlapping points will be averaged, otherwise they
        will be left in the final data file.

        :param pre_cut: number of points to cut at the start of the distribution
        :param post_cut: number of points to cut at the end of the distribution
        :param average_overlap: if True, the overlapping points will be averaged
        :param q_summing: if True, constant-Q binning will be used
        :param bck_in_q: if True, and constant-Q binning is used, the background will be estimated
                         along constant-Q lines rather than along TOF/pixel boundaries.
    """
    # Call the reduction using the template
    qz_mid, refl, d_refl, meta_data = template.process_from_template_ws(ws, template_file,
                                                                        q_summing=q_summing,
                                                                        tof_weighted=q_summing,
                                                                        bck_in_q=bck_in_q, info=True)

    # Save partial results
    coll = output.RunCollection()
    npts = len(qz_mid)
    coll.add(qz_mid[pre_cut:npts-post_cut], refl[pre_cut:npts-post_cut],
             d_refl[pre_cut:npts-post_cut], meta_data=meta_data)
    coll.save_ascii(os.path.join(output_dir, 'REFL_%s_%s_%s_partial.txt' % (meta_data['sequence_id'],
                                                                            meta_data['sequence_number'],
                                                                            meta_data['run_number'])),
                    meta_as_json=True)

    # Assemble partial results into a single R(q)
    seq_list, run_list = assemble_results(meta_data['sequence_id'], output_dir, average_overlap)

    # Save template
    write_template(seq_list, run_list, template_file, output_dir)

    # Return the sequence identifier
    return run_list[0]


def assemble_results(first_run, output_dir, average_overlap=False):
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
            coll.add_from_file(os.path.join(output_dir, item))

    coll.save_ascii(os.path.join(output_dir, 'REFL_%s_combined_data_auto.txt' % first_run))

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
    with open(os.path.join(output_dir, 'REFL_%s_auto_template.xml' % run_list[0]), 'w') as fd:
        fd.write(xml_str)


def reduce_fixed_two_theta(ws, template_file, output_dir, pre_cut=1, post_cut=1, average_overlap=False,
           q_summing=False, bck_in_q=False, peak_width=10, offset_from_first=True):
    """
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
    ths_value = ws.getRun()['ths'].value[0]

    # Read template so we can load the direct beam run
    template_data = template.read_template(template_file, sequence_number)

    # Load normalization run
    ws_db = mtd_api.LoadEventNexus("REF_L_%s" % template_data.norm_file)
    tthd_value = ws.getRun()['tthd'].value[0]

    # Look for parameters that might have been determined earlier for this measurement
    options_file = os.path.join(output_dir, 'REFL_%s_options.json' % sequence_id)
    if offset_from_first and sequence_number > 1 and os.path.isfile(options_file):
        with open(options_file, 'r') as fd:
            options = json.load(fd)
        pixel_offset = options['pixel_offset']
        tthd_calibration = options['tthd_db']
        twotheta = 2*ths_value + options['twotheta_offset']
    else:
        # Fit direct beam position
        x_min=template_data.norm_peak_range[0]-25
        x_max=template_data.norm_peak_range[1]+25
        tof, _x, _y = peak_finding.process_data(ws_db, summed=True, tof_step=200)
        peak_center = np.argmax(_y)
        db_center, db_width, _ = peak_finding.fit_signal_flat_bck(_x, _y, x_min=x_min, x_max=x_max, center=peak_center)
        print("    DB center: %g\t Width: %g from [%g %g]" % (db_center, db_width,
                                                              template_data.norm_peak_range[0],
                                                              template_data.norm_peak_range[1]))

        # Fit the reflected beam position
        x_min=template_data.data_peak_range[0]-peak_width
        x_max=template_data.data_peak_range[1]+peak_width
        tof, _x, _y = peak_finding.process_data(ws, summed=True, tof_step=200)
        peak_center = np.argmax(_y)
        sc_center, sc_width, _ = peak_finding.fit_signal_flat_bck(_x, _y, x_min=x_min, x_max=x_max, center=peak_center)
        pixel_offset = sc_center - (template_data.data_peak_range[1] + template_data.data_peak_range[0])/2.0
        print("    SC center: %g\t Width: %g" % (sc_center, sc_width))

        pixel_width = float(ws.getInstrument().getNumberParameter("pixel-width")[0]) / 1000.0
        sample_det_distance = event_reduction.EventReflectivity.DEFAULT_4B_SAMPLE_DET_DISTANCE
        twotheta = np.arctan((db_center-sc_center)*pixel_width / sample_det_distance) / 2.0 * 180 / np.pi

        # Store the tthd of the direct beam and account for the fact that it may be
        # different from our reflected beam for this calibration data.
        # This will allow us to be compatible with both fixed and moving detector arm.
        tthd_db = ws_db.getRun()['tthd'].value[0]
        twotheta = twotheta + tthd_value - tthd_db

        # If this is the first angle, keep the value for later
        options = dict(twotheta_offset = twotheta - 2*ths_value,
                       pixel_offset = pixel_offset,
                       tthd_db = tthd_db, tthd_sc = tthd_value)
        with open(options_file, 'w') as fp:
            json.dump(options, fp)
    print("    Two-theta = %g" % twotheta)

    # Modify the template with the fitted results
    print("    Template peak: [%g %g]" % (template_data.data_peak_range[0],
                                          template_data.data_peak_range[1]))

    template_data.data_peak_range = np.rint(np.asarray(template_data.data_peak_range) + pixel_offset).astype(int)
    template_data.background_roi = np.rint(np.asarray(template_data.background_roi) + pixel_offset).astype(int)
    print("    New peak:      [%g %g]" % (template_data.data_peak_range[0],
                                          template_data.data_peak_range[1]))
    print("    New bck:       [%g %g]" % (template_data.background_roi[0],
                                          template_data.background_roi[1]))

    # Call the reduction using the template
    qz_mid, refl, d_refl, meta_data = template.process_from_template_ws(ws, template_data,
                                                                        q_summing=q_summing,
                                                                        tof_weighted=q_summing,
                                                                        bck_in_q=bck_in_q, info=True,
                                                                        theta_value = twotheta/2.0,
                                                                        ws_db = ws_db)

    # Save partial results
    coll = output.RunCollection()
    npts = len(qz_mid)
    coll.add(qz_mid[pre_cut:npts-post_cut], refl[pre_cut:npts-post_cut],
             d_refl[pre_cut:npts-post_cut], meta_data=meta_data)
    coll.save_ascii(os.path.join(output_dir, 'REFL_%s_%s_%s_partial.txt' % (meta_data['sequence_id'],
                                                                            meta_data['sequence_number'],
                                                                            meta_data['run_number'])),
                    meta_as_json=True)

    # Assemble partial results into a single R(q)
    seq_list, run_list = assemble_results(meta_data['sequence_id'], output_dir, average_overlap)

    # Save template
    write_template(seq_list, run_list, template_file, output_dir)

    # Return the sequence identifier
    return run_list[0]

