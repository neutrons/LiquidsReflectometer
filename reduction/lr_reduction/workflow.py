"""
    Autoreduction process for the Liquids Reflectometer
"""
import sys
import os
import numpy as np

import mantid
import mantid.simpleapi as mtd_api

from . import template
from . import reduction_template_reader
from . import output


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
    idx = np.fabs(refl) > 0
    npts = len(qz_mid[idx])
    coll.add(qz_mid[idx][pre_cut:npts-post_cut], refl[idx][pre_cut:npts-post_cut],
             d_refl[idx][pre_cut:npts-post_cut], meta_data=meta_data)
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
            if not len(toks) == 5:
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
