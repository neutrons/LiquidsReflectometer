import os
import pytest
import warnings

import numpy as np

import mantid
import mantid.simpleapi as mtd_api
mantid.kernel.config.setLogLevel(3)

from lr_reduction import template
from lr_reduction import event_reduction
from lr_reduction import workflow


def no_test_full_reduction():
    """
        Test the fill reduction chain
    """
    template_path = 'data/template.xml'

    # Number of data points to cut at the beginning and end of each run
    pre_cut = 1
    post_cut = 1

    qz_all = []
    refl_all = []
    d_refl_all = []
    first_run = None

    for run_number in range(198409, 198417):
        ws_sc = mtd_api.Load("REF_L_%s" % run_number)
        qz_mid, refl, d_refl = template.process_from_template_ws(ws_sc, template_path)

        if first_run is None:
            first_run = run_number
            resolution = event_reduction.compute_resolution(ws_sc)

        idx = np.fabs(refl) > 0
        qz_mid = qz_mid[idx][pre_cut:-post_cut]
        refl = refl[idx][pre_cut:-post_cut]
        d_refl = d_refl[idx][pre_cut:-post_cut]

        for i in range(len(qz_mid)-1, -1, -1):
            qz_all.append(qz_mid[i])
            refl_all.append(refl[i])
            d_refl_all.append(d_refl[i])

    qz_all = np.asarray(qz_all)
    refl_all = np.asarray(refl_all)
    d_refl_all = np.asarray(d_refl_all)
    idx = np.argsort(qz_all)

    qz_all = np.take_along_axis(qz_all, idx, axis=None)
    refl_all = np.take_along_axis(refl_all, idx, axis=None)
    d_refl_all = np.take_along_axis(d_refl_all, idx, axis=None)

    assert(resolution == 0.02785205863936946)
    ref_data = np.loadtxt('data/reference_rq.txt').T
    assert(len(ref_data[1]) == len(refl_all))
    assert(np.sum(ref_data[1]-refl_all) == 0)


def test_reduce_workflow():
    ws = mtd_api.Load("REF_L_198409")
    template_path = 'data/template.xml'
    workflow.reduce(ws, template_path, output_dir='/tmp')

    
    
def test_template_processing():
    """
        Test that we can read and write a template
    """
    template_path = 'data/template.xml'
    ws = mtd_api.Load("REF_L_198409")
    workflow.process_template(ws, template_path, '/tmp')
    template_data = template.read_template('/tmp/REF_L_198409_auto_template.xml', 1)
    assert(template_data.data_files[0] == 198409)
