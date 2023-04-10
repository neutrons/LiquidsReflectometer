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


def test_full_reduction():
    """
        Test the fill reduction chain
    """
    template_path = 'data/template.xml'

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

        for i in range(len(qz_mid)):
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
    assert(np.fabs(np.sum(ref_data[1]-refl_all)) < 1e-10)


def test_reduce_workflow():
    template_path = 'data/template.xml'
    output_dir = '/tmp'
    for i in range(198409, 198417):
        ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir,
                        average_overlap=False)

    reduced_path = 'data/reference_rq.txt'
    if os.path.isfile(reduced_path):
        _data = np.loadtxt(reduced_path).T

    reduced_path = os.path.join(output_dir, 'REFL_198409_combined_data_auto.txt')
    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(3):
        assert(np.fabs(np.sum(_data[i]-_refl[i])) < 1e-10)

    # The reference was computed with a constant dq/q but our approach recalculates
    # it for each run, so we expect a small discrepancy within 1%.
    assert(np.sum((_data[3]-_refl[3])/_refl[3])/len(_refl[3]) < 0.01)
