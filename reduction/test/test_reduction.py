import os

import mantid
import mantid.simpleapi as mtd_api
import numpy as np
mtd_api.config.appendDataSearchDir('tests/data/liquidsreflectometer-data/nexus')
mtd_api.config["default.facility"] = "SNS"
mtd_api.config["default.instrument"] = "REF_L"

mantid.kernel.config.setLogLevel(0)

from lr_reduction import event_reduction, template, workflow


def test_full_reduction():
    """
        Test the fill reduction chain
    """
    template_path = 'data/template.xml'
    
    #print(os.listdir('tests/data/liquidsreflectometer-data/nexus/'))
    qz_all = []
    refl_all = []
    d_refl_all = []
    first_run = None
    for run_number in range(198409, 198417):
        ws_sc = mtd_api.Load("REF_L_%s.nxs.h5" % run_number)
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

    assert(np.fabs(resolution - 0.02785205863936946) < 1e-5)
    ref_data = np.loadtxt('data/reference_rq.txt').T
    assert(len(ref_data[1]) == len(refl_all))
    assert(np.fabs(np.sum(ref_data[1]-refl_all)) < 1e-10)


def test_reduce_workflow():
    template_path = 'data/template.xml'
    output_dir = '/tmp'
    reduced_path = os.path.join(output_dir, 'REFL_198409_combined_data_auto.txt')
    if os.path.isfile(reduced_path):
        os.remove(reduced_path)

    for i in range(198409, 198417):
        ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir,
                        average_overlap=False)

    reference_path = 'data/reference_rq.txt'
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(3):
        assert(np.fabs(np.sum(_data[i]-_refl[i])) < 1e-10)

    # The reference was computed with a constant dq/q but our approach recalculates
    # it for each run, so we expect a small discrepancy within 1%.
    assert(np.sum((_data[3]-_refl[3])/_refl[3])/len(_refl[3]) < 0.01)


def test_reduce_workflow_201282():
    """
        Test to reproduce autoreduction output
    """
    template_path = 'data/template_201282.xml'
    output_dir = '/tmp'
    reduced_path = os.path.join(output_dir, 'REFL_201282_combined_data_auto.txt')
    os.remove(reduced_path)

    for i in range(201282, 201289):
        ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir,
                        average_overlap=False)

    reference_path = 'data/reference_rq_201282.txt'
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(3):
        assert(np.fabs(np.sum(_data[i]-_refl[i])) < 1e-10)

    # The reference was computed with a constant dq/q but our approach recalculates
    # it for each run, so we expect a small discrepancy within 1%.
    assert(np.sum((_data[3]-_refl[3])/_refl[3])/len(_refl[3]) < 0.01)


def test_background_subtraction():
    """
        Test with background subtraction off for the data and on for the normalization
    """
    template_path = 'data/template_short_nobck.xml'
    output_dir = '/tmp'
    reduced_path = os.path.join(output_dir, 'REFL_198382_combined_data_auto.txt')
    os.remove(reduced_path)

    for i in range(198388, 198390):
        ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir,
                        average_overlap=False)

    reference_path = 'data/reference_short_nobck.txt'
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(3):
        assert(np.fabs(np.sum(_data[i]-_refl[i])) < 1e-10)

    # The reference was computed with a constant dq/q but our approach recalculates
    # it for each run, so we expect a small discrepancy within 1%.
    assert(np.sum((_data[3]-_refl[3])/_refl[3])/len(_refl[3]) < 0.01)
