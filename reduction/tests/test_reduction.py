# standard imports
import os
from pathlib import Path

# third-party imports
import mantid
import mantid.simpleapi as mtd_api
import numpy as np
import pytest

# lr_reduction imports
from lr_reduction import event_reduction, template, workflow
from lr_reduction.utils import amend_config

mtd_api.config["default.facility"] = "SNS"
mtd_api.config["default.instrument"] = "REF_L"
mantid.kernel.config.setLogLevel(3)


def cleanup_partial_files(output_dir, runs):
    """
    Clean up reduced files left behind after reduction
    """
    for i, r in enumerate(runs):
        reduced_path = os.path.join(output_dir, "REFL_%s_%d_%s_partial.txt" % (runs[0], i + 1, r))
        if os.path.isfile(reduced_path):
            os.remove(reduced_path)


def test_info(nexus_dir):
    """
    Test utility functions to get basic info
    """
    with amend_config(data_dir=nexus_dir):
        ws_sc = mtd_api.Load("REF_L_198409")
    wl_min, wl_max = event_reduction.get_wl_range(ws_sc)
    assert wl_min == 13.7
    assert wl_max == 16.3


def test_attenuation(nexus_dir):
    """
    Test attenuation calculation can complete
    """
    with amend_config(data_dir=nexus_dir):
        ws_sc = mtd_api.Load("REF_L_198409")
    event_reduction.process_attenuation(ws_sc, 0.005)


def test_q_summing(nexus_dir):
    """
    Test Q summing process
    """
    template_path = "data/template.xml"
    template.read_template(template_path, 7)
    with amend_config(data_dir=nexus_dir):
        ws_sc = mtd_api.Load("REF_L_%s" % 198415)
    qz_mid0, refl0, _, meta_data = template.process_from_template_ws(ws_sc, template_path, info=True)

    assert np.fabs(meta_data["dq_over_q"] - 0.02759) < 1e-3

    # Now try with Q summing, which should have similar results
    qz_mid, refl, _, meta_data = template.process_from_template_ws(ws_sc, template_path, tof_weighted=True, info=True, q_summing=True)

    assert np.fabs(meta_data["dq_over_q"] - 0.009354) < 1e-5

    # Note that TOF weighted may have a slightly different range, so here we skip
    # the extra point.
    assert len(qz_mid0) == len(qz_mid[1:])
    assert np.fabs(np.mean(refl[1:] - refl0)) < 1e-6

    # Cleanup
    output_dir = "data/"
    cleanup_partial_files(output_dir, range(198409, 198417))


def test_full_reduction(nexus_dir):
    """
    Test the full reduction chain
    """
    template_path = "data/template.xml"
    qz_all = []
    refl_all = []
    d_refl_all = []
    first_run = None
    for run_number in range(198409, 198417):
        with amend_config(data_dir=nexus_dir):
            ws_sc = mtd_api.Load("REF_L_%s" % run_number)
        qz_mid, refl, d_refl = template.process_from_template_ws(ws_sc, template_path)

        if first_run is None:
            first_run = run_number
            resolution = event_reduction.compute_resolution(ws_sc)
            wavelength, d_lambda = event_reduction.compute_wavelength_resolution(ws_sc)

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

    assert np.fabs(resolution - 0.02785205863936946) < 1e-5
    assert len(wavelength) == len(ws_sc.readX(1))
    assert len(d_lambda) == len(ws_sc.readY(1))
    assert np.all(np.array(wavelength) > 0)
    ref_data = np.loadtxt("data/reference_rq.txt").T
    assert len(ref_data[1]) == len(refl_all)
    assert np.fabs(np.sum(ref_data[1] - refl_all)) < 1e-10

    # Cleanup
    output_dir = "data/"
    cleanup_partial_files(output_dir, range(198409, 198417))


def test_reduce_workflow(nexus_dir):
    template_path = "data/template.xml"
    output_dir = "data/"
    reduced_path = os.path.join(output_dir, "REFL_198409_combined_data_auto.txt")
    if os.path.isfile(reduced_path):
        os.remove(reduced_path)

    for i in range(198409, 198417):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir, average_overlap=False)

    reference_path = "data/reference_rq.txt"
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(3):
        assert np.fabs(np.sum(_data[i] - _refl[i])) < 1e-10

    # The reference was computed with a constant dq/q but our approach recalculates
    # it for each run, so we expect a small discrepancy within 1%.
    assert np.sum((_data[3] - _refl[3]) / _refl[3]) / len(_refl[3]) < 0.01

    # Cleanup
    cleanup_partial_files(output_dir, range(198409, 198417))


def test_reduce_functional_bck(nexus_dir, template_dir):
    os.chdir(Path(template_dir).parent)
    template_path = "data/template_fbck.xml"
    output_dir = "data/"
    reduced_path = os.path.join(output_dir, "REFL_198409_combined_data_auto.txt")
    if os.path.isfile(reduced_path):
        os.remove(reduced_path)

    for i in range(198409, 198417):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)

        sequence_number = ws.getRun().getProperty("sequence_number").value[0]
        template_data = template.read_template(template_path, sequence_number)
        template_data.two_backgrounds = True

        workflow.reduce(ws, template_data, output_dir=output_dir, average_overlap=False)

    reference_path = "data/reference_fbck.txt"
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(2):
        assert np.fabs(np.sum(_data[i] - _refl[i])) < 1e-9

    # Error bars from fit might be different
    assert np.fabs(np.sum(_data[2] - _refl[2])) < 1e-8

    # The reference was computed with a constant dq/q but our approach recalculates
    # it for each run, so we expect a small discrepancy within 1%.
    assert np.sum((_data[3] - _refl[3]) / _refl[3]) / len(_refl[3]) < 0.01

    # Cleanup
    cleanup_partial_files(output_dir, range(198409, 198417))


def test_compute_wavelength_resolution_n_spectra():
    """
    Call compute wavelength resolution method with a workspace
    of multiple spectra
    """
    ws = mtd_api.CreateSampleWorkspace(WorkspaceType="Event")

    with pytest.raises(ValueError):
        _, _ = event_reduction.compute_wavelength_resolution(ws)


def test_reduce_bck_option_mismatch(nexus_dir):
    """
    Ask for functional background but pass by a background range with
    only a single region. This will revert to simple averaging over the range.
    """
    template_path = "data/template.xml"
    output_dir = "data/"
    reduced_path = os.path.join(output_dir, "REFL_198409_combined_data_auto.txt")
    if os.path.isfile(reduced_path):
        os.remove(reduced_path)

    for i in range(198409, 198417):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)
        sequence_number = ws.getRun().getProperty("sequence_number").value[0]
        template_data = template.read_template(template_path, sequence_number)
        template_data.background_roi = template_data.background_roi[:2]
        template_data.two_backgrounds = True
        workflow.reduce(ws, template_data, output_dir=output_dir, average_overlap=False)

    reference_path = "data/reference_rq.txt"
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(3):
        assert np.fabs(np.sum(_data[i] - _refl[i])) < 1e-10

    # The reference was computed with a constant dq/q but our approach recalculates
    # it for each run, so we expect a small discrepancy within 1%.
    assert np.sum((_data[3] - _refl[3]) / _refl[3]) / len(_refl[3]) < 0.01

    # Cleanup
    cleanup_partial_files(output_dir, range(198409, 198417))


def test_reduce_workflow_with_overlap_avg(nexus_dir):
    """
    Test the complete working, but this time we average the point in the
    overlap regions.
    """
    template_path = "data/template.xml"
    output_dir = "data/"
    reduced_path = os.path.join(output_dir, "REFL_198409_combined_data_auto.txt")
    if os.path.isfile(reduced_path):
        os.remove(reduced_path)

    for i in range(198409, 198417):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir, average_overlap=True)

    reference_path = "data/reference_rq_avg.txt"
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(3):
        assert np.fabs(np.sum(_data[i] - _refl[i])) < 1e-10

    # The reference was computed with a constant dq/q but our approach recalculates
    # it for each run, so we expect a small discrepancy within 1%.
    assert np.sum((_data[3] - _refl[3]) / _refl[3]) / len(_refl[3]) < 0.01

    # Cleanup
    cleanup_partial_files(output_dir, range(198409, 198417))


def test_quick_reduce(nexus_dir):
    """
    Test the quick reduction workflow
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_201284")
        ws_db = mtd_api.Load("REF_L_201045")

    _refl = workflow.reduce_explorer(ws, ws_db, center_pixel=145, db_center_pixel=145)
    reference_path = "data/reference_r201284_quick.txt"
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    for i in range(3):
        assert np.fabs(np.sum(_data[i] - _refl[i])) < 1e-10


def test_reduce_workflow_201282(nexus_dir):
    """
    Test to reproduce autoreduction output
    """
    template_path = "data/template_201282.xml"
    output_dir = "data/"
    reduced_path = os.path.join(output_dir, "REFL_201282_combined_data_auto.txt")
    if os.path.isfile(reduced_path):
        os.remove(reduced_path)

    for i in range(201282, 201289):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir, average_overlap=False)

    reference_path = "data/reference_rq_201282.txt"
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(3):
        assert np.fabs(np.sum(_data[i] - _refl[i])) < 1e-10

    # The reference was computed with a constant dq/q but our approach recalculates
    # it for each run, so we expect a small discrepancy within 1%.
    assert np.sum((_data[3] - _refl[3]) / _refl[3]) / len(_refl[3]) < 0.01


def test_background_subtraction(nexus_dir):
    """
    Test with background subtraction off for the data and on for the normalization
    """
    template_path = "data/template_short_nobck.xml"
    output_dir = "data/"
    reduced_path = os.path.join(output_dir, "REFL_198382_combined_data_auto.txt")
    if os.path.isfile(reduced_path):
        os.remove(reduced_path)

    for i in range(198388, 198390):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir, average_overlap=False)

    reference_path = "data/reference_short_nobck.txt"
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(3):
        assert np.fabs(np.sum(_data[i] - _refl[i])) < 1e-10

    # The reference was computed with a constant dq/q but our approach recalculates
    # it for each run, so we expect a small discrepancy within 1%.
    assert np.sum((_data[3] - _refl[3]) / _refl[3]) / len(_refl[3]) < 0.01

    cleanup_partial_files(output_dir, range(198382, 198390))
