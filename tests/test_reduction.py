# standard imports
import os
from pathlib import Path
from unittest.mock import patch

# third-party imports
import mantid
import mantid.simpleapi as mtd_api
import numpy as np
import pytest

# lr_reduction imports
from lr_reduction import event_reduction, template, workflow
from lr_reduction.scaling_factors.calculate import StitchingType
from lr_reduction.template import read_template
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


def test_q_summing(template_dir, nexus_dir):
    """
    Test Q summing process
    """
    template_path = os.path.join(template_dir, "template.xml")
    template.read_template(template_path, 7)
    with amend_config(data_dir=nexus_dir):
        ws_sc = mtd_api.Load("REF_L_%s" % 198415)
    qz_mid0, refl0, _, meta_data = template.process_from_template_ws(ws_sc, template_path, info=True)

    assert np.fabs(meta_data["dq_over_q"] - 0.02261) < 1e-3

    # Now try with Q summing, which should have similar results
    qz_mid, refl, _, meta_data = template.process_from_template_ws(
        ws_sc, template_path, tof_weighted=True, info=True, q_summing=True
    )

    assert np.fabs(meta_data["dq_over_q"] - 0.009354) < 1e-5

    # Note that TOF weighted may have a slightly different range, so here we skip
    # the extra point.
    assert len(qz_mid0) == len(qz_mid[1:])

    # plt.plot(qz_mid,refl, 'o', label='q summing and tof weighted')
    # plt.plot(qz_mid0, refl0, 'x', label='no q summing')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend()
    # plt.show()

    # Q-summing for lowest angle gives odd edge effects at low q. Will need addressing
    # as part of future regridding but not typically used for lowest angle and check
    # agreement excluding these lowest points for now. (Graph shows change.)
    assert np.fabs(np.mean(refl[3:] - refl0[2:])) < 1e-6

    # Cleanup
    output_dir = "data/"
    cleanup_partial_files(output_dir, range(198409, 198417))


@pytest.mark.parametrize(
    "template_file, q_summing, expected_q_summing, tof_weighted",
    [
        ("template.xml", None, False, False),
        ("template_with_const_q_true.xml", None, True, False),
        ("template.xml", True, True, True),
        ("template.xml", False, False, True),
    ],
)
def test_q_summing_as_option(template_dir, nexus_dir, template_file, q_summing, expected_q_summing, tof_weighted):
    """
    Test Q summing with and without supplying q_summing option
    """
    template_path = os.path.join(template_dir, template_file)
    template.read_template(template_path, 7)
    with amend_config(data_dir=nexus_dir):
        ws_sc = mtd_api.Load("REF_L_%s" % 198415)

    with patch.object(event_reduction.EventReflectivity, "specular") as mock_specular:
        mock_specular.return_value = (np.array([1, 2, 3]), np.array([0.1, 0.2]), np.array([0.01, 0.02]))

        qz_mid, refl, _, meta_data = template.process_from_template_ws(
            ws_sc, template_path, tof_weighted=tof_weighted, info=True, q_summing=q_summing
        )

        mock_specular.assert_called_once()
        call_kwargs = mock_specular.call_args[1]
        assert call_kwargs["q_summing"] == expected_q_summing

    output_dir = "data/"
    cleanup_partial_files(output_dir, range(198409, 198417))


def test_full_reduction(template_dir, nexus_dir):
    """
    Test the full reduction chain
    """
    template_path = os.path.join(template_dir, "template.xml")
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

    assert np.fabs(resolution - 0.022751) < 1e-5
    ref_data = np.loadtxt("data/reference_rq.txt").T

    # Optional plotting for checking tests:
    # import matplotlib.pyplot as plt
    # plt.plot(qz_all,refl_all, 'o', label='New reduction')
    # plt.plot(ref_data[0], ref_data[1], 'x', label='Prior data')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend()
    # plt.show()

    assert len(ref_data[1]) == len(refl_all)
    relative_difference = (ref_data[1] - refl_all) / ref_data[1]
    average_relative_difference = np.fabs(np.sum(relative_difference) / len(refl_all))
    assert average_relative_difference < 0.05

    # Cleanup
    output_dir = "data/"
    cleanup_partial_files(output_dir, range(198409, 198417))


def test_reduce_workflow(template_dir, nexus_dir, tmp_path):
    template_path = os.path.join(template_dir, "template.xml")
    output_dir = tmp_path
    reduced_path = os.path.join(output_dir, "REFL_198409_combined_data_auto.txt")

    for i in range(198409, 198417):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir, average_overlap=False)

    reference_path = os.path.join(template_dir, "reference_rq.txt")
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(3):
        fractional_differences = (_data[i] - _refl[i]) / _data[i]
        average_fractional_difference = np.fabs(np.sum(fractional_differences) / len(_refl[i]))
        assert average_fractional_difference < 0.07

def test_reduce_workflow_with_stitching_automatic_average(template_dir, nexus_dir, tmp_path):
    """
    Test the complete working, but this time we average the point in the
    overlap regions.
    """
    template_path = os.path.join(template_dir, "template_stitching_automatic_average.xml")
    output_dir = tmp_path

    for i in range(198409, 198417):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir, average_overlap=True)

    template_out_path = os.path.join(output_dir, "REF_L_198409_auto_template.xml")

    assert os.path.isfile(template_out_path)

    template_out = read_template(template_out_path, 198409)

    assert template_out.stitching_reflectivity_scale_factor == 1.0
    assert template_out.stitching_configuration.type == StitchingType.AUTOMATIC_AVERAGE
    assert template_out.stitching_configuration.scale_factor_qmin == 0.01
    assert template_out.stitching_configuration.scale_factor_qmax == 0.03
    assert template_out.stitching_configuration.normalize_first_angle is False
    # TODO: Add to this test once save/loading is complete in ewm13786


def test_reduce_functional_bck(nexus_dir, template_dir, tmp_path):
    os.chdir(Path(template_dir).parent)
    template_path = os.path.join(template_dir, "template_fbck.xml")
    output_dir = tmp_path
    reduced_path = os.path.join(output_dir, "REFL_198409_combined_data_auto.txt")

    for i in range(198409, 198417):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)

        sequence_number = ws.getRun().getProperty("sequence_number").value[0]
        template_data = template.read_template(template_path, sequence_number)
        template_data.two_backgrounds = True

        workflow.reduce(ws, template_data, output_dir=output_dir, average_overlap=False)

    reference_path = os.path.join(template_dir, "reference_fbck.txt")
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    # Optional plotting for checking tests:
    # plt.errorbar(_refl[0], _refl[1], _refl[2], marker='o', label='New reduction')
    # plt.errorbar(_data[0], _data[1],_data[2], marker='x', label='Prior data')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend()
    # plt.show()

    for i in range(3):
        fractional_differences = (_data[i] - _refl[i]) / _data[i]
        average_fractional_difference = np.fabs(np.sum(fractional_differences) / len(_refl[i]))
        assert average_fractional_difference < 0.07


def test_compute_wavelength_resolution_n_spectra():
    """
    Call compute wavelength resolution method with a workspace
    of multiple spectra
    """
    ws = mtd_api.CreateSampleWorkspace(WorkspaceType="Event")

    with pytest.raises(ValueError):
        _, _ = event_reduction.compute_wavelength_resolution(ws)


def test_reduce_bck_option_mismatch(template_dir, nexus_dir, tmp_path):
    """
    Ask for functional background but pass by a background range with
    only a single region. This will revert to simple averaging over the range.
    """
    template_path = os.path.join(template_dir, "template.xml")
    output_dir = tmp_path
    reduced_path = os.path.join(output_dir, "REFL_198409_combined_data_auto.txt")

    for i in range(198409, 198417):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)
        sequence_number = ws.getRun().getProperty("sequence_number").value[0]
        template_data = template.read_template(template_path, sequence_number)
        template_data.background_roi = template_data.background_roi[:2]
        template_data.two_backgrounds = True
        workflow.reduce(ws, template_data, output_dir=output_dir, average_overlap=False)

    reference_path = os.path.join(template_dir, "reference_rq.txt")
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    for i in range(3):
        fractional_differences = (_data[i] - _refl[i]) / _data[i]
        average_fractional_difference = np.fabs(np.sum(fractional_differences) / len(_refl[i]))
        assert average_fractional_difference < 0.07


def test_reduce_workflow_with_overlap_avg(template_dir, nexus_dir, tmp_path):
    """
    Test the complete working, but this time we average the point in the
    overlap regions.
    """
    template_path = os.path.join(template_dir, "template.xml")
    output_dir = tmp_path
    reduced_path = os.path.join(output_dir, "REFL_198409_combined_data_auto.txt")

    for i in range(198409, 198417):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir, average_overlap=True)

    reference_path = os.path.join(template_dir, "reference_rq_avg.txt")
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    # Optional plotting for checking tests:
    # import matplotlib.pyplot as plt
    # plt.errorbar(_refl[0], _refl[1], _refl[2], marker='o', label='New reduction')
    # plt.errorbar(_data[0], _data[1],_data[2], marker='x', label='Prior data')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend()
    # plt.show()

    for i in range(3):
        fractional_differences = (_data[i] - _refl[i]) / _data[i]
        average_fractional_difference = np.fabs(np.sum(fractional_differences) / len(_refl[i]))
        assert average_fractional_difference < 0.07


def test_quick_reduce(nexus_dir, datarepo_dir):
    """
    Test the quick reduction workflow
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_201284")
        ws_db = mtd_api.Load("REF_L_201045")
    _refl = workflow.reduce_explorer(ws, ws_db, center_pixel=145, db_center_pixel=145)
    reference_path = os.path.join(datarepo_dir, "reference_r201284_quick.txt")
    _data = np.loadtxt(reference_path).T
    # Optional plotting for checking tests:
    # plt.errorbar(_refl[0], _refl[1], _refl[2], marker='o', label='New reduction')
    # plt.errorbar(_data[0], _data[1],_data[2], marker='x', label='Prior data')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend()
    # plt.show()
    for i in range(3):
        assert np.fabs(np.sum(_data[i] - _refl[i])) < 1e-5


def test_reduce_workflow_201282(template_dir, nexus_dir):
    """
    Test to reproduce autoreduction output
    """
    template_path = os.path.join(template_dir, "template_201282.xml")
    output_dir = "data/"
    reduced_path = os.path.join(output_dir, "REFL_201282_combined_data_auto.txt")
    if os.path.isfile(reduced_path):
        os.remove(reduced_path)

    for i in range(201282, 201289):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir, average_overlap=False)

    reference_path = os.path.join(template_dir, "reference_rq_201282.txt")
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    # Optional plotting for checking tests:
    # plt.errorbar(_refl[0], _refl[1], _refl[2], marker='o', label='New reduction')
    # plt.errorbar(_data[0], _data[1],_data[2], marker='x', label='Prior data')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend()
    # plt.show()

    for i in range(3):
        fractional_differences = (_data[i] - _refl[i]) / _refl[i]
        average_fractional_difference = np.fabs(np.sum(fractional_differences) / len(_refl[i]))
        assert average_fractional_difference < 0.07


def test_background_subtraction(template_dir, nexus_dir):
    """
    Test with background subtraction off for the data and on for the normalization
    """
    template_path = os.path.join(template_dir, "template_short_nobck.xml")
    output_dir = "data/"
    reduced_path = os.path.join(output_dir, "REFL_198382_combined_data_auto.txt")
    if os.path.isfile(reduced_path):
        os.remove(reduced_path)

    for i in range(198388, 198390):
        with amend_config(data_dir=nexus_dir):
            ws = mtd_api.Load("REF_L_%s" % i)
        workflow.reduce(ws, template_path, output_dir=output_dir, average_overlap=False)

    reference_path = os.path.join(template_dir, "reference_short_nobck.txt")
    if os.path.isfile(reference_path):
        _data = np.loadtxt(reference_path).T

    if os.path.isfile(reduced_path):
        _refl = np.loadtxt(reduced_path).T

    # Optional plotting for checking tests:
    # plt.errorbar(_refl[0], _refl[1], _refl[2], marker='o', label='New reduction')
    # plt.errorbar(_data[0], _data[1],_data[2], marker='x', label='Prior data')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend()
    # plt.show()

    for i in range(3):
        fractional_differences = (_data[i] - _refl[i]) / _refl[i]
        average_fractional_difference = np.fabs(np.sum(fractional_differences) / len(_refl[i]))
        assert average_fractional_difference < 0.02

    cleanup_partial_files(output_dir, range(198382, 198390))
