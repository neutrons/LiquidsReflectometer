import unittest
import os
import pytest

import mantid
import mantid.simpleapi as mtd_api
import numpy as np
from mantid import config

mantid.kernel.config.setLogLevel(3)

from reduction.lr_reduction import event_reduction, template, workflow


class ReductionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.getcwd().endswith("LiquidsReflectometer"):
            os.chdir("tests")

        cwd = os.getcwd()
        dirs = [26010, 26776, 28662, 29196, 31279]
        for dir_num in dirs:
            config.appendDataSearchDir(str(os.path.join(cwd, f"data/liquidsreflectometer-data/SNS/REF_L/IPTS-{dir_num}/nexus")))

    @pytest.mark.datarepo()
    def test_full_reduction(self):
        """
        Test the fill reduction chain
        """
        template_path = "data/template.xml"
        qz_all = []
        refl_all = []
        d_refl_all = []
        first_run = None
        for run_number in range(198409, 198417):
            ws_sc = mtd_api.Load(f"REF_L_{run_number}")
            qz_mid, refl, d_refl = template.process_from_template_ws(ws_sc, template_path)

            if first_run is None:
                first_run = run_number
                resolution = event_reduction.compute_resolution(ws_sc)
                print(resolution)
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

        # assert(resolution == 0.02785205863936946)
        ref_data = np.loadtxt("data/reference_rq.txt").T
        assert len(ref_data[1]) == len(refl_all)
        assert np.fabs(np.sum(ref_data[1] - refl_all)) < 1e-10

    @pytest.mark.datarepo()
    def test_reduce_workflow(self):
        template_path = "data/template.xml"
        output_dir = "results"
        reduced_path = os.path.join(output_dir, "REFL_198409_combined_data_auto.txt")
        if os.path.exists(reduced_path):
            os.remove(reduced_path)

        for i in range(198409, 198417):
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

    @pytest.mark.datarepo()
    def test_reduce_workflow_201282(self):
        """
        Test to reproduce autoreduction output
        """
        template_path = "data/template_201282.xml"
        output_dir = "results"
        reduced_path = os.path.join(output_dir, "REFL_201282_combined_data_auto.txt")
        if os.path.exists(reduced_path):
            os.remove(reduced_path)

        for i in range(201282, 201289):
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

    @pytest.mark.datarepo()
    def test_background_subtraction(self):
        """
        Test with background subtraction off for the data and on for the normalization
        """
        template_path = "data/template_short_nobck.xml"
        output_dir = "results"
        reduced_path = os.path.join(output_dir, "REFL_198382_combined_data_auto.txt")
        if os.path.exists(reduced_path):
            os.remove(reduced_path)

        for i in range(198388, 198390):
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
