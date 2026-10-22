import os

import mantid.simpleapi as mtd_api
import numpy as np

mtd_api.config["default.facility"] = "SNS"
mtd_api.config["default.instrument"] = "REF_L"

from lr_reduction.scaling_factors import workflow as sf_workflow
from lr_reduction.utils import amend_config


def check_results(data_file, reference):
    """
    Check scaling factor file output against reference
    """
    # Read data and skip header
    with open(data_file, "r") as fd:
        _cfg_data = fd.readlines()
        cfg_data = []
        for l in _cfg_data:
            if not l.startswith("#"):
                cfg_data.append(l)

    with open(reference, "r") as fd:
        _cfg_ref = fd.readlines()
        cfg_ref = []
        for l in _cfg_ref:
            if not l.startswith("#"):
                cfg_ref.append(l)

    for i in range(len(cfg_ref)):
        # Newly generated data
        toks = cfg_data[i].split(" ")
        for t in toks:
            kv = t.split("=")

        # Reference data
        toks = cfg_ref[i].split(" ")
        for t in toks:
            kv_ref = t.split("=")

        for j in range(len(kv_ref)):
            v_calc = float(kv[1])
            v_ref = float(kv_ref[1])
            delta = np.fabs((v_ref - v_calc) / v_ref)
            assert delta < 0.02


def test_compute_sf(nexus_dir):
    """
    Test the computation of scaling factors
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_197912")

    output_dir = "/tmp"

    # We are passing the first run of the set. For the autoreduction,
    # we would be missing runs from the complete set so we will want to
    # wait for the whole set to be acquired.
    output = sf_workflow.process_scaling_factors(ws, output_dir, use_deadtime=False, wait=True, postfix="_test")
    assert output is False

    output_cfg = os.path.join(output_dir, "sf_197912_Si_test.cfg")
    if os.path.isfile(output_cfg):
        os.remove(output_cfg)

    output = sf_workflow.process_scaling_factors(ws, output_dir, use_deadtime=False, wait=False, postfix="_test")
    assert output is True

    check_results(output_cfg, "data/sf_197912_Si_auto.cfg")


def test_compute_sf_with_deadtime(nexus_dir):
    """
    Test the computation of scaling factors
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_197912")

    output_dir = "/tmp"

    output_cfg = os.path.join(output_dir, "sf_197912_Si_test_dt.cfg")
    if os.path.isfile(output_cfg):
        os.remove(output_cfg)

    output = sf_workflow.process_scaling_factors(ws, output_dir, use_deadtime=True, wait=False, postfix="_test_dt")
    assert output is True

    check_results(output_cfg, "data/sf_197912_Si_dt_par_42_200.cfg")


def test_compute_sf_with_deadtime_tof_300(nexus_dir):
    """
    Test the computation of scaling factors
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_197912")

    output_dir = "/tmp"

    output_cfg = os.path.join(output_dir, "sf_197912_Si_test_dt.cfg")
    if os.path.isfile(output_cfg):
        os.remove(output_cfg)

    output = sf_workflow.process_scaling_factors(
        ws,
        output_dir,
        use_deadtime=True,
        deadtime=4.6,
        deadtime_tof_step=300,
        paralyzable=False,
        wait=False,
        postfix="_test_dt",
    )
    assert output is True

    check_results(output_cfg, "data/sf_197912_Si_dt_par_46_300.cfg")


def test_compute_sf_with_deadtime_tof_200(nexus_dir):
    """
    Test the computation of scaling factors
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_197912")

    output_dir = "/tmp"

    output_cfg = os.path.join(output_dir, "sf_197912_Si_test_dt.cfg")
    if os.path.isfile(output_cfg):
        os.remove(output_cfg)

    output = sf_workflow.process_scaling_factors(
        ws,
        output_dir,
        use_deadtime=True,
        deadtime=4.6,
        deadtime_tof_step=200,
        paralyzable=False,
        wait=False,
        postfix="_test_dt",
    )
    assert output is True

    check_results(output_cfg, "data/sf_197912_Si_dt_par_46_200.cfg")


def test_compute_sf_with_deadtime_tof_200_sort(nexus_dir):
    """
    Test the computation of scaling factors
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_197912")

    output_dir = "/tmp"

    output_cfg = os.path.join(output_dir, "sf_197912_Si_test_dt.cfg")
    if os.path.isfile(output_cfg):
        os.remove(output_cfg)

    output = sf_workflow.process_scaling_factors(
        ws,
        output_dir,
        order_by_runs=False,
        use_deadtime=True,
        deadtime=4.6,
        deadtime_tof_step=200,
        paralyzable=False,
        wait=False,
        postfix="_test_dt",
    )
    assert output is True

    check_results(output_cfg, "data/sf_197912_Si_dt_par_46_200.cfg")
