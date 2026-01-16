import os

import pytest

from lr_autoreduce.reduce_REF_L import autoreduce
from lr_reduction import template
from lr_reduction.utils import amend_config


@pytest.mark.datarepo
def test_autoreduce_reflected(nexus_dir, template_dir, tmp_path):
    """Test REF_L autoreduction for reflected data"""
    events_file = "REF_L_201282.nxs.h5"
    output_dir = tmp_path / "output"
    template_file = os.path.join(template_dir, "template_201282.xml")
    template_data = template.read_template(template_file, 1)
    template_data.scaling_factor_flag = False

    with amend_config(data_dir=nexus_dir):
        autoreduce(events_file, output_dir, template_data, publish=False)

    # Check that output files are created
    assert (output_dir / "REF_L_201282.html").exists()
    assert (output_dir / "REFL_201282_1_201282_partial.txt").exists()
    assert (output_dir / "REFL_201282_combined_data_auto.txt").exists()


@pytest.mark.datarepo
def test_autoreduce_direct(nexus_dir, template_dir, tmp_path):
    """Test REF_L autoreduction for direct beam data"""
    events_file = "REF_L_201043.nxs.h5"
    output_dir = tmp_path / "output"
    template_file = os.path.join(template_dir, "template_201282.xml")
    template_data = template.read_template(template_file, 1)

    with amend_config(data_dir=nexus_dir):
        autoreduce(events_file, output_dir, template_data, publish=False)

    # Check that the HTML report is created
    assert (output_dir / "REF_L_201043.html").exists()
