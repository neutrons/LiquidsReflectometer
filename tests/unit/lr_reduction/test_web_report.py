import os

import numpy as np
import pytest
from mantid.kernel import amend_config
from mantid.simpleapi import LoadEventNexus

from lr_reduction import template
from lr_reduction.data_info import DataType
from lr_reduction.event_reduction import EventReflectivity
from lr_reduction.mantid_utils import SampleLogValues
from lr_reduction.web_report import (
    assemble_report,
    generate_report_plots,
    generate_report_section_reduction_parameters,
    generate_report_section_run_meta_data,
    generate_report_sections,
)


@pytest.fixture(scope="module")
def workspace_sc(nexus_dir):
    """Fixture to create and return a Mantid workspace (LoadEventNexus) for testing."""
    with amend_config(data_dir=nexus_dir):
        ws = LoadEventNexus("REF_L_201288")
    return ws


@pytest.fixture(scope="module")
def workspace_db(nexus_dir):
    """Fixture to load a Mantid workspace for direct beam data."""
    with amend_config(data_dir=nexus_dir):
        ws_db = LoadEventNexus("REF_L_198399")
    return ws_db


@pytest.fixture(scope="module")
def template_data(workspace_sc, template_dir):
    """Fixture to create template data for testing."""
    sequence_number = workspace_sc.getRun().getProperty("sequence_number").value[0]
    template_path = os.path.join(template_dir, "template.xml")
    template_data = template.read_template(template_path, sequence_number)
    return template_data


@pytest.fixture(scope="module")
def meta_data(workspace_sc, workspace_db, template_data):
    """Fixture to create metadata dictionary for testing."""
    peak = template_data.data_peak_range
    peak_center = (peak[0] + peak[1]) / 2.0
    sample_logs = SampleLogValues(workspace_sc)
    theta = sample_logs["ths"] * np.pi / 180.0  # radians

    event_refl = EventReflectivity(workspace_sc,
                                   workspace_db,
                                   peak,
                                   template_data.background_roi,
                                   template_data.norm_peak_range,
                                   template_data.norm_background_roi,
                                   peak_center,
                                   template_data.data_x_range,
                                   template_data.norm_x_range,
                                   theta=theta)

    meta_data = event_refl.to_dict()
    meta_data["tof_weighted"] = False
    meta_data["bck_in_q"] = False
    meta_data["theta_offset"] = template_data.angle_offset
    return meta_data


def test_generate_report_section_reduction_parameters(workspace_sc, template_data, meta_data):
    report = generate_report_section_reduction_parameters(workspace_sc, template_data, meta_data)
    assert len(report) == 942


def test_generate_report_plots_reflected_beam(workspace_sc, template_data):
    html_plots = generate_report_plots(workspace_sc, template_data, DataType.REFLECTED_BEAM)
    assert len(html_plots) == 5
    assert None not in html_plots


def test_generate_report_plots_direct_beam(workspace_db, template_data):
    html_plots = generate_report_plots(workspace_db, template_data, DataType.DIRECT_BEAM)
    assert len(html_plots) == 5
    assert None not in html_plots


def test_generate_report_section_run_meta_data(workspace_sc):
    html_meta_data = generate_report_section_run_meta_data(workspace_sc)
    assert len(html_meta_data) == 132


def test_generate_report_sections(workspace_sc, template_data, meta_data):
    report_sections = generate_report_sections(workspace_sc, template_data, meta_data)
    assert report_sections.run_meta_data is not None
    assert report_sections.reduction_parameters is not None
    assert report_sections.plots is not None


def test_generate_report_section_direct_beam(workspace_db, template_data):
    report_sections = generate_report_sections(workspace_db, template_data)
    assert report_sections.run_meta_data is not None
    assert report_sections.reduction_parameters is not None
    assert report_sections.plots is not None


def test_assemble_report(workspace_sc, template_data, meta_data):
    report_sections = generate_report_sections(workspace_sc, template_data, meta_data)
    report = assemble_report(None, report_sections)
    assert "<div>" in report
    assert "</div>" in report
