import os

import numpy as np
import pytest
from mantid.kernel import amend_config
from mantid.simpleapi import LoadEventNexus

from lr_reduction import template
from lr_reduction.event_reduction import EventReflectivity
from lr_reduction.mantid_utils import SampleLogValues
from lr_reduction.web_report import (
    ReportGenerator,
    generate_event_count_info,
    generate_plots,
    generate_web_report,
    html_wrapper,
)


@pytest.fixture(scope="module")
def workspace_sc(nexus_dir):
    """Fixture to create a Report object for testing."""
    with amend_config(data_dir=nexus_dir):
        ws = LoadEventNexus("REF_L_201288")
    return ws


@pytest.fixture(scope="module")
def workspace_db(nexus_dir):
    """Fixture to create a Report object for testing."""
    with amend_config(data_dir=nexus_dir):
        ws_db = LoadEventNexus("REF_L_201051")
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
    """Fixture to create a Report object for testing."""
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


def test_generate_web_report(workspace_sc, template_data, meta_data):
    report = generate_web_report(workspace_sc, template_data, meta_data)
    with open("/home/u5z/output.html", "w") as file:
        file.write(report)
    # assert len(report) == 2


def test_generate_plots(workspace_sc, template_data):
    html_plots = generate_plots(workspace_sc, template_data)
    plot_html = "<table style='width:100%'>\n"
    plot_html += "<tr>\n"
    for plot in html_plots:
        if plot is not None:
            plot_html += "<td>%s</td>\n" % plot
    plot_html += "</tr>\n"
    plot_html += "</table>\n"
    with open("/home/u5z/plots_output.html", "w") as file:
        file.write(html_wrapper(plot_html))


def test_generate_event_count_info(workspace_sc):
    info_html = generate_event_count_info(workspace_sc)
    with open("/home/u5z/info_output.html", "w") as file:
        file.write(html_wrapper(info_html))


def test_generate_report(workspace_sc, template_data, meta_data):
    generator = ReportGenerator(workspace_sc, template_data, meta_data)
    report = generator.generate()
    assert report
