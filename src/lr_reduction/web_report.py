"""
Report class used to populate the web monitor
"""

import math
import time
from typing import NamedTuple, Optional, Union

import numpy as np
import plotly.graph_objs as go
import plotly.offline as pyo
import requests
from mantid.simpleapi import Integration, Rebin, RefRoi, SumSpectra, Transpose, logger
from plot_publisher import publish_plot

from lr_reduction import template
from lr_reduction.data_info import DataType
from lr_reduction.mantid_utils import SampleLogValues
from lr_reduction.reduction_template_reader import ReductionParameters
from lr_reduction.types import MantidWorkspace

XY_PLOT_ZOOM_X_RANGE = [25, 225]
XY_PLOT_ZOOM_Y_RANGE = [100, 200]
YTOF_PLOT_ZOOM_X_RANGE = [100, 200]


class ReportSections(NamedTuple):
    """
    HTML fragments generated from a reduction run.

    Attributes
    ----------
    run_meta_data
        HTML <div> with run-level metadata.
    plots
        List of HTML <div>:s containing plots of the data run.
    reduction_parameters
        HTML <div> describing the reduction parameters used.
    """
    run_meta_data: str
    plots: list[str]
    reduction_parameters: str


def html_wrapper(report: Union[str, None]) -> str:
    """Wraps a report (set of <div> elements) in a complete HTML document

    Adds the javascript engine (PlotLy.js) address, HTML head, and body tags.

    Parameters
    ----------
    report : str
        The HTML content to be wrapped. This should contain on or more <div> elements and possibly
        summary <table> elements.

    Returns
    -------
    str
        A complete HTML document as a string, with the provided report content embedded within the body.
    """
    js_version = pyo.get_plotlyjs_version()
    url = f"https://cdn.plot.ly/plotly-{js_version}.js"
    try:
        response = requests.head(url, timeout=5)
        if response.status_code != 200:
            raise requests.RequestException(f"CDN returned status {response.status_code}")
    except requests.RequestException:
        logger.error(f"Plotly.js version {js_version} not found, using version 3.0.0 instead")
        url = "https://cdn.plot.ly/plotly-3.0.0.js"

    prefix = f"""<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Plotly Chart</title>
        <script src="{url}"></script>
    </head>
    <body>

    """
    suffix = """

    </body>
    </html>
    """
    report = "" if report is None else report
    return prefix + report + suffix  # allow for report being `None`


def save_report(html_report: str, report_file: str):
    """Save report to a local file

    Parameters
    ----------
    html_report : str
        HTML report containing <div> and <table> elements with information and plots from data
        reduction of one run.
    report_file : str
        File path where the report will be saved as an HTML file.
    """
    logger.notice(f"Saving report to: {report_file}")
    with open(report_file, "w", encoding="utf-8") as f:
        f.write(html_wrapper(html_report))


def upload_report(html_report: str, run_number: Union[str, int]) -> Optional[requests.Response]:
    r"""Upload report to the live data server

    Parameters
    ----------
    html_report: str
        HTML report containing <div> and <table> elements with information and plots from data
        reduction of one run.
    run_number: str, int
        Run number (e.g. '123435').

    Returns
    -------
    Optional[requests.Response]
        `Response` object returned by the livedata server, or `None` if the function is unable to do find the
        library to generate the `request.post()`
    """
    logger.notice(f"Uploading report for run {run_number} to livedata server")
    return publish_plot("REF_L", run_number, files={"file": html_report})


def assemble_report(reflectivity_plot_div: str | None = None, report_sections: ReportSections | None = None) -> str:
    """Assemble report HTML snippets into on final HTML report

    Parameters
    ----------
    reflectivity_plot_div: str | None
        HTML snippet for the reflectivity plot
    report_sections: ReportSections | None
        HTML snippets for the report sections

    Returns
    -------
    str
        HTML report
    """
    logger.notice("Assembling web report")
    plot_html = "<div></div>"

    if reflectivity_plot_div is not None:
        plot_html += "<div>%s</div>\n" % reflectivity_plot_div

    if report_sections is not None:
        plot_html += "<div>%s</div>\n" % report_sections.run_meta_data
        plot_html += "<div>%s</div>\n" % report_sections.reduction_parameters
        plot_html += "<table style='width:100%'>\n"
        plots = [p for p in report_sections.plots if p is not None]
        # arrange plots in table with two plots per row
        for i in range(0, len(plots), 2):
            plot_html += "<tr>\n"
            plot_html += "<td>%s</td>\n" % plots[i]
            if i + 1 < len(plots):
                plot_html += "<td>%s</td>\n" % plots[i + 1]
            else:
                plot_html += "<td></td>\n"  # optional filler cell
            plot_html += "</tr>\n"
        plot_html += "</table>\n"
        plot_html += "<hr>\n"

    return plot_html


def generate_report_sections(
    workspace: MantidWorkspace,
    template_file: str | ReductionParameters | None,
    meta_data: dict = None,
) -> ReportSections:
    """
    Generate diagnostics plots and report metadata from a reduction workspace.

    Parameters
    ----------
    workspace: MantidWorkspace
        Mantid workspace
    template_file: str | ReductionParameters | None
        Template file path or pre-parsed template data, or `None` for a direct beam
    meta_data : dict
        Metadata to embed in the report, or `None` for a direct beam

    Returns
    -------
    ReportSections
        Container holding HTML fragments for the reduction report
        (parameters, plots, and run metadata).
    """
    sample_logs = SampleLogValues(workspace)
    data_type = DataType.from_workspace(workspace)
    logger.notice(f"  - Data type: {data_type.name}")

    try:
        sequence_number = sample_logs["sequence_number"]
    except Exception:  # noqa: BLE001
        sequence_number = 1

    try:
        number_events = workspace.getNumberEvents()
    except Exception:  # noqa: BLE001
        number_events = 0

    plots = []

    logger.notice(f"  - generating report [{number_events}]")

    if data_type == DataType.REFLECTED_BEAM:
        # Read template if needed
        if isinstance(template_file, str):
            template_data = template.read_template(template_file, sequence_number)
        else:
            template_data = template_file
        report = generate_report_section_reduction_parameters(workspace, template_data, meta_data)
    elif data_type == DataType.DIRECT_BEAM:
        report = generate_report_section_direct_beam_parameters(workspace)
        template_data = None
    else:
        logger.error("Invalid data type for report: %s", data_type.name)
        return ReportSections("", [], "")

    run_meta_data = generate_report_section_run_meta_data(workspace)

    try:
        plots = generate_report_plots(workspace, template_data, data_type)
    except Exception as e:  # noqa: BLE001
        logger.notice(f"Could not generate plots: {e}")
        logger.error("Could not generate plots", exc_info=True)

    logger.notice(f"  - report: {len(report)} {len(plots)}")

    return ReportSections(run_meta_data, plots, report)


def generate_report_section_run_meta_data(workspace: MantidWorkspace) -> str:
    """Generate an HTML table containing run information

    Parameters
    ----------
    workspace : mantid.api.Workspace
        Event workspace

    Returns
    -------
    str

    Notes
    -----
    This is useful for live reduction. For a finished run, WebMon fetches this from ONCat.
    """
    logger.notice("  - generating run meta data")

    meta = "<p>\n<table style='width:80%'>"
    meta += "<tr><td># events:</td><td>%s</td></tr>" % workspace.getNumberEvents()

    try:
        p_charge = SampleLogValues(workspace)["gd_prtn_chrg"]
        meta += "<tr><td>p-charge [uAh]:</td><td>%6.4g</td></tr>" % p_charge
    except KeyError:
        meta += "<tr><td>p-charge [uAh]:</td><td>N/A</td></tr>"
    meta += "</table>\n<p>\n"
    return meta


def generate_report_section_reduction_parameters(workspace: MantidWorkspace, template_data: ReductionParameters, meta_data: dict) -> str:
    """Generate HTML report from a reduced workspace and template data

    Parameters
    ----------
    workspace: MantidWorkspace
        Reflected beam workspace
    template_data
        Reduction parameters
    meta_data
        Reduction metadata

    Returns
    -------
    str
        Reduction configuration in the form of an HTML table
    """
    logger.notice("  - generating reduction parameters")

    sample_logs = SampleLogValues(workspace)
    direct_beam = template_data.norm_file
    two_backgrounds = meta_data["use_functional_bck"]

    meta = "<table style='width:80%'>"
    meta += "<tr><td>Run:</td><td><b>%s</b></td><td><b>Direct beam: %s</b></td></tr>" % (
        int(sample_logs["run_number"]),
        direct_beam,
    )
    meta += "<tr><td>Q-binning:</td><td>%s</td><td>-</td></tr>" % meta_data["q_summing"]
    if meta_data["q_summing"]:
        meta += "<tr><td>Specular peak:</td><td>%g</td><td>-</td></tr>" % (
            meta_data["specular_pixel"],
        )
    meta += "<tr><td>Peak range:</td><td>%s - %s</td><td>%s - %s</td></tr>" % (
        template_data.data_peak_range[0],
        template_data.data_peak_range[1],
        template_data.norm_peak_range[0],
        template_data.norm_peak_range[1],
    )
    meta += "<tr><td>Two backgrounds:</td><td>%s</td><td>-</td></tr>" % two_backgrounds
    meta += "<tr><td>Background:</td><td>%s - %s</td><td>%s - %s</td></tr>" % (
        template_data.background_roi[0],
        template_data.background_roi[1],
        template_data.norm_background_roi[0],
        template_data.norm_background_roi[1],
    )
    if two_backgrounds:
        meta += "<tr><td>Background2:</td><td>%s - %s</td><td>-</td></tr>" % (
            template_data.background_roi[2],
            template_data.background_roi[3],
        )
    meta += "<tr><td>X range:</td><td>%s - %s</td><td>%s - %s</td></tr>" % (
        template_data.data_x_range[0],
        template_data.data_x_range[1],
        template_data.norm_x_range[0],
        template_data.norm_x_range[1],
    )
    meta += "<tr><td>Sequence:</td><td>%s: %s/%s</td></tr>" % (
        sample_logs["sequence_id"],
        sample_logs["sequence_number"],
        sample_logs["sequence_total"],
    )

    meta += "<tr><td>Stitching Type:</td><td>%s</td></tr>" % (template_data.stitching_configuration.type.value)
    meta += "<tr><td>Stitching Q range:</td><td>%s - %s</td></tr>" % (
        template_data.stitching_configuration.scale_factor_qmin,
        template_data.stitching_configuration.scale_factor_qmax,
    )
    meta += "<tr><td>Stitching normalize first angle:</td><td>%s</td></tr>" % (
        template_data.stitching_configuration.normalize_first_angle
    )

    meta += "<tr><td>Report time:</td><td>%s</td></tr>" % time.ctime()
    meta += "</table>\n"

    meta += "<table style='width:100%'>"
    meta += "<tr><th>Wavelength</th><th>Q</th><th>Thi</th><th>Ths</th><th>Offset</th><th>Theta used</th></tr>"  # noqa E501
    meta += "<tr><td>%6.4g - %6.4g</td><td>%6.4g - %6.4g</td><td>%6.4g</td><td>%6.4g</td><td>%6.4g</td><td>%6.4g</td></tr>\n" % (
        meta_data["wl_min"],
        meta_data["wl_max"],
        meta_data["q_min"],
        meta_data["q_max"],
        sample_logs["thi"],
        sample_logs["ths"],
        template_data.angle_offset,
        meta_data["theta"] * 180.0 / np.pi,
    )
    meta += "</table>\n"
    return meta


def generate_report_section_direct_beam_parameters(workspace: MantidWorkspace) -> str:
    """Generate HTML report from a direct beam workspace

    Parameters
    ----------
    workspace: MantidWorkspace
        Direct beam workspace

    Returns
    -------
    str
        Reduction configuration in the form of an HTML table
    """
    logger.notice("  - generating direct beam parameters")

    sample_logs = SampleLogValues(workspace)

    meta = "<table style='width:80%'>"
    meta += "<tr><td>Run:</td><td><b>%s</b> (direct beam)</td></tr>" % (
        int(sample_logs["run_number"]),
    )
    meta += "<tr><td>Sequence:</td><td>%s: %s/%s</td></tr>" % (
        sample_logs["sequence_id"],
        sample_logs["sequence_number"],
        sample_logs["sequence_total"],
    )
    meta += "<tr><td>Report time:</td><td>%s</td></tr>" % time.ctime()
    meta += "</table>\n"

    return meta


def generate_report_plots(workspace: MantidWorkspace, template_data: ReductionParameters, data_type: DataType) -> list[str]:
    """
    Generate diagnostic plots from the event workspace

    The generated plots are:

    - X-Y plot
    - Y-TOF plot
    - Counts per Y pixel
    - Counts per X pixel
    - TOF distribution

    Parameters
    ----------
    workspace: MantidWorkspace
        Workspace for reflected or direct beam run
    template_data: ReductionParameters
        Reduction parameters
    data_type: DataType
        Data type (reflected or direct beam)

    Returns
    -------
    list[str]
        List of HTML snippets for the diagnostic plots
    """
    n_x = int(workspace.getInstrument().getNumberParameter("number-of-x-pixels")[0])
    n_y = int(workspace.getInstrument().getNumberParameter("number-of-y-pixels")[0])

    if data_type == DataType.DIRECT_BEAM:
        scatt_peak = None
        scatt_low_res = None
        scatt_bck = None
        tof_range_ms = None
        tof_zoom_range = None
    else:
        scatt_peak = template_data.data_peak_range
        scatt_low_res = template_data.data_x_range
        scatt_bck = template_data.background_roi
        tof_range_ms = [t / 1000.0 for t in template_data.tof_range]
        # convert to ms and add margins
        tof_zoom_range = [template_data.tof_range[0]/1000.0 - 5.0, template_data.tof_range[1]/1000.0 + 5.0]

    # X-Y plot
    xy_plot = None
    try:
        integrated = Integration(workspace)
        signal = integrated.extractY()
        z = np.reshape(signal, (n_x, n_y))
        xy_plot = _plot2d(
            z=z.T,
            x=list(range(n_x)),
            y=list(range(n_y)),
            x_range=scatt_low_res,
            y_range=scatt_peak,
            y_bck_range=scatt_bck,
            x_zoom_range=XY_PLOT_ZOOM_X_RANGE,
            y_zoom_range=XY_PLOT_ZOOM_Y_RANGE,
        )
    except Exception:  # noqa E722
        logger.warning("  - Could not generate XY plot")
        xy_plot = _plotText("Could not generate XY plot")

    logger.notice("  - generating X-TOF plot")
    # Y-TOF plot
    y_tof_plot = None
    try:
        tof_min = workspace.getTofMin()
        tof_max = workspace.getTofMax()
        workspace = Rebin(workspace, params="%s, 50, %s" % (tof_min, tof_max))
        # algorithm RefRoi sums up the intensities in a region of interest on a 2D detector
        # returns a MatrixWorkspace
        direct_summed = RefRoi(
            InputWorkspace=workspace,
            NXPixel=n_x,
            NYPixel=n_y,
            ConvertToQ=False,
            IntegrateY=False,
            OutputWorkspace="direct_summed",
        )
        signal = np.transpose(direct_summed.extractY())
        tof_axis = direct_summed.extractX()[0] / 1000.0
        tof_axis = (tof_axis[:-1] + tof_axis[1:]) / 2.0  # average TOF values

        y_tof_plot = _plot2d(
            z=signal,
            y=tof_axis,
            x=list(range(signal.shape[1])),
            x_range=scatt_peak,
            x_bck_range=scatt_bck,
            y_range=tof_range_ms,
            x_label="Y pixel",
            y_label="TOF (ms)",
            swap_axes=True,
            x_zoom_range=YTOF_PLOT_ZOOM_X_RANGE,
            y_zoom_range=tof_zoom_range,
        )
    except Exception:  # noqa E722
        logger.warning("  - Could not generate X-TOF plot")
        y_tof_plot = _plotText("Could not generate X-TOF plot")

    logger.notice("  - generating Y count distribution")
    # Count per Y pixel
    peak_pixels = None
    try:
        direct_summed = RefRoi(
            InputWorkspace=workspace,
            IntegrateY=False,
            NXPixel=n_x,
            NYPixel=n_y,
            ConvertToQ=False,
            OutputWorkspace="direct_summed",
        )
        integrated = Integration(direct_summed)
        integrated = Transpose(integrated)
        signal_y = integrated.readY(0)
        signal_x = np.arange(len(signal_y))
        peak_pixels = _plot1d(
            signal_x,
            signal_y,
            x_range=scatt_peak,
            bck_range=scatt_bck,  # TODO: handle two backgrounds
            x_label="Y pixel",
            y_label="Counts",
        )
    except Exception:  # noqa E722
        logger.warning("  - Could not generate Y count distribution")
        peak_pixels = _plotText(
            "Could not generate Y count distribution"
        )

    logger.notice("  - generating X count distribution")
    # Count per X pixel
    low_res_profile = None
    try:
        direct_summed = RefRoi(
            InputWorkspace=workspace,
            NXPixel=n_x,
            NYPixel=n_y,
            ConvertToQ=False,
            OutputWorkspace="direct_summed",
        )
        integrated = Integration(direct_summed)
        integrated = Transpose(integrated)
        signal_y = integrated.readY(0)
        signal_x = list(range(len(signal_y)))
        low_res_profile = _plot1d(
            signal_x,
            signal_y,
            y_log=False,
            x_range=scatt_low_res,
            x_label="X pixel",
            y_label="Counts",
        )
    except Exception:  # noqa E722
        logger.warning("  - Could not generate X count distribution")
        low_res_profile = _plotText("Could not generate X count distribution")

    # TOF distribution
    tof_dist = None
    try:
        workspace = SumSpectra(workspace)
        signal_x = workspace.readX(0) / 1000.0
        signal_y = workspace.readY(0)
        tof_dist = _plot1d(
            signal_x,
            signal_y,
            y_log=False,
            x_range=tof_range_ms,
            x_label="TOF (ms)",
            y_label="Counts",
        )
    except Exception:  # noqa E722
        logger.warning("  - Could not generate TOF distribution")
        tof_dist = _plotText(
            "Could not generate TOF distribution"
        )

    return [xy_plot, y_tof_plot, peak_pixels, low_res_profile, tof_dist]


def _plot2d(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray | list,
    x_range: list = None,
    y_range: list = None,
    x_label: str = "X pixel",
    y_label: str = "Y pixel",
    title: str = "",
    x_bck_range: list = None,
    y_bck_range: list = None,
    swap_axes: bool = False,
    x_zoom_range: list = None,
    y_zoom_range: list = None,
):
    """
    Generate a simple 2D plot as an HTML snippet containing a Plotly graph embedded within a web page

    Parameters
    ----------
    x
        x-axis values
    y
        y-axis values
    z
        z-axis counts
    x_range
        x-axis range for the plot
    y_range
        y-axis range for the plot
    x_label
        label for the x-axis
    y_label
        label for the y-axis
    title
        the title of the plot
    x_bck_range
        x-axis background range for the plot
    y_bck_range
        y-axis background range for the plot
    swap_axes
        whether to swap the x and y axes

    Returns
    -------
    str
    """
    colorscale = [
        [0, "rgb(0,0,131)"],
        [0.125, "rgb(0,60,170)"],
        [0.375, "rgb(5,255,255)"],
        [0.625, "rgb(255,255,0)"],
        [0.875, "rgb(250,0,0)"],
        [1, "rgb(128,0,0)"],
    ]

    if swap_axes:
        x, y = y, x
        z = z.T
        x_range, y_range = y_range, x_range
        x_bck_range, y_bck_range = y_bck_range, x_bck_range
        x_label, y_label = y_label, x_label
        x_zoom_range, y_zoom_range = y_zoom_range, x_zoom_range

    # Eliminate items in array Z that are not finite and below a certain threshold
    x_grid, y_grid = np.meshgrid(x, y)
    x_flat, y_flat, z_flat = x_grid.flatten(), y_grid.flatten(), z.flatten()
    threshold = 0.01 * np.max(z_flat)
    mask = np.isfinite(z_flat) & (z_flat > threshold)  # Keep only significant values (exclude NaN, -inf, inf)
    x_sparse, y_sparse, z_sparse = x_flat[mask], y_flat[mask], z_flat[mask]

    # Round the remaining values to a certain number of decimal places, for instance 0.003455245 to 0.0034.
    # This will later save disk space when writing the figure to file
    def leading_decimal_places(x: float):
        """Calculate the number of leading decimal places for a number between 0 and 1.

        Returns 0 if x is outside the range (0, 1).
        """
        if x <= 0 or x >= 1:
            return 0  # Safe default for out-of-range values
        return abs(math.floor(math.log10(x)))

    z_sparse = np.round(z_sparse, 1 + leading_decimal_places(threshold))

    plotly_objects = []

    heatmap = go.Heatmap(
        x=x_sparse,
        y=y_sparse,
        z=z_sparse,
        autocolorscale=False,
        showscale=False,
        hoverinfo="x+y+z",
        colorscale=colorscale,
    )
    plotly_objects.append(heatmap)

    x_range_color = "rgba(152, 0, 0, .8)"
    y_range_color = "rgba(0, 128, 0, 1)"
    if swap_axes:
        x_range_color = "rgba(0, 128, 0, 1)"
        y_range_color = "rgba(152, 0, 0, .8)"

    # Set the color scale limits
    if x_range is not None:
        x_left = go.Scatter(
            name="",
            x=[x_range[0], x_range[0]],
            y=[min(y), max(y)],
            marker=dict(
                color=x_range_color,
            ),
        )
        x_right = go.Scatter(
            name="",
            x=[x_range[1], x_range[1]],
            y=[min(y), max(y)],
            marker=dict(
                color=x_range_color,
            ),
        )
        plotly_objects.append(x_left)
        plotly_objects.append(x_right)

    if x_bck_range is not None:
        x_left = go.Scatter(
            name="",
            x=[x_bck_range[0], x_bck_range[0]],
            y=[min(y), max(y)],
            marker=dict(
                color="rgba(152, 152, 152, .8)",
            ),
        )
        x_right = go.Scatter(
            name="",
            x=[x_bck_range[1], x_bck_range[1]],
            y=[min(y), max(y)],
            marker=dict(
                color="rgba(152, 152, 152, .8)",
            ),
        )
        plotly_objects.append(x_left)
        plotly_objects.append(x_right)

    if y_range is not None:
        y_left = go.Scatter(
            name="",
            y=[y_range[0], y_range[0]],
            x=[min(x), max(x)],
            marker=dict(
                color=y_range_color,
            ),
        )
        y_right = go.Scatter(
            name="",
            y=[y_range[1], y_range[1]],
            x=[min(x), max(x)],
            marker=dict(
                color=y_range_color,
            ),
        )
        plotly_objects.append(y_left)
        plotly_objects.append(y_right)

    if y_bck_range is not None:
        y_left = go.Scatter(
            name="",
            y=[y_bck_range[0], y_bck_range[0]],
            x=[min(x), max(x)],
            marker=dict(
                color="rgba(152, 152, 152, .8)",
            ),
        )
        y_right = go.Scatter(
            name="",
            y=[y_bck_range[1], y_bck_range[1]],
            x=[min(x), max(x)],
            marker=dict(
                color="rgba(152, 152, 152, .8)",
            ),
        )
        plotly_objects.append(y_left)
        plotly_objects.append(y_right)

    x_layout = dict(
        title=x_label,
        zeroline=False,
        exponentformat="power",
        showexponent="all",
        showgrid=True,
        showline=True,
        mirror="all",
        ticks="inside",
        range=x_zoom_range,
    )
    y_layout = dict(
        title=y_label,
        zeroline=False,
        exponentformat="power",
        showexponent="all",
        showgrid=True,
        showline=True,
        mirror="all",
        ticks="inside",
        range=y_zoom_range,
    )
    layout = go.Layout(
        title=title,
        showlegend=False,
        autosize=True,
        width=300,
        height=300,
        margin=dict(t=40, b=40, l=40, r=20),
        hovermode="closest",
        bargap=0,
        xaxis=x_layout,
        yaxis=y_layout,
    )
    fig = go.Figure(data=plotly_objects, layout=layout)

    div_html = pyo.plot(fig, output_type="div", include_plotlyjs=False, show_link=False)

    return div_html


def _plot1d(
    x, y, x_range=None, y_range=None, x_label="", y_label="Counts", title="", bck_range=None, x_log=False, y_log=True
) -> str:
    """Generate a simple 1D plot as an HTML snippet containing a Plotly graph embedded within a web page

    Parameters
    ----------
    x: np.ndarray
        x-axis values
    y: np.ndarray
    x_label: str
        x-axis label
    y_label: str
        y-axis label
    title: str
        plot title
    bck_range: np.ndarray, list
        array of length 2 to specify a background region in x
    x_log: bool
        if true, the x-axis will be a log scale
    y_log: bool
        if true, the y-axis will be a log scale

    Returns
    -------
    str
    """
    data = [go.Scatter(name="", x=x, y=y)]
    range_min_default = 0.001

    if x_range is not None:
        y_pos = [v for v in y if v > 0]
        min_y = min(y_pos) if y_pos else range_min_default
        x_left = go.Scatter(
            name="",
            x=[x_range[0], x_range[0]],
            y=[min_y, max(y)],
            marker=dict(
                color="rgba(152, 0, 0, .8)",
            ),
        )
        x_right = go.Scatter(
            name="",
            x=[x_range[1], x_range[1]],
            y=[min_y, max(y)],
            marker=dict(
                color="rgba(152, 0, 0, .8)",
            ),
        )
        data.append(x_left)
        data.append(x_right)

    if y_range is not None:
        x_pos = [v for v in x if v > 0]
        min_x = min(x_pos) if x_pos else range_min_default
        y_left = go.Scatter(
            name="",
            y=[y_range[0], y_range[0]],
            x=[min_x, max(x)],
            marker=dict(
                color="rgba(0, 128, 0, 1)",
            ),
        )
        y_right = go.Scatter(
            name="",
            y=[y_range[1], y_range[1]],
            x=[min_x, max(x)],
            marker=dict(
                color="rgba(0, 128, 0, 1)",
            ),
        )
        data.append(y_left)
        data.append(y_right)

    if bck_range is not None:
        y_pos = [v for v in y if v > 0]
        min_y = min(y_pos) if y_pos else range_min_default
        x_left = go.Scatter(
            name="",
            x=[bck_range[0], bck_range[0]],
            y=[min_y, max(y)],
            marker=dict(
                color="rgba(152, 152, 152, .8)",
            ),
        )
        x_right = go.Scatter(
            name="",
            x=[bck_range[1], bck_range[1]],
            y=[min_y, max(y)],
            marker=dict(
                color="rgba(152, 152, 152, .8)",
            ),
        )
        data.append(x_left)
        data.append(x_right)

    x_layout = dict(
        title=x_label,
        zeroline=False,
        exponentformat="power",
        showexponent="all",
        showgrid=True,
        showline=True,
        mirror="all",
        ticks="inside",
    )
    if x_log:
        x_layout["type"] = "log"

    y_layout = dict(
        title=y_label,
        zeroline=False,
        exponentformat="power",
        showexponent="all",
        showgrid=True,
        showline=True,
        mirror="all",
        ticks="inside",
    )
    if y_log:
        y_layout["type"] = "log"

    layout = go.Layout(
        title=title,
        showlegend=False,
        autosize=True,
        width=300,
        height=300,
        margin=dict(t=40, b=40, l=40, r=20),
        hovermode="closest",
        bargap=0,
        xaxis=x_layout,
        yaxis=y_layout,
    )

    fig = go.Figure(data=data, layout=layout)
    return pyo.plot(fig, output_type="div", include_plotlyjs=False, show_link=False)


def _plotText(text, title=""):
    """
    Displays an informative message as a plot

    Parameters
    ----------
    text : str
        Text to display
    title : str
        Title of the text area

    Returns
    -------
    str
    """

    layout = go.Layout(
        annotations=[
            dict(
                text=text,
                font=dict(size=13, color="red"),
                xref="paper",  # `paper` sets relative coordinates
                yref="paper",
                align="center",
                x=0.5,
                y=0.5,
                showarrow=False,
            )
        ],
        autosize=True,
        title=title,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        showlegend=False,
    )

    fig = go.Figure(layout=layout)
    plot = pyo.plot(fig, output_type="div", include_plotlyjs=False, show_link=False)

    return plot
