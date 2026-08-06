"""ORSO format I/O module."""

import datetime
from pathlib import Path

import numpy as np
from orsopy import fileio
from orsopy.fileio import (
    DataSource,
    Experiment,
    InstrumentSettings,
    Measurement,
    Orso,
    Person,
    Reduction,
    Sample,
    Software,
    Value,
)

from lr_reduction import __version__ as lr_reduction_version
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult
from lr_reduction.types import MantidWorkspace
from lr_reduction.utils.logging import get_logger

logger = get_logger(__name__)


def read_config(filepath: str) -> ReductionConfig:
    """Read an ORSO file and return a ReductionConfig populated with metadata from the ORSO file.

    Placeholder: real metadata extraction from the ORSO header's recorded configuration
    (§3.3.8.1) is not yet implemented; this fabricates a minimal valid single-run config from
    the run number embedded in the filename.
    """
    ...


def read_single_run(filepath: str) -> ReductionResult:
    """Read an ORSO file and return a ReductionResult populated with data from the ORSO file."""
    ...


def read_partials(partial_dir: str, sequence_id: int) -> list[ReductionResult]:
    """Discover and load partial ORSO files for a given sequence ID."""
    ...


# TODO: Update this function to extract relevant Mantid logs and return the ORSO metadata fields accordingly.
def _get_mantid_logs(workspace: MantidWorkspace) -> dict[str, str]:
    """Extract Mantid logs from a Mantid workspace and return them as a dictionary."""
    logs = {}
    for log in workspace.getRun().getLogData():
        logs[log.name()] = log.value()
    return logs


# TODO:
#   - Double check comment/metadata and fields/values/units with CISes
#   - Parse Results, ReductionConfig, and Nexus logs for info to populate ORSO fields
#   - Desired filename format:       ???
#   - Desired name and contact info: ???
#   - Determine source for metadata:
#       - experiment title
#       - sample name, description, and category
def write(
    results: ReductionResult | CombinedReductionResult,
    output_dir: str | Path = ".",
    title: str | None = None,
) -> str:
    """Write reduction results to an ORSO format file and return the path of the written file."""

    logger.info(f"Writing ORSO reduced data for {type(results).__name__} to {output_dir}")

    header = fileio.orso.Orso.empty()

    user = Person(
        name="TODO: PLACEHOLDER USER",
        affiliation="Oak Ridge National Laboratory",
        contact="TODO: PLACEHOLDER EMAIL",
    )

    if title is None:
        title = (
            f"REF_L_{results.reduction_config.runs[0].run_number}_{results.reduction_config.runs[0].sequence_number}"
        )

    experiment = Experiment(
        title=title,
        instrument="Liquids Reflectometer",
        start_date=datetime.datetime.today(),
        probe="neutron",
        facility="ORNL/SNS",  # ORNL?
        # proposalID = "PLACEHOLDER PROPOSAL",
        # doi = "PLACEHOLDER DOI"
    )

    sample = Sample(
        name="PLACEHOLDER SAMPLE",
        description="PLACEHOLDER DESCRIPTION",
        category="PLACEHOLDER CATEGORY",
    )

    # TODO: use Value/ValueRange instead when we start implementing
    #       see https://www.reflectometry.org/orsopy/orsopy.fileio.html#orsopy.fileio.Value
    incident_angle = Value(90.0)  # or [0.0, 90.0]?
    wavelength = Value(1.0)  # or [0.0, 1.0]?
    instrument_settings = InstrumentSettings(
        incident_angle=incident_angle,
        wavelength=wavelength,
        configuration="PLACEHOLDER CONFIGURATION",
    )

    files = ["PLACEHOLDER FILES"]
    measurement = Measurement(
        instrument_settings=instrument_settings,
        data_files=files,  # type: ignore
    )

    data_source = DataSource(owner=user, experiment=experiment, sample=sample, measurement=measurement)

    # Interpreted units are ["1/angstrom", "1/nm", "1", "1/s", None]
    q_column = fileio.base.Column(name="Q", unit="1/angstrom", physical_quantity="momentum transfer")
    r_column = fileio.base.Column(name="R", unit=None, physical_quantity="reflectivity")
    dr_column = fileio.base.ErrorColumn(error_of="dR", error_type="uncertainty", value_is="sigma")
    dq_column = fileio.base.ErrorColumn(error_of="dQ", error_type="resolution", value_is="sigma")
    header.columns = [q_column, r_column, dr_column, dq_column]

    # TODO: placeholder data for now, will be replaced with actual data from results
    q = np.array([0.01, 0.02, 0.03])
    r = np.array([0.1, 0.2, 0.3])
    dr = np.array([0.001, 0.002, 0.003])
    dq = q * 0.02

    # TODO: Parse reduction config to populate comment field in Reduction object.
    metadata = {
        "wl_min": 13.05504715524623,
        "wl_max": 16.39347348797596,
        "q_min": 0.008027891007325118,
        "q_max": 0.010080776946106923,
        "theta": 0.010472986057602888,
        "start_time": "2023-01-16T10:47:39.627455667",
        "experiment": "IPTS-29196",
        "run_number": "201282",
        "run_title": "Expt 8 Cu-B BF4 noEtOH Full OCV 1-201282-1.",
        "norm_run": 201043,
        "time": "Fri Feb 20 11:31:18 2026",
        "dq0": 0,
        "dq_over_q": 0.022671471200417816,
        "sequence_number": 1,
        "sequence_id": 201282,
        "q_summing": False,
        "specular_pixel": 141.0,
        "use_functional_bck": False,
        "scaling_factors": {"a": 1, "err_a": 0, "b": 0, "err_b": 0},
        "tof_weighted": False,
        "bck_in_q": False,
        "theta_offset": 0,
    }
    reduction = Reduction(
        software=Software(name="lr_reduction", version=lr_reduction_version),
        comment=metadata,
    )

    orso_class = Orso(
        data_source,
        reduction=reduction,
        columns=header.columns,
        # comment=f"metadata={metadata}",
    )

    dataset = fileio.orso.OrsoDataset(info=orso_class, data=np.array([q, r, dr, dq]).T)

    if isinstance(results, ReductionResult):
        fn = Path(output_dir) / f"REF_L_{results.sequence_id}_{results.sequence_number}.ort"
    else:  # CombinedReductionResult
        fn = Path(output_dir) / f"REF_L_{results.sequence_id}.ort"
    fileio.orso.save_orso(datasets=[dataset], fname=str(fn))
