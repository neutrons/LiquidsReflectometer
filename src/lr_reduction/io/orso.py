"""ORSO format I/O module."""

# TODO:
#   - load_partials(partial_dir: Path, sequence_id: int) -> list[SingleRunResult] discovering REF_L_{sequence_id}_*_partial.ort,
#   returning an empty list for now.

import datetime

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
    ValueRange,
)

from lr_reduction import __version__ as lr_reduction_version
from lr_reduction.utils.logging import get_logger

logger = get_logger(__name__)

header = fileio.orso.Orso.empty()


def read(filepath: str):
    """Read an ORSO file and return a ReductionConfig populated with data from the ORSO file."""
    orso = fileio.orso.load_orso(filepath)
    logger.info(f"ORSO file read from {filepath}.")
    print(f"{orso=}")
    # TODO: parse header/comments, etc. into ReductionConfig
    return


def read_partials(partial_dir: str, sequence_id: int):
    """Discover and load partial ORSO files for a given sequence ID."""
    fn = f"{partial_dir}/REF_L_{sequence_id}_*_partial.ort"
    logger.info(f"Discovering partial ORSO files matching {fn}...")
    # TODO: implement discovery and loading of partial ORSO files into SingleRunResult
    return []  # return an empty list for now


# TODO: double check comment/metadata and fields/values/units with CISes
# TODO: parse Results (and ReductionConfig?) for info to populate ORSO fields
def write(sequence_id: str, sequence_number: str = ""):
    """Write an ORSO file for a given sequence ID and optional sequence number."""
    header = fileio.orso.Orso.empty()

    user = Person(name="PLACEHOLDER USER", affiliation="ORNL", contact="PLACEHOLDER EMAIL")

    experiment = Experiment(
        title="PLACEHOLDER TITLE",  # where does partial meta title come from? use that?
        instrument="Liquids Reflectometer",
        start_date=datetime.datetime.today(),
        probe="neutron",
        facility="SNS",  # ORNL?
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
    incident_angle = 90.0  # or [0.0, 90.0]?
    wavelength = 1.0  # or [0.0, 1.0]?
    instrument_settings = InstrumentSettings(
        incident_angle=incident_angle,
        wavelength=wavelength,
        # polarization=,
        configuration="PLACEHOLDER CONFIGURATION",
    )

    files = ["PLACEHOLDER FILES"]
    measurement = Measurement(
        instrument_settings=instrument_settings,
        data_files=files,
    )

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

    data_source = DataSource(owner=user, experiment=experiment, sample=sample, measurement=measurement)

    # Interpreted units are ["1/angstrom", "1/nm", "1", "1/s", None]
    q_column = fileio.base.Column(name="Q", unit="1/angstrom", physical_quantity="momentum transfer")
    r_column = fileio.base.Column(name="R", unit=None, physical_quantity="reflectivity")
    dr_column = fileio.base.ErrorColumn(error_of="dR", error_type="uncertainty", value_is="sigma")
    dq_column = fileio.base.ErrorColumn(error_of="dQ", error_type="resolution", value_is="sigma")
    header.columns = [q_column, r_column, dr_column, dq_column]

    # We can also make some data so that this code example will write something out
    q = np.array([0.01, 0.02, 0.03])
    R = np.array([0.1, 0.2, 0.3])
    dR = np.array([0.001, 0.002, 0.003])
    dq = q * 0.02

    reduction = Reduction(
        software=Software(name="lr_reduction", version=lr_reduction_version),
    )

    orso_class = Orso(data_source, reduction=reduction, columns=header.columns, comment=metadata)

    dataset = fileio.orso.OrsoDataset(info=orso_class, data=np.array([q, R, dR, dq]).T)

    fn = f"{sequence_id}_{sequence_number}.ort" if sequence_number else f"{sequence_id}.ort"
    fileio.orso.save_orso(datasets=[dataset], fname=fn)
