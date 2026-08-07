from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
import numpy.typing as npt
from orsopy.fileio import Software

from lr_reduction.models.config import ReductionConfig


@dataclass
class ReflectivityCurve:
    """The reduced R(Q) data carried by both result types.

    Columns are stored as named arrays; the on-disk column order (Q, R, dR, dQ,
    then T, L, dT, dL when present) is produced by :meth:`to_array`, so producers
    never assemble positional columns themselves.

    Attributes
    ----------
    q
        Momentum transfer [1/Angstrom]
    r
        Reflectivity
    dr
        Reflectivity uncertainty
    dq
        Q resolution [1/Angstrom]
    theta
        Incidence angle [deg]; the optional 8-column T column
    lam
        Wavelength [Angstrom]; the optional 8-column L column
    dtheta
        Incidence angle uncertainty [deg]
    dlam
        Wavelength uncertainty [Angstrom]
    """

    q: npt.NDArray[np.float64]
    r: npt.NDArray[np.float64]
    dr: npt.NDArray[np.float64]
    dq: npt.NDArray[np.float64]
    theta: npt.NDArray[np.float64] | None = None
    lam: npt.NDArray[np.float64] | None = None
    dtheta: npt.NDArray[np.float64] | None = None
    dlam: npt.NDArray[np.float64] | None = None

    def __post_init__(self):
        optional_columns = [self.theta, self.lam, self.dtheta, self.dlam]
        if any(column is not None for column in optional_columns) and not all(
            column is not None for column in optional_columns
        ):
            raise ValueError("theta, lam, dtheta, dlam must be provided together")

        columns = [
            self.q,
            self.r,
            self.dr,
            self.dq,
            *(column for column in optional_columns if column is not None),
        ]
        if any(column.ndim != 1 for column in columns):
            raise ValueError("reflectivity curve columns must be one-dimensional")
        if len({len(column) for column in columns}) != 1:
            raise ValueError("reflectivity curve columns must have equal lengths")

    @property
    def is_eight_column(self) -> bool:
        return self.theta is not None

    def to_array(self) -> npt.NDArray[np.float64]:
        """Return the curve as an (N, 4) or (N, 8) array in on-disk column order."""
        columns = [self.q, self.r, self.dr, self.dq]
        if self.is_eight_column:
            columns += [self.theta, self.lam, self.dtheta, self.dlam]
        return np.column_stack(columns)

    @classmethod
    def empty(cls) -> "ReflectivityCurve":
        """Return a zero-point curve, used by placeholder results until Op 2/Op 3 exist."""
        return cls(q=np.empty(0), r=np.empty(0), dr=np.empty(0), dq=np.empty(0))


@dataclass
class ReductionResult:
    """Result of reducing a single run: one ORSO partial.

    Attributes
    ----------
    curve
        The reduced R(Q) data; not necessarily Q-sorted
    run_numbers
        Constituent run number(s); more than one when the run is a sum of
        source runs, one data_files entry each
    sequence_id
        Identifies the reflectivity curve this run belongs to
    sequence_number
        1-based position of the run within the sequence
    reduction_config
        The one-run configuration that reproduces this partial when embedded
        in its ORSO header
    lr_reduction_info
        lr_reduction package information
    mantid_info
        mantid package information
    reduction_timestamp
        orsopy.fileio.reduction.timestamp
    """

    curve: ReflectivityCurve
    run_numbers: list[int]
    sequence_id: int
    sequence_number: int
    reduction_config: ReductionConfig
    lr_reduction_info: Software
    mantid_info: Software
    reduction_timestamp: datetime


@dataclass
class CombinedReductionResult:
    """Result of assembling a sequence: one combined ORSO file.

    Not a subclass of ReductionResult: a combined curve must not be accepted
    where a single-run partial is expected (e.g. as sequence-assembly input).

    Attributes
    ----------
    curve
        The assembled curve, Q-monotone
    partials
        The single-run results the curve was assembled from, whether reduced
        in-memory or loaded from on-disk partials; supplies the data_files
        and additional_files entries of the combined output
    scale_factors
        Per-angle stitching scale factors computed by assembly, recorded in
        the ORSO header
    reduction_config
        The full sequence configuration that reproduces this curve when
        embedded in the ORSO header
    html_report
        Diagnostic HTML report for the sequence
    lr_reduction_info
        lr_reduction package information
    mantid_info
        mantid package information
    reduction_timestamp
        orsopy.fileio.reduction.timestamp
    """

    curve: ReflectivityCurve
    reduction_config: ReductionConfig
    lr_reduction_info: Software
    mantid_info: Software
    reduction_timestamp: datetime
    partials: list[ReductionResult] = field(default_factory=list)
    scale_factors: list[float] = field(default_factory=list)
    html_report: str | None = None
