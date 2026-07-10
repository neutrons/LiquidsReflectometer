from dataclasses import dataclass
from datetime import datetime
from typing import Any, List, Optional, TypeAlias

import numpy as np
import numpy.typing as npt
from mantid.api import Workspace
from mantid.dataobjects import EventWorkspace
from orsopy.fileio import Software

MantidWorkspace = str | Workspace

"""
The assembled direct beam produced by Op 1 (§2.2)
Consumed by single-run reduction (§3.2.1, §4.1.1)
"""
CompositeDirectBeam: TypeAlias = EventWorkspace

@dataclass
class SingleRunResult:
    """
    The result of single-run reduction (§2.3) for one angle, and the hand-off to in-memory sequence assembly (§2.4.1.b).
    One run's reduced output plus the metadataneeded to write an ORSO partial (§4.2)
    (§11.3)
    """

    # (N, 4) with columns: [Q (Å - 1), R(Q) (dimensionless), δR(Q) (dimensionless), δQ (Å⁻1)]
    # (N, 8) if `eight_col` is enabled, with the following columns appended:
    #   [T (Å), L (Å), δT (Å), δL (Å)]
    #   as defined by the legacy 8-column output format (§4.2.2) and frozen by §7.2.
    # The per-run output is not required to be Q-sorted; ordering is the responsibility of sequence assembly (§2.4.1).
    # The column order shall match the on-disk partial file (§11.4) exactly.
    reduction_output: npt.NDArray[np.float64]

    sequence_id: str                            # the sequence id of the run
    sequence_number: int                        # the sequence number of the run
    lr_reduction_info: Software   # lr_reduction package information
    mantid_info: Software         # mantid package information
    reduction_timestamp: datetime               # orsopy.fileio.reduction.timestamp
    experiment_id: str                          # the SNS IPTS / proposal number (ORSO data_source.experiment.proposalID
    configuration_yaml: dict[str, Any]          # the complete configuration used for the reduction so
                                                # the reduced data is reproducible and self-describing

@dataclass
class SequenceResult:
    """
    SequenceResult (§4.3). The assembled curve: concatenated Q-monotone arrays plus the
    per-angle stitching scale factors (§4.3.3), written as the combined ORSO file.
    """
    sequence_id: str
    reflected_runs: List[str]                       # all reflected runs stitched into the curve
                                                    # (orsopy.data_source.measurement.data_files)
    assembly_result: npt.NDArray[np.float64]        # the merged R(Q) curve
    stitching_scale_factors: List[float]            # the per-angle scale factors used to stitch
                                                    # the reflected runs together (§4.3.3)
    html_report: str

    directed_beam_runs: Optional[List[str]] = None  # orsopy.data_source.measurement.additional_files
