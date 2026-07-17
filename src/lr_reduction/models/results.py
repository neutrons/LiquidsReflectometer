from dataclasses import dataclass
from datetime import datetime

import numpy as np
import numpy.typing as npt
from orsopy.fileio import Software


@dataclass
class ReductionResult:
    """Result of reducing data from a single run.

    Attributes
    ----------
    TBD
    """
    reflectivity: npt.NDArray[np.float64]
    reduction_config: ReductionConfig
    html_report: str | None
    # Reduction metadata
    lr_reduction_info: Software         # lr_reduction package information
    mantid_info: Software               # mantid package information
    reduction_timestamp: datetime       # orsopy.fileio.reduction.timestamp

@dataclass
class CombinedReductionResult(ReductionBase):
    """Result of assembling reduced data from multiple runs.

    Attributes
    ----------
    TBD
    """
    assembly_config: AssemblyConfig
