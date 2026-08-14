from datetime import datetime

from mantid import __version__ as mantid_version
from orsopy.fileio import Software

from lr_reduction import __version__ as lr_reduction_version
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult, ReflectivityCurve
from lr_reduction.operations.interfaces import OperationInterface
from lr_reduction.utils import get_logger

logger = get_logger(__name__)

_LR_REDUCTION_SOFTWARE = Software(name="lr_reduction", version=lr_reduction_version)
_MANTID_SOFTWARE = Software(name="mantid", version=mantid_version)


class CombineResultsOperation(OperationInterface[list[ReductionResult], ReductionConfig, CombinedReductionResult]):
    """
    A class for assembling reflectivity data from multiple runs.

    Performs scaling and stitching of multiple partial reflectivity runs into a single composite reflectivity dataset.
    """

    def __init__(self, data: list[ReductionResult], config: ReductionConfig):
        super().__init__(data, config)

    def validate_input(self) -> None:
        """
        Validate the input data and configuration for the reflectivity assembly.
        """
        # TODO: Implement validation logic here
        pass

    def process(self) -> CombinedReductionResult:
        """
        Perform the reflectivity assembly operations and return the result.
        """
        # TODO: Implement processing logic here
        return CombinedReductionResult(
            curve=ReflectivityCurve.empty(),
            reduction_config=self.config,
            lr_reduction_info=_LR_REDUCTION_SOFTWARE,
            mantid_info=_MANTID_SOFTWARE,
            reduction_timestamp=datetime.now(),
            partials=[],
            scale_factors=[],
            html_report=None,
        )

    def cleanup(self) -> None:
        """
        Perform any necessary cleanup after processing.
        """
        # TODO: Implement cleanup logic here
        pass
