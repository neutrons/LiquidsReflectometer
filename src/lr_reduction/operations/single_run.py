from datetime import datetime

from mantid import __version__ as mantid_version
from orsopy.fileio import Software

from lr_reduction import __version__ as lr_reduction_version
from lr_reduction.models.config import ReflectedRunConfig
from lr_reduction.models.results import ReductionResult, ReflectivityCurve
from lr_reduction.models.run_data import RunData
from lr_reduction.operations.interfaces import OperationInterface
from lr_reduction.utils import get_logger

logger = get_logger(__name__)

_LR_REDUCTION_SOFTWARE = Software(name="lr_reduction", version=lr_reduction_version)
_MANTID_SOFTWARE = Software(name="mantid", version=mantid_version)


class SingleRunReductionOperation(OperationInterface[RunData, ReflectedRunConfig, ReductionResult]):
    """
    A single run reduction operation.
    """

    def __init__(self, data: RunData, config: ReflectedRunConfig):
        super().__init__(data, config)

    def validate_input(self) -> None:
        """
        Validate the input data and configuration for the single run reduction.
        """
        # TODO: Implement validation logic here
        ...

    def process(self) -> ReductionResult:
        """
        Perform the single run reduction processing and return the result.
        """
        # TODO: Implement the actual processing logic here. For now, we return a placeholder result.
        return ReductionResult(
            curve=ReflectivityCurve.empty(),
            run_numbers=[],
            sequence_id=0,
            sequence_number=0,
            reduction_config=self.config,
            lr_reduction_info=_LR_REDUCTION_SOFTWARE,
            mantid_info=_MANTID_SOFTWARE,
            reduction_timestamp=datetime.now(),
        )

    def cleanup(self) -> None:
        """
        Perform any necessary cleanup after processing.
        """
        # TODO: Implement cleanup logic here
        ...
