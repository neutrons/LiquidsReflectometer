from datetime import datetime

from mantid import __version__ as mantid_version
from orsopy.fileio import Software

from lr_reduction import __version__ as lr_reduction_version
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import ReductionResult, ReflectivityCurve
from lr_reduction.models.run_data import RunData
from lr_reduction.operations.interfaces import OperationInterface
from lr_reduction.types import ID, CompositeDirectBeam
from lr_reduction.utils import get_logger

logger = get_logger(__name__)

_LR_REDUCTION_SOFTWARE = Software(name="lr_reduction", version=lr_reduction_version)
_MANTID_SOFTWARE = Software(name="mantid", version=mantid_version)


class SingleRunReductionOperation(OperationInterface[RunData, ReductionConfig, ReductionResult]):
    """
    A single run reduction operation.

    `config` is the full, multi-run reduction configuration; `sequence_number` identifies which
    run within it this operation reduces (`config.runs[sequence_number]`). `ReductionConfig.for_run`
    could scope `config` down to a self-contained one-run config for ORSO-header embedding, but for
    now the full config is threaded through as-is and `sequence_number` is used directly instead.
    """

    def __init__(
        self,
        data: RunData,
        config: ReductionConfig,
        composite_direct_beam: CompositeDirectBeam,
        sequence_number: ID,
    ):
        super().__init__(data, config)
        self.composite_direct_beam = composite_direct_beam
        self.sequence_number = sequence_number

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
        run = self.config.runs[self.sequence_number]
        return ReductionResult(
            curve=ReflectivityCurve.empty(),
            run_numbers=run.resolved_source_runs,
            sequence_id=0,
            sequence_number=self.sequence_number,
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
