from lr_reduction.models.config import AssemblyConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult
from lr_reduction.operations.interfaces import OperationInterface


class CombineResults(OperationInterface[list[ReductionResult], AssemblyConfig, CombinedReductionResult]):
    """
    A class for assembling reflectivity data from multiple runs.

    Performs scaling and stitching of multiple partial reflectivity runs into a single composite reflectivity dataset.
    """

    def __init__(self, data: list[ReductionResult], config: AssemblyConfig):
        super().__init__(data, config)

    def validate_input(self) -> None:
        """
        Validate the input data and configuration for the reflectivity assembly.
        """
        # Implement validation logic here
        pass

    def process(self) -> CombinedReductionResult:
        """
        Perform the reflectivity assembly operations and return the result.
        """
        # Implement processing logic here
        ...

    def cleanup(self) -> None:
        """
        Perform any necessary cleanup after processing.
        """
        # Implement cleanup logic here
        pass
