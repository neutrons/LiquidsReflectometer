from lr_reduction.models.config import ReflectedRunConfig
from lr_reduction.models.results import SingleRunResult
from lr_reduction.models.run_data import RunData
from lr_reduction.operations.interfaces import OperationInterface


class SingleRunReduction(OperationInterface[RunData, ReflectedRunConfig, SingleRunResult]):
    """
    A single run reduction operation.
    """

    def __init__(self, data: RunData, config: ReflectedRunConfig):
        super().__init__(data, config)

    def validate_input(self) -> None:
        """
        Validate the input data and configuration for the single run reduction.
        """
        # Implement validation logic here
        pass

    def process(self) -> SingleRunResult:
        """
        Perform the single run reduction processing and return the result.
        """
        # Implement processing logic here
        ...

    def cleanup(self) -> None:
        """
        Perform any necessary cleanup after processing.
        """
        # Implement cleanup logic here
        pass
