from lr_reduction.models.config import DirectBeamConfig
from lr_reduction.models.run_data import RunData
from lr_reduction.operations.interfaces import OperationInterface
from lr_reduction.types import CompositeDirectBeam


class DirectBeamComposition(OperationInterface[list[RunData], DirectBeamConfig, CompositeDirectBeam]):
    """
    """

    def __init__(self, data: list[RunData], config: DirectBeamConfig):
        super().__init__(data, config)

    def validate_input(self) -> None:
        """
        Validate the input data and configuration for the single run reduction.
        """
        # Implement validation logic here
        pass

    def process(self) -> CompositeDirectBeam:
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
