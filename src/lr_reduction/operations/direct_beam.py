from lr_reduction.models.config import DirectBeamConfig
from lr_reduction.models.run_data import RunData
from lr_reduction.operations.interfaces import OperationInterface
from lr_reduction.types import CompositeDirectBeam


class DirectBeamCompositionOperation(OperationInterface[list[RunData], DirectBeamConfig, CompositeDirectBeam]):
    """Combine multiple direct beam runs into a single composite direct beam result."""

    def __init__(self, data: list[RunData], config: DirectBeamConfig):
        super().__init__(data, config)

    def validate_input(self) -> None:
        """
        Validate the input data and configuration for the direct beam composition.
        """
        # Implement validation logic here
        pass

    def process(self) -> CompositeDirectBeam:
        """
        Perform the direct beam composition operations and return the result.
        """
        # Implement processing logic here
        ...

    def cleanup(self) -> None:
        """
        Perform any necessary cleanup after processing.
        """
        # Implement cleanup logic here
        pass
