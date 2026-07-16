from lr_reduction.models.results import SequenceResult, SingleRunResult
from lr_reduction.processing.interfaces import OperationInterface


# TODO: Change config type to AssemblyConfig when it is implemented
class ReflectivityAssembly(OperationInterface[list[SingleRunResult], None, SequenceResult]):
    """
    A class for assembling reflectivity data from multiple runs.

    Performs scaling and stitching of multiple partial reflectivity runs into a single composite reflectivity dataset.
    """

    def __init__(self, data: list[SingleRunResult], config: None):
        super().__init__(data, config)

    def validate_input(self) -> None:
        """
        Validate the input data and configuration for the reflectivity assembly.
        """
        # Implement validation logic here
        pass

    def process(self) -> SequenceResult:
        """
        Perform the reflectivity assembly processing and return the result.
        """
        # Implement processing logic here
        ...

    def cleanup(self) -> None:
        """
        Perform any necessary cleanup after processing.
        """
        # Implement cleanup logic here
        pass
