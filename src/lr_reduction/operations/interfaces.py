from abc import ABC, abstractmethod


class OperationInterface[T, S, U](ABC):
    """Abstract base class for operations that can be executed with a configuration and data.

    Type parameters:
        T: The type of the input data.
        S: The type of the configuration.
        U: The type of the result produced by the operation.
    """

    def __init__(self, data: T, config: S) -> None:
        self.data = data
        self.config = config

    @abstractmethod
    def validate_input(self) -> None:
        """Validate the input data before processing."""
        pass

    @abstractmethod
    def process(self) -> U:
        """Perform the main processing logic and return the result."""
        pass

    @abstractmethod
    def cleanup(self) -> None:
        """Perform any necessary cleanup after processing."""
        pass

    def execute(self) -> U:
        """Execute the operation and return the result."""
        self.validate_input()
        result = self.process()
        self.cleanup()
        return result
