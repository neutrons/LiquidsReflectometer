from abc import ABC, abstractmethod
from typing import Generic, TypeVar

T = TypeVar("T")  # Input data type
S = TypeVar("S")  # Configuration type
U = TypeVar("U")  # Output type


class OperationInterface(ABC, Generic[T, S, U]):
    """Abstract base class for operations that can be executed with a configuration and data."""

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
