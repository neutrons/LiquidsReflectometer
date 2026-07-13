from abc import ABC, abstractmethod
from typing import Generic, TypeVar

InT = TypeVar("InT")
OutT = TypeVar("OutT", None)


class OperationInterface(ABC, Generic[InT, OutT]):
    """Abstract base class for operations that can be executed with a configuration and data."""

    def __init__(self, data: InT) -> None:
        self.data = data
        self.result: OutT | None = None

    @abstractmethod
    def validate_input(self) -> None:
        """Validate the input data before processing."""
        pass

    @abstractmethod
    def process(self) -> OutT:
        """Perform the main processing logic and return the result."""
        pass

    @abstractmethod
    def cleanup(self) -> None:
        """Perform any necessary cleanup after processing."""
        pass

    def execute(self) -> OutT:
        """Execute the operation and return the result."""
        self.validate_input()
        self.result = self.process()
        self.cleanup()
        return self.result
