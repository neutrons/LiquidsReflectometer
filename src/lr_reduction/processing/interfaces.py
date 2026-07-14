from abc import ABC, abstractmethod
from typing import Generic, TypeVar

DataT = TypeVar("DataT")
ConfigT = TypeVar("ConfigT")
OutT = TypeVar("OutT")


class OperationInterface(ABC, Generic[DataT, ConfigT, OutT]):
    """Abstract base class for operations that can be executed with a configuration and data."""

    def __init__(self, data: DataT, config: ConfigT) -> None:
        self.data = data
        self.config = config

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
        result = self.process()
        self.cleanup()
        return result
