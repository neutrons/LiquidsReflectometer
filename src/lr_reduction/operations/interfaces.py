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
        """Release what `process` created and did not return; called after `process` succeeds or fails."""
        pass

    def execute(self) -> U:
        """Validate the input, process it, and clean up, also when processing fails.

        `cleanup` runs in a `finally`, so it must not raise: an exception from it would replace
        the one `process` raised.
        """
        self.validate_input()
        try:
            return self.process()
        finally:
            self.cleanup()
