from __future__ import annotations

from abc import ABC, abstractmethod

from lr_reduction.models.config import ReductionConfig


class Entrypoint[T, S](ABC):
    """
    Template Method skeleton shared by every reduction invocation surface.

    A fixed :meth:`execute` runs the steps in a guaranteed order; each concrete
    entrypoint overrides only the steps that differ.

    Type parameters:
        T: The type of the input data this entrypoint operates on.
        S: The type of the result this entrypoint produces.

    Usage:

        class MyEntrypoint(Entrypoint[MyInput, MyResult]):
            def load_configuration(self) -> ReductionConfig:
                ...

            def load_data(self, config: ReductionConfig) -> MyInput:
                ...

            def call_operations(self, data: MyInput, config: ReductionConfig) -> MyResult:
                ...

            def save_output(self, result: MyResult) -> None:
                ...
    """

    @abstractmethod
    def load_configuration(self) -> ReductionConfig:
        """Resolve and load this entrypoint's reduction configuration."""

    @abstractmethod
    def load_data(self, config: ReductionConfig) -> T:
        """Obtain the run data this entrypoint operates on, given its configuration."""

    @abstractmethod
    def call_operations(self, data: T, config: ReductionConfig) -> S:
        """Drive the operations and return their result. No I/O."""

    @abstractmethod
    def save_output(self, result: S) -> None:
        """Persist (and/or publish) `result`."""

    def execute(self) -> S:
        """Invariant skeleton: configuration, then data, then operations, then output."""
        config = self.load_configuration()
        data = self.load_data(config)
        result = self.call_operations(data, config)
        self.save_output(result)
        return result
