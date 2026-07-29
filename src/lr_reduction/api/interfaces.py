from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Generic, TypeVar

from lr_reduction.models.config import ReductionConfig

T = TypeVar("T")


class Entrypoint(ABC, Generic[T]):
    """
    Template Method skeleton shared by every reduction invocation surface (§6.4, §8).

    A fixed :meth:`execute` runs the steps in a guaranteed order; each concrete
    entrypoint overrides only the steps that differ.
    """

    @abstractmethod
    def load_configuration(self) -> ReductionConfig:
        """Resolve and load this entrypoint's configuration (§3.3.8.1)."""

    @abstractmethod
    def load_data(self, config: ReductionConfig) -> Any:
        """Obtain the run data this entrypoint operates on, given its configuration."""

    @abstractmethod
    def call_operations(self, config: ReductionConfig, data: Any) -> T:
        """Drive the operations of §1.2 and return their result. No I/O."""

    @abstractmethod
    def save_output(self, result: T) -> None:
        """Persist (and/or publish) *result*."""

    def execute(self) -> T:
        """Invariant skeleton: configuration, then data, then operations, then output."""
        config = self.load_configuration()
        data = self.load_data(config)
        result = self.call_operations(config, data)
        self.save_output(result)
        return result
