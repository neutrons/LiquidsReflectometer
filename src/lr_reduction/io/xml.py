"""Legacy XML I/O module. Reads legacy RefRed XML templates as configuration input."""

from lr_reduction.models.config import ReductionConfig
from lr_reduction.utils import get_logger

logger = get_logger(__name__)


def read_config(filepath: str) -> ReductionConfig:
    """Read a legacy XML template and convert it to a ReductionConfig."""
    ...
