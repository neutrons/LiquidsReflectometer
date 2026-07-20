"""Legacy XML I/O module.

Reads legacy RefRed XML templates and converts them to a `ReductionConfig`, kept here for
backward compatibility during the coexistence period (this branch is `@deprecated`, see
`io/config_loader.py`). Per §3.3.10 the new workflow itself shall not read XML — the sanctioned
path from a legacy template into the new workflow is a standalone conversion utility (§3.3.10.2:
each legacy single direct-beam run is wrapped as a one-run composite direct beam).
"""

from lr_reduction.models.config import ReductionConfig
from lr_reduction.utils import get_logger

logger = get_logger(__name__)


def read_config(filepath: str) -> ReductionConfig:
    """Read a legacy XML template and convert it to a ReductionConfig."""
    logger.info(f"Reading legacy XML configuration from {filepath}")
    # Placeholder implementation; map legacy_reduction_template_reader.from_xml() output
    # (a list of ReductionParameters) onto ReflectedRunConfig/DirectBeamConfig entries
    ...
