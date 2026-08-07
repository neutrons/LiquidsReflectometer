from datetime import datetime

import numpy as np
from orsopy.fileio import Software

from lr_reduction.io.orso import write
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import ReductionResult


def test_read_orso(template_dir): ...


def test_read_partials(template_dir): ...


# TODO: read output file and check contents
def test_write_orso(tmp_path):
    """Test writing an ORSO file from a ReductionResult."""
    config = ReductionConfig()
    reduction_result = ReductionResult(
        reflectivity=np.array([0.1, 0.2, 0.3]),
        reduction_config=config,
        html_report=None,
        lr_reduction_info=Software(name="lr_reduction", version="1.0.0"),
        mantid_info=Software(name="mantid", version="1.0.0"),
        reduction_timestamp=datetime.fromisoformat("2024-01-01T12:00:00"),
    )
    result = write(results=reduction_result, output_dir=tmp_path, title="Test ORSO Output")
    assert result.exists()
