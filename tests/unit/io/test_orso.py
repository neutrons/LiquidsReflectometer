import os
from datetime import datetime

import numpy as np
from orsopy.fileio import Software

from lr_reduction.io.orso import read, read_partials, write
from lr_reduction.types import SingleRunResult


def test_read_orso(template_dir):
    filepath = f"{template_dir}/REF_L_111111.ort"
    config = read(filepath)
    assert config.sequence_id == "PLACEHOLDER SEQUENCE ID"
    assert len(config.direct_beams) == 1
    assert len(config.reflected_runs) == 1


# TODO: add a test for read_partials that actually reads a partial file and checks the contents
def test_read_partials(template_dir):
    sequence_id = 111111
    partials = read_partials(template_dir, sequence_id)
    assert isinstance(partials, list)


# TODO: read output file and check contents
def test_write_orso(template_dir):
    # Create a dummy SingleRunResult
    single_run_result = SingleRunResult(
        reduction_output=np.array([[1.0, 2.0, 3.0, 4.0]]),
        sequence_id="test_sequence",
        sequence_number=1,
        lr_reduction_info=Software(name="lr_reduction", version="1.0"),
        mantid_info=Software(name="mantid", version="6.0"),
        reduction_timestamp=datetime.now(),
        experiment_id="IPTS-12345",
        configuration_yaml={"param": "value"},
    )

    # Write the ORSO file
    write(single_run_result, output_dir=template_dir)

    # Check that the ORSO file was created
    orso_file = f"{template_dir}/REF_L_test_sequence_1.ort"
    assert os.path.exists(orso_file)

    # Clean up the created ORSO file
    os.remove(orso_file)
