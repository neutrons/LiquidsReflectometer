"""
Meta-data information for LR reduction
"""

from enum import IntEnum

import numpy as np

from lr_reduction.mantid_utils import SampleLogValues
from lr_reduction.typing import MantidWorkspace


class DataType(IntEnum):
    """
    Enum to represent the typical types of data for the Liquids Reflectometer

    Attributes:
        UNKNOWN (int): Represents unknown data TYPE, usually because of low neutron count in the associated run.
        REFLECTED_BEAM (int): Represents reflected beam data.
        DIRECT_BEAM (int): Represents direct beam data.
    """

    UNKNOWN = -1
    REFLECTED_BEAM = 0
    DIRECT_BEAM = 1

    @classmethod
    def from_workspace(cls, input_workspace: MantidWorkspace):
        """
        Determine the data type from the given workspace.

        Currently, the 'data_type' property is not set automatically in
        the DAS, therefore, use the incoming angle and scattering angle
        to determine the data type.
        """
        sample_logs = SampleLogValues(input_workspace)
        try:
            tthd = sample_logs["tthd"]
            ths = sample_logs["ths"]
            if np.fabs(tthd) < 0.001 and np.fabs(ths) < 0.001:
                value = cls.DIRECT_BEAM
            else:
                value = cls.REFLECTED_BEAM
        except Exception:  # noqa E722
            value = cls.REFLECTED_BEAM  # missing logs, assume reflected beam
        return value

    def __str__(self):
        return self.name
