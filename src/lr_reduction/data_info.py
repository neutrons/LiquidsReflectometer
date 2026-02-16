"""
Meta-data information for LR reduction
"""

from enum import IntEnum

import numpy as np
from mantid.simpleapi import logger

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
        value = cls.REFLECTED_BEAM
        try:
            coordinate_system = CoordinateSystem.from_workspace(input_workspace)
            # Determine whether this is a direct beam based on the geometry
            if coordinate_system == CoordinateSystem.EARTH_CENTERED:
                thi = sample_logs["thi"]
                tthd = sample_logs["tthd"]
                if np.isclose(thi, tthd, atol=0.01):
                    value = cls.DIRECT_BEAM
            elif coordinate_system == CoordinateSystem.BEAM_CENTERED:
                # Beam-centered coordinate system
                ths = sample_logs["ths"]
                tthd = sample_logs["tthd"]
                if np.fabs(tthd) < 0.001 and np.fabs(ths) < 0.001:
                    value = cls.DIRECT_BEAM
            else:
                logger.warning("Unknown coordinate system; assuming reflected beam")
        except RuntimeError as e:
            logger.warning(f"Missing sample log {e}, assuming reflected beam")
        return value

    def __str__(self):
        return self.name


class CoordinateSystem(IntEnum):
    """
    Enum to represent the coordinate system used in the experiment.

    Attributes:
        UNKNOWN (int): Represents unknown coordinate system.
        EARTH_CENTERED (int): Represents the Earth-centered coordinate system.
        BEAM_CENTERED (int): Represents the Beam-centered coordinate system.
    """

    UNKNOWN = -1
    EARTH_CENTERED = 0
    BEAM_CENTERED = 1

    @classmethod
    def from_workspace(cls, input_workspace: MantidWorkspace):
        """
        Determine the coordinate system from the given workspace.
        """
        sample_logs = SampleLogValues(input_workspace)
        value = cls.UNKNOWN
        try:
            # This is the new log for coordinate system mode
            if "BL4B:CS:Mode:Coordinates" in sample_logs:
                if sample_logs["BL4B:CS:Mode:Coordinates"] == 0:
                    value = cls.EARTH_CENTERED
                else:
                    value = cls.BEAM_CENTERED
            # Fallback to older method using operating mode
            elif "BL4B:CS:ExpPl:OperatingMode" in sample_logs:
                if sample_logs["BL4B:CS:ExpPl:OperatingMode"] == "Free Liquid":
                    value = cls.EARTH_CENTERED
                else:
                    value = cls.BEAM_CENTERED
        except RuntimeError as e:
            logger.warning(f"Missing sample log {e}, unable to determine coordinate system")
        return value

    def __str__(self):
        return self.name
