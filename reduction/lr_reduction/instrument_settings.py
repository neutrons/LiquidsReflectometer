from dataclasses import dataclass


@dataclass
class InstrumentSettings:
    """
    Dataclass to store the instrument geometry parameters.

    Attributes
    ----------
    source_detector_distance : float - Distance from the source to the detector in meters
    sample_detector_distance : float - Distance from the sample to the detector in meters
    num_x_pixels : int - Number of pixels in the x direction
    num_y_pixels : int - Number of pixels in the y direction
    pixel_width : float - Width of the pixels in mm
    xi_reference : float - xi reference value in mm
    s1_sample_distance : float - Distance from the sample to the s1 aperture in mm
    """

    source_detector_distance: float = None
    sample_detector_distance: float = None
    num_x_pixels: int = None
    num_y_pixels: int = None
    pixel_width: float = None
    xi_reference: float = None
    s1_sample_distance: float = None
