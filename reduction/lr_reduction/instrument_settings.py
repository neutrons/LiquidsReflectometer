from dataclasses import dataclass


@dataclass
class InstrumentSettings:
    """
    Dataclass to store the instrument geometry parameters.
    Default values determined from settings.json

    Attributes
    ----------
    source_detector_distance : float - Distance from the source to the detector in meters
    sample_detector_distance : float - Distance from the sample to the detector in meters
    num_x_pixels : int - Number of pixels in the x direction
    num_y_pixels : int - Number of pixels in the y direction
    pixel_width : float - Width of the pixels in mm
    xi_reference : float - Zero-position of slit relative to sample, in mm
    s1_sample_distance : float - Distance from the sample to the s1 aperture in mm
    """

    apply_instrument_settings: bool = False
    source_detector_distance: float = 15.75
    sample_detector_distance: float = 1.83
    num_x_pixels: int = 256
    num_y_pixels: int = 304
    pixel_width: float = 0.70
    xi_reference: float = 445
    s1_sample_distance: float = 1.485
