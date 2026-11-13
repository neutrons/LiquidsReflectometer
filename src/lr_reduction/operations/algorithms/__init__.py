"""Mantid algorithms of the reduction, each registered by `lr_reduction.utils.registration.register` on import.

Importing this package registers all of them, so they can be called through `mantid.simpleapi`.
"""

from lr_reduction.operations.algorithms import (
    apply_dead_time_correction,  # noqa: F401
    calculate_dead_time_correction,  # noqa: F401
    cd_attenuation_correction,  # noqa: F401
    convert_tof_to_wavelength,  # noqa: F401
    merge_direct_beam_segments,  # noqa: F401
    roi_extraction,  # noqa: F401
    slit_scan_correction,  # noqa: F401
    trim_direct_beam_tof,  # noqa: F401
)
