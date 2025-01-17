import numpy as np
from lmfit.models import LinearModel


def find_ranges_without_overlap(r1, r2):
    """
    Returns the part of r1 that does not contain r2
    When summing pixels for reflectivity, include the full range,
    which means that for a range [a, b], b is included.
    The range that we return must always exclude the pixels
    included in r2.

    Parameters
    ----------
    r1 : list
        Range of pixels to consider
    r2 : list
        Range of pixels to exclude

    Returns
    -------
    list
        List of ranges that do not overlap with r2
    """
    x1, x2 = r1
    x3, x4 = r2

    # r2 completely outside r1
    if x4 < x1 or x3 > x2:
        return [r1]  # r1 is the range without r2

    # r2 is entirely within r1
    if x1 <= x3 and x2 >= x4:
        return [[x1, x3 - 1], [x4 + 1, x2]]  # ranges before and after r2

    # r2 overlaps r1 from the right
    if x1 <= x3:
        return [[x1, x3 - 1]]  # range before r2

    # r2 overlaps r1 from the left
    if x2 >= x4:
        return [[x4 + 1, x2]]  # range after r2

    return []  # no range without r2


def functional_background(
    ws, event_reflectivity, peak, bck, low_res, normalize_to_single_pixel=False, q_bins=None, wl_dist=None, wl_bins=None, q_summing=False
):
    """
    Estimate background using a linear function over a background range that may include the specular peak.
    In the case where the peak is included in the background range, the peak is excluded from the background.

    Parameters
    ----------
    ws : Mantid workspace
        Workspace containing the data
    event_reflectivity : EventReflectivity
        EventReflectivity object
    peak : list
        Range of pixels that define the peak
    bck : list
        Range of pixels that define the background. It contains 4 pixels, defining up to two ranges.
    low_res : list
        Range in the x direction on the detector
    normalize_to_single_pixel : bool
        If True, the background is normalized to the number of pixels used to integrate the signal
    q_bins : numpy.ndarray
        Array of Q bins
    wl_dist : numpy.ndarray
        Wavelength distribution for the case where we use weighted events for normatization
    wl_bins : numpy.ndarray
        Array of wavelength bins for the case where we use weighted events for normatization
    q_summing : bool
        If True, sum the counts in Q bins

    Returns
    -------
    numpy.ndarray
        Reflectivity background
    numpy.ndarray
        Reflectivity background error
    """
    charge = ws.getRun().getProtonCharge()
    # For each background range, exclue the peak
    bck_ranges = find_ranges_without_overlap([bck[0], bck[1]], peak)
    bck_ranges.extend(find_ranges_without_overlap([bck[2], bck[3]], peak))
    print("Functional background: %s" % bck_ranges)

    # Iterate through the ranges and gather the background points
    bck_counts = []
    d_bck_counts = []
    pixels = []
    for r in bck_ranges:
        # This condition takes care of rejecting empty ranges,
        # including the case where the second range was left to [0, 0],
        # which has been the default before implementing this more flexible
        # approach.
        if not r[0] == r[1]:
            _b, _d_b = event_reflectivity._reflectivity(
                ws,
                peak_position=0,
                q_bins=q_bins,
                peak=r,
                low_res=low_res,
                theta=event_reflectivity.theta,
                q_summing=q_summing,
                wl_dist=wl_dist,
                wl_bins=wl_bins,
                sum_pixels=False,
            )
            bck_counts.append(_b)
            d_bck_counts.append(_d_b)
            pixels.extend(list(range(r[0], r[1] + 1)))

    # Put all those points together
    _bck = np.vstack(bck_counts)
    _d_bck = np.vstack(d_bck_counts)
    pixels = np.asarray(pixels)

    # Loop over the Q or TOF bins and fit the background
    refl_bck = np.zeros(_bck.shape[1])
    d_refl_bck = np.zeros(_bck.shape[1])

    for i in range(_bck.shape[1]):
        # Use average signal for background estimate
        _estimate = np.mean(_bck[:, i])
        linear = LinearModel()
        pars = linear.make_params(slope=0, intercept=_estimate)

        weights = 1 / _d_bck[:, i]
        # Here we have counts normalized by proton charge, so if we want to
        # assign an error of 1 on the counts, it should be 1/charge.
        weights[_bck[:, i] == 0] = charge

        fit = linear.fit(_bck[:, i], pars, method="leastsq", x=pixels, weights=weights)

        slope = fit.params["slope"].value
        intercept = fit.params["intercept"].value
        d_slope = np.sqrt(fit.covar[0][0])
        d_intercept = np.sqrt(fit.covar[1][1])

        # Compute background under the peak
        total_bck = 0
        total_err = 0
        for k in range(peak[0], peak[1] + 1):
            total_bck += intercept + k * slope
            total_err += d_intercept**2 + k**2 * d_slope**2

        _pixel_area = peak[1] - peak[0] + 1.0

        refl_bck[i] = (slope * (peak[1] + peak[0] + 1) + 2 * intercept) * _pixel_area / 2
        d_refl_bck[i] = (
            np.sqrt(d_slope**2 * (peak[1] + peak[0] + 1) ** 2 + 4 * d_intercept**2 + 4 * (peak[1] + peak[0] + 1) * fit.covar[0][1])
            * _pixel_area
            / 2
        )

        # In case we neen the background per pixel as opposed to the total sum under the peak
        if normalize_to_single_pixel:
            _pixel_area = peak[1] - peak[0] + 1.0
            refl_bck /= _pixel_area
            d_refl_bck /= _pixel_area

    return refl_bck, d_refl_bck


def side_background(
    ws, event_reflectivity, peak, bck, low_res, normalize_to_single_pixel=False, q_bins=None, wl_dist=None, wl_bins=None, q_summing=False
):
    """
    Original background substration done using two pixels defining the
    area next to the specular peak that are considered background.

    Parameters
    ----------
    ws : Mantid workspace
        Workspace containing the data
    event_reflectivity : EventReflectivity
        EventReflectivity object
    peak : list
        Range of pixels that define the peak
    bck : list
        Range of pixels that define the background
    low_res : list
        Range in the x direction on the detector
    normalize_to_single_pixel : bool
        If True, the background is normalized to the number of pixels used to integrate the signal
    q_bins : numpy.ndarray
        Array of Q bins
    wl_dist : numpy.ndarray
        Wavelength distribution for the case where we use weighted events for normatization
    wl_bins : numpy.ndarray
        Array of wavelength bins for the case where we use weighted events for normatization
    q_summing : bool
        If True, sum the counts in Q bins

    Returns
    -------
    numpy.ndarray
        Reflectivity background
    numpy.ndarray
        Reflectivity background error
    """
    q_bins = event_reflectivity.q_bins if q_bins is None else q_bins

    # Background on the left of the peak only. We allow the user to overlap the peak
    # on the right, but only use the part left of the peak.
    if bck[0] < peak[0] - 1 and bck[1] < peak[1] + 1:
        right_side = min(bck[1], peak[0] - 1)
        _left = [bck[0], right_side]
        print("Left side background: [%s, %s]" % (_left[0], _left[1]))
        refl_bck, d_refl_bck = event_reflectivity._roi_integration(
            ws, peak=_left, low_res=low_res, q_bins=q_bins, wl_dist=wl_dist, wl_bins=wl_bins, q_summing=q_summing
        )
    # Background on the right of the peak only. We allow the user to overlap the peak
    # on the left, but only use the part right of the peak.
    elif bck[0] > peak[0] - 1 and bck[1] > peak[1] + 1:
        left_side = max(bck[0], peak[1] + 1)
        _right = [left_side, bck[1]]
        print("Right side background: [%s, %s]" % (_right[0], _right[1]))
        refl_bck, d_refl_bck = event_reflectivity._roi_integration(
            ws, peak=_right, low_res=low_res, q_bins=q_bins, wl_dist=wl_dist, wl_bins=wl_bins, q_summing=q_summing
        )
    # Background on both sides
    elif bck[0] < peak[0] - 1 and bck[1] > peak[1] + 1:
        _left = [bck[0], peak[0] - 1]
        refl_bck, d_refl_bck = event_reflectivity._roi_integration(
            ws, peak=_left, low_res=low_res, q_bins=q_bins, wl_dist=wl_dist, wl_bins=wl_bins, q_summing=q_summing
        )
        _right = [peak[1] + 1, bck[1]]
        _refl_bck, _d_refl_bck = event_reflectivity._roi_integration(
            ws, peak=_right, low_res=low_res, q_bins=q_bins, wl_dist=wl_dist, wl_bins=wl_bins, q_summing=q_summing
        )
        print("Background on both sides: [%s %s] [%s %s]" % (_left[0], _left[1], _right[0], _right[1]))

        refl_bck = (refl_bck + _refl_bck) / 2.0
        d_refl_bck = np.sqrt(d_refl_bck**2 + _d_refl_bck**2) / 2.0
    else:
        print("Invalid background: [%s %s]" % (bck[0], bck[1]))
        refl_bck = np.zeros(q_bins.shape[0] - 1)
        d_refl_bck = refl_bck

    # At this point we have integrated the region of interest and obtain the average per
    # pixel, so unless that's what we want we need to multiply by the number of pixels
    # used to integrate the signal.
    if not normalize_to_single_pixel:
        _pixel_area = peak[1] - peak[0] + 1.0
        refl_bck *= _pixel_area
        d_refl_bck *= _pixel_area

    return refl_bck, d_refl_bck
