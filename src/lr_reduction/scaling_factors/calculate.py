from dataclasses import dataclass
from typing import Literal

import numpy as np
from mantid import mtd
from mantid.simpleapi import CreateWorkspace, Fit, ReplaceSpecialValues


@dataclass
class ReducedData:
    """Class for reduced data."""

    q: np.ndarray
    r: np.ndarray
    err: np.ndarray
    # Optional temporary arrays for calculations
    temp_r: np.ndarray | None = None
    temp_err: np.ndarray | None = None


@dataclass
class OverlapInfo:
    """Class for overlap axis information."""

    min_x: float
    max_x: float
    index_min_in_left: int
    index_max_in_right: int
    no_overlap: bool


def scaling_factor_critical_edge(q_min: float, q_max: float, data: list[ReducedData]) -> float:
    """Calculate scaling factor of critical edge for normalization."""
    values = np.zeros(shape=0, dtype=float)
    errors = np.zeros(shape=0, dtype=float)

    for reduced_data in data:
        low_bound = reduced_data.q >= q_min
        high_bound = reduced_data.q <= q_max
        indices = np.argwhere(low_bound & high_bound).T[0]

        _values_i = reduced_data.r[indices]
        values = np.concatenate((values, _values_i))

        _errors_i = reduced_data.err[indices]
        errors = np.concatenate((errors, _errors_i))

    if len(values) > 1:
        sf = 1 / np.average(values, weights=1 / errors**2).item()
    else:
        sf = 1.0

    return sf


class OverlapScalingFactor:
    """
    Container for methods to calculate scaling factor for overlapping regions
    between two reduced data sets.

    Attributes
    ----------
    left_data : ReducedData
        The left reduced data set.
    right_data : ReducedData
        The right reduced data set.
    sf_auto : float | None
        An optional automatic scaling factor to apply to the left data set.
    """

    def __init__(
        self,
        left_data: ReducedData,
        right_data: ReducedData,
        sf_auto: float = 1.0,
    ):
        self.left_data = left_data
        self.right_data = right_data
        self.sf_auto = sf_auto

    def get_scaling_factor(self) -> float:
        """
        Main function to get the scaling factor for the overlapping region
        between the two data sets.
        """
        left_set = self.apply_sf(self.left_data)
        right_set = self.right_data
        left_x_axis = left_set.q
        right_x_axis = right_set.q
        overlap = self.calculate_axis_overlap(left_x_axis, right_x_axis)

        if overlap.no_overlap:
            sf = 1.0
        else:
            [a_left, b_left] = self.fit_data(left_set, overlap.index_min_in_left, data_type="left")
            [a_right, b_right] = self.fit_data(right_set, overlap.index_max_in_right, data_type="right")

            nbr_points = 10
            fit_range_to_use = self.get_fitting_overlap_range(overlap.min_x, overlap.max_x, nbr_points)

            sf = self.scale_to_apply_for_best_overlap(fit_range_to_use, a_left, b_left, a_right, b_right)

        return sf

    def apply_sf(self, data: ReducedData) -> ReducedData:
        """
        Apply the auto scaling factor to the data
        """
        r = data.r * self.sf_auto
        err = data.err * self.sf_auto
        data.temp_r = r
        data.temp_err = err
        return data

    def get_fitting_overlap_range(self, min_x, max_x, nbr_points):
        step = (float(max_x) - float(min_x)) / float(nbr_points)
        _fit_range = np.arange(min_x, max_x + step, step)
        return _fit_range

    def calculate_axis_overlap(self, left_axis: np.ndarray, right_axis: np.ndarray) -> OverlapInfo:
        """
        Calculate the overlap region of the two axis
        """
        overlap = OverlapInfo(
            min_x=-1,
            max_x=-1,
            index_min_in_left=0,
            index_max_in_right=0,
            no_overlap=True,
        )

        if left_axis[-1] <= right_axis[0]:  # no overlap
            return overlap

        overlap.no_overlap = False
        overlap.min_x = right_axis[0]
        overlap.max_x = left_axis[-1]
        overlap.index_min_in_left = self.find_nearest(left_axis, overlap.min_x)
        overlap.index_max_in_right = self.find_nearest(right_axis, overlap.max_x)

        return overlap

    def fit_data(
        self, data: ReducedData, threshold_index: int, data_type: Literal["left", "right"] = "right"
    ) -> list[float]:
        """
        will fit the data with linear fitting y=ax + b
        """
        if data_type == "left":
            assert data.temp_r is not None and data.temp_err is not None, "Temporary data arrays should not be None."
            x_axis = data.q[threshold_index:]
            y_axis = data.temp_r[threshold_index:]
            e_axis = data.temp_err[threshold_index:]
        else:
            x_axis = data.q[: threshold_index + 1]
            y_axis = data.r[: threshold_index + 1]
            e_axis = data.err[: threshold_index + 1]

        data_to_fit = CreateWorkspace(DataX=x_axis, DataY=y_axis, DataE=e_axis, Nspec=1)

        data_to_fit = ReplaceSpecialValues(
            InputWorkspace=data_to_fit, NaNValue=0, NaNError=0, InfinityValue=0, InfinityError=0
        )

        Fit(InputWorkspace=data_to_fit, Function="name=UserFunction, Formula=a+b*x, a=1, b=2", Output="res")

        res = mtd["res_Parameters"]

        b = res.cell(0, 1)
        a = res.cell(1, 1)

        return [a, b]

    def find_nearest(self, array: np.ndarray, value: float) -> int:
        idx = (np.abs(array - value)).argmin().item()
        return idx

    def scale_to_apply_for_best_overlap(
        self, fit_range_to_use: np.ndarray, a_left: float, b_left: float, a_right: float, b_right: float
    ) -> float:
        """
        This function will use the same overlap region and will determine the scaling to apply to
        the second fit to get the best match
        """
        left_mean = self.calculate_mean_over_range(fit_range_to_use, a_left, b_left)
        right_mean = self.calculate_mean_over_range(fit_range_to_use, a_right, b_right)
        _sf = right_mean / left_mean
        return _sf

    def calculate_mean_over_range(self, range_to_use: np.ndarray, a: float, b: float) -> float:
        """
        Calculate the average value of the function over the given range
        """
        sz_range = range_to_use.size
        _sum = 0
        for i in range(sz_range):
            _value = a * range_to_use[i] + b
            _sum += _value
        _mean = float(_sum) / float(sz_range)
        return _mean
