from dataclasses import dataclass
from enum import Enum
from typing import Literal

import numpy as np
from mantid import mtd
from mantid.simpleapi import CreateWorkspace, Fit, ReplaceSpecialValues, logger


class StitchingType(Enum):
    NONE = "None"
    AUTOMATIC_AVERAGE = "AutomaticAverage"

    @classmethod
    def from_value(cls, value: str):
        if value is None:
            return cls.NONE
        for item in cls:
            if item.value.casefold() == value.casefold():
                return item
        # reached if value cannot be matched
        raise ValueError(f"Invalid StitchingType value: {value}")

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

        new_values = reduced_data.r[indices]
        values = np.concatenate((values, new_values))

        new_errors = reduced_data.err[indices]
        errors = np.concatenate((errors, new_errors))
    if len(values) > 1 and np.all(np.isfinite(errors)) and not np.any(errors == 0):
        sf = 1 / np.average(values, weights=1 / errors**2).item()
    else:
        logger.warning("Insufficient or invalid error data; setting scaling factor to 1.0.")
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

            sf = self.scale_factor_for_overlap_region(fit_range_to_use, a_left, b_left, a_right, b_right)

        return sf

    def apply_sf(self, data: ReducedData) -> ReducedData:
        """Apply the auto scaling factor to a data set."""
        r = data.r * self.sf_auto
        err = data.err * self.sf_auto
        data.temp_r = r
        data.temp_err = err
        return data

    def get_fitting_overlap_range(self, min_x: float, max_x: float, nbr_points: int) -> np.ndarray:
        """Get the range of values for fitting within the overlap region."""
        step = (float(max_x) - float(min_x)) / float(nbr_points)
        fit_range = np.arange(min_x, max_x + step, step)
        return fit_range

    def calculate_axis_overlap(self, left_axis: np.ndarray, right_axis: np.ndarray) -> OverlapInfo:
        """Calculate the overlap region between two 1D axes."""
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
        """Perform a linear least-squares fit of the reflectivity data using `y = a * x + b`.

        The region of data to be fitted is selected based on `threshold_index` and `data_type`:
        * If `data_type == "right"` (default), the fit is performed on the data from
            the beginning of the arrays up to and including `threshold_index`.
        * If `data_type == "left"`, the fit is performed on the data from
            `threshold_index` to the end of the arrays, using the temporary arrays
            `temp_r` and `temp_err` stored in `data`.

        Parameters
        ----------
        data : ReducedData
            The reduced data to be fitted.
        threshold_index : int
            Index in the `q` axis that defines the boundary of the region to fit.
            For `data_type="right"`, data up to and including this index are used.
            For `data_type="left"`, data from this index to the end are used.
        data_type : {"left", "right"}, optional
            Selects which side of `threshold_index` is fitted. `"right"` fits the
            low-Q side (start to `threshold_index`); `"left"` fits the high-Q side
            (`threshold_index` to end) using `temp_r` and `temp_err`.

        Returns
        -------
        list[float]
            A two-element list `[a, b]` containing the fitted linear parameters for
            the model `y = a * x + b`, where `a` is the slope and `b` is the intercept.
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
        """Find the index of the array element nearest to the given value."""
        idx = (np.abs(array - value)).argmin().item()
        return idx

    def scale_factor_for_overlap_region(
        self, fit_range_to_use: np.ndarray, a_left: float, b_left: float, a_right: float, b_right: float
    ) -> float:
        """Calculate the scaling factor to apply for the overlap region between two fits."""
        left_mean = self.calculate_mean_over_range(fit_range_to_use, a_left, b_left)
        right_mean = self.calculate_mean_over_range(fit_range_to_use, a_right, b_right)
        if np.isclose(left_mean, 0.0):
            logger.warning("Left mean value is zero; setting scaling factor to 1.0 to avoid division by zero.")
            sf = 1.0
        else:
            sf = right_mean / left_mean
        return sf

    def calculate_mean_over_range(self, range_to_use: np.ndarray, a: float, b: float) -> float:
        """Calculate the average value of the function over the given range."""
        mean_value = float(np.mean(a * range_to_use + b))
        return mean_value
