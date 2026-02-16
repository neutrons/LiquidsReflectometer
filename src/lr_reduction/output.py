"""
Write R(q) output
"""

import json

import numpy as np
from plot_publisher import plot1d

from lr_reduction.scaling_factors.calculate import (
    OverlapScalingFactor,
    ReducedData,
    StitchingConfiguration,
    StitchingType,
    scaling_factor_critical_edge,
)

from . import __version__ as VERSION


class RunCollection:
    """
    A collection of runs to assemble into a single R(Q)
    """

    def __init__(self, average_overlap=False, stitching_configuration: StitchingConfiguration | None = None):
        self.collection = []
        self.stitching_configuration = stitching_configuration if stitching_configuration else StitchingConfiguration(StitchingType.NONE)
        self.stitching_reflectivity_scale_factors = []
        self.average_overlap = average_overlap
        self.qz_all = []
        self.refl_all = []
        self.d_refl_all = []
        self.d_qz_all = []

    def add(self, q, r, dr, meta_data, dq=None):
        """
        Add a partial R(q) to the collection

        Parameters
        ----------
        q : array
            Q values
        r : array
            R values
        dr : array
            Error in R values
        meta_data : dict
            Meta data for the run
        dq : array, optional
            Q resolution
        """
        if dq is None:
            resolution = meta_data["dq_over_q"]
            dq = resolution * q
        self.collection.append(dict(q=q, r=r, dr=dr, dq=dq, info=meta_data))
        self.stitching_reflectivity_scale_factors.append(1.0)

    def calculate_scale_factors(self):
        """
        Calculate scale factors for each run in the collection
        """

        if self.stitching_configuration.type == StitchingType.NONE:
            return
        elif self.stitching_configuration.type == StitchingType.AUTOMATIC_AVERAGE:

            # Convert collection to type used in scaling factor calculation
            # Sorted by increasing q with references to original indices
            sorted_indices, sorted_collection_reduced_data = zip(
                *sorted(
                    ((i, ReducedData(run['q'], run['r'], run['dr'])) for i, run in enumerate(self.collection)),
                    key=lambda t: t[1].q[0]
                )
            )

            # Track a cumulative scale factor to properly scale subsequent runs
            cumulative_scale_factor = 1.0

            # Apply critical edge scaling if enabled
            if self.stitching_configuration.normalize_first_angle:
                ce = scaling_factor_critical_edge(self.stitching_configuration.scale_factor_qmin, self.stitching_configuration.scale_factor_qmax, sorted_collection_reduced_data)
                self.stitching_reflectivity_scale_factors[sorted_indices[0]] = ce
                cumulative_scale_factor = ce

            # Apply scaling for remaining runs
            if len(sorted_collection_reduced_data) > 1:
                for i in range(1, len(sorted_indices)):
                    overlap_sf_calculator = OverlapScalingFactor(
                        left_data=sorted_collection_reduced_data[i - 1],
                        right_data=sorted_collection_reduced_data[i]
                    )
                    sf = overlap_sf_calculator.get_scaling_factor()
                    cumulative_scale_factor *= sf
                    self.stitching_reflectivity_scale_factors[sorted_indices[i]] = cumulative_scale_factor

    def merge(self):
        """
        Merge the collection of runs
        """
        qz_all = []
        refl_all = []
        d_refl_all = []
        d_qz_all = []

        for idx, item in enumerate(self.collection):
            for i in range(len(item["q"])):
                qz_all.append(item["q"][i])
                refl_all.append(item["r"][i] * self.stitching_reflectivity_scale_factors[idx])
                d_refl_all.append(item["dr"][i] * self.stitching_reflectivity_scale_factors[idx])
                d_qz_all.append(item["dq"][i])

        qz_all = np.asarray(qz_all)
        refl_all = np.asarray(refl_all)
        d_refl_all = np.asarray(d_refl_all)
        d_qz_all = np.asarray(d_qz_all)
        idx = np.argsort(qz_all)

        self.qz_all = np.take_along_axis(qz_all, idx, axis=None)
        self.refl_all = np.take_along_axis(refl_all, idx, axis=None)
        self.d_refl_all = np.take_along_axis(d_refl_all, idx, axis=None)
        self.d_qz_all = np.take_along_axis(d_qz_all, idx, axis=None)

        if self.average_overlap:
            # New full list of points
            qz_all = []
            refl_all = []
            d_refl_all = []
            d_qz_all = []

            # Average information for groups of points
            qz = self.qz_all[0]
            total = self.refl_all[0]
            err2 = self.d_refl_all[0] ** 2
            dq = self.d_qz_all[0]
            npts = 1.0

            for i in range(1, len(self.qz_all)):
                if (self.qz_all[i] - qz) / qz > 0.000001:
                    # Store the previous point
                    qz_all.append(qz)
                    refl_all.append(total / npts)
                    d_refl_all.append(np.sqrt(err2) / npts)
                    d_qz_all.append(dq)

                    # Start a new point
                    qz = self.qz_all[i]
                    total = self.refl_all[i]
                    err2 = self.d_refl_all[i] ** 2
                    dq = self.d_qz_all[i]
                    npts = 1.0
                else:
                    total += self.refl_all[i]
                    err2 += self.d_refl_all[i] ** 2
                    npts += 1.0

            # Store the last point
            qz_all.append(qz)
            refl_all.append(total / npts)
            d_refl_all.append(np.sqrt(err2) / npts)
            d_qz_all.append(dq)

            self.qz_all = np.asarray(qz_all)
            self.refl_all = np.asarray(refl_all)
            self.d_refl_all = np.asarray(d_refl_all)
            self.d_qz_all = np.asarray(d_qz_all)

    def save_ascii(self, file_path, meta_as_json=False):
        """
        Save R(Q) in ASCII format.
        This function merges the data before saving. It writes metadata and R(Q) data
        to the specified file in ASCII format. The metadata includes experiment details,
        reduction version, run title, start time, reduction time, and other optional
        parameters. The R(Q) data includes Q, R, dR, and dQ values.

        Parameters
        ----------
        file_path : str
            The path to the file where the ASCII data will be saved.
        meta_as_json : bool, optional
            If True, metadata will be written in JSON format. Default is False.
        """
        self.calculate_scale_factors()
        self.merge()

        with open(file_path, "w") as fd:
            # Write meta data
            initial_entry_written = False
            for i, item in enumerate(self.collection):
                _meta = item["info"]
                if not initial_entry_written:
                    fd.write("# Experiment %s Run %s\n" % (_meta["experiment"], _meta["run_number"]))
                    fd.write("# Reduction %s\n" % VERSION)
                    fd.write("# Run title: %s\n" % _meta["run_title"])
                    fd.write("# Run start time: %s\n" % _meta["start_time"])
                    fd.write("# Reduction time: %s\n" % _meta["time"])
                    if "q_summing" in _meta:
                        fd.write("# Q summing: %s\n" % _meta["q_summing"])
                    if "tof_weighted" in _meta:
                        fd.write("# TOF weighted: %s\n" % _meta["tof_weighted"])
                    if "bck_in_q" in _meta:
                        fd.write("# Bck in Q: %s\n" % _meta["bck_in_q"])
                    if "theta_offset" in _meta:
                        fd.write("# Theta offset: %s\n" % _meta["theta_offset"])
                    fd.write("# Stitching type: %s\n" % self.stitching_configuration.type.value)
                    if self.stitching_configuration.type != StitchingType.NONE:
                        fd.write("# Scale factor q min: %s\n" % self.stitching_configuration.scale_factor_qmin)
                        fd.write("# Scale factor q max: %s\n" % self.stitching_configuration.scale_factor_qmax)
                    if meta_as_json:
                        fd.write("# Meta:%s\n" % json.dumps(_meta))
                    fd.write("# DataRun   NormRun   TwoTheta(deg)  LambdaMin(A)   ")
                    fd.write("LambdaMax(A) Qmin(1/A)    Qmax(1/A)    SF_A         SF_B          SF\n")
                    fd.write("")
                if "scaling_factors" in _meta:
                    a = _meta["scaling_factors"]["a"]
                    b = _meta["scaling_factors"]["b"]
                else:
                    a = 1
                    b = 0
                two_theta = _meta["theta"] * 360 / np.pi
                value_list = (
                    _meta["run_number"],
                    _meta["norm_run"],
                    two_theta,
                    _meta["wl_min"],
                    _meta["wl_max"],
                    _meta["q_min"],
                    _meta["q_max"],
                    a,
                    b,
                    self.stitching_reflectivity_scale_factors[i],
                )
                fd.write("# %-9s %-9s %-14.6g %-14.6g %-12.6g %-12.6s %-12.6s %-12.6s %-12.6s %-12.6s\n" % value_list)
                initial_entry_written = True

            # Write R(q)
            fd.write("# %-21s %-21s %-21s %-21s\n" % ("Q [1/Angstrom]", "R", "dR", "dQ [FWHM]"))
            fd.writelines(
                "%20.16f  %20.16f  %20.16f  %20.16f\n"
                % (self.qz_all[i], self.refl_all[i], self.d_refl_all[i], self.d_qz_all[i])
                for i in range(len(self.qz_all))
            )

    def add_from_file(self, file_path):
        """
        Read a partial result file and add it to the collection

        Parameters
        ----------
        file_path : str
            The path to the file to be read
        """
        _q, _r, _dr, _dq, _meta = read_file(file_path)
        self.add(_q, _r, _dr, _meta, dq=_dq)

    def plot(self):
        """
        Plot the combined reflectivity curve for a collection of runs

        Returns
        -------
        str
            HTML div containing the combined reflectivity curve plot
        """
        refl_curves = []
        run_names = []

        for i, item in enumerate(self.collection):
            refl_curves.append([item["q"], item["r"], item["dr"], item["dq"]])
            run_names.append(f"Run: {item['info']['run_number']}  SF: {self.stitching_reflectivity_scale_factors[i]:.3f}")

        # run_number parameter is only used when publish=True
        return plot1d(run_number="dummy_run", data_list=refl_curves, data_names=run_names,
                      instrument='REF_L',
                      x_title=u"Q (1/A)", x_log=True,
                      y_title="Reflectivity", y_log=True, show_dx=False, publish=False)


def read_file(file_path):
    """
    Read a data file and extract meta data

    Parameters
    ----------
    file_path : str
        The path to the file to be read
    """
    _meta = dict()
    with open(file_path, "r") as fd:
        for l in fd:
            if l.startswith("# Meta:"):
                _meta = json.loads(l[len("# Meta:") : -1])
    try:
        _q, _r, _dr, _dq = np.loadtxt(file_path).T
    except:
        print("Could not read file. It may have no points")
        _q = _r = _dr = _dq = []
    return _q, _r, _dr, _dq, _meta
