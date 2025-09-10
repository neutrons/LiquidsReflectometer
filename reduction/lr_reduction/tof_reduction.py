"""
Implementation of reflectivity reduction for time-of-flight data using histograms of time-of-flight values.
"""

import numpy as np

from lr_reduction.event_reduction import EventReflectivity


class TOFReduction(EventReflectivity):
    """
    Traditional implementation of reflectivity reduction for time-of-flight data.
    As opposed to the event-based reduction, this implementation is based on
    histograms of time-of-flight values.
    """

    def _reflectivity(
        self,
        ws,
        peak_position,  # noqa: ARG002
        peak,
        low_res,
        sum_pixels=True,
        **__,
    ):
        """
        Overloads the _reflectivity method from EventReflectivity.
        This method histograms counts into time-of-flight bins in the selected region
        of interest of the detector.

        TODO: This method should be renamed to something more descriptive.

        Parameters
        ----------
        ws : Mantid workspace
            The workspace containing the data
        peak_position : int
            The position of the peak, unused in this implementation
        peak : list
            The range of pixels that define the peak
        low_res : list
            The range in the transverse direction on the detector
        sum_pixels : bool
            If True, the counts are summed over the pixels in the peak range

        Returns
        -------
        tuple
            A tuple containing:
            - refl : numpy.ndarray
                The reflectivity values
            - d_refl_sq : numpy.ndarray
                The uncertainties in the reflectivity values
        """
        charge = ws.getRun().getProtonCharge()
        shape = len(self.tof_bins) - 1 if sum_pixels else ((peak[1] - peak[0] + 1), len(self.tof_bins) - 1)

        refl = np.zeros(shape)
        d_refl_sq = np.zeros(shape)

        for i in range(low_res[0], int(low_res[1] + 1)):
            for j in range(peak[0], int(peak[1] + 1)):
                if self.instrument == self.INSTRUMENT_4A:
                    pixel = j * self.n_y + i
                else:
                    pixel = i * self.n_y + j
                evt_list = ws.getSpectrum(pixel)
                if evt_list.getNumberEvents() == 0:
                    continue

                tofs = evt_list.getTofs()
                if self.use_emission_time:
                    tofs = self.emission_time_correction(ws, tofs)

                event_weights = evt_list.getWeights()

                _counts, _ = np.histogram(tofs, bins=self.tof_bins, weights=event_weights)
                _err, _ = np.histogram(tofs, bins=self.tof_bins, weights=event_weights**2)

                if sum_pixels:
                    refl += _counts
                    d_refl_sq += _err
                else:
                    refl[j - peak[0]] += _counts
                    d_refl_sq[j - peak[0]] += _err

        d_refl_sq = np.sqrt(d_refl_sq) / charge
        refl /= charge

        return refl, d_refl_sq

    def specular(self, q_summing=False, normalize=True, number_of_bins=100):
        """
        Compute specular reflectivity.

        For constant-Q binning, it's preferred to use tof_weighted=True.

        Parameters
        ----------
        q_summing : bool
            Turns on constant-Q binning
        normalize : bool
            If True, and tof_weighted is False, normalization will be skipped
        number_of_bins : int
            Number of bins for the time-of-flight axis

        Returns
        -------
        tuple
            A tuple containing:
            - q_bins : numpy.ndarray
                The Q bin boundaries
            - refl : numpy.ndarray
                The reflectivity values
            - d_refl_sq : numpy.ndarray
                The uncertainties in the reflectivity values
        """

        # TODO: make this a parameter
        self.n_tof_bins = number_of_bins

        tof_bins = np.linspace(self.tof_range[0], self.tof_range[1], self.n_tof_bins + 1)
        self.tof_bins = tof_bins

        # Scattering data
        refl, d_refl = self._reflectivity(
            self._ws_sc,
            peak_position=self.specular_pixel,
            peak=self.signal_peak,
            low_res=self.signal_low_res,
            theta=self.theta,
            sum_pixels=not q_summing,
        )
        # Remove background
        if self.signal_bck is not None:
            refl_bck, d_refl_bck = self.bck_subtraction(normalize_to_single_pixel=True, q_summing=q_summing)
            refl -= refl_bck
            d_refl = np.sqrt(d_refl**2 + d_refl_bck**2)

        if normalize:
            norm, d_norm = self._reflectivity(self._ws_db, peak_position=0, peak=self.norm_peak, low_res=self.norm_low_res, sum_pixels=True)

            # Direct beam background could be added here. The effect will be negligible.
            if self.norm_bck is not None:
                norm_bck, d_norm_bck = self.norm_bck_subtraction()
                norm -= norm_bck
                d_norm = np.sqrt(d_norm**2 + d_norm_bck**2)

            refl = refl / norm
            d_refl = np.sqrt(d_refl**2 / norm**2 + refl**2 * d_norm**2 / norm**4)

        # Convert to Q
        wl_list = self.tof_bins / self.constant
        d_theta = self.gravity_correction(self._ws_sc, wl_list)

        if q_summing:
            refl, d_refl = self.constant_q(wl_list, refl, d_refl, d_theta)
            self.refl = refl
            self.d_refl = d_refl
        else:
            self.q_bins = np.flip(4.0 * np.pi / wl_list * np.sin(self.theta - d_theta))
            self.refl = np.flip(refl)
            self.d_refl = np.flip(d_refl)

        # Remove leading zeros
        r = np.trim_zeros(self.refl, "f")
        trim = len(self.refl) - len(r)
        self.refl = self.refl[trim:]
        self.d_refl = self.d_refl[trim:]
        self.q_bins = self.q_bins[trim:]

        return self.q_bins, self.refl, self.d_refl

    def constant_q(self, wl_values, signal, signal_err, d_theta):
        """
        Compute reflectivity using constant-Q binning

        Parameters
        ----------
        wl_values : numpy.ndarray
            The wavelength values
        signal : numpy.ndarray
            The normalized counts per pixel and wavelength
        signal_err : numpy.ndarray
            The uncertainties in the counts
        d_theta : numpy.ndarray
            Theta offset due to gravity correction

        Returns
        -------
        tuple
            A tuple containing:
            - refl : numpy.ndarray
                The reflectivity values
            - refl_err : numpy.ndarray
                The uncertainties in the reflectivity values
        """
        pixels = np.arange(self.signal_peak[0], self.signal_peak[1] + 1)

        _pixel_width = self.pixel_width

        x_distance = _pixel_width * (pixels - self.specular_pixel)
        delta_theta_f = np.arctan(x_distance / self.det_distance) / 2.0

        # Sign will depend on reflect up or down
        ths_value = self._ws_sc.getRun()["ths"].value[-1]
        delta_theta_f *= np.sign(ths_value)

        # Calculate Qz for each pixel
        LL, TT = np.meshgrid(wl_values, delta_theta_f)
        dTT, _ = np.meshgrid(d_theta, delta_theta_f)
        qz = 4 * np.pi / LL * np.sin(self.theta - dTT + TT)
        qz = qz.T

        # We use np.digitize to bin. The output of digitize is a bin
        # number for each entry, starting at 1. The first bin (0) is
        # for underflow entries, and the last bin is for overflow entries.
        n_q_values = len(self.q_bins)
        refl = np.zeros(n_q_values - 1)
        refl_err = np.zeros(n_q_values - 1)
        signal_n = np.zeros(n_q_values - 1)

        # Number of wavelength bins
        n_wl = qz.shape[0] - 1
        # Number of pixels
        n_pix = qz.shape[1]

        for tof in range(n_wl):
            # When digitizing, the first and last bins are out-of-bounds values
            z_inds = np.digitize(qz[tof], self.q_bins)

            # Move the indices so the valid bin numbers start at zero,
            # since this is how we are going to address the output array
            z_inds -= 1

            for ix in range(n_pix):
                if z_inds[ix] < n_q_values - 1 and z_inds[ix] >= 0 and not np.isnan(signal[ix][tof]):
                    refl[z_inds[ix]] += signal[ix][tof]
                    refl_err[z_inds[ix]] += signal_err[ix][tof] ** 2
                    signal_n[z_inds[ix]] += 1.0

        signal_n = np.where(signal_n > 0, signal_n, 1)
        refl = float(signal.shape[0]) * refl / signal_n
        refl_err = float(signal.shape[0]) * np.sqrt(refl_err) / signal_n

        # The following is for information purposes
        x0 = _pixel_width * (self.specular_pixel - self.signal_peak[0])
        x1 = _pixel_width * (self.specular_pixel - self.signal_peak[1])
        delta_theta_f0 = np.arctan(x0 / self.det_distance) / 2.0
        delta_theta_f1 = np.arctan(x1 / self.det_distance) / 2.0

        qz_max = 4.0 * np.pi / self.tof_range[1] * self.constant * np.fabs(np.sin(self.theta + delta_theta_f0))
        qz_min = 4.0 * np.pi / self.tof_range[1] * self.constant * np.fabs(np.sin(self.theta + delta_theta_f1))
        mid_point = (qz_max + qz_min) / 2.0

        print("Qz range: ", mid_point)
        self.summing_threshold = mid_point

        return refl, refl_err
