import numpy as np

from lr_reduction.event_reduction import EventReflectivity


class TOFReduction(EventReflectivity):

    def __init__(self, **kwargs):
        super(TOFReduction, self).__init__(**kwargs)

    def _reflectivity(
        self, ws, peak_position, peak, low_res, theta, wl_bins=None, sum_pixels=True, **kwargs
    ):
        """
        Unused parameters:
        - q_bins


        """
        print(peak_position, peak, low_res, theta, wl_bins, sum_pixels, kwargs)
        charge = ws.getRun().getProtonCharge()
        self.n_tof_bins = 100

        shape = len(wl_bins) - 1 if sum_pixels else ((peak[1] - peak[0] + 1), len(wl_bins) - 1)
        tof_bins = np.linspace(self.tof_range[0], self.tof_range[1], self.n_tof_bins + 1)

        refl = np.zeros(shape)
        d_refl_sq = np.zeros(shape)
        #counts = np.zeros(shape)

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

                _counts, _ = np.histogram(tofs, bins=tof_bins, weights=event_weights)
                if sum_pixels:
                    refl += _counts
                else:
                    refl[j - peak[0]] += _counts

            d_refl_sq = np.sqrt(np.fabs(refl)) / charge
            refl /= charge

        return refl, d_refl_sq

    def specular(self, q_summing=False, normalize=True):
        """
        Compute specular reflectivity.

        For constant-Q binning, it's preferred to use tof_weighted=True.

        Parameters
        ----------
        q_summing : bool
            Turns on constant-Q binning
        normalize : bool
            If True, and tof_weighted is False, normalization will be skipped

        Returns
        -------
        q_bins
            The Q bin boundaries
        refl
            The reflectivity values
        d_refl
            The uncertainties in the reflectivity values
        """
        # Scattering data
        refl, d_refl = self._reflectivity(
            self._ws_sc,
            peak_position=self.specular_pixel,
            peak=self.signal_peak,
            low_res=self.signal_low_res,
            theta=self.theta,
        )
        # Remove background
        if False and self.signal_bck is not None:
            refl_bck, d_refl_bck = self.bck_subtraction()
            refl -= refl_bck
            d_refl = np.sqrt(d_refl**2 + d_refl_bck**2)

        if normalize:
            # Normalize by monitor
            self._normalize(refl, d_refl)
            norm, d_norm = self._reflectivity(
                self._ws_db, peak_position=0, peak=self.norm_peak, low_res=self.norm_low_res, theta=self.theta, q_summing=False
            )

            # Direct beam background could be added here. The effect will be negligible.
            if self.norm_bck is not None:
                norm_bck, d_norm_bck = self.norm_bck_subtraction()
                norm -= norm_bck
                d_norm = np.sqrt(d_norm**2 + d_norm_bck**2)
            db_bins = norm > 0

            refl[db_bins] = refl[db_bins] / norm[db_bins]
            d_refl[db_bins] = np.sqrt(
                d_refl[db_bins] ** 2 / norm[db_bins] ** 2 + refl[db_bins] ** 2 * d_norm[db_bins] ** 2 / norm[db_bins] ** 4
            )

            # Hold on to normalization to be able to diagnose issues later
            self.norm = norm[db_bins]
            self.d_norm = d_norm[db_bins]

            # Clean up points where we have no direct beam
            zero_db = [not v for v in db_bins]
            refl[zero_db] = 0
            d_refl[zero_db] = 0

        # Convert to Q

        self.refl = refl
        self.d_refl = d_refl



        # Remove leading zeros
        r = np.trim_zeros(self.refl, "f")
        trim = len(self.refl) - len(r)
        self.refl = self.refl[trim:]
        self.d_refl = self.d_refl[trim:]
        self.q_bins = self.q_bins[trim:]



        # Compute Q resolution
        #self.dq_over_q = compute_resolution(self._ws_sc, theta=self.theta, q_summing=q_summing)
        self.q_summing = q_summing

        return self.q_bins, self.refl, self.d_refl
