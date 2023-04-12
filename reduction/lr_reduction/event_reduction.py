"""
    Event based reduction for the Liquids Reflectometer
"""
import sys
import time

import mantid
import mantid.simpleapi as api

import numpy as np


def get_wl_range(ws):
    """
        Determine TOF range from the data
        :param workspace ws: workspace to work with
    """
    run_object = ws.getRun()

    wl = run_object.getProperty('LambdaRequest').value[0]
    chopper_speed = run_object.getProperty('SpeedRequest1').value[0]

    # Cut the edges by using a width of 2.6 A
    wl_min = (wl - 1.3 * 60.0 / chopper_speed)
    wl_max = (wl + 1.3 * 60.0 / chopper_speed)

    return [wl_min, wl_max]

def get_q_binning(q_min=0.001, q_max=0.15, q_step=-0.02):
    """
        Determine Q binning
    """
    if q_step > 0:
        n_steps = int((q_max-q_min)/q_step)
        return q_min + np.asarray([q_step * i for i in range(n_steps)])
    else:
        _step = 1.0+np.abs(q_step)
        n_steps = int(np.log(q_max/q_min)/np.log(_step))
        return q_min * np.asarray([_step**i for i in range(n_steps)])

def quicknxs_scale(theta, peak, low_res, norm_peak, norm_low_res):
    """
        Scaling factor to multiply by to be compatible with QuickNXS 1.0.
    """
    quicknxs_scale = (float(norm_peak[1])-float(norm_peak[0])) * (float(norm_low_res[1])-float(norm_low_res[0]))
    quicknxs_scale /= (float(peak[1])-float(peak[0])) * (float(low_res[1])-float(low_res[0]))
    _scale = 0.005 / np.fabs(np.sin(theta)) if theta > 0.0002 else 1.0
    quicknxs_scale *= _scale
    return quicknxs_scale


class EventReflectivity(object):
    """
        Event based reflectivity calculation.
        List of items to be taken care of outside this class:
          - Edge points cropping
          - Angle offset
          - Putting runs together in one R(q) curve
          - Scaling factors
    """
    QX_VS_QZ = 0
    KZI_VS_KZF = 1
    DELTA_KZ_VS_QZ = 3
    THETAF_VS_WL = 4
    INSTRUMENT_4A = 0
    INSTRUMENT_4B = 1
    DEFAULT_4B_SAMPLE_DET_DISTANCE = 1.83
    DEFAULT_4B_SOURCE_DET_DISTANCE = 15.75

    def __init__(self, scattering_workspace, direct_workspace,
                 signal_peak, signal_bck, norm_peak, norm_bck,
                 specular_pixel, signal_low_res, norm_low_res,
                 q_min=None, q_step=-0.02, q_max=None,
                 tof_range=None, theta=1.0, instrument=None):
        """
            Pixel ranges include the min and max pixels.

            :param scattering_workspace: Mantid workspace containing the reflected data
            :param direct_workspace: Mantid workspace containing the direct beam data [if None, normalization won't be applied]
            :param signal_peak: pixel min and max for the specular peak
            :param signal_bck: pixel range of the background [if None, the background won't be subtracted]
            :param norm_peak: pixel range of the direct beam peak
            :param norm_bck: direct background subtraction is not used [deprecated]
            :param specular_pixel: pixel of the specular peak
            :param signal_low_res: pixel range of the specular peak out of the scattering plane
            :param norm_low_res: pixel range of the direct beam out of the scattering plane
            :param q_min: value of lowest q point
            :param q_step: step size in Q. Enter a negative value to get a log scale
            :param q_min: value of largest q point
            :param tof_range: TOF range,or None
            :param theta: theta scattering angle in radians
        """
        if instrument in [self.INSTRUMENT_4A, self.INSTRUMENT_4B]:
            self.instrument = instrument
        else:
            self.instrument = self.INSTRUMENT_4B
        self.signal_peak = signal_peak
        self.signal_bck = signal_bck
        self.norm_peak = norm_peak
        self.norm_bck = norm_bck
        self.signal_low_res = signal_low_res
        self.norm_low_res = norm_low_res
        self.specular_pixel = specular_pixel
        self.q_min = q_min
        self.q_max = q_max
        self.q_step = q_step
        self.tof_range = tof_range
        self.theta = theta
        self._offspec_x_bins = None
        self._offspec_z_bins = None
        self.summing_threshold = None

        # Process workspaces
        if self.tof_range is not None:
            self._ws_sc = api.CropWorkspace(InputWorkspace=scattering_workspace,
                                            XMin=tof_range[0], XMax=tof_range[1],
                                            OutputWorkspace='_'+str(scattering_workspace))
            if direct_workspace is not None:
                self._ws_db = api.CropWorkspace(InputWorkspace=direct_workspace,
                                                XMin=tof_range[0], XMax=tof_range[1],
                                                OutputWorkspace='_'+str(direct_workspace))
            else:
                self._ws_db = None
        else:
            self._ws_sc = scattering_workspace
            self._ws_db = direct_workspace

        # Extract meta data
        self.extract_meta_data()

    def extract_meta_data(self):
        """
            Extract meta data from the data file.
        """
        # Set up basic data
        self.n_x = int(self._ws_sc.getInstrument().getNumberParameter("number-of-x-pixels")[0])
        self.n_y = int(self._ws_sc.getInstrument().getNumberParameter("number-of-y-pixels")[0])

        self.pixel_width = float(self._ws_sc.getInstrument().getNumberParameter("pixel-width")[0]) / 1000.0

        if self.instrument == self.INSTRUMENT_4B:
            self.extract_meta_data_4B()
        else:
            self.extract_meta_data_4A()

        h = 6.626e-34  # m^2 kg s^-1
        m = 1.675e-27  # kg
        self.constant = 1e-4 * m * self.source_detector_distance / h

        if self.tof_range is None:
            self.wl_range = get_wl_range(self._ws_sc)
        else:
            self.wl_range = [self.tof_range[0] / self.constant, self.tof_range[1] /  self.constant]

        if self.q_min is None:
            self.q_min = 4.0*np.pi/self.wl_range[1] * np.fabs(np.sin(self.theta))
        if self.q_max is None:
            self.q_max = 4.0*np.pi/self.wl_range[0] * np.fabs(np.sin(self.theta))

        # Q binning to use
        self.q_bins = get_q_binning(self.q_min, self.q_max, self.q_step)

        # Catch options that can be turned off
        if self.signal_low_res is None:
            self.signal_low_res = [1, self.n_x-1]
        if self.norm_low_res is None:
            self.norm_low_res = [1, self.n_x-1]

    def extract_meta_data_4A(self):
        """
            4A-specific meta data
        """
        run_object = self._ws_sc.getRun()
        self.det_distance = run_object['SampleDetDis'].getStatistics().mean
        source_sample_distance = run_object['ModeratorSamDis'].getStatistics().mean
        if not run_object['SampleDetDis'].units in ['m', 'meter']:
            self.det_distance /= 1000.0
        if not run_object['ModeratorSamDis'].units in ['m', 'meter']:
            source_sample_distance /= 1000.0
        self.source_detector_distance = source_sample_distance + self.det_distance

    def extract_meta_data_4B(self):
        """
            4B-specific meta data

            Distance from source to sample was 13.63 meters prior to the source
            to detector distance being determined with Bragg edges to be 15.75 m. 
        """
        if self._ws_sc.getInstrument().hasParameter("sample-det-distance"):
            self.det_distance = self._ws_sc.getInstrument().getNumberParameter("sample-det-distance")[0]
        else:
            self.det_distance = self.DEFAULT_4B_SAMPLE_DET_DISTANCE

        if self._ws_sc.getInstrument().hasParameter("source-det-distance"):
            self.source_detector_distance = self._ws_sc.getInstrument().getNumberParameter("source-det-distance")[0]
        else:
            self.source_detector_distance = self.DEFAULT_4B_SOURCE_DET_DISTANCE

    def __repr__(self):
        output = "sample-det: %s\n" % self.det_distance
        output += "pixel: %s\n" % self.pixel_width
        output += "WL: %s %s\n" % (self.wl_range[0], self.wl_range[1])
        output += "Q: %s %s\n" % (self.q_min, self.q_max)
        output += "Theta = %s" % self.theta
        return output

    def to_dict(self):
        """
            Returns meta-data to be used/stored.
        """
        if self._ws_sc.getRun().hasProperty("start_time"):
            start_time = self._ws_sc.getRun().getProperty("start_time").value
        else:
            start_time = 'live'
        experiment = self._ws_sc.getRun().getProperty("experiment_identifier").value
        run_number = self._ws_sc.getRun().getProperty("run_number").value
        sequence_number = int(self._ws_sc.getRun().getProperty("sequence_number").value[0])
        sequence_id = int(self._ws_sc.getRun().getProperty("sequence_id").value[0])
        run_title = self._ws_sc.getTitle()

        if self._ws_db:
            norm_run = self._ws_db.getRunNumber()
        else:
            norm_run = 0

        dq0 = 0
        dq_over_q = compute_resolution(self._ws_sc)
        return dict(wl_min=self.wl_range[0], wl_max=self.wl_range[1],
                    q_min=self.q_min, q_max=self.q_max, theta=self.theta,
                    start_time=start_time, experiment=experiment, run_number=run_number,
                    run_title=run_title, norm_run=norm_run, time=time.ctime(),
                    dq0=dq0, dq_over_q=dq_over_q, sequence_number=sequence_number,
                    sequence_id=sequence_id)

    def specular(self, q_summing=False, tof_weighted=False, bck_in_q=False,
                 clean=False, normalize=True):
        """
            Compute specular reflectivity.

            For constant-Q binning, it's preferred to use tof_weighted=True.
            
            :param q_summing: turns on constant-Q binning
            :param tof_weighted: if True, binning will be done by weighting each event to the DB distribution
            :param bck_in_q: if True, the background will be estimated in Q space using the constant-Q binning approach
            :param clean: if True, and Q summing is True, then leading artifact will be removed
            :param normalize: if True, and tof_weighted is False, normalization will be skipped
        """
        if tof_weighted:
            self.specular_weighted(q_summing=q_summing, bck_in_q=bck_in_q)
        else:
            self.specular_unweighted(q_summing=q_summing, normalize=normalize)

        # Remove leading zeros
        r = np.trim_zeros(self.refl, 'f')
        trim = len(self.refl) - len(r)
        self.refl = self.refl[trim:]
        self.d_refl = self.d_refl[trim:]
        self.q_bins = self.q_bins[trim:]

        # Remove leading artifact from the wavelength coverage
        # Remember that q_bins is longer than refl by 1 because
        # it contains bin boundaries
        if clean and self.summing_threshold:
            print("Summing threshold: %g" % self.summing_threshold)
            idx = self.q_bins > self.summing_threshold
            self.refl = self.refl[idx[:-1]]
            self.d_refl = self.d_refl[idx[:-1]]
            self.q_bins = self.q_bins[idx]

        return self.q_bins, self.refl, self.d_refl

    def specular_unweighted(self, q_summing=False, normalize=True):
        """
            Simple specular reflectivity calculation. This is the same approach as the
            original LR reduction, which sums up pixels without constant-Q binning.
            The original approach bins in TOF, then rebins the final results after
            transformation to Q. This approach bins directly to Q.
        """
        # Scattering data
        refl, d_refl = self._reflectivity(self._ws_sc, peak_position=self.specular_pixel,
                                          peak=self.signal_peak, low_res=self.signal_low_res,
                                          theta=self.theta, q_summing=q_summing)

        # Remove background
        if self.signal_bck is not None:
            refl_bck, d_refl_bck = self.bck_subtraction()
            refl -= refl_bck
            d_refl = np.sqrt(d_refl**2 + d_refl_bck**2)

        if normalize:
            # Since this approach does not involve constant-Q binning, and since the transform
            # from TOF to Q is the same across the TOF range for a given run (constant angle),
            # we can bin the DB according to the same transform instead of binning and dividing in TOF.
            # This is mathematically equivalent and convenient in terms of abstraction for later
            # use for the constant-Q calculation elsewhere in the code.
            norm, d_norm = self._reflectivity(self._ws_db, peak_position=0,
                                              peak=self.norm_peak, low_res=self.norm_low_res,
                                              theta=self.theta, q_summing=False)

            # Direct beam background could be added here. The effect will be negligible.
            if self.norm_bck is not None:
                norm_bck, d_norm_bck = self.norm_bck_subtraction()
                norm -= norm_bck
                d_norm = np.sqrt(d_norm**2 + d_norm_bck**2)
            db_bins = norm>0

            refl[db_bins] = refl[db_bins]/norm[db_bins]
            d_refl[db_bins] = np.sqrt(d_refl[db_bins]**2 / norm[db_bins]**2 + refl[db_bins]**2 * d_norm[db_bins]**2 / norm[db_bins]**4)

            # Clean up points where we have no direct beam
            zero_db = [not v for v in db_bins]
            refl[zero_db] = 0
            d_refl[zero_db] = 0

        self.refl = refl
        self.d_refl = d_refl
        return self.q_bins, refl, d_refl

    def specular_weighted(self, q_summing=True, bck_in_q=False):
        """
            Compute reflectivity by weighting each event by flux.
            This allows for summing in Q and to estimate the background in either Q
            or pixels next to the peak.
        """
        # Event weights for normalization
        db_charge = self._ws_db.getRun().getProtonCharge()
        wl_events = self._get_events(self._ws_db, self.norm_peak, self.norm_low_res)
        wl_dist, wl_bins = np.histogram(wl_events, bins=60)
        wl_dist = wl_dist/db_charge/(wl_bins[1]-wl_bins[0])
        wl_middle = [(wl_bins[i+1]+wl_bins[i])/2.0 for i in range(len(wl_bins)-1)]

        refl, d_refl = self._reflectivity(self._ws_sc, peak_position=self.specular_pixel,
                                          peak=self.signal_peak, low_res=self.signal_low_res,
                                          theta=self.theta, q_summing=q_summing, wl_dist=wl_dist, wl_bins=wl_middle)

        if self.signal_bck is not None:
            refl_bck, d_refl_bck = self.bck_subtraction(wl_dist=wl_dist, wl_bins=wl_middle,
                                                        q_summing=bck_in_q)
            refl -= refl_bck
            d_refl = np.sqrt(d_refl**2 + d_refl_bck**2)

        self.refl = refl
        self.d_refl = d_refl
        return self.q_bins, refl, d_refl

    def _roi_integration(self, ws, peak, low_res, q_bins=None, wl_dist=None, wl_bins=None, q_summing=False):
        """
            Integrate a region of interest and normalize by the number of included pixels.

            The options are the same as for the reflectivity calculation.
            If wl_dist and wl_bins are supplied, the events will be weighted by flux.
            If q_summing is True, the angle of each neutron will be recalculated according to
            their position on the detector and place in the proper Q bin.
        """
        q_bins = self.q_bins if q_bins is None else q_bins
        refl_bck, d_refl_bck = self._reflectivity(ws, peak_position=0, q_bins=q_bins,
                                                  peak=peak, low_res=low_res,
                                                  theta=self.theta, q_summing=q_summing,
                                                  wl_dist=wl_dist, wl_bins=wl_bins)

        _pixel_area = (peak[1]-peak[0]+1.0)
        refl_bck /= _pixel_area
        d_refl_bck /= _pixel_area
        return refl_bck, d_refl_bck

    def _bck_subtraction(self, ws, peak, bck, low_res, normalize_to_single_pixel=False,
                         q_bins=None, wl_dist=None, wl_bins=None, q_summing=False):
        """
            Abstracted out background subtraction process.

            The options are the same as for the reflectivity calculation.
            If wl_dist and wl_bins are supplied, the events will be weighted by flux.
            If q_summing is True, the angle of each neutron will be recalculated according to
            their position on the detector and place in the proper Q bin.
        """
        q_bins = self.q_bins if q_bins is None else q_bins

        # Background on the left of the peak only. We allow the user to overlap the peak on the right,
        # but only use the part left of the peak.
        if bck[0] < peak[0]-1 and bck[1] < peak[1]+1:
            right_side = min(bck[1], peak[0]-1)
            _left = [bck[0], right_side]
            print("Left side background: [%s, %s]" % (_left[0], _left[1]))
            refl_bck, d_refl_bck = self._roi_integration(ws, peak=_left, low_res=low_res,
                                                         q_bins=q_bins, wl_dist=wl_dist,
                                                         wl_bins=wl_bins, q_summing=q_summing)
        # Background on the right of the peak only. We allow the user to overlap the peak on the left,
        # but only use the part right of the peak.
        elif bck[0] > peak[0]-1 and bck[1] > peak[1]+1:
            left_side = max(bck[0], peak[1]+1)
            _right = [left_side, bck[1]]
            print("Right side background: [%s, %s]" % (_right[0], _right[1]))
            refl_bck, d_refl_bck = self._roi_integration(ws, peak=_right, low_res=low_res,
                                                         q_bins=q_bins, wl_dist=wl_dist,
                                                         wl_bins=wl_bins, q_summing=q_summing)
        # Background on both sides
        elif bck[0] < peak[0]-1 and bck[1] > peak[1]+1:
            _left = [bck[0], peak[0]-1]
            refl_bck, d_refl_bck = self._roi_integration(ws, peak=_left, low_res=low_res,
                                                         q_bins=q_bins, wl_dist=wl_dist,
                                                         wl_bins=wl_bins, q_summing=q_summing)
            _right = [peak[1]+1, bck[1]]
            _refl_bck, _d_refl_bck = self._roi_integration(ws, peak=_right, low_res=low_res,
                                                           q_bins=q_bins, wl_dist=wl_dist,
                                                           wl_bins=wl_bins, q_summing=q_summing)
            print("Background on both sides: [%s %s] [%s %s]" % (_left[0], _left[1], _right[0], _right[1]))

            refl_bck = (refl_bck + _refl_bck)/2.0
            d_refl_bck = np.sqrt(d_refl_bck**2 + _d_refl_bck**2)/2.0
        else:
            print("Invalid background: [%s %s]" % (bck[0], bck[1]))
            refl_bck = np.zeros(q_bins.shape[0]-1)
            d_refl_bck = refl_bck

        # At this point we have integrated the region of interest and obtain the average per
        # pixel, so unless that's what we want we need to multiply by the number of pixels
        # used to integrate the signal.
        if not normalize_to_single_pixel:
            _pixel_area = peak[1] - peak[0]+1.0
            refl_bck *= _pixel_area
            d_refl_bck *= _pixel_area

        return refl_bck, d_refl_bck

    def bck_subtraction(self, normalize_to_single_pixel=False, q_bins=None, wl_dist=None, wl_bins=None,
                        q_summing=False):
        """
            Higher-level call for background subtraction. Hides the ranges needed to define the ROI.
        """
        return self._bck_subtraction(self._ws_sc, self.signal_peak, self.signal_bck, self.signal_low_res,
                                     normalize_to_single_pixel=normalize_to_single_pixel, q_bins=q_bins,
                                     wl_dist=wl_dist, wl_bins=wl_bins, q_summing=q_summing)

    def norm_bck_subtraction(self):
        """
            Higher-level call for background subtraction for the normalization run.
        """
        return self._bck_subtraction(self._ws_db, self.norm_peak, self.norm_bck, self.norm_low_res,
                                     normalize_to_single_pixel=False)

    def slice(self, x_min=0.002, x_max=0.004, x_bins=None, z_bins=None,
              refl=None, d_refl=None, normalize=False):
        """
            Retrieve a slice from the off-specular data.
        """
        x_bins = self._offspec_x_bins if x_bins is None else x_bins
        z_bins = self._offspec_z_bins if z_bins is None else z_bins
        refl = self._offspec_refl if refl is None else refl
        d_refl = self._offspec_d_refl if d_refl is None else d_refl

        i_min = len(x_bins[x_bins<x_min])
        i_max = len(x_bins[x_bins<x_max])

        _spec = np.sum(refl[i_min:i_max], axis=0)
        _d_spec = np.sum( (d_refl[i_min:i_max])**2, axis=0)
        _d_spec = np.sqrt(_d_spec)
        if normalize:
            _spec /= (i_max-i_min)
            _d_spec /= (i_max-i_min)

        return z_bins, _spec, _d_spec

    def _reflectivity(self, ws, peak_position, peak, low_res, theta, q_bins=None, q_summing=False, wl_dist=None, wl_bins=None):
        """
            Assumes that the input workspace is normalized by proton charge.
        """
        charge = ws.getRun().getProtonCharge()
        _q_bins = self.q_bins if q_bins is None else q_bins

        refl = np.zeros(len(_q_bins)-1)
        d_refl_sq = np.zeros(len(_q_bins)-1)
        counts = np.zeros(len(_q_bins)-1)
        _pixel_width = self.pixel_width if q_summing else 0.0

        for i in range(low_res[0], int(low_res[1]+1)):
            for j in range(peak[0], int(peak[1]+1)):
                if self.instrument == self.INSTRUMENT_4A:
                    pixel = j * self.n_y + i
                else:
                    pixel = i * self.n_y + j
                evt_list = ws.getSpectrum(pixel)
                if evt_list.getNumberEvents() == 0:
                    continue

                wl_list = evt_list.getTofs() / self.constant

                # Gravity correction
                d_theta = self.gravity_correction(ws, wl_list)

                # Sign will depend on reflect up or down
                x_distance = _pixel_width * (j - peak_position)
                delta_theta_f = np.arctan(x_distance / self.det_distance) / 2.0
                qz=4.0*np.pi/wl_list * np.sin(theta + delta_theta_f - d_theta)
                qz = np.fabs(qz)

                if wl_dist is not None and wl_bins is not None:
                    wl_weights = 1.0/np.interp(wl_list, wl_bins, wl_dist, np.inf, np.inf)
                    hist_weigths = wl_weights * qz / wl_list
                    _counts, _ = np.histogram(qz, bins=_q_bins, weights=hist_weigths)
                    refl += _counts
                    _counts, _ = np.histogram(qz, bins=_q_bins)
                    counts += _counts
                else:
                    _counts, _ = np.histogram(qz, bins=_q_bins)
                    refl += _counts

        # The following is for information purposes
        if q_summing:
            x0 = _pixel_width * (peak_position - peak[0])
            x1 = _pixel_width * (peak_position - peak[1])
            delta_theta_f0 = np.arctan(x0 / self.det_distance) / 2.0
            delta_theta_f1 = np.arctan(x1 / self.det_distance) / 2.0

            qz_max = 4.0*np.pi/self.tof_range[1]*self.constant * np.fabs(np.sin(theta + delta_theta_f0))
            qz_min = 4.0*np.pi/self.tof_range[1]*self.constant * np.fabs(np.sin(theta + delta_theta_f1))
            mid_point = (qz_max + qz_min)/2.0
            print("Qz range: ", qz_min, mid_point, qz_max)
            self.summing_threshold = mid_point

        if wl_dist is not None and wl_bins is not None:
            bin_size = _q_bins[1:] - _q_bins[:-1]
            non_zero = counts > 0
            d_refl_sq[non_zero] =  refl[non_zero] / np.sqrt(counts[non_zero]) / charge / bin_size[non_zero]
            refl[non_zero] = refl[non_zero] / charge / bin_size[non_zero]
        else:
            d_refl_sq = np.sqrt(np.fabs(refl)) / charge
            refl /= charge 

        return refl, d_refl_sq

    def _get_events(self, ws, peak, low_res):
        """
            Return an array of wavelengths for a given workspace.
        """
        wl_events = np.asarray([])

        for i in range(low_res[0], int(low_res[1]+1)):
            for j in range(peak[0], int(peak[1]+1)):
                if self.instrument == self.INSTRUMENT_4A:
                    pixel = j * self.n_y + i
                else:
                    pixel = i * self.n_y + j
                evt_list = ws.getSpectrum(pixel)
                wl_list = evt_list.getTofs() / self.constant
                wl_events = np.concatenate((wl_events, wl_list))

        return wl_events

    def off_specular(self, x_axis=None, x_min=-0.015, x_max=0.015, x_npts=50,
                     z_min=None, z_max=None, z_npts=-120, bck_in_q=None):
        """
            Compute off-specular
            :param x_axis: Axis selection
            :param x_min: Min value on x-axis
            :param x_max: Max value on x-axis
            :param x_npts: Number of points in x (negative will produce a log scale)
            :param z_min: Min value on z-axis (if none, default Qz will be used)
            :param z_max: Max value on z-axis (if none, default Qz will be used)
            :param z_npts: Number of points in z (negative will produce a log scale)
        """
        # Z axis binning
        qz_bins = self.q_bins
        if z_min is not None and z_max is not None:
            if z_npts < 0:
                qz_bins = np.logspace(np.log10(z_min), np.log10(z_max), num=np.abs(z_npts))
            else:
                qz_bins = np.linspace(z_min, z_max, num=z_npts)

        # X axis binning
        if x_npts > 0:
            qx_bins = np.linspace(x_min, x_max, num=x_npts)
        else:
            qx_bins = np.logspace(np.log10(x_min), np.log10(x_max), num=np.abs(x_npts))

        wl_events = self._get_events(self._ws_db, self.norm_peak, self.norm_low_res)
        wl_dist, wl_bins = np.histogram(wl_events, bins=60)
        wl_middle = [(wl_bins[i+1]+wl_bins[i])/2.0 for i in range(len(wl_bins)-1)]

        _refl, _d_refl = self._off_specular(self._ws_sc, wl_dist, wl_middle, qx_bins, qz_bins,
                                            self.specular_pixel, self.theta, x_axis=x_axis)
        db_charge = self._ws_db.getRun().getProtonCharge()
        _refl *= db_charge * (wl_bins[1]-wl_bins[0])
        _d_refl *= db_charge * (wl_bins[1]-wl_bins[0])

        # Background
        if self.signal_bck:
            if bck_in_q is None:
                print("Not implemented")
            else:
                _, refl_bck, d_refl_bck = self.slice(bck_in_q[0], bck_in_q[1],
                                                     x_bins=qx_bins, z_bins=qz_bins,
                                                     refl=_refl, d_refl=_d_refl,
                                                     normalize=True)
                _refl -= refl_bck
                _d_refl = np.sqrt(_d_refl**2 + d_refl_bck**2)

        self._offspec_x_bins = qx_bins
        self._offspec_z_bins = qz_bins
        self._offspec_refl = _refl
        self._offspec_d_refl = _d_refl

        return qx_bins, qz_bins, _refl, _d_refl

    def _off_specular(self, ws, wl_dist, wl_bins, x_bins, z_bins, peak_position, theta, x_axis=None):
        charge = ws.getRun().getProtonCharge()
        refl = np.zeros([len(x_bins)-1, len(z_bins)-1])
        counts = np.zeros([len(x_bins)-1, len(z_bins)-1])

        for j in range(0, self.n_x):
            wl_list = np.asarray([])
            for i in range(self.signal_low_res[0], int(self.signal_low_res[1]+1)):
                if self.instrument == self.INSTRUMENT_4A:
                    pixel = j * self.n_y + i
                else:
                    pixel = i * self.n_y + j
                evt_list = ws.getSpectrum(pixel)
                wl_events = evt_list.getTofs() / self.constant
                wl_list = np.concatenate((wl_events, wl_list))

            k = 2.0 * np.pi / wl_list
            wl_weights = 1.0/np.interp(wl_list, wl_bins, wl_dist, np.inf, np.inf)

            #TODO: Sign with depend on reflect up or down
            x_distance = float(j-peak_position) * self.pixel_width
            delta_theta_f = np.arctan(x_distance / self.det_distance)
            theta_f = theta + delta_theta_f

            qz = k * (np.sin(theta_f) + np.sin(theta))
            qz = np.fabs(qz)
            qx = k * (np.cos(theta_f) - np.cos(theta))
            ki_z = k * np.sin(theta)
            kf_z = k * np.sin(theta_f)

            _x = qx
            _z = qz
            if x_axis == EventReflectivity.DELTA_KZ_VS_QZ:
                _x = (ki_z - kf_z)
            elif x_axis == EventReflectivity.KZI_VS_KZF:
                _x = ki_z
                _z = kf_z
            elif x_axis == EventReflectivity.THETAF_VS_WL:
                _x = wl_list
                _z = np.ones(len(wl_list))*theta_f

            histo_weigths = wl_weights * _z / wl_list
            _counts, _, _ = np.histogram2d(_x, _z, bins=[x_bins, z_bins], weights=histo_weigths)
            refl += _counts
            _counts, _, _ = np.histogram2d(_x, _z, bins=[x_bins, z_bins])
            counts += _counts

        bin_size = z_bins[1] - z_bins[0]
        d_refl_sq =  refl / np.sqrt(counts) / charge / bin_size
        refl /= charge * bin_size

        return refl, d_refl_sq

    def gravity_correction(self, ws, wl_list):
        """
            Gravity correction for each event
        """
        # Xi reference would be the position of xi if the si slit were to be positioned
        # at the sample. The distance from the sample to si is then xi_reference - xi.
        xi_reference = 445
        if ws.getInstrument().hasParameter("xi-reference"):
            xi_reference = ws.getInstrument().getNumberParameter("xi-reference")[0]

        # Distance between the s1 and the sample
        s1_sample_distance = 1485
        if ws.getInstrument().hasParameter("s1-sample-distance"):
            s1_sample_distance = ws.getInstrument().getNumberParameter("s1-sample-distance")[0]*1000

        xi = 310
        if ws.getInstrument().hasParameter("BL4B:Mot:xi.RBV"):
            xi = abs(ws.getRun().getProperty("BL4B:Mot:xi.RBV").value[0])

        sample_si_distance = xi_reference - xi
        slit_distance = s1_sample_distance - sample_si_distance

        # Angle of the incident beam on a horizontal sample
        theta_in=-4.0

        # Calculation from the ILL paper. This works for inclined beams.
        # Calculated theta is the angle on the sample

        g = 9.8067                        # m/s^2
        h = 6.6260715e-34                 # Js=kg m^2/s
        mn = 1.67492749804e-27            # kg

        v = h/(mn*wl_list*1e-10)
        k = g/(2*v**2)

        # Define the sample position as x=0, y=0. increasing x is towards moderator
        xs=0

        # positions of slits
        x1 = sample_si_distance / 1000
        x2 = (sample_si_distance + slit_distance) / 1000

        #height of slits determined by incident theta, y=0 is the sample height
        y1=x1*np.tan(theta_in*np.pi/180)
        y2=x2*np.tan(theta_in*np.pi/180)

        # This is the location of the top of the parabola
        x0=(y1-y2+k*(x1**2-x2**2))/(2*k*(x1-x2))

        # Angle is arctan(dy/dx) at sample
        theta_sample = np.arctan(2*k*(x0-xs)) * 180/np.pi

        return (theta_sample-theta_in) * np.pi / 180.0


def compute_resolution(ws, default_dq=0.027):
    """
        Compute the Q resolution from the meta data.
    """
    # We can't compute the resolution if the value of xi is not in the logs.
    # Since it was not always logged, check for it here.
    if not ws.getRun().hasProperty("BL4B:Mot:xi.RBV"):
        # For old data, the resolution can't be computed, so use
        # the standard value for the instrument
        print("Could not find BL4B:Mot:xi.RBV: using supplied dQ/Q")
        return default_dq

    # Xi reference would be the position of xi if the si slit were to be positioned
    # at the sample. The distance from the sample to si is then xi_reference - xi.
    xi_reference = 445
    if ws.getInstrument().hasParameter("xi-reference"):
        xi_reference = ws.getInstrument().getNumberParameter("xi-reference")[0]

    # Distance between the s1 and the sample
    s1_sample_distance = 1485
    if ws.getInstrument().hasParameter("s1-sample-distance"):
        s1_sample_distance = ws.getInstrument().getNumberParameter("s1-sample-distance")[0]*1000

    s1h = abs(ws.getRun().getProperty("S1VHeight").value[0])
    ths = abs(ws.getRun().getProperty("ths").value[0])
    xi = abs(ws.getRun().getProperty("BL4B:Mot:xi.RBV").value[0])
    sample_si_distance = xi_reference - xi
    slit_distance = s1_sample_distance - sample_si_distance
    dq_over_q = s1h / slit_distance * 180 / 3.1416 / ths
    return dq_over_q
