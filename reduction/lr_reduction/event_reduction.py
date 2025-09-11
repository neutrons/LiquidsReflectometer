"""
Event based reduction for the Liquids Reflectometer
"""

import datetime
import json
import os
import time
from typing import Optional

import mantid.simpleapi as api
import numpy as np
from scipy.optimize import brentq

from lr_reduction.instrument_settings import InstrumentSettings
from lr_reduction.utils import mantid_algorithm_exec

from . import background, dead_time_correction

PLANCK_CONSTANT = 6.626e-34  # m^2 kg s^-1
NEUTRON_MASS = 1.675e-27  # kg

# Attenuators: PV name and thickness in cm
CD_ATTENUATORS = [
    ["BL4B:Actuator:50MRb", 0.0058],
    ["BL4B:Actuator:100MRb", 0.0122],
    ["BL4B:Actuator:200MRb", 0.0244],  # Uncalibrated
    ["BL4B:Actuator:400MRb", 0.0488],  # Uncalibrated
]


def get_wl_range(ws):
    """
    Determine TOF range from the data

    Parameters
    ----------
    ws
        Mantid workspace to work with

    Returns
    -------
      list
        [min, max] wavelength range
    """
    run_object = ws.getRun()

    wl = run_object.getProperty("LambdaRequest").value[0]
    chopper_speed = run_object.getProperty("SpeedRequest1").value[0]

    # Cut the edges by using a width of 2.6 A
    wl_min = wl - 1.3 * 60.0 / chopper_speed
    wl_max = wl + 1.3 * 60.0 / chopper_speed

    return [wl_min, wl_max]


def get_q_binning(q_min=0.001, q_max=0.15, q_step=-0.02):
    """
    Determine Q binning.

    This function calculates the binning for Q values based on the provided minimum, maximum, and step values.
    If the step value is positive, it generates a linear binning. If the step value is negative, it generates
    a logarithmic binning.
    Corrected with +1 on number of steps.

    Parameters
    ----------
    q_min : float
        The minimum Q value.
    q_max : float
        The maximum Q value.
    q_step : float
        The step size for Q binning. If positive, linear binning is used.
        If negative, logarithmic binning is used.

    Returns
    -------
    numpy.ndarray
        A numpy array of Q values based on the specified binning.
    """
    if q_step > 0:
        n_steps = int((q_max - q_min) / q_step) + 1  # the +1 is to include q_max
        return q_min + np.asarray([q_step * i for i in range(n_steps)])
    else:
        _step = 1.0 + np.abs(q_step)
        n_steps = int(np.log(q_max / q_min) / np.log(_step)) + 1  # the +1 is to include q_max
        return q_min * np.asarray([_step**i for i in range(n_steps)])


def get_attenuation_info(ws):
    """
    Retrieve information about attenuation from a Mantid workspace.
    This function calculates the total thickness of all attenuators that are
    in the path of the beam by summing up the thicknesses of the attenuators
    specified in the global variable `CD_ATTENUATORS`.

    Parameters
    ----------
    ws
        Mantid workspace from which to retrieve the attenuation information.

    Returns
    -------
    float
        The total thickness of the attenuators in the path of the beam.
    """
    run_info = ws.getRun()
    attenuator_thickness = 0

    # Sum up all the attenuators that are in the path of the beam
    for att in CD_ATTENUATORS:
        if att[0] in run_info and run_info[att[0]].value[0] == 0:
            attenuator_thickness += att[1]

    return attenuator_thickness


def read_settings(ws) -> InstrumentSettings:
    """
    Read settings file and return values for the given timestamp

    Parameters
    ----------
    ws: Mantid workspace

    Returns
    -------
    settings: InstrumentSettings
    """
    settings_dict = dict()
    package_dir, _ = os.path.split(__file__)

    t = ws.getRun()["start_time"].value.split("T")[0]
    timestamp = datetime.date.fromisoformat(t)

    with open(os.path.join(package_dir, "settings.json"), "r") as fd:
        data = json.load(fd)
        for key in data.keys():
            chosen_value = None
            delta_time = None
            for item in data[key]:
                valid_from = datetime.date.fromisoformat(item["from"])
                delta = valid_from - timestamp
                if delta_time is None or (delta.total_seconds() < 0 and delta > delta_time):
                    delta_time = delta
                    chosen_value = item["value"]
            settings_dict[key] = chosen_value
    settings = InstrumentSettings(
        apply_instrument_settings=False,
        source_detector_distance=settings_dict["source-det-distance"],
        sample_detector_distance=settings_dict["sample-det-distance"],
        num_x_pixels=settings_dict["number-of-x-pixels"],
        num_y_pixels=settings_dict["number-of-y-pixels"],
        pixel_width=settings_dict["pixel-width"],
        xi_reference=settings_dict["xi-reference"],
        s1_sample_distance=settings_dict["s1-sample-distance"],
        wavelength_resolution_function=settings_dict["wavelength-resolution-function"],
    )
    return settings


def process_attenuation(ws, thickness=0):
    """
    Correct for absorption by assigning weight to each neutron event

    Parameters
    ----------
    ws
        Mantid workspace to correct
    thickness: float
        Attenuator thickness in cm (default is 0).

    Returns
    -------
        Mantid workspace
        Corrected Mantid workspace
    """
    settings = read_settings(ws)
    if settings.source_detector_distance is not None:
        SDD = settings.source_detector_distance
    else:
        SDD = EventReflectivity.DEFAULT_4B_SOURCE_DET_DISTANCE

    constant = 1e-4 * NEUTRON_MASS * SDD / PLANCK_CONSTANT

    package_dir, _ = os.path.split(__file__)
    mu_abs = np.loadtxt(os.path.join(package_dir, "Cd-abs-factors.txt")).T
    wl_model = mu_abs[0]

    # Turn model into a histogram
    wl_step = wl_model[-1] - wl_model[-2]
    final_wl = wl_model[-1] + wl_step
    wl_model = np.append(
        wl_model,
        [
            final_wl,
        ],
    )
    mu_model = mu_abs[1]
    tof_model = constant * wl_model
    transmission = 1 / np.exp(-mu_model * thickness)
    transmission_ws = api.CreateWorkspace(OutputWorkspace="transmission",
                                          DataX=tof_model, DataY=transmission, UnitX="TOF", NSpec=1)

    ws = api.Multiply(ws, transmission_ws, OutputWorkspace=str(ws))
    return ws


def get_dead_time_correction(ws, template_data):
    """
    Compute dead time correction to be applied to the reflectivity curve.
    The method will also try to load the error events from each of the
    data files to ensure that we properly estimate the dead time correction.

    Parameters
    ----------
    ws
        Workspace with raw data to compute correction for
    template_data : reduction_template_reader.ReductionParameters
        Reduction parameters

    Returns
    -------
    mantid.api.Workspace
        Workspace with dead time correction to apply
    """
    tof_min = ws.getTofMin()
    tof_max = ws.getTofMax()

    run_number = ws.getRun().getProperty("run_number").value
    error_ws = api.LoadErrorEventsNexus("REF_L_%s" % run_number)
    corr_ws = mantid_algorithm_exec(
        dead_time_correction.SingleReadoutDeadTimeCorrection,
        InputWorkspace=ws,
        InputErrorEventsWorkspace=error_ws,
        Paralyzable=template_data.paralyzable,
        DeadTime=template_data.dead_time_value,
        UseDeadTimeThreshold=template_data.use_dead_time_threshold,
        DeadTimeThreshold=template_data.dead_time_threshold,
        TOFStep=template_data.dead_time_tof_step,
        TOFRange=[tof_min, tof_max],
        OutputWorkspace="corr",
    )
    corr_ws = api.Rebin(corr_ws, [tof_min, 10, tof_max])
    return corr_ws


def apply_dead_time_correction(ws, template_data):
    """
    Apply dead time correction, and ensure that it is done only once
    per workspace.

    Parameters
    ----------
    ws
        Workspace with raw data to compute correction for
    template_data : reduction_template_reader.ReductionParameters
        Reduction parameters

    Returns
    -------
    mantid.api.Workspace
        Workspace with dead time correction applied
    """
    if "dead_time_applied" not in ws.getRun():
        corr_ws = get_dead_time_correction(ws, template_data)
        ws = api.Multiply(ws, corr_ws, OutputWorkspace=str(ws))
        api.AddSampleLog(Workspace=ws, LogName="dead_time_applied", LogText="1", LogType="Number")
    return ws


class EventReflectivity:
    """
    Data reduction for the Liquids Reflectometer.
    List of items to be taken care of outside this class:

    - Edge points cropping
    - Angle offset
    - Putting runs together in one R(q) curve
    - Scaling factors

    Pixel ranges include the min and max pixels.

    Parameters
    ----------
    scattering_workspace
        Mantid workspace containing the reflected data
    direct_workspace
        Mantid workspace containing the direct beam data [if None, normalization won't be applied]
    signal_peak : list
        Pixel min and max for the specular peak
    signal_bck : list
        Pixel range of the background [if None, the background won't be subtracted]
    norm_peak : list
        Pixel range of the direct beam peak
    norm_bck : list
        Direct background subtraction is not used [deprecated]
    specular_pixel : float
        Pixel of the specular peak
    signal_low_res : list
        Pixel range of the specular peak out of the scattering plane
    norm_low_res : list
        Pixel range of the direct beam out of the scattering plane
    q_min : float
        Value of lowest q point
    q_step : float
        Step size in Q. Enter a negative value to get a log scale
    q_min : float
        Value of largest q point
    tof_range : list, None
        TOF range,or None
    theta : float
        Theta scattering angle in radians
    dead_time : float
        If not zero, dead time correction will be used
    paralyzable : bool
        If True, the dead time calculation will use the paralyzable approach
    dead_time_value : float
        value of the dead time in microsecond
    dead_time_tof_step : float
        TOF bin size in microsecond
    use_dead_time_threshold : bool
        If True, use a correction of 0 for TOF bins requiring corrections greater than ``DeadTimeThreshold``
    dead_time_threshold : float [optional]
        If ``UseDeadTimeThreshold`` is True, this is the upper limit for dead-time correction ratios
    use_emmission_time : bool
        If True, the emission time delay will be computed
    """

    QX_VS_QZ = 0
    KZI_VS_KZF = 1
    DELTA_KZ_VS_QZ = 3
    THETAF_VS_WL = 4
    INSTRUMENT_4A = 0
    INSTRUMENT_4B = 1
    DEFAULT_4B_SAMPLE_DET_DISTANCE = 1.83
    DEFAULT_4B_SOURCE_DET_DISTANCE = 15.75

    def __init__(
        self,
        scattering_workspace,
        direct_workspace,
        signal_peak,
        signal_bck,
        norm_peak,
        norm_bck,
        specular_pixel,
        signal_low_res,
        norm_low_res,
        q_min=None,
        q_step=-0.02,
        q_max=None,
        tof_range=None,
        theta=1.0,
        instrument=None,
        functional_background=False,
        dead_time=False,
        paralyzable=True,
        dead_time_value=4.2,
        dead_time_tof_step=100,
        use_dead_time_threshold=False,
        dead_time_threshold: Optional[float] = 1.5,
        instrument_settings: InstrumentSettings = None,
        use_emission_time=True,
    ):
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
        self.q_summing = False
        self.dq_over_q = 0
        self.dead_time = dead_time
        self.paralyzable = paralyzable
        self.dead_time_value = dead_time_value
        self.dead_time_tof_step = dead_time_tof_step
        self.use_dead_time_threshold = use_dead_time_threshold
        self.dead_time_threshold = dead_time_threshold
        self.instrument_settings = instrument_settings
        self.use_emission_time = use_emission_time

        # Turn on functional background estimation
        self.use_functional_bck = functional_background

        # Process workspaces
        if self.tof_range is not None:
            self._ws_sc = api.CropWorkspace(
                InputWorkspace=scattering_workspace, XMin=tof_range[0], XMax=tof_range[1],
                OutputWorkspace="_" + str(scattering_workspace)
            )
            if direct_workspace is not None:
                self._ws_db = api.CropWorkspace(
                    InputWorkspace=direct_workspace, XMin=tof_range[0], XMax=tof_range[1],
                    OutputWorkspace="_" + str(direct_workspace)
                )
            else:
                self._ws_db = None
        else:
            self._ws_sc = scattering_workspace
            self._ws_db = direct_workspace

        # Extract meta data
        self.extract_meta_data()

    def extract_meta_data(self):
        """
        Extract meta data from the loaded data file.
        """
        # Get instrument parameters
        if (self.instrument_settings is None) or (not self.instrument_settings.apply_instrument_settings):
            settings = read_settings(self._ws_sc)
        else:
            settings = self.instrument_settings

        self.n_x = settings.num_x_pixels
        self.n_y = settings.num_y_pixels
        self.pixel_width = settings.pixel_width / 1000.0

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
            self.wl_range = [self.tof_range[0] / self.constant, self.tof_range[1] / self.constant]

        # q_min and q_max are the boundaries for the final Q binning
        # We also hold on to the true Q range covered by the measurement
        self.q_min_meas = 4.0 * np.pi / self.wl_range[1] * np.fabs(np.sin(self.theta))
        self.q_max_meas = 4.0 * np.pi / self.wl_range[0] * np.fabs(np.sin(self.theta))

        # Template defaults specify q_min but not q_max. TODO: Need to check if it needs changing.
        if self.q_min is None:
            self.q_min = self.q_min_meas
        if self.q_max is None:
            self.q_max = self.q_max_meas

        # Q binning to use
        self.q_bins = get_q_binning(self.q_min, self.q_max, self.q_step)

        # Catch options that can be turned off
        if self.signal_low_res is None:
            self.signal_low_res = [1, self.n_x - 1]
        if self.norm_low_res is None:
            self.norm_low_res = [1, self.n_x - 1]

    def extract_meta_data_4A(self):
        """
        4A-specific meta data
        """
        run_object = self._ws_sc.getRun()
        self.sample_detector_distance = run_object["SampleDetDis"].getStatistics().mean
        source_sample_distance = run_object["ModeratorSamDis"].getStatistics().mean
        if run_object["SampleDetDis"].units not in ["m", "meter"]:
            self.sample_detector_distance /= 1000.0
        if run_object["ModeratorSamDis"].units not in ["m", "meter"]:
            source_sample_distance /= 1000.0
        self.source_detector_distance = source_sample_distance + self.sample_detector_distance

    def extract_meta_data_4B(self):
        """
        4B-specific meta data

        Distance from source to sample was 13.63 meters prior to the source
        to detector distance being determined with Bragg edges to be 15.75 m.
        """
        if (self.instrument_settings is None) or (not self.instrument_settings.apply_instrument_settings):
            settings = read_settings(self._ws_sc)
        else:
            settings = self.instrument_settings

        self.sample_detector_distance = settings.sample_detector_distance

        # Check that we have the needed meta data for the emission delay calculation
        if self.use_emission_time:
            moderator_available = "BL4B:Chop:Skf2:ChopperModerator" in self._ws_sc.getRun()
            if not moderator_available:
                print("Moderator information unavailable: skipping emission time calculation")
                self.use_emission_time = False

        # Get the source-detector distance
        if self.use_emission_time and not settings.apply_instrument_settings:
            self.source_detector_distance = self._ws_sc.getRun().getProperty("BL4B:Det:TH:DlyDet:BasePath").value[0]
        else:
            self.source_detector_distance = settings.source_detector_distance

    def __repr__(self):
        """
        Generate a string representation of the reduction settings.

        Returns
        -------
        str
            String representation of the reduction settings
        """
        output = "Reduction settings:\n"
        output += f"    source-det: {self.source_detector_distance}\n"
        output += f"    sample-det: {self.sample_detector_distance}\n"
        output += f"    pixel-width: {self.pixel_width}\n"
        output += f"    pixel-dimensions: {self.n_x} x {self.n_y}\n"
        output += f"    WL: {self.wl_range[0]} {self.wl_range[1]}\n"
        output += f"    Q: {self.q_min_meas} {self.q_max_meas}\n"
        output += f"    Theta = {self.theta * 180 / np.pi}\n"
        output += f"    Emission delay = {self.use_emission_time}"
        return output

    def to_dict(self):
        """
        Returns meta-data to be used/stored.

        Returns
        -------
        dict
            Dictionary with meta-data
        """
        if self._ws_sc.getRun().hasProperty("start_time"):
            start_time = self._ws_sc.getRun().getProperty("start_time").value
        else:
            start_time = "live"
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
        return dict(
            wl_min=self.wl_range[0],
            wl_max=self.wl_range[1],
            q_min=self.q_min_meas,
            q_max=self.q_max_meas,
            theta=self.theta,
            start_time=start_time,
            experiment=experiment,
            run_number=run_number,
            run_title=run_title,
            norm_run=norm_run,
            time=time.ctime(),
            dq0=dq0,
            dq_over_q=self.dq_over_q,
            sequence_number=sequence_number,
            sequence_id=sequence_id,
            q_summing=self.q_summing,
        )

    def specular(self, q_summing=False, tof_weighted=False, bck_in_q=False, clean=False, normalize=True):
        """
        Compute specular reflectivity.

        For constant-Q binning, it's preferred to use tof_weighted=True.

        Parameters
        ----------
        q_summing : bool
            Turns on constant-Q binning
        tof_weighted : bool
            If True, binning will be done by weighting each event to the DB distribution
        bck_in_q : bool
            If True, the background will be estimated in Q space using the constant-Q binning approach
        clean : bool
            If True, and Q summing is True, then leading artifact will be removed
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
        if tof_weighted:
            self.specular_weighted(q_summing=q_summing, bck_in_q=bck_in_q)
        else:
            self.specular_unweighted(q_summing=q_summing, normalize=normalize)

        # Remove leading zeros
        r = np.trim_zeros(self.refl, "f")
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

        # Compute Q resolution
        self.dq_over_q = compute_resolution(self._ws_sc, theta=self.theta, q_summing=q_summing)
        # TODO: integrate wavelength component to self.dq_over_q
        self.q_summing = q_summing

        return self.q_bins, self.refl, self.d_refl

    def specular_unweighted(self, q_summing=False, normalize=True):
        """
        Simple specular reflectivity calculation. This is the same approach as the
        original LR reduction, which sums up pixels without constant-Q binning.
        The original approach bins in TOF, then rebins the final results after
        transformation to Q. This approach bins directly to Q.

        Parameters
        ----------
        q_summing : bool
            If True, sum the data in Q-space.
        normalize : bool
            If True, normalize the reflectivity by the direct beam.

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
            q_summing=q_summing,
        )

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
            direct_beam, d_direct_beam = self._reflectivity(
                self._ws_db, peak_position=0, peak=self.norm_peak, low_res=self.norm_low_res,
                theta=self.theta, q_summing=False
            )

            # Direct beam background could be added here. The effect will be negligible.
            if self.norm_bck is not None:
                direct_beam_bck, d_direct_beam_bck = self.norm_bck_subtraction()
                direct_beam -= direct_beam_bck
                d_direct_beam = np.sqrt(d_direct_beam**2 + d_direct_beam_bck**2)
            db_bins = direct_beam > 0

            # Error propagation for ratio without dividing by refl (robust if some items of refl are zero)
            refl_normalized = np.zeros_like(refl)  # array ofzeros with same shape as refl
            d_refl_normalized = np.zeros_like(d_refl)
            r = refl[db_bins]
            r_error = d_refl[db_bins]
            db = direct_beam[db_bins]
            db_error = d_direct_beam[db_bins]
            refl_normalized[db_bins] = r / db
            d_refl_normalized[db_bins] = np.sqrt((r_error**2) / (db**2) + (r**2) * (db_error**2) / (db**4))

            # rename to match prior code and avoid issues
            refl[db_bins] = refl_normalized[db_bins]
            d_refl[db_bins] = d_refl_normalized[db_bins]

            # Hold on to normalization to be able to diagnose issues later
            self.norm = direct_beam[db_bins]
            self.d_norm = d_direct_beam[db_bins]

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

        Parameters
        ----------
        q_summing : bool
            If True, sum the data in Q-space.
        bck_in_q : bool
            If True, subtract background along Q lines.

        Returns
        -------
        q_bins
            The Q bin boundaries
        refl
            The reflectivity values
        d_refl
            The uncertainties in the reflectivity values
        """
        # Event weights for normalization
        ## this gets the DB events and produces an array of weights and wavelength values
        ## think this could be replaced by a pre-processed transmission array...i.e. need the wl_dist and wl_bins.
        db_charge = self._ws_db.getRun().getProtonCharge()
        wl_events, wl_weights = self._get_events(self._ws_db, self.norm_peak, self.norm_low_res)
        wl_dist, wl_bins = np.histogram(wl_events, bins=100, weights=wl_weights)
        _bin_width = wl_bins[1:] - wl_bins[:-1]
        wl_dist = wl_dist / db_charge / _bin_width
        wl_middle = [(wl_bins[i + 1] + wl_bins[i]) / 2.0 for i in range(len(wl_bins) - 1)]

        wl_var, _ = np.histogram(wl_events, bins=wl_bins, weights=wl_weights**2)
        wl_std = np.sqrt(wl_var) / (db_charge * _bin_width)

        refl, d_refl = self._reflectivity(
            self._ws_sc,
            peak_position=self.specular_pixel,
            peak=self.signal_peak,
            low_res=self.signal_low_res,
            theta=self.theta,
            q_summing=q_summing,
            wl_dist=wl_dist,
            wl_bins=wl_middle,
            wl_error=wl_std
        )

        if self.signal_bck is not None:
            refl_bck, d_refl_bck = self.bck_subtraction(wl_dist=wl_dist, wl_bins=wl_middle,
                                                        q_summing=bck_in_q, wl_std=wl_std)
            refl -= refl_bck
            d_refl = np.sqrt(d_refl**2 + d_refl_bck**2)

        # Does there need to be a clean-up of non-zeros like in the unweighted version?

        self.refl = refl
        self.d_refl = d_refl
        return self.q_bins, refl, d_refl

    def _roi_integration(self, ws, peak, low_res, q_bins=None, wl_dist=None, wl_bins=None,
                         wl_std=None, q_summing=False):
        """
        Integrate a region of interest and normalize by the number of included pixels.

        The options are the same as for the reflectivity calculation.
        If wl_dist and wl_bins are supplied, the events will be weighted by flux.
        If q_summing is True, the angle of each neutron will be recalculated according to
        their position on the detector and place in the proper Q bin.
        """
        q_bins = self.q_bins if q_bins is None else q_bins
        refl_bck, d_refl_bck = self._reflectivity(
            ws,
            peak_position=0,
            q_bins=q_bins,
            peak=peak,
            low_res=low_res,
            theta=self.theta,
            q_summing=q_summing,
            wl_dist=wl_dist,
            wl_bins=wl_bins,
            wl_error=wl_std
        )

        _pixel_area = peak[1] - peak[0] + 1.0
        refl_bck /= _pixel_area
        d_refl_bck /= _pixel_area
        return refl_bck, d_refl_bck

    def bck_subtraction(self, normalize_to_single_pixel=False, q_bins=None,
                        wl_dist=None, wl_bins=None, wl_std=None, q_summing=False):
        """
        Perform background subtraction on the signal.
        This method provides a higher-level call for background subtraction,
        hiding the ranges needed to define the Region of Interest (ROI).

        Parameters
        ----------
        normalize_to_single_pixel : bool
            If True, normalize the background to a single pixel.
        q_bins
            array of bins for the momentum transfer (q) values.
        wl_dist
            Array of wavelength (wl) values.
        wl_bins
            Array of bins for the wavelength (wl) values.
        wl_std
            Array of errors for wavelength normalization.
        q_summing : bool
            If True, sum the q values.

        Returns
        -------
        mantid.api.Workspace
            The workspace with the background subtracted.
        """
        # Sanity check
        if len(self.signal_bck) == 2 and self.use_functional_bck:
            msg = "Background range incompatible with functional background: "
            msg += "switching to averaging"
            print(msg)
            self.use_functional_bck = False

        if self.use_functional_bck:
            return background.functional_background(
                self._ws_sc,
                self,
                self.signal_peak,
                self.signal_bck,
                self.signal_low_res,
                normalize_to_single_pixel=normalize_to_single_pixel,
                q_bins=q_bins,
                wl_dist=wl_dist,
                wl_bins=wl_bins,
                wl_std=wl_std,
                q_summing=q_summing
            )
        else:
            return background.side_background(
                self._ws_sc,
                self,
                self.signal_peak,
                self.signal_bck,
                self.signal_low_res,
                normalize_to_single_pixel=normalize_to_single_pixel,
                q_bins=q_bins,
                wl_dist=wl_dist,
                wl_bins=wl_bins,
                wl_std=wl_std,
                q_summing=q_summing
            )

    def norm_bck_subtraction(self):
        """
        Higher-level call for background subtraction for the normalization run.
        """
        return background.side_background(
            self._ws_db, self, self.norm_peak, self.norm_bck, self.norm_low_res, normalize_to_single_pixel=False
        )

    def slice(self, x_min=0.002, x_max=0.004, x_bins=None, z_bins=None, refl=None, d_refl=None, normalize=False):
        """
        Retrieve a slice from the off-specular data.
        """
        x_bins = self._offspec_x_bins if x_bins is None else x_bins
        z_bins = self._offspec_z_bins if z_bins is None else z_bins
        refl = self._offspec_refl if refl is None else refl
        d_refl = self._offspec_d_refl if d_refl is None else d_refl

        i_min = len(x_bins[x_bins < x_min])
        i_max = len(x_bins[x_bins < x_max])

        _spec = np.sum(refl[i_min:i_max], axis=0)
        _d_spec = np.sum((d_refl[i_min:i_max]) ** 2, axis=0)
        _d_spec = np.sqrt(_d_spec)
        if normalize:
            _spec /= i_max - i_min
            _d_spec /= i_max - i_min

        return z_bins, _spec, _d_spec

    def _reflectivity(
        self, ws, peak_position, peak, low_res, theta, q_bins=None, q_summing=False,
        wl_dist=None, wl_bins=None, wl_error=None, sum_pixels=True
    ):
        """
        Assumes that the input workspace is not normalized by proton charge.

        In the current workflow, for specular_unweighted this is applied separately to the reflectivity
        run and the normalisation run.
        For specular_weighted, the wl_dist and wl_bins are used to provide the normalisation instead.
        sum_pixels and q_bins are not specified and used with defaults
        low_res is the x_range of interest from template
        peak is data_peak_range from template (i.e. specular ROI)

        It returns reflectivity and error.
        """
        charge = ws.getRun().getProtonCharge()
        ## in standard workflow, q_bins not specified and calculated using the get_q_binning algorithm,
        # with ranges and steps taken from the template.
        _q_bins = self.q_bins if q_bins is None else q_bins

        shape = len(_q_bins) - 1 if sum_pixels else ((peak[1] - peak[0] + 1), len(_q_bins) - 1)
        refl = np.zeros(shape)
        d_refl = np.zeros(shape)
        counts = np.zeros(shape)
        keep_weights = np.zeros(shape)
        # pixel_width only used for q_summing:
        _pixel_width = self.pixel_width if q_summing else 0.0

        ## this first section is very similar to _get_events.
        # TODO: look at refactoring. Linked to tof vs lam vs q too and desired outputs.
        for i in range(low_res[0], int(low_res[1] + 1)):
            for j in range(peak[0], int(peak[1] + 1)):
                # set pixel mapping based on orientation of detector.
                if self.instrument == self.INSTRUMENT_4A:
                    pixel = j * self.n_y + i
                else:
                    pixel = i * self.n_y + j
                # collect the events
                evt_list = ws.getSpectrum(pixel)
                if evt_list.getNumberEvents() == 0:
                    continue

                # collect tof values and apply emission time correction based on flag from template
                tofs = evt_list.getTofs()
                if self.use_emission_time:
                    tofs = self.emission_time_correction(ws, tofs=tofs)
                # convert tof values to wavelength.
                wl_list = tofs / self.constant

                # Gravity correction
                d_theta = self.gravity_correction(ws, wl_list)
                # collect weighted events
                event_weights = evt_list.getWeights()

                # Calculate per-spectum offset in theta for q-summing.
                # For q_summing=True, peak_position is fit in template.py using the peak_finding.py functions.
                # For q_summing=False, _pixel_width=0 and no offset is applied.
                x_distance = _pixel_width * (j - peak_position)
                delta_theta_f = np.arctan(x_distance / self.sample_detector_distance) / 2.0
                # Sign will depend on reflect up or down
                ths_value = ws.getRun()["ths"].value[-1]
                delta_theta_f *= np.sign(ths_value)

                # convert wavelengths into qz. This could be separated to enable output in lam and q.
                qz = 4.0 * np.pi / wl_list * np.sin(theta + delta_theta_f - d_theta)
                # Remove unfeasible negative values:
                valid_mask = qz > 0
                qz = qz[valid_mask]
                wl_list = wl_list[valid_mask]
                event_weights = event_weights[valid_mask]

                # this workflow is used for specular_weighted:
                # matches the bins (in wavelength) of the norm run and refl run.
                # Then converts into q and applies the event weighting.
                if wl_dist is not None and wl_bins is not None:
                    if wl_error is None:
                        wl_error = np.zeros(len(wl_dist))
                        print("No error provided for the normalization weighting")
                    # Computer wavelength weights
                    wl_weights = 1.0 / np.interp(wl_list, wl_bins, wl_dist, np.inf, np.inf)
                    hist_weights = wl_weights * qz / wl_list

                    # Relative error in weights from Poisson statistics.
                    # ***** This needs looking at as wl_dist is already processed and not just the raw counts
                    wl_interp_std = np.interp(wl_list, wl_bins, wl_error, left=0, right=0)
                    wl_interp_dist = np.interp(wl_list, wl_bins, wl_dist, left=np.inf, right=np.inf)
                    rel_err_weights = wl_interp_std / wl_interp_dist

                    #Relative error in wl_list also from Poisson statistics
                    rel_err_wl_list = np.sqrt(wl_list) / wl_list

                    #Combine relative errors
                    rel_err_hist_weights = np.sqrt(rel_err_weights**2 + rel_err_wl_list**2)

                    # Absolute error in histogram weights
                    d_hist_weights = hist_weights * rel_err_hist_weights

                    # Apply event weights
                    hist_weights *= event_weights
                    d_hist_weights *= event_weights

                    # Histogram the data and errors
                    _counts, _ = np.histogram(qz, bins=_q_bins, weights=hist_weights)
                    _norm, _ = np.histogram(qz, bins=_q_bins)
                    _for_errors, _ = np.histogram(qz, bins=_q_bins, weights=d_hist_weights**2)
                    if sum_pixels:
                        refl += _counts
                        counts += _norm
                        keep_weights += _for_errors
                    else:
                        refl[j - peak[0]] += _counts
                        counts[j - peak[0]] += _norm
                        keep_weights[j - peak[0]] += _for_errors
                # this workflow is used for specular_unweighted:
                else:
                    _counts, _ = np.histogram(qz, bins=_q_bins, weights=event_weights)
                    if sum_pixels:
                        refl += _counts
                    else:
                        refl[j - peak[0]] += _counts

        # The following is for information purposes
        # This also sets the summing_threshold for cutting off low Q points in the specular function.
        if q_summing:
            x0 = _pixel_width * (peak_position - peak[0])
            x1 = _pixel_width * (peak_position - peak[1])
            delta_theta_f0 = np.arctan(x0 / self.sample_detector_distance) / 2.0
            delta_theta_f1 = np.arctan(x1 / self.sample_detector_distance) / 2.0

            qz_max = 4.0 * np.pi / self.tof_range[1] * self.constant * np.fabs(np.sin(theta + delta_theta_f0))
            qz_min = 4.0 * np.pi / self.tof_range[1] * self.constant * np.fabs(np.sin(theta + delta_theta_f1))
            mid_point = (qz_max + qz_min) / 2.0
            print("Qz range: ", qz_min, mid_point, qz_max)
            self.summing_threshold = mid_point

        # for the specular weighted workflow, normalises to charge and bin_size
        if wl_dist is not None and wl_bins is not None:
            bin_size = _q_bins[1:] - _q_bins[:-1]
            #non_zero = counts > 0
            # Deal with the case where we don't sum all the bins
            if not sum_pixels:
                bin_size = np.tile(bin_size, [counts.shape[0], 1])

            # Normalize the reflection and error
            refl_norm = refl / charge / bin_size
            d_refl = np.sqrt(keep_weights) / charge / bin_size

            ## Rename back to prior to avoid errors:
            refl = refl_norm

        else:
            d_refl = np.sqrt(np.fabs(refl)) / charge
            refl /= charge

        return refl, d_refl

    def _get_events(self, ws, peak, low_res):
        """
        Return an array of wavelengths for a given workspace.
        """
        wl_events = np.asarray([])
        wl_weights = np.asarray([])

        for i in range(low_res[0], int(low_res[1] + 1)):
            for j in range(peak[0], int(peak[1] + 1)):
                if self.instrument == self.INSTRUMENT_4A:
                    pixel = j * self.n_y + i
                else:
                    pixel = i * self.n_y + j
                evt_list = ws.getSpectrum(pixel)
                tofs = evt_list.getTofs()
                # Correct for emission time as needed
                if self.use_emission_time:
                    tofs = self.emission_time_correction(ws, tofs=tofs)
                wl_list = tofs / self.constant
                wl_events = np.concatenate((wl_events, wl_list))
                weights = evt_list.getWeights()
                wl_weights = np.concatenate((wl_weights, weights))
        return wl_events, wl_weights

    def off_specular(self, x_axis=None, x_min=-0.015, x_max=0.015, x_npts=50, z_min=None, z_max=None,
                      z_npts=-120, bck_in_q=None):
        """
        Compute off-specular

        Parameters
        ----------
        x_axis : int
            Axis selection from QX_VS_QZ, KZI_VS_KZF, DELTA_KZ_VS_QZ
        x_min : float
            Min value on x-axis
        x_max : float
            Max value on x-axis
        x_npts : int
            Number of points in x (negative will produce a log scale)
        z_min : float
            Min value on z-axis (if none, default Qz will be used)
        z_max : float
            Max value on z-axis (if none, default Qz will be used)
        z_npts : int
            Number of points in z (negative will produce a log scale)
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

        wl_events, wl_weights = self._get_events(self._ws_db, self.norm_peak, self.norm_low_res)
        wl_dist, wl_bins = np.histogram(wl_events, bins=100, weights=wl_weights)
        wl_middle = [(wl_bins[i + 1] + wl_bins[i]) / 2.0 for i in range(len(wl_bins) - 1)]

        _refl, _d_refl = self._off_specular(
            self._ws_sc, wl_dist, wl_middle, qx_bins, qz_bins, self.specular_pixel, self.theta, x_axis=x_axis
        )
        db_charge = self._ws_db.getRun().getProtonCharge()
        _refl *= db_charge * (wl_bins[1] - wl_bins[0])
        _d_refl *= db_charge * (wl_bins[1] - wl_bins[0])

        # Background
        if self.signal_bck:
            if bck_in_q is None:
                print("Not implemented")
            else:
                _, refl_bck, d_refl_bck = self.slice(
                    bck_in_q[0], bck_in_q[1], x_bins=qx_bins, z_bins=qz_bins, refl=_refl, d_refl=_d_refl, normalize=True
                )
                _refl -= refl_bck
                _d_refl = np.sqrt(_d_refl**2 + d_refl_bck**2)

        self._offspec_x_bins = qx_bins
        self._offspec_z_bins = qz_bins
        self._offspec_refl = _refl
        self._offspec_d_refl = _d_refl

        return qx_bins, qz_bins, _refl, _d_refl

    def _off_specular(self, ws, wl_dist, wl_bins, x_bins, z_bins, peak_position, theta, x_axis=None):
        charge = ws.getRun().getProtonCharge()
        refl = np.zeros([len(x_bins) - 1, len(z_bins) - 1])
        counts = np.zeros([len(x_bins) - 1, len(z_bins) - 1])

        for j in range(0, self.n_x):
            wl_list = np.asarray([])
            for i in range(self.signal_low_res[0], int(self.signal_low_res[1] + 1)):
                if self.instrument == self.INSTRUMENT_4A:
                    pixel = j * self.n_y + i
                else:
                    pixel = i * self.n_y + j
                evt_list = ws.getSpectrum(pixel)
                wl_events = evt_list.getTofs() / self.constant
                wl_list = np.concatenate((wl_events, wl_list))

            k = 2.0 * np.pi / wl_list
            wl_weights = 1.0 / np.interp(wl_list, wl_bins, wl_dist, np.inf, np.inf)

            x_distance = float(j - peak_position) * self.pixel_width
            delta_theta_f = np.arctan(x_distance / self.sample_detector_distance)
            # Sign will depend on reflect up or down
            ths_value = ws.getRun()["ths"].value[-1]
            delta_theta_f *= np.sign(ths_value)


            theta_f = theta + delta_theta_f

            qz = k * (np.sin(theta_f) + np.sin(theta))
            qz = np.fabs(qz)
            qx = k * (np.cos(theta_f) - np.cos(theta))
            ki_z = k * np.sin(theta)
            kf_z = k * np.sin(theta_f)

            _x = qx
            _z = qz
            if x_axis == EventReflectivity.DELTA_KZ_VS_QZ:
                _x = ki_z - kf_z
            elif x_axis == EventReflectivity.KZI_VS_KZF:
                _x = ki_z
                _z = kf_z
            elif x_axis == EventReflectivity.THETAF_VS_WL:
                _x = wl_list
                _z = np.ones(len(wl_list)) * theta_f

            histo_weigths = wl_weights * _z / wl_list
            _counts, _, _ = np.histogram2d(_x, _z, bins=[x_bins, z_bins], weights=histo_weigths)
            refl += _counts
            _counts, _, _ = np.histogram2d(_x, _z, bins=[x_bins, z_bins])
            counts += _counts

        bin_size = z_bins[1] - z_bins[0]
        d_refl_sq = refl / np.sqrt(counts) / charge / bin_size
        refl /= charge * bin_size

        return refl, d_refl_sq

    def emission_time_correction(self, ws, tofs):
        """
        Correct TOF for emission time delay in the moderator.

        Parameters
        ----------
        ws : mantid.api.Workspace
            Mantid workspace to extract correction meta-data from
        tofs : numpy.ndarray
            Array of uncorrected TOF values

        Returns
        -------
        numpy.ndarray
            Array of corrected TOF values
        """

        mt_run = ws.getRun()
        use_emission_delay = False
        if "BL4B:Chop:Skf2:ChopperModerator" in mt_run:
            moderator_calc = mt_run.getProperty("BL4B:Chop:Skf2:ChopperModerator").value[0]
            t_mult = mt_run.getProperty("BL4B:Chop:Skf2:ChopperMultiplier").value[0]
            t_off = mt_run.getProperty("BL4B:Chop:Skf2:ChopperOffset").value[0]
            use_emission_delay = moderator_calc == 1

        # Parameters from nexus logs are for correction in lambda. The /self.constant converts the tofs to lam.
        if use_emission_delay:
            tofs -= t_off + t_mult * tofs / self.constant
        return tofs

    def gravity_correction(self, ws, wl_list):
        """
        Gravity correction for each event
        Think this works on an array of wavelengths so could work for non-event list too.

        Parameters
        ----------
        ws : mantid.api.Workspace
            Mantid workspace to extract correction meta-data from.
        wl_list : numpy.ndarray
            Array of wavelengths for each event.

        Returns
        -------
        numpy.ndarray
            Array of gravity-corrected theta values for each event, in radians.
        """
        # Xi reference would be the position of xi if the si slit were to be positioned
        # at the sample. The distance from the sample to si is then xi_reference - xi.
        xi_reference = 445
        if ws.getInstrument().hasParameter("xi-reference"):
            xi_reference = ws.getInstrument().getNumberParameter("xi-reference")[0]

        # Distance between the s1 and the sample
        s1_sample_distance = 1485
        if ws.getInstrument().hasParameter("s1-sample-distance"):
            s1_sample_distance = ws.getInstrument().getNumberParameter("s1-sample-distance")[0] * 1000

        xi = 310
        if ws.getInstrument().hasParameter("BL4B:Mot:xi.RBV"):
            xi = abs(ws.getRun().getProperty("BL4B:Mot:xi.RBV").value[0])

        sample_si_distance = xi_reference - xi
        slit_distance = s1_sample_distance - sample_si_distance

        # Angle of the incident beam on a horizontal sample
        # TODO: this will need to be updated with logged value.
        theta_in = -4.0

        # Calculation from the ILL paper. This works for inclined beams.
        # Calculated theta is the angle on the sample

        g = 9.8067  # m/s^2
        h = 6.6260715e-34  # Js=kg m^2/s
        mn = 1.67492749804e-27  # kg

        v = h / (mn * wl_list * 1e-10)
        k = g / (2 * v**2)

        # Define the sample position as x=0, y=0. increasing x is towards moderator
        xs = 0

        # positions of slits
        x1 = sample_si_distance / 1000
        x2 = (sample_si_distance + slit_distance) / 1000

        # height of slits determined by incident theta, y=0 is the sample height
        y1 = x1 * np.tan(theta_in * np.pi / 180)
        y2 = x2 * np.tan(theta_in * np.pi / 180)

        # This is the location of the top of the parabola
        x0 = (y1 - y2 + k * (x1**2 - x2**2)) / (2 * k * (x1 - x2))

        # Angle is arctan(dy/dx) at sample
        theta_sample = np.arctan(2 * k * (x0 - xs)) * 180 / np.pi

        return (theta_sample - theta_in) * np.pi / 180.0


def compute_resolution(ws, default_dq=0.027, theta=None, q_summing=False):
    """
    Compute the Q resolution from the meta data.
    Corrected to include both slits.

    Parameters
    ----------
    ws : mantid.api.Workspace
        Mantid workspace to extract correction meta-data from.
    theta : float
        Scattering angle in radians
    q_summing : bool
        If True, the pixel size will be used for the resolution

    Returns
    -------
    float
        The dQ/Q resolution (FWHM)
    """
    settings = read_settings(ws)

    if theta is None:
        theta = abs(ws.getRun().getProperty("ths").value[0]) * np.pi / 180.0

    if q_summing:
        # Q resolution is reported as FWHM, so here we consider this to be
        # related to the pixel width
        sdd = 1830
        if settings.sample_detector_distance:
            sdd = settings.sample_detector_distance * 1000
        pixel_size = 0.7
        if settings.pixel_width:
            pixel_size = settings.pixel_width

        # All angles here in radians, assuming small angles
        dq_over_q = np.arcsin(pixel_size / sdd) / theta
        print("Q summing: %g" % dq_over_q)
        return dq_over_q

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
    if settings.s1_sample_distance is not None:
        s1_sample_distance = settings.s1_sample_distance * 1000

    # Get slit gap openings.
    s1h = abs(ws.getRun().getProperty("S1VHeight").value[0])
    sih = abs(ws.getRun().getProperty("SiVHeight").value[0])
    xi = abs(ws.getRun().getProperty("BL4B:Mot:xi.RBV").value[0])
    sample_si_distance = xi_reference - xi
    slit_distance = s1_sample_distance - sample_si_distance
    dq_over_q = (s1h + sih) * 0.5 / slit_distance / theta
    return dq_over_q

## New function for resolution, ready for testing.
## Check choices of returned values once ready to implement.
## Needs update for q-summing.
def trapezoidal_distribution_params(ws, Theta_deg=None, FootPrint=None, SlitRatio=None):
    """
    Calculate trapezoidal parameters L, l from beam/slit geometry and
    compute equivalent normal sigma for 68% probability interval to
    enable calculation of angular part of dq_over_q.
    Option to provide footprint and slitratio rather than reading from file.

    Parameters:
    -----------
    ws : mantid.api.Workspace
        Mantid workspace to extract correction meta-data from.
    Theta_deg : float
        Incident beam angle in degrees.
    FootPrint : float (optional)
        Beam footprint length on sample [mm].
    SlitRatio : float (optional)
        Ratio of slit openings S1/Si.

    Returns:
    --------
    L_bottom : float
        Half-width of trapezoidal support [degrees].
    l_top : float
        Half-width of trapezoidal flat-top [degrees].
    sigma_equiv : float
        Equivalent normal sigma for central 68% interval [degrees].
    dth_over_th : float
        Angular part of dq/q
    """

    settings = read_settings(ws)

    # Xi reference would be the position of xi if the si slit were to be positioned
    # at the sample. The distance from the sample to si is then xi_reference - xi.
    xi_ref = 445
    if ws.getInstrument().hasParameter("xi-reference"):
        xi_ref = ws.getInstrument().getNumberParameter("xi-reference")[0]

    # Distance between the s1 and the sample
    dS1Samp = 1485
    if settings.s1_sample_distance is not None:
        dS1Samp = settings.s1_sample_distance * 1000

    xi_pos = abs(ws.getRun().getProperty("BL4B:Mot:xi.RBV").value[0])

    # Derived distances
    dSiSamp = xi_ref - xi_pos
    dS1Si = dS1Samp - dSiSamp

    if Theta_deg is None:
        Theta_deg = abs(ws.getRun().getProperty("ths").value[0])

    # Slit openings - can be taken from a FootPrint and SlitRatio if provided, or the logs
    if FootPrint is not None and SlitRatio is not None:
        SiY = FootPrint * np.sin(np.radians(Theta_deg)) / (1 + (dSiSamp/dS1Si)*(1 + SlitRatio))
        S1Y = SiY * SlitRatio
    else:
        S1Y = abs(ws.getRun().getProperty("S1VHeight").value[0])
        SiY = abs(ws.getRun().getProperty("SiVHeight").value[0])

    # Trapezoidal half-width angles (deg)
    L_bottom = np.degrees(np.arctan((S1Y + SiY) / (2 * dS1Si)))  # full half-width
    l_top = np.degrees(np.arctan((S1Y - SiY) / (2 * dS1Si)))  # flat-top half-width

    sigma_equiv = _find_sigma_68(L_bottom, l_top)
    dth_over_th = sigma_equiv / Theta_deg

    return L_bottom, l_top, sigma_equiv, dth_over_th

# Analytic trapezoidal CDF for use in trapezoial_distribution_params
def _trapezoidal_cdf_analytic(x, L_, l_):
    h = 1 / (L_ + l_)
    cdf = np.zeros_like(x)
    left_tail = x < -L_
    left_slope = (x >= -L_) & (x < -l_)
    flat_top = (x >= -l_) & (x < l_)
    right_slope = (x >= l_) & (x <= L_)
    right_tail = x > L_

    cdf[left_tail] = 0
    cdf[left_slope] = (h / (L_ - l_)) * 0.5 * (x[left_slope] + L_)**2
    F_neg_l = (h / (L_ - l_)) * 0.5 * (L_ - l_)**2
    cdf[flat_top] = F_neg_l + h * (x[flat_top] + l_)
    cdf[right_slope] = 1 - (h / (L_ - l_)) * 0.5 * (L_ - x[right_slope])**2
    cdf[right_tail] = 1
    return cdf

# Find equivalent normal sigma for 68% central interval for use in trapezoidal_distribution_params
def _find_sigma_68(L_, l_, target_prob=0.68):
    def func(x):
        return _trapezoidal_cdf_analytic(np.array([x]), L_, l_)[0] - \
                _trapezoidal_cdf_analytic(np.array([-x]), L_, l_)[0] - target_prob
    return brentq(func, 0, L_)

## Fix the resolution to include si - Done
## Add new function for correction for the shape correction - Done.
## Add to new function the q-summing part
## Add new function for the lambda correction (Jose?)
## Fix the error propagation in the reflectivity - Needs testing
## Fix the qz set to 0
## Fix the pixel angle using the top/bottom
## Function to create DB
## Function to process from the DB run number or pre-processed
## Fix the q-summing saving into the template and being read back in
## Fix q-bin error and auto trimming points

def compute_wavelength_resolution(ws):
    """
    Compute the wavelength resolution from the meta data.

    Parameters
    ----------
    ws : mantid.api.Workspace
        Mantid workspace to extract correction meta-data from.

    Returns
    -------
    tuple of np.ndarray
        (wavelength, d_lambda):
            wavelength: the fitted wavelength values
            d_lambda: the difference between wavelength and the fit

    Raises
    ------
    ValueError : if ws does not have exactly one spectrum
    """

    if ws.spectrumInfo().size() != 1:
        raise ValueError("Workspace must have only one spectrum")

    settings = read_settings(ws)

    out = api.EvaluateFunction(Function=settings.wavelength_resolution_function,
                         InputWorkspace=ws,
                         OutputWorkspace='out')

    wavelength = out.readX(1)
    d_lambda = out.readY(1)

    return np.array(wavelength), np.array(d_lambda)
