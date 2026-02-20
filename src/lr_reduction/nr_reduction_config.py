from pathlib import Path

class NRReductionConfig:
    """
    Configuration class for NR reduction parameters.
    This defines setable parameters for the workflow including defaults for when they are not provided.
    These can be read and set from a template file.
    """
    
    def __init__(self, method='meanTheta'):
        """
        Initialize reduction configuration.
        Method allows different conversion routes from lambda to q to be selected based on the experimental setup.
        It will allow expansion to other methods as developments occur.
        
        Parameters
        ----------
        method : str
            'constantQ', 'meanTheta', or 'constantTOF'
        """
        self.method = method.lower()
        if self.method not in ['constantq', 'meantheta', 'constanttof']:
            raise ValueError(f"Unknown method: {method}. Use 'constantQ', 'meanTheta', or 'constantTOF'")
        
        # Data configuration - arrays for the data and settings to be processed. These should be of equal length
        self.DBname = []    # This assumes a pre-processed DB file at the moment.
        self.RBnum = []
        self.RB_Ymin = []
        self.RB_Ymax = []
        self.Sname = "reduction_output" # TODO: change the default here to align with defaults from refred
        self.experiment_id = "" # IPTS number
        self.data_x_range = [50,200]

        # Path configuration - Defaults assume IPTS specified and saved into that folder.
        self._Spath_override = None
        self._NEXUSpathRB_override = None
        self._DBpath_override = None
        self._BINpath_override = None
        self.DTCsubname = '_DTC'
        self.BINsubname = '_DTC'
        self.errBINsubname = '_err_DTC'

        # Processing flags
        self.Normalize = False  # Scales the reflectivity to 1 based on a critical edge region (defined by Qnorm)
        self.AutoScale = False  # Toggle to automatically scale between angle settings calculated in the overlap region
        self.useCalcTheta = False   # Toggle to use a fitted specular peak position (relative to the DB position) for theta which overwrites the THS/THI values
        self.plotON = True  # Toggle for plots to show during the reduction steps. Turn off for batch processing etc.
        self.plotQ4 = False # Toggle for the NR plots to be RQ4 vs RQ.
        
        # Background configuration
        self.BkgROI = []    # Need to work out a default...and explain the format here.
        self.useBS = []  # Toggle to use background subtraction per angle
        
        # Q-space configuration #TODO: make better defaults here!!
        # These should be specified within the templates but enable ranges for lambda and q. #TODO: check if zeros are removed?
        self.qmin = 0.001
        self.qmax = 0.5
        self.dqbin = 0.005
        self.Qline_threshold = 1.0  # Fraction of q-line required to be within the bin to be included in the non-constantTOF mode.
        
        self.LambdaMin = None    #If None this will be calculated from the chopper ranges. If supplied should be an array.
        self.LambdaMax = None
        self.tof_max = [] # Need a better way to autoset this...
        self.tof_min = []
        self.tof_bin = 50   # Check defaults and how these go through.

        # Theta/Lambda configuration
        self.ThetaShift = []  # Per-angle theta shift in degrees
        self.ScaleFactor = []  # Per-angle scale factor
        self.Qnorm = 0.015  # Q threshold for normalization
        
        # Instrument geometry - generally read from the instrument settings (if None) but can be overwritten.
        self.mmpix = None  # pixel size in mm
        self.dSampDet = None # sample-to-detector distance
        self.ny = None  # number of vertical pixels
        self.nx = None  # number of vertical pixels
        self.dMod = None  # moderator-to-detector distance
        self.xi_ref = None  # xi=0 reference distance
        self.dS1Samp = None  # slit S1 to sample distance

        self.IncidentTheta = 4.0  # degrees Angle of the beamline relative to earth. Positive is downwards. Will become PV.
        self.emission_coefficients = None  # TOF emission time coefficients
 
        # Dead-time parameters
        self.dead_time = 4.2
        self.dead_time_tof_step = 50
        # TODO: add the dead_time threshold behaviour.

        # Detector resolution parameters
        self.DetResFn = 'rectangular'  # 'rectangular' or 'gaussian'
        self.DetSigma = 0.8  # detector resolution sigma

        # Parameters for peak fitting
        self.peak_pad = 1 # number of extra pixels to include outside the background ranges fit the peak fit range
        self.peak_type = 'supergauss'   # function for the fit. Current options: 'gauss' or 'supergauss'

    # Path configuration - Defaults assume IPTS specified and saved into that folder.
    # #TODO: test the defaults loading part...!        
    @property
    def base_path(self) -> Path:
        return Path("/SNS/REF_L") / self.experiment_id
    
    @property
    def Spath(self) -> Path:
        if self._Spath_override is not None:
            return Path(self._Spath_override)
        return self.base_path / "reduced"
    @property
    def BINpath(self) -> Path:
        if self._BINpath_override is not None:
            return Path(self._BINpath_override)
        return self.base_path / "reduced"
    @property
    def NEXUSpathRB(self) -> Path:
        if self._NEXUSpathRB_override is not None:
            return Path(self._NEXUSpathRB_override)
        return self.base_path / "nexus"
    @property
    def DBpath(self) -> Path:
        if self._DBpath_override is not None:
            return Path(self._DBpath_override)
        return self.base_path / "reduced"
    @Spath.setter
    def Spath(self, value):
        self._Spath_override = value

    @NEXUSpathRB.setter
    def NEXUSpathRB(self, value):
        self._NEXUSpathRB_override = value

    @DBpath.setter
    def DBpath(self, value):
        self._DBpath_override = value
    @BINpath.setter
    def BINpath(self, value):
        self._BINpath_override = value
