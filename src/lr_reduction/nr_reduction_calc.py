"""
Main calculation file for reflectivity reduction. It allows for different q conversion methods.
NR_Reduction.reduce() is main function to call from a script to process multiple runs from an array.
NR_Reduction._reduce_single_run() is the function to call if only process a single run and want to 
bypass the combining (e.g. for the autoreduction process).
For now, it assumes a pre-processed DB file. 
"""
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter1d, uniform_filter1d
import os
import binary_processing as BP
import nr_tools as tools
import datetime
import json


class NR_Reduction:
    """
    Main calculation file for reflectivity reduction. It allows for different q conversion methods.
    NR_Reduction.reduce() is main function to call from a script to process multiple runs from an array.
    NR_Reduction._reduce_single_run() is the function to call if only process a single run and want to 
    bypass the combining (e.g. for the autoreduction process).
    For now, it assumes a pre-processed DB file. 
    """
    
    def __init__(self, config):
        """
        Initialize
        
        Parameters
        ----------
        config : NRReductionConfig
            Configuration object
        """
        self.config = config
        self._validate_config()
        
    def _validate_config(self):
        """
        Validate configuration consistency
        Must include RBnum and DBname in arrays of equal length.
        For optional arrays the defaults are set if not explicitly provided.
        """
        n_settings = len(self.config.RBnum)
        
        if not self.config.DBname:
            raise ValueError("DBname must be set")
        if n_settings == 0:
            raise ValueError("RBnum must be set")
        if len(self.config.DBname) != n_settings:
            raise ValueError(f"DBname length {len(self.config.DBname)} != RBnum length {n_settings}")
        if len(self.config.RB_Ymin) != n_settings or len(self.config.RB_Ymax) != n_settings:
            raise ValueError(f"RB_Ymin/RB_Ymax length must match RBnum length {n_settings}")
        if self.config.LambdaMin is not None and len(self.config.LambdaMin) != n_settings:
            raise ValueError(f"If supplied, LambdaMin expects list equal to number of runs, length {n_settings}")
        if self.config.LambdaMax is not None and len(self.config.LambdaMax) != n_settings:
            raise ValueError(f"If supplied, LambdaMax expects list equal to number of runs, length {n_settings}")
        
        # Set defaults for optional arrays
        if not self.config.ThetaShift:
            self.config.ThetaShift = [0] * n_settings
        if not self.config.useBS:
            self.config.useBS = [1] * n_settings
        if not self.config.ScaleFactor:
            self.config.ScaleFactor = [1] * n_settings
        if not self.config.tof_min:
            self.config.tof_min = [0] * n_settings  
        if not self.config.tof_max:
            self.config.tof_max = [50000] * n_settings 
            
    def reduce(self, save=True, save_all=False, plot=True):
        """
        Perform the reduction of all angle settings and combine into an output.
        
        save: (optional) save out the final combined Q, R, dR, dQ file
        save_all: (optiona) save out the combined and individual settings
        plot: (optional) show the NR plot on completion

        Returns
        -------
        dict
            Dictionary containing Q, R, dR, dQ arrays
        """
        Q, R, dR, dQ = [], [], [], []
        
        # TODO: Add handling for summed run files.
        for i, rb_num in enumerate(self.config.RBnum):
            print(f'--------------------------------------------')
            result = self._reduce_single_run(i, rb_num)

            # TODO: add better autoscaling options.
            if self.config.AutoScale and i != 0:
                y1=R[i-1][np.where(Q[i-1] >= min(result['q']))]
                y2=result['r'][np.where(result['q'] <= max(Q[i-1]))]
                e1=dR[i-1][np.where(Q[i-1] >= min(result['q']))]
                e2=result['dr'][np.where(result['q'] <= max(Q[i-1]))]   
        
                scale, sigma_scale = tools.weighted_mean(y1,y2,e1,e2)
                result['r'] = result['r'] * scale
                result['dr'] = result['dr'] * scale
                print('Scaling factor: ', np.round(scale,3))

            Q.append(result['q'])
            R.append(result['r'])
            dR.append(result['dr'])
            dQ.append(result['dq'])
            if save_all:
                # save out individual parts
                # TODO: Need to fix the saving logic for multiple runs!!
                self.save_results(result, sname=f"{self.config.Sname}_{i}")
        
        # Combine results for all settings
        Q_combined = np.concatenate(Q)
        R_combined = np.concatenate(R)
        dR_combined = np.concatenate(dR)
        dQ_combined = np.concatenate(dQ)
        
        # Sort by Q for combined data
        idx = np.argsort(Q_combined)
        combine_results = {'Q': Q_combined[idx], 'R': R_combined[idx], 'dR': dR_combined[idx], 'dQ': dQ_combined[idx]}

        if save or save_all:    #TODO: fix the saving parts here this is messy!
            self.save_results(combine_results)
        # TODO: Decide whether to keep in here or have as a separate part after the reduciton....?
        if plot:
            for i, rb_num in enumerate(self.config.RBnum):
                if self.config.plotQ4: 
                    plt.errorbar(Q[i],R[i]*Q[i]**4, yerr=dR[i]*Q[i]**4, xerr=dQ[i],  fmt='o', markersize=1)
                    plt.ylabel(r'$R \cdot Q^4$', fontsize=14)
                    plt.xscale('linear')
                else: 
                    plt.errorbar(Q[i],R[i], yerr=dR[i], xerr=dQ[i],  fmt='o', markersize=1)
                    plt.ylabel('R')
                    plt.xscale('log')
            plt.title('NR data: '+str(rb_num), fontsize=16)
            plt.yscale('log')   # Do we need a toggle on this?
            Angstrom = '\u212B'
            plt.xlabel('Q [1/'+Angstrom+']', fontsize=14)
            plt.show()

        return {
            'Q': Q_combined[idx],
            'R': R_combined[idx],
            'dR': dR_combined[idx],
            'dQ': dQ_combined[idx],
            'Q_per_run': Q,
            'R_per_run': R,
            'dR_per_run': dR,
            'dQ_per_run': dQ,
        }
    
    def _load_binary_from_disk(self, rb_num, i):
        """
        This is not used in the main workflow as want to recalculate based on x-ranges. In future if store binary as 
        both x and y this might be useful.
        """
        # Define expected binary file paths
        nRB_bin = self.config.BINpath / f"{rb_num}{self.config.BINsubname}.npy"
        nRBE_bin = self.config.BINpath / f"{rb_num}{self.config.errBINsubname}.npy"
        nRB_tof = self.config.BINpath / f"{rb_num}{self.config.DTCsubname}.dat"
        
        # Check if all binary files exist
        files_exist = (
            os.path.isfile(nRB_bin) and 
            os.path.isfile(nRBE_bin) and 
            os.path.isfile(nRB_tof)
        )
    
        if files_exist:
            # Load from cached binary files
            print(f"Loading cached binary files for run {rb_num}")
            y_tof_corr = np.load(nRB_bin)
            error_array_corr = np.load(nRBE_bin)
            tof_crc = np.loadtxt(nRB_tof, unpack=True)
            tof_array = tof_crc[0] if isinstance(tof_crc, np.ndarray) and len(tof_crc.shape) > 1 else tof_crc
            log_values = {} # TODO: Better handling here, possible read of header? For now this is not in the workflow.
            return tof_array, y_tof_corr, error_array_corr, log_values
        
        else:
            # Binary files don't exist, compute from Nexus file
            print(f"Computing binary data from Nexus file for run {rb_num}...")
            tof_array, y_tof_corr, error_array_corr, log_values = self._make_binary_files(rb_num, self.config.tof_min[i], self.config.tof_max[i])
            return tof_array, y_tof_corr, error_array_corr, log_values

    def _make_binary_files(self, rb_num, tof_min=0, tof_max=50000):
        """
        Make binary files.
        
        If computes from nexus, returns the TOF array, Y vs TOF data, and error array without saving to disk.
        This allows parameters like lowres to be updated without using stale cached files.
        
        Parameters
        ----------
        rb_num : int
            Run block number
        
        Returns
        -------
        tuple
            (tof_array, y_tof_corr, error_array_corr) - Arrays ready for reduction
        
        Raises
        ------
        FileNotFoundError
            If source Nexus file not found
        RuntimeError
            If binary processing fails
        """
                
        # Get Nexus file path
        nNxRB = self.config.NEXUSpathRB / f"REF_L_{rb_num}.nxs.h5"
        
        if not os.path.isfile(nNxRB):
            raise FileNotFoundError(f"Nexus file not found: {nNxRB}")
        
        lowres = self.config.data_x_range
        
        # Convert to binary
        try:
            tof_array, y_tof_corr, error_array_corr, log_values = BP.convert_to_binary(
                nNxRB,
                lowres=lowres,
                collapse_x=True,
                tofbin=self.config.tof_bin,
                tofmax=tof_max,
                tofmin=tof_min,
                deadtime=self.config.dead_time,
                tof_step=self.config.dead_time_tof_step,
                n_y=304,    # TODO: Work out how to add here when settings not created yet.
                n_x=256
            )
            
            print(f"Binary data computed successfully for run {rb_num}")
            return tof_array, y_tof_corr, error_array_corr, log_values
            
        except Exception as e:
            raise RuntimeError(f"Failed to compute binary data for run {rb_num}: {str(e)}")
    
    def _load_and_extract_lambda(self, i, rb_num):
        """
        Load binary TOF data and convert to lambda space including the emission time correction.
        Return arrays for run and direct beam with matching bins.
        Return additional metadata for later calculations.
        
        Parameters
        ----------
        i : int
            Run index
        rb_num : int
            Run number
            
        Returns
        -------
        tuple
            (lDB, iDB, eDB, DBpixel, DBtthd, q, ypix, RB, RBE, LAMBDA, LambdaBinSize,
             ThCen)
        """
        # Get binary data recompute to ensure correct x-ranges etc.
        tRB, nRB, nRBE, log_values = self._make_binary_files(rb_num, self.config.tof_min[i], self.config.tof_max[i])
        self.log_values = log_values
        # Read in the instrument settings file from the json. # TODO: A little more logic should be added to mimic prior setup.
        settings = self.read_settings()
        self.settings = self.apply_config_overrides(settings)
        self.settings['si_sample_distance'] = self.settings['xi_reference'] - self.log_values['xi']
        self.settings['interslit_distance'] = self.settings['s1_sample_distance'] - self.settings['si_sample_distance']    

        # Calculate lam range if not provided - # TODO: Add smarter calculation.
        if self.config.LambdaMin is None or self.config.LambdaMax is None:
            lam_range = tools.get_lam_range(self.log_values["lam_request"], self.log_values["frequency"])
        if self.config.LambdaMin is None:
            self.config.LambdaMinUse = lam_range[0]
        else:
            self.config.LambdaMinUse = self.config.LambdaMin[i] # TODO: check a better way to do this. Might need different lam per settings.
        if self.config.LambdaMax is None:
            self.config.LambdaMaxUse = lam_range[1]
        else:
            self.config.LambdaMaxUse = self.config.LamdbaMax[i]

        # Flip the arrays to give detector pixel ascending.
        RB = np.flipud(nRB)
        RBE = np.flipud(nRBE)
        
        # Read in pre-processed direct beam file and extract header information from pre-process step.
        lDB, iDB, eDB, metaDB = tools.load_db_file(self.config.DBpath, self.config.DBname[i])
        DBpixel = float(metaDB['db_pixel']) # TODO: Add handling for if this doesn't exist.
        DBtthd = float(metaDB['tthd'])
        # add these to the logs
        self.log_values['DBpixel'] = DBpixel
        self.log_values['DBtthd'] = DBtthd

        # Adjust theta using theta shift if supplied.
        # #TODO: follow through to check abs(). Need sign for gravity etc.
        theta_motor, mode = self._choose_theta_log()
        ThCen = theta_motor + self.config.ThetaShift[i]
        self.log_values['ThCen'] = ThCen
        print(f"Theta, TTHD = {np.round(ThCen, 3)}, {np.round(self.log_values['tthd'], 3)}")
        
        # Create Q vector
        q = tools.log_qvector(self.config.qmin, self.config.qmax, self.config.dqbin)
        # Y pixel array
        ypix = np.linspace(self.settings['num_y_pixels'] - 1, 0, self.settings['num_y_pixels'])

        # Convert TOF to Lambda including the emission time correction.
        if self.config.emission_coefficients is not None:
            t0 = self.config.emission_coefficients
        else:
            t0 = self.log_values['emission_coefficients']
        LAMBDA = 3956 * (tRB) / self.settings['source_detector_distance']
        LAMBDA = 3956 * (tRB - (t0[0] + t0[1] * LAMBDA)) / self.settings['source_detector_distance']
        LambdaBinSize = abs(LAMBDA[1] - LAMBDA[0])
        
        # Trim to desired lambda range
        mask = ((LAMBDA >= self.config.LambdaMinUse) & (LAMBDA <= self.config.LambdaMaxUse))
        LAMBDA = LAMBDA[mask]
        RB = RB[:, mask]
        RBE = RBE[:, mask]
        
        # Interpolate DB lambda to match new RB binning
        iDB = np.interp(LAMBDA, lDB, iDB)
        eDB = np.interp(LAMBDA, lDB, eDB)
        
        # TODO: check what is stored and whether this is the best place.
        self.q = q
        
        return lDB, iDB, eDB, DBpixel, DBtthd, q, ypix, RB, RBE, LAMBDA, LambdaBinSize, ThCen, mode # TODO: should q be returned and a self.q?

    def _fit_and_calculate_theta(self, i, ypix, RB, DBpixel, DBtthd, ThCen):
        """
        Fit specular peak position on detector and calculate corresponding theta
        
        Parameters
        ----------
        i : int
            Run index
        ypix : np.ndarray
            Y pixel array
        RB : np.ndarray
            Run data (2D)
        DBpixel : float
            Direct beam pixel position
        DBtthd : float
            Direct beam tthd angle
        ThCen : float
            Theta value from logs (i.e. starting theta)
            
        Returns
        -------
        tuple
            (Yfit, Ifit, RBpixel, Icalc, ThCen, tthd)
        """
        # Find pixel range to fit based on run Ymin/max values and the background parameters.
        sorted_bkgdROI = self._background_roi_sorter(self.config.BkgROI[i], self.config.RB_Ymin[i], self.config.RB_Ymax[i])
        p1 = tools.find_pixel_index(ypix, sorted_bkgdROI[1]) + int(self.config.peak_pad)
        p0 = tools.find_pixel_index(ypix, sorted_bkgdROI[2]) - int(self.config.peak_pad)
        # TODO: Check the logic on p1 and p0 here? Should it be different to this? I haven't quite matched Erik's...

        Ydata = ypix[p0:p1]
        # Collapse intensity per pixel
        Idata = np.sum(RB[p0:p1, :], axis=1)
        # Fit the peak and extract parameters
        par, fit, bkg = tools.fit_peak(Ydata, Idata, peaktype=self.config.peak_type, bkgtype='linear')
        Yfit = fit[:, 0]
        Ifit = fit[:, 1]
        
        # Calculate theta from fit pixel position and comparison to expected peak position based on DB pixel/tthd positions.
        RBpixel = par[1]
        dPix = RBpixel - DBpixel
        dMM = dPix * self.settings['pixel_width']
        ThetaCalc = np.arcsin(dMM / self.settings['sample_detector_distance']) * 180 / np.pi
        ThetaCalc = ThetaCalc + (self.log_values['tthd'] - DBtthd) / 2
        print(f'Calculated theta: {np.round(ThetaCalc, 3)}, dTheta: {np.round(ThetaCalc - ThCen, 3)}')
        
        # Use calculated value if flag
        if self.config.useCalcTheta:
            ThCen = ThetaCalc
            self.log_values['ThCen'] = ThCen
        else:
            RBpixel = DBpixel
        
        # Calculate expected beam profile on detector using logs
        Icalc = tools.calc_beam_on_detector(Yfit, RBpixel, self.log_values['siY'], self.log_values['s1Y'],
                                            self.settings['interslit_distance'], self.settings['si_sample_distance'], 
                                         self.settings['sample_detector_distance'], self.settings['pixel_width'],
                                         self.config.DetSigma, self.config.DetResFn)
        Icalc = (Icalc * par[0]) + bkg
        
        # Plot beam profile if requested - compares to calculated profile from instrument geometry.
        if self.config.plotON:
            fig, ax = plt.subplots()
            ax.plot(Ydata, Idata, 'ok')
            ax.plot(Yfit, Ifit, '-r', label=f'{self.config.peak_type} fit')
            ax.plot(Yfit, bkg, '--r', label='background')
            ax.plot(Yfit, Icalc, '-b', label='Calculated')
            ax.plot([self.config.RB_Ymin[i], self.config.RB_Ymin[i]], [min(Idata), max(Idata)],
                    '--k', linewidth=1, label='Data ROI')
            ax.plot([self.config.RB_Ymax[i], self.config.RB_Ymax[i]], [min(Idata), max(Idata)],
                    '--k', linewidth=1)
            ax.set_title('Vertical beam distribution', fontsize=16)
            ax.set_xlabel('Detector Pixel', fontsize=14)
            ax.legend()
            plt.show()
        
        return Yfit, Ifit, RBpixel, Icalc, ThCen, bkg

    def _calculate_jacobian(self, LAMBDA, Thv, include_dqdtheta=True):
        """Calculate Jacobian determinant for Q-space transformation.
        
        Args:
            LAMBDA: Wavelength array
            Thv: Theta values (can be 1D or 2D)
            include_dqdtheta: If True, use full 2D Jacobian; if False, use only dqdlambda
        
        Returns:
            J: Jacobian magnitude
        """
        dqdtheta = (4 * np.pi / LAMBDA) * np.cos(np.radians(Thv)) * np.pi / 180
        dqdlambda = -(4 * np.pi / LAMBDA**2) * np.sin(np.radians(Thv))
        
        if include_dqdtheta:
            jacobian = abs(np.sqrt(dqdtheta**2 + dqdlambda**2))
        else:
            jacobian = abs(dqdlambda)
        
        return jacobian
    
    def _pixel_to_q(self, Thv, LAMBDA, LambdaBinSize, dTheta, dLambda, ThetaBinSize=0):
        """Calculate Q values and resolution for a given theta and wavelength.
        
        Args:
            Thv: Theta value
            LAMBDA: Wavelength value
            LambdaBinSize: Bin size in wavelength
            dTheta: Theta sigma values
            dLambda: Wavelength sigma values
            ThetaBinSize: Optional inclusion of theta bin size for non-constantTOF
        
        Returns:
            qcen, qlo, qhi, dq_val
        """
        qcen = 4 * np.pi * np.sin(np.radians(Thv)) / LAMBDA
        qlo = 4 * np.pi * np.sin(np.radians(Thv - ThetaBinSize / 2)) / (LAMBDA + LambdaBinSize / 2)
        qhi = 4 * np.pi * np.sin(np.radians(Thv + ThetaBinSize / 2)) / (LAMBDA - LambdaBinSize / 2)       
        dq_val = qcen * np.sqrt((dLambda / LAMBDA)**2 + (dTheta / Thv)**2)
        
        return qcen, qlo, qhi, dq_val
    
    def _calculate_theta_and_bins(self, ymmRB, ThCen):
        """Calculate theta mapping based on reduction method (constantQ or meanTheta).
        
        Args:
            ymmRB: Y pixel positions in mm relative to RB center
            ThCen: Center theta angle to be used in calculation
        
        Returns:
            Theta: Theta array for each pixel
            ThetaBinSize: Bin sizes for array
            dTheta: Theta sigma
        """
        # calculate the mean theta from slits
        print(self.config.method)
        if self.config.method == 'meantheta':
            # Calculate mean theta from slits
            MeanTheta, dTheta = self.calc_mean_theta(ymmRB, self.log_values['s1Y'], self.log_values['siY'], 
                                                    self.settings['sample_detector_distance'], self.settings['si_sample_distance'],
                                                    self.settings['interslit_distance'], self.config.DetSigma,
                                                    self.config.DetResFn, self.config.plotON)
            
            MeanTheta = MeanTheta + ThCen
            Theta = MeanTheta
        elif self.config.method == 'constantq':
            # Constant Q-line: fixed theta based on geometry
            Theta = ThCen + np.arctan(ymmRB / self.settings['sample_detector_distance']) * 180 / np.pi

            # TODO: Come back to here!!
            sigma_pix = self.settings['pixel_width'] / np.sqrt(12)   #convert tophat to equivalent gaussian
            sigma_y = np.sqrt(self.config.DetSigma**2 + sigma_pix**2)
            dTheta_val = np.degrees(np.arctan((sigma_y/self.settings['sample_detector_distance'])))
            dTheta = np.full(len(Theta), dTheta_val)
        else:
            raise ValueError(f"Theta calculation only defined for config.method 'constantQ' or 'meanTheta'")
            
        # Store theta bins for next calculation.
        ThetaBinSize = abs(np.diff(Theta))
        lastBinSize = ThetaBinSize[-1]
        ThetaBinSize = np.concatenate([ThetaBinSize, [lastBinSize]])
        
        return Theta, ThetaBinSize, dTheta
    
    def _calc_detector_convolution(self,m, b, ResFn, sigma, plotON=True):
        
        verts = [
                tools.intersect(m[0], b[0], m[1], b[1]),
                tools.intersect(m[1], b[1], m[2], b[2]),
                tools.intersect(m[2], b[2], m[3], b[3]),
                tools.intersect(m[3], b[3], m[0], b[0]),
            ]
            
        verts = np.array(verts)
        Det_corners = np.vstack([verts, verts[0]])
            
        #pad the y axis to prepare for convolution
        if ResFn == 'rectangular':
            pad = sigma * np.sqrt(12) 
            width = sigma * np.sqrt(12)

        if ResFn == 'gaussian':
            pad = sigma * 4.0 

            
        # define the theta and Y vectors, Y is padded to allow convolution
        nt=1000
        ny=1000
        tstep = ((max(verts[:,0]))-(min(verts[:,0])))/nt
        ystep = ((max(verts[:,1])+pad)-(min(verts[:,1])-pad))/ny
        tvec = np.linspace(min(verts[:,0]), max(verts[:,0])+tstep, nt)
        yvec = np.linspace(min(verts[:,1])-pad, max(verts[:,1])+pad+ystep, ny)

        # create a smeared array
        smear=np.zeros((ny,nt))
            
        #loop through theta
        for i in range(nt):
            #calculate the hi and lo Y values in this theta slice
            hi = np.minimum(tvec[i]*m[3] + b[3], tvec[i]*m[2] + b[2])
            lo = np.maximum(tvec[i]*m[0] + b[0], tvec[i]*m[1] + b[1])
            
            # make a tophat function
            p = np.where((yvec >= lo) & (yvec <= hi))
            d = yvec*0
            d[p] = d[p]+1

            # smear it and store it
            if ResFn == 'rectangular': 
                smear[:,i] = uniform_filter1d(d, size=int(width/ystep)) 
            if ResFn == 'gaussian': 
                smear[:,i] = gaussian_filter1d(d, sigma = sigma/ystep)     
        
        # mask out empty rows
        mask = np.any(smear != 0, axis = 1)
        smear = smear[mask]
        yvec = yvec[mask]
            
        if plotON:
            fig, ax = plt.subplots()
            ax.imshow(smear, origin='lower', aspect='auto', extent=[min(tvec),max(tvec),min(yvec),max(yvec)], cmap='magma')
            ax.plot(Det_corners[:,0], Det_corners[:,1], '--k', linewidth=1)
            ax.set_xlabel('dTheta [deg]', fontsize=14)
            ax.set_ylabel('Detector position [mm]', fontsize=14)
            ax.set_title('Angular beam distribution', fontsize=16)
            plt.show()
        return yvec, smear, tvec

    def calc_mean_theta(self, yvalues, si_H, s1_H, d_sam_det, d_si_sam, d_s1_si, sigma, ResFn, plotON):
        # Calculate the Mean and Sigma of theta angles as a function of Y position on the detector
        # create a fine grid of Y positions, convolve with the detector resolution
        # calculate the mean and sigma for provided Y values from the convolved array
        
        # define positive angles and up and negative angles as down
        a, y, m, b = tools.calc_beam_geometry_from_slits(si_H, s1_H, d_s1_si, d_si_sam, d_sam_det, radians=False)

        # convolve theta-y polygon with detector resolution
        yvec, smear, tvec = self._calc_detector_convolution(m, b, ResFn, sigma, plotON)

        # Compute mean and sigma analytically
        mean_theta = []
        sigma_theta = []

        for Y in yvalues:
            pos = np.where(abs(Y-yvec) == min(abs(Y-yvec)))
            
            # simple line cut through convolved array
            mu = np.sum(smear[pos,:] * tvec) / np.sum(smear[pos,:])
            var = np.sum(smear[pos,:] * (tvec - mu)**2) / np.sum(smear[pos,:])
            sig = np.sqrt(var)
            
            mean_theta.append(mu)
            sigma_theta.append(sig)
            
        mean_theta = np.array(mean_theta)
        sigma_theta = np.array(sigma_theta)

        return mean_theta, sigma_theta


    def plot_theta_lam(self, LAMBDA, Theta, Rarr, iDB, qlines=True, set_ylims=None, set_xlims=None): #TODO: decide if this should be in the class?
            '''
            Generate 2D plot of theta-lambda with option to overlay lines of constant q from calculation
            
            LAMBDA: array of lambda values
            Theta: corresponding array of theta values
            Rarr: 
            iDB: 
            qlines: Option to overlay q lines
            set_ylims: Option to specify y limits, otherwise uses extent of Theta array. (e.g. [0.3,0.8])
            set_xlims: Option to specify x limits, otherwise uses LambdaMin and LambdaMax from config. (e.g. [2.2, 9])
            '''
            fig, ax = plt.subplots()
            RN = Rarr * iDB
            positives = RN[RN > 0]
            vmin = np.percentile(positives, 25)
            cmap = 'magma'
            ax.pcolormesh(LAMBDA, abs(Theta), RN, norm=LogNorm(vmin=vmin, vmax=RN.max()),
                          shading='auto', cmap=cmap)
            
            if qlines:
                dqbin_plot = self.config.dqbin * 20
                qvs = tools.log_qvector(self.config.qmin, self.config.qmax, dqbin_plot)
                for qv in qvs:
                    qline = np.degrees(np.arcsin(qv * LAMBDA / 4 / np.pi))
                    ax.plot(LAMBDA, qline, '--g', linewidth=1.5)
            
            if set_ylims:
                ax.set_ylim(set_ylims[0],set_ylims[1])
            else:
                ax.set_ylim(min(abs(Theta)), max(abs(Theta)))
            if set_xlims:
                ax.set_xlim(set_xlims[0],set_xlims[1])
            else:
                ax.set_xlim(self.config.LambdaMinUse, self.config.LambdaMaxUse)

            ax.set_xscale('log')
            ax.set_title('Q integration lines', fontsize=16)
            ax.set_xlabel('Lambda [Å]', fontsize=14)
            ax.set_ylabel('Theta [deg]', fontsize=14)
            plt.show()
    
    def _background_roi_sorter(self, BkgROI, ymin, ymax):
        # Wrapper to handle the background region when set by the template.
        # This has 2 zeroes in the array when only 1 background region is selected.

        # Find number of zero values
        sorted_offsets = np.sort(BkgROI)
        num_zeros = np.sum(sorted_offsets == 0)
        if num_zeros == 0:
            return sorted_offsets
        elif num_zeros == 2:
            # assume adjacent to ROI and replace the zeros with the ymin and ymax
            sorted_offsets[0] = ymin
            sorted_offsets[1] = ymax
            resorted_offsets = np.sort(sorted_offsets)
            return resorted_offsets

    def background_subtract(self, LAMBDA, R, E, R_crop, E_crop, y_roi, ypix,
                            ymin, ymax, BkgROI, ploton):
        """
        Calculate and apply background subtraction.
        This is handled in lambda space so uses a different function to lr_reduction.
        Background range supplied as a single array of size 4 to be start and stop
        values either side of the specular peak. There will be a separate helper function 
        for handling when this is adjacent to the specular (i.e. provides zeros in the template)
        
        LAMBDA: Description
        R: Description
        E: Description
        ypix: Description
        ymin: Description
        ymax: Description
        BkgROI_offsets: Description
        ploton: Description
        """

        background_idx = self._background_roi_sorter(BkgROI, ymin, ymax)

        # masks
        mask_b1 = (ypix >= background_idx[0]) & (ypix <= background_idx[1])
        mask_b2 = (ypix >= background_idx[2]) & (ypix <= background_idx[3])

        y1 = 0.5 * (background_idx[0] + background_idx[1])
        y2 = 0.5 * (background_idx[2] + background_idx[3])

        # mean background and errors
        b1 = np.mean(R[mask_b1, :], axis=0)
        b2 = np.mean(R[mask_b2, :], axis=0)

        eb1 = np.sqrt(np.sum(E[mask_b1, :]**2, axis=0)) / np.sum(mask_b1)
        eb2 = np.sqrt(np.sum(E[mask_b2, :]**2, axis=0)) / np.sum(mask_b2)

        # line parameters per TOF
        a = (b2 - b1) / (y2 - y1)
        c = b1 - a * y1

        # uncertainties
        var_bkg = (eb1**2 + eb2**2) / 2
        
        if ploton:
            ll=np.where((ypix > min(background_idx[0]-5,background_idx[1]-5)) & (ypix < max(background_idx[2]+5,background_idx[3]+5)))
            fig, ax = plt.subplots()
            # TODO: Need to look at the plot settings here. Think the vmin, vmax should allow for varying signals to stay visible.
            log_data = np.log(R[ll[0],:]+0.00001)
            ax.imshow(log_data, vmin=np.percentile(log_data,0.5), vmax=np.percentile(log_data, 99.5), aspect='auto', extent=[LAMBDA.min(), LAMBDA.max(),
                                                        min(background_idx[0]-5,background_idx[1]-5),
                                                        max(background_idx[2]+5,background_idx[3]+5)], cmap='magma')
            ax.axhline(y=ymin, color='green', linestyle='--')
            ax.axhline(y=ymax, color='green', linestyle='--')
            ax.axhline(y=background_idx[0], color='red', linestyle='--', linewidth=1)
            ax.axhline(y=background_idx[1], color='red', linestyle='--', linewidth=1)
            ax.axhline(y=background_idx[2], color='red', linestyle='--', linewidth=1)
            ax.axhline(y=background_idx[3], color='red', linestyle='--', linewidth=1) 
            ax.set_title('Background subtraction ROIs', fontsize=16)
            ax.set_xlabel('Lambda [Å]', fontsize=14)
            ax.set_ylabel('Detector Pixel', fontsize=14)
            plt.show()
            
        # subtract background
        for i in range(R.shape[1]):
            bkg = a[i] * y_roi + c[i]
            
            R_crop[:, i] -= bkg
            E_crop[:, i] = np.sqrt(E_crop[:, i]**2 + var_bkg[i])

        return R_crop, E_crop

    def _choose_theta_log(self):
        """
        Look-up of earth vs beam centered from log values and choose thi or ths respectively.
        Used for initial theta value and flag on gravity correction.

        return:
            theta: value of thi or ths
            mode: 0 for "earth", 1 for "beam"
        """
        # Determine if earth or beam centered
        if "coordinates" in self.log_values:
            mode = self.log_values["coordinates"] # This is one that shows earth vs beam center 0=earth; 1=beam
        else:
            if self.log_values["op_mode"] == "Free Liquid":
                mode = 0 # i.e. earth
            else:
                mode = 1 # i.e. beam

        # Pick correct thi/ths log for angle
        if mode == 0:
            theta = self.log_values["thi"]
        else:
            theta = self.log_values["ths"]

        return theta, mode

    def _reduce_single_run(self, i, rb_num, save=True):
        """
        Reduce a single run setting using the pre-defined config.
        
        Parameters
        ----------
        i : int
            Run index within the set to be combined (e.g. for angle-variable settings)
        rb_num : int
            Run number
            
        Returns
        -------
        dict
            Reduced data for this run setting of q, r, dr, dq
        """
        # Load and extract data to lambda space
        # Note: q and lDB are available as self.q, not needed in unpacking
        _, iDB, eDB, _, _, _, ypix, RB, RBE, LAMBDA, LambdaBinSize, _, mode = self._load_and_extract_lambda(i, rb_num)
        ## NOTE: Is the q value needed from the output above?! Track through the self.q...?

        # Crop to y-pixel ROI
        mask_roi = (ypix >= self.config.RB_Ymin[i]) & (ypix <= self.config.RB_Ymax[i])

        R_mask = RB[mask_roi, :].copy()
        E_mask = RBE[mask_roi, :].copy()
        y_roi = ypix[mask_roi]

        # Apply background subtraction
        if self.config.useBS[i]:
            Rarr, REarr = self.background_subtract(LAMBDA, RB, RBE, R_mask, E_mask, y_roi, ypix, self.config.RB_Ymin[i],
                                             self.config.RB_Ymax[i], self.config.BkgROI[i], self.config.plotON)
        else:
            Rarr = R_mask
            REarr = E_mask
        #TODO: pull the plotting out from the BS part?!

        # For constantTOF, use 1D TOF binning
        if self.config.method == "constanttof":
            iRB = np.sum(Rarr, axis=0)
            eRB = np.sqrt(np.sum(REarr**2, axis=0))
            # TODO: check the zero removal part!!
            # Remove zeros
            mask = (iDB != 0) & (iRB != 0)
            LAMBDA = LAMBDA[mask]
            iDB = iDB[mask]
            eDB = eDB[mask]
            iRB = iRB[mask]
            eRB = eRB[mask]
            # for continuity, probably needs clearing up:
            Rarr = iRB
            REarr = eRB
        
        # Normalize by incident spectrum & propagate error
        R0 = Rarr.copy()       
        Rarr, REarr = tools.divide_propagate_error(R0, REarr, iDB, eDB)
        
        # Remove NaNs - #TODO: check if this is still correct...
        nan_idx = ~np.isfinite(Rarr)
        Rarr[nan_idx] = 0
        REarr[nan_idx] = 0

        # Fit peak and calculate theta #TODO: check the call and whether need all of these?
        _, _, RBpixel, _, _, _ = self._fit_and_calculate_theta(
            i, ypix, RB, self.log_values['DBpixel'], self.log_values['DBtthd'], self.log_values['ThCen'])

        if self.config.method != "constanttof":
            # Calculate theta for each pixel over ROI
            ypixRB = ypix[(ypix >= self.config.RB_Ymin[i]) & (ypix <= self.config.RB_Ymax[i])] - RBpixel
            ymmRB = ypixRB * self.settings['pixel_width']
            Theta, ThetaBinSize, dTheta = self._calculate_theta_and_bins(ymmRB, self.log_values['ThCen'])

            # Plot 2D lambda vs theta data including overlay of q-summing lines
            if self.config.plotON:
                self.plot_theta_lam(LAMBDA, Theta, Rarr, iDB, qlines=True)
        
            # Remove truncated Q-lines based on qline threshold
            lmin = 4 * np.pi * np.sin(np.radians(min(abs(Theta)))) / self.q #TODO: track through self.q vs qvals...!!
            lmax = 4 * np.pi * np.sin(np.radians(max(abs(Theta)))) / self.q
            Qline_fraction = (np.minimum(self.config.LambdaMaxUse, lmax) - 
                            np.maximum(self.config.LambdaMinUse, lmin)) / (lmax - lmin)
        
            mask = (Qline_fraction >= self.config.Qline_threshold)
            q_vals = self.q[mask]
            Qline_fraction = Qline_fraction[mask]
        else:
            # Make arrays matching sizes for constantTOF to have common loops below
            # Probably a neater way to do this!!
            dummy_array = np.zeros([1, Rarr.shape[0]])
            dummy_array[0] = Rarr
            Rarr = dummy_array
            dummy_array = np.zeros([1, REarr.shape[0]])
            dummy_array[0] = REarr
            REarr = dummy_array
            Theta = np.array([self.log_values['ThCen']])
            ThetaBinSize = np.array([LAMBDA[1] - LAMBDA[0]])
            dTheta = np.array([tools.dTheta_Sigma(self.log_values['siY'], self.log_values['s1Y'], self.settings['interslit_distance'])])
            q_vals = self.q
            Qline_fraction = np.ones(len(q_vals)) # check this!!
        
        # Lambda resolution
        dLambda = tools.dLambda_Sigma(LAMBDA)
        
        # Gravity correction
        # For this function expects e.g. 4.0 for downward angle
        if mode == 1:
            IncTheta = self.config.IncidentTheta - self.log_values["thi"]
        else:
            IncTheta = self.log_values["thi"]
        ThetaGC = tools.gravity_correct(LAMBDA, IncTheta, self.settings['si_sample_distance'], self.settings['interslit_distance'])
        
        # Transform to Q-space
        ## something is wrong with q_vals shape. Might need to reassign self.q to q_vals...?
        #r = self.q * 0
        #dr = self.q * 0
        #dq = self.q * 0
        #Jsum = self.q * 0
        r = q_vals * 0
        dr = q_vals * 0
        dq = q_vals * 0
        Jsum = q_vals * 0

        for T in range(Rarr.shape[0]):
            # Apply gravity correction
            Thv = abs(Theta[T] + ThetaGC)
            
            if self.config.method == "constanttof":
                # Jacobian determinant
                J = self._calculate_jacobian(LAMBDA, Thv, include_dqdtheta=False)
                theta_bin = 0 # TODO: check this logic is correct!!
            else:
                # Jacobian determinant
                J = self._calculate_jacobian(LAMBDA, Thv, include_dqdtheta=True)
                theta_bin = ThetaBinSize[T]
            
            # Loop through lambda to map to Q space
            for L in range(Rarr.shape[1]):
                # Calculate Q values
                qcen, qlo, qhi, dqval = self._pixel_to_q(
                    Thv[L], LAMBDA[L], LambdaBinSize, dTheta[T],
                    dLambda[L], ThetaBinSize=theta_bin)
                # TODO: Comparison of this and the lambda expression.
                if (qhi <= max(q_vals)) & (qlo >= min(q_vals)):
                    # find q bins overlapped by the pixel
                    idx = np.flatnonzero((q_vals >= qlo) & (q_vals <= qhi))
                    # if no overlap, put into nearest q bin
                    if idx.size == 0:
                        idx = np.array([np.argmin(np.abs(qcen - q_vals))])
                    # spread evenly over all overlapped bins. #TODO: Future, investigate options here.
                    wt = J[L] / idx.size
                    r[idx] = r[idx] + Rarr[T, L] * wt
                    dr[idx] = dr[idx] + (REarr[T, L] * wt)**2
                    dq[idx] = dqval
                    Jsum[idx] = Jsum[idx] + wt

        FAC = Jsum / Rarr.shape[0]
        
        # Remove NaNs and zeros and keep region within qline fraction #TODO: work out whether the qline_fraction part is needed both here and above.
        mask = (np.isfinite(r) & (FAC != 0))
        q_vals, r, dr, dq, FAC, Qline_fraction = (x[mask] for x in (q_vals, r, dr, dq, FAC, Qline_fraction))
        
        r = r / FAC * Qline_fraction
        dr = np.sqrt(dr) / FAC * Qline_fraction
        
        # Apply scale factor if specified.
        r = r * self.config.ScaleFactor[i]
        dr = dr * self.config.ScaleFactor[i]

        # Normalize to critical edge region if flag applied.
        if self.config.Normalize and i == 0:
            NormV = np.mean(r[np.where(q_vals <= self.config.Qnorm)])
            r = r / NormV
            dr = dr / NormV
            print(f'Normalization factor: {np.round(1/NormV, 3)}')

        return {'q': q_vals, 'r': r, 'dr': dr, 'dq': dq}

    def save_results(self, results, sname = None):
        """
        Save results as .dat file with header
        
        Parameters
        ----------
        results : dict
            Results from reduce() method
        """
        array = np.column_stack((results['Q'], results['R'], results['dR'], results['dQ']))
        
        # TODO: Sort out the header to include extra information...
        head = (
            f"NR_runs = {self.config.RBnum}\n"
            f"DB = {self.config.DBname}\n"
            f"Method = {self.config.method}\n"
            f"Normalize = {self.config.Normalize}\n"
            f"Autoscale = {self.config.AutoScale}\n"
            f"columns = Q, R, dR, dQ\n"
            f"{'---' * 20}"
        )
        if not sname:
            output_file = self.config.Spath / f"{self.config.Sname}.dat"
        else:
            output_file = self.config.Spath / f"{sname}.dat"
        np.savetxt(output_file,
                  array, header=head, delimiter='\t')
        print(f"Saved combined result to {output_file}")

    # TODO: This should align with prior workflow but needs more work to link into the logic of flags
    #  for using the flag of overwriting instrumnet settings. Is essentially the same as event_reduction
    #  version but has a different logic for checking the date.
    def read_settings(self):
        """
        Read settings file and return values for the given timestamp

        Returns
        -------
        settings
        """
        settings_dict = dict()
        package_dir, _ = os.path.split(__file__)

        timestamp = datetime.datetime.fromisoformat(self.log_values['start_time']).date()

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
        key_map = {
            'source_detector_distance': "source-det-distance",
            'sample_detector_distance': "sample-det-distance",
            'num_x_pixels': "number-of-x-pixels",
            'num_y_pixels': "number-of-y-pixels",
            'pixel_width': "pixel-width",
            'xi_reference': "xi-reference",
            's1_sample_distance': "s1-sample-distance",
            'wavelength_resolution_function': "wavelength-resolution-function",
        }

        settings_output = {
            new_key: settings_dict[old_key]
            for new_key, old_key in key_map.items()
        }

        settings_output['sample_detector_distance'] *= 1000 # Code here expects these in mm.
        settings_output['source_detector_distance'] *= 1000
        settings_output['s1_sample_distance'] *= 1000

        return settings_output
    
    def apply_config_overrides(self, settings: dict) -> dict:
        """
        Apply overrides to the instrument settings if not None in the config
        and returns updated dictionary.
        
        """
        overrides = {
            'xi_reference': getattr(self.config, 'xi_ref', None),
            'num_x_pixels': getattr(self.config, 'nx', None),
            'num_y_pixels': getattr(self.config, 'ny', None),
            'source_detector_distance': getattr(self.config, 'dMod', None),
            'sample_detector_distance': getattr(self.config, 'dSampDet', None),
            'pixel_width': getattr(self.config, 'mmpix', None),
            's1_sample_distance': getattr(self.config, 'dS1Samp', None),
        }

        for key, value in overrides.items():
            if value is not None:
                settings[key] = value

        return settings
