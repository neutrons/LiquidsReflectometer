import os
from pathlib import Path

import lr_reduction.binary_processing as BP
import lr_reduction.nr_tools as tools
import numpy as np
from matplotlib import pyplot as plt


class Direct_Beam:

    def __init__(self, experiment_id=None, NEXUSpath=None, savepath=None):
        # if experiment_id is provided the IPTS-based convention will be used. Otherwise callers
        # should supply explicit NEXUSpath/savepath when constructing the
        # object or when calling create_db.
        self.savepath = None
        self.NEXUSpath = None
        # experiment identifier (e.g. 'IPTS-36119' or '36119'). If provided to
        # __init__ or create_db the code will use the IPTS convention to build
        # paths when explicit paths are not supplied.
        self.experiment_id = None
        # apply constructor overrides: explicit NEXUSpath/savepath or experiment_id
        if NEXUSpath is not None:
            self.NEXUSpath = NEXUSpath
        if savepath is not None:
            self.savepath = savepath
        if experiment_id:
            self.experiment_id = self._normalize_experiment_id(experiment_id)
            # override the paths using convention when experiment_id provided
            self.NEXUSpath = f"/SNS/REF_L/{self.experiment_id}/nexus/"
            self.savepath = f"/SNS/REF_L/{self.experiment_id}/shared/transmission/"
        #self.Cd_foils = [57.5, 126.5, 123.0]      #microns for 50, 100A, 100B Cd foils
        self.Cd = [57.5,
            126.5,
            126.5+123.0,
            2*126.5+2*123.0] #microns Cd for each attenuator
        self.Chop2_cut_fn = [2.077, -16818.0]        # linear fit to chopper cut time
        self.Icut = 1e-10                 # threshold cut off any data below Icut (noisy)
        self.DTCcut = 1.25                # threshold to cut off any data above DTCcut
        self.DTCcut_config1 = 1.5         # threshold to cut off any data above DTCcut_config1 for the lowest lambda run
        self.CutOffset = 1                # offset to adjust the chopper cut position
        self.dMod=15500                   # moderator-to-detector distance
        self.t0=[0.114,0.0295]            # for a linear fit to emission time as a function of lambda: t0[0]+t0[1]*Lambda=emission time
        self.y_ROI = [130,170]            # y pixels to include in the direct beam spectrum
        self.low_res = [75,190]           # x pixels to include in the direct beam spectrum
        self.n_y = 304
        self.n_x = 256

    def _normalize_experiment_id(self, expid):
        """Normalize experiment id to form 'IPTS-XXXX'. Accepts numeric or 'IPTS-XXXX' strings."""
        if expid is None:
            return None
        s = str(expid).strip()
        if s.upper().startswith('IPTS-'):
            return s
        # raw numeric id
        if s.isdigit():
            return f'IPTS-{s}'
        # try to extract digits
        import re
        m = re.search(r'(\d{4,})', s)
        if m:
            return f'IPTS-{m.group(1)}'
        return s

    def create_db(self, run_list, save_name, plot=True, mu_file = None, flip_atten=False, return_traces=False, experiment_id=None):
        """
        Create a direct beam spectrum from the given run list including an attenuation correction for Cd foils (if used).
        Saves the file to the specified location, with header information of the DB pixel position.

        Parameters:
        run_list: list of run numbers to include in the direct beam spectrum
        save_name: name of the output file to save the direct beam spectrum to
        plot: whether to plot the direct beam spectrum
        mu_file: name of the file containing the Cd linear attenuation coefficient data including path. If None will read from settings file.
        flip_atten: whether to flip the attenuator values (earlier runs had an issue in the log files)
        return_traces: whether to return the individual traces for each run in addition to the combined spectrum
        experiment_id: optional experiment ID to use for determining default paths (overrides instance-level ID)

        Returns:
        lam_out: array of wavelength values for the direct beam spectrum
        int_out: array of intensity values for the direct beam spectrum
        err_out: array of error values for the direct beam spectrum
        """

        # loop over the Cd spectra measurements
        LAM = []
        INT = []
        ERR = []
        traces = []
        db_pixel = []
        Cd_values = []

        if plot:
            plt.figure()
        # determine which NEXUS and save paths to use for this create_db call.
        # precedence: explicit experiment_id argument -> instance experiment_id set in __init__ -> instance NEXUSpath/savepath
        nexus_base = None
        save_base = None
        if experiment_id:
            expid = self._normalize_experiment_id(experiment_id)
            nexus_base = f"/SNS/REF_L/{expid}/nexus/"
            save_base = f"/SNS/REF_L/{expid}/shared/transmission/"
        elif getattr(self, 'experiment_id', None):
            nexus_base = f"/SNS/REF_L/{self.experiment_id}/nexus/"
            save_base = f"/SNS/REF_L/{self.experiment_id}/shared/transmission/"
        else:
            nexus_base = self.NEXUSpath
            save_base = self.savepath

        # validate we have usable paths
        if not nexus_base or not save_base:
            raise ValueError('NEXUSpath/savepath not set: provide experiment_id or supply NEXUSpath/savepath when constructing Direct_Beam or call create_db with experiment_id')

        for i, run in enumerate(run_list):
            # get header info from Nexus: atten and chop2 phase
            fname = os.path.join(nexus_base, f'REF_L_{run}.nxs.h5')

            tof_array, y_tof_corr, error_array_corr, log_values, DTC_corr = BP.convert_to_binary(fname, self.low_res, collapse_x = True, tofbin=50, tofmax=100000, tofmin=0, deadtime=4.2, tof_step=100)
            T = tof_array * 1000
            DTC = DTC_corr
            y_tof_corr = np.flipud(y_tof_corr) # flip the y-axis to match the orientation of the detector
            error_array_corr = np.flipud(error_array_corr)

            # Find the DB pixel
            ypix = np.linspace(self.n_y-1, 0, self.n_y)
            mask = (ypix >= self.y_ROI[0]+10) & (ypix < self.y_ROI[1]+10)
            y_db = ypix[mask]
            i_db = np.sum(y_tof_corr[mask,:], axis=1)
            par, fit, bkg = tools.fit_peak(y_db, i_db, peaktype="supergauss", bkgtype="linear")
            db_pixel.append(par[1])

            y_crop = y_tof_corr[self.y_ROI[0]:self.y_ROI[1], :]
            I = np.sum(y_crop, axis=0)
            E = np.sum(error_array_corr[self.y_ROI[0]:self.y_ROI[1], :], axis=0)

            # TODO: can this be outside the loop, but need the time.
            if mu_file is None:
                # Read from settings.json if not provided
                settings = tools.read_settings(time=log_values['start_time'])
                mu_file = settings['cd-attenuator-correction-file']
                print(f'Using Cd attenuation file: {mu_file}')

            # read in Cd linear attenuation coefficient data
            L_ENDF, mu_ENDF = np.loadtxt(mu_file, unpack=True, skiprows=1)

            # Need to split some parts out into separate functions if the logic is correct.
            Cd_thickness = self._extract_cd_values(log_values, flip_atten)
            print(f'Run {run}: Cd thickness = {Cd_thickness:.5f} cm')
            Cd_values.append(Cd_thickness)

            if run == run_list[-1]: low_tag = True
            else: low_tag = False
            T, I, E, DTC = self._trim_and_chop(T, I, E, DTC, log_values["chop2_PD"], lowest = low_tag)
            T = T * 1e-3 # convert to ms
            #print(f'Run {run}: After trimming and chopping, {len(T)} points remain')

            # TODO: use logic from nr_reduction_calc
            # convert to Lambda
            L=3956*T/self.dMod
            L=3956*(T-(self.t0[0]+self.t0[1]*L))/self.dMod

            mu = np.interp(L, L_ENDF, mu_ENDF)
            trans = np.exp(-mu*Cd_thickness)

            # Need to interpret the scale multiplier and apply the correction
            # This allows the program to handle slit-scan DBs
            try:
                n = log_values["scale_multiplier"]
                print(f'Using {n} as scale multiplier')
            except Exception as e:
                print(f'Using 1.0 as scale multiplier b/c {e} not in {repr(log_values)}')
                n = 1.0

            I_trans = I / trans
            if plot:
                plt.plot(L, I_trans, 'o', markersize=1)
            # store trace for optional return
            traces.append((L, I_trans, run))
            LAM.extend(L)
            INT.extend(I_trans * n * n)
            ERR.extend((E * n * n)/trans)

        LAM = np.array(LAM)
        INT = np.array(INT)
        ERR = np.array(ERR)

        lam_out, int_out, err_out = self._lam_error_sort(LAM, INT, ERR)

        if plot:
            plt.plot(lam_out,int_out, '-k')
            plt.ylim(min(int_out), max(int_out)*1.5)
            plt.ylabel('I')
            plt.xlabel('Lambda [A]')
            plt.yscale('log')
            plt.show()

        db_pixel_mean, db_pixel_sigma, tthd_mean, tthd_sigma = self._calc_average_db_pixel(db_pixel, log_values["tthd"])

        header = (
            f"Processed DB spectrum from spectrum maker\n"
            f"Cd attenuator = {Cd_values}\n"
            f"db_runs = {run_list}\n"
            f"db_pixel = {db_pixel_mean}\n"
            f"db_pixel_err = {db_pixel_sigma}\n"
            f"tthd = {tthd_mean}\n"
            f"tthd_err = {tthd_sigma}\n"
            "columns = lambda intensity error"
        )

        array = np.column_stack((lam_out,int_out,err_out))
        # check if folder exists and create if not
        Path(save_base).mkdir(parents=True, exist_ok=True)
        np.savetxt(Path(save_base, save_name), array, header = header, delimiter='\t')
        if return_traces:
            return lam_out, int_out, err_out, traces
        return lam_out, int_out, err_out

    def _calc_average_db_pixel(self, db_pixel_list, tthd):
        """
        Wrapper to calculate the average DB pixel position and its standard deviation,
        as well as the average tthd and its standard deviation for use in header.
        """
        MeanPos = np.round(np.mean(db_pixel_list),2)
        SigmaPos = np.round(np.std(db_pixel_list),2)
        MeanTTHD = np.round(np.mean(tthd),2)
        SigmaTTHD = np.round(np.std(tthd),2)
        print(f'Average DB pixel: {MeanPos} +/- {SigmaPos}')
        return MeanPos, SigmaPos, MeanTTHD, SigmaTTHD

    def _lam_error_sort(self, LAM, INT, ERR):
        """
        Wrapper to sort the lambda arrays and handle duplicates in overlap regions, alongside errors.
        """
        idx = np.argsort(LAM)
        LAM = LAM[idx]
        INT = INT[idx]
        ERR = ERR[idx]

        unique_lam = np.unique(LAM)

        lam_out = []
        int_out = []
        err_out = []

        for lam in unique_lam:
            mask = (LAM == lam)
            i = INT[mask]
            e = ERR[mask]

            w = 1.0 / e**2

            i_bar = np.sum(w * i) / np.sum(w)
            e_bar = np.sqrt(1.0 / np.sum(w))

            lam_out.append(lam)
            int_out.append(i_bar)
            err_out.append(e_bar)

        return lam_out, int_out, err_out



    def _get_chop2_cut(self, chop2_phase):
        return (chop2_phase * self.Chop2_cut_fn[0] + self.Chop2_cut_fn[1])/1000 + self.CutOffset


    def _extract_cd_values(self, log_values, flip_atten=False):
        Atten = log_values['Atten']

        if flip_atten:
            Atten = np.floor(1 - Atten)

        Cd_thickness = (self.Cd[0] * Atten[0] +
                        self.Cd[1] * Atten[1] +
                        self.Cd[2] * Atten[2] +
                        self.Cd[3] * Atten[3])* 1e-4  #cm

        return Cd_thickness

    def _trim_and_chop(self, T, I, E, DTC, chop2_phase, lowest = False):

        # trim off points that drop below the intensity threshold
        p = np.where(I >= self.Icut)
        T=T[p]
        I=I[p]
        E=E[p]
        DTC=DTC[p]

        #trim based on chop2 phase
        TCut = self._get_chop2_cut(chop2_phase)
        TCut = TCut * 1000
        p = np.where(T >= TCut)
        T=T[p]
        I=I[p]
        E=E[p]
        DTC=DTC[p]

        # trim off parts above the DTC threshold
        if lowest: DTC_threshold = self.DTCcut_config1
        else: DTC_threshold = self.DTCcut

        above = np.where(DTC > DTC_threshold)[0]
        if above.size == 0:
            p = np.where(DTC < DTC_threshold)[0]
        else:
            p0 = above[-1]
            p = np.where(DTC[p0+1:] < DTC_threshold)[0] + p0 +1

        T=T[p]
        I=I[p]
        E=E[p]
        DTC=DTC[p]

        return T, I, E, DTC

def plot_db(filename):
    data = np.loadtxt(filename, unpack=True, skiprows=1)
    LAM = data[0]
    INT = data[1]
    ERR = data[2]

    plt.errorbar(LAM, INT, yerr=ERR, fmt='o', markersize=1)
    plt.yscale('log')
    plt.xlabel('Lambda [A]')
    plt.ylabel('Intensity')
    plt.title('Direct beam spectrum')
    plt.show()
