import binary_processing as BP
from pathlib import Path
import os
import numpy as np
from matplotlib import pyplot as plt
import nr_tools as tools

class Direct_Beam:

    def __init__(self):
        # TODO: these need better integration with the rest of the config etc.
        self.MUpath = '/SNS/REF_L/IPTS-36119/shared/Cd_transmission/'
        self.savepath = '/SNS/REF_L/IPTS-36119/shared/Cd_transmission/'
        self.NEXUSpath = '/SNS/REF_L/IPTS-34978/nexus/' # TODO this should have similar path setup to the nr_reduction_config.py
        self.Cd_foils = [57.5, 126.5, 123.0]      #microns for 50, 100A, 100B Cd foils
        self.Cd = [self.Cd_foils[0], 
            self.Cd_foils[1], 
            self.Cd_foils[1]+self.Cd_foils[2], 
            2*self.Cd_foils[1]+2*self.Cd_foils[2]] #microns Cd for each attenuator
        self.Chop2_cut_fn = [2.077, -16.8180]        # linear fit to chopper cut time
        self.Icut = 1e-10                 # cut off any data below Icut (noisy)
        self.DTCcut = 1.25                # cut off any data above DTCcut (artifact-y)
        self.DTCcut_config1 = 1.5         # cut off any data above DTCcut_config1 for the first run
        self.CutOffset = 1                # cut off TOF below the chopper cut (chop-y)
        self.dMod=15500                   # moderator-to-detector distance
        self.t0=[0.114,0.0295]            # for a linear fit to emission time as a function of lambda: t0[0]+t0[1]*Lambda=emission time
        self.y_ROI = [130,170]            # y pixels to include in the direct beam spectrum
        self.low_res = [75,190]           # low res TOF range to include in the direct beam spectrum
        self.n_y = 304
        self.n_x = 256

    def create_db(self, run_list, save_name, plot=True, mu_file = 'Cd_mu_2025.dat', flip_atten=False, return_traces=False):
        # read in Cd linear attenuation coefficient data obtained from ENDf
        L_ENDF, mu_ENDF = np.loadtxt(self.MUpath+mu_file, unpack=True, skiprows=1)

        # loop over the Cd spectra measurements
        LAM = []
        INT = []
        ERR = []
        traces = []
        db_pixel = []
        Cd_values = []

        if plot:
            plt.figure()
        for i, run in enumerate(run_list):
            # get header info from Nexus: atten and chop2 phase
            fname = self.NEXUSpath + 'REF_L_'+str(run)+'.nxs.h5'

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

            # Need to split some parts out into separate functions if the logic is correct.
            Cd_thickness = self._extract_cd_values(log_values, flip_atten)

            #print some info for debugging
            print(f'Run {run}: Cd thickness = {Cd_thickness:.2f} cm')
            Cd_values.append(Cd_thickness)

            ## This is the wrong shape. The I and E are 2 dimensional. Need to fix the output used above.
            if run == run_list[-1]: low_tag = True
            else: low_tag = False
            T, I, E, DTC = self._trim_and_chop(T, I, E, DTC, log_values["chop2_PD"], lowest = low_tag)

            T = T * 1e-3 # convert to ms for easier handling
            #print some info for debugging
            print(f'Run {run}: After trimming and chopping, {len(T)} points remain')
    
            # TODO: use logic from nr_reduction_calc
            # convert to Lambda
            L=3956*T/self.dMod
            L=3956*(T-(self.t0[0]+self.t0[1]*L))/self.dMod
            
            mu = np.interp(L, L_ENDF, mu_ENDF)
            trans = np.exp(-mu*Cd_thickness)

            I_trans = I / trans
            if plot:
                plt.plot(L, I_trans, 'o', markersize=1)
            # store trace for optional return
            traces.append((L, I_trans, run))
            LAM.extend(L)
            INT.extend(I_trans)
            ERR.extend(E/trans)

        LAM = np.array(LAM)
        INT = np.array(INT)
        ERR = np.array(ERR)

        lam_out, int_out, err_out = self._lam_error_sort(LAM, INT, ERR)

        #LAM = lam_out
        #INT = int_out
        #ERR = err_out

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
        np.savetxt(self.savepath+save_name, array, header = header, delimiter='\t')
        if return_traces:
            return lam_out, int_out, err_out, traces
        return lam_out, int_out, err_out

    def _calc_average_db_pixel(self, db_pixel_list, tthd):
        MeanPos = np.round(np.mean(db_pixel_list),2)
        SigmaPos = np.round(np.std(db_pixel_list),2)
        MeanTTHD = np.round(np.mean(tthd),2)
        SigmaTTHD = np.round(np.std(tthd),2)
        print(f'Average DB pixel: {MeanPos} +/- {SigmaPos}')
        return MeanPos, SigmaPos, MeanTTHD, SigmaTTHD

    def _lam_error_sort(self, LAM, INT, ERR):
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

    def plot_db(self, filename):
        data = np.loadtxt(self.savepath+filename, unpack=True, skiprows=1)
        LAM = data[0]
        INT = data[1]
        ERR = data[2]

        plt.errorbar(LAM, INT, yerr=ERR, fmt='o', markersize=1)
        plt.yscale('log')
        plt.xlabel('Lambda [A]')
        plt.ylabel('Intensity')
        plt.title('Direct beam spectrum')
        plt.show()

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
        if lowest:
            DTC_threshold = self.DTCcut_config1
        else:
            DTC_threshold = self.DTCcut

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