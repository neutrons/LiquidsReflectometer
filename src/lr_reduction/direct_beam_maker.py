import binary_processing as BP
from pathlib import Path
import os
import numpy as np
from matplotlib import pyplot as plt

class Direct_Beam:

    def __init__(self):
        # TODO: these need better integration with the rest of the config etc.
        self.MUpath = '/SNS/REF_L/IPTS-36119/shared/Cd_transmission/'
        self.savepath = '/SNS/REF_L/IPTS-36119/shared/Cd_transmission/'
        self.NEXUSpath = '/SNS/REF_L/IPTS-34978/nexus/' # TODO this should have similar path setup to the nr_reduction_config.py
        self.Cd_foils = [57.5, 126.0, 122.5]      #microns for 50, 100A, 100B Cd foils
        self.Cd = [self.Cd_foils[0], 
            self.Cd_foils[1], 
            self.Cd_foils[1]+self.Cd_foils[2], 
            2*self.Cd_foils[1]+2*self.Cd_foils[2]] #microns Cd for each attenuator
        self.Chop2_cut_fn = [2.077, -16818.0]        # linear fit to chopper cut time
        self.Icut = 1e-10                 # cut off any data below Icut (noisy)
        self.DTCcut = 1.25                # cut off any data above DTCcut (artifact-y)
        self.CutOffset = 1                # cut off TOF below the chopper cut (chop-y)
        self.dMod=15500                               # moderator-to-detector distance
        self.t0=[0.114,0.0295]                        # for a linear fit to emission time as a function of lambda: t0[0]+t0[1]*Lambda=emission time


    def create_db(self, run_list, save_name, low_res, plot=True, mu_file = 'Cd_mu_2025.dat', flip_atten=False):
        # read in Cd linear attenuation coefficient data obtained from ENDf
        L_ENDF, mu_ENDF = np.loadtxt(self.MUpath+mu_file, unpack=True, skiprows=1)


        # loop over the Cd spectra measurements
        LAM = []
        INT = []
        ERR = []

        for i, run in enumerate(run_list):
            # get header info from Nexus: atten and chop2 phase
            fname = self.NEXUSpath + 'REF_L_'+str(run)+'.nxs.h5'

            e_offset, event_id, error_event_offset, pcharge, cPc, log_values = BP.load_and_extract(fname)
            tof_array, DTC, error_counts = BP.get_deadtime_correction(error_event_offset, e_offset, cPc, 
                                                                      tofbin=50, tofmax=100000, tofmin=0, deadtime=4.2, tof_step=100)
            # Convert to other nomenclature for testing
            tDTC = tof_array
            DTC = DTC
            tof_array, y_tof, error_array = BP.get_y_tof(tof_array, event_id, e_offset, low_res, pcharge, n_y=304, n_x=256)
            T = tof_array
            I = y_tof
            E = error_array

            # Need to split some parts out into separate functions if the logic is correct.
            Cd_thickness = self._extract_cd_values(log_values, flip_atten)

            #print some info for debugging
            print(f'Run {run}: Cd thickness = {Cd_thickness:.2f} cm')

            T, I, E, DTC = self._trim_and_chop(T, I, E, DTC, log_values["chop2_PD"])

            #print some info for debugging
            print(f'Run {run}: After trimming and chopping, {len(T)} points remain')
    
            # TODO: use logic from nr_reduction_calc
            # convert to Lambda
            L=3956*T/self.dMod
            L=3956*(T-(self.t0[0]+self.t0[1]*L))/self.dMod
            
            mu = np.interp(L, L_ENDF, mu_ENDF)
            trans = np.exp(-mu*Cd_thickness)

            if plot:
                plt.plot(L, I/trans, 'o', markersize=1)
            LAM.extend(L)
            INT.extend(I/trans)
            ERR.extend(E/trans)
    
        LAM = np.array(LAM)
        INT = np.array(INT)
        ERR = np.array(ERR)

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

        LAM = lam_out
        INT = int_out
        ERR = err_out

        if plot:
            plt.plot(LAM,INT, '-k')
            plt.ylim(min(INT), max(INT)*1.5)
            plt.ylabel('I')
            plt.xlabel('Lambda [A]')
            plt.yscale('log')

        array = np.column_stack((LAM,INT,ERR))
        np.savetxt(self.savepath+save_name, array, header = 'Lambda\tIntensity\tError', delimiter='\t')

    def _extract_cd_values(self, log_values, flip_atten=False):
        Atten = np.array([log_values['Att1'],log_values['Att2'],log_values['Att3'],log_values['Att4']])
    
        if flip_atten:
            Atten = np.floor(1 - Atten)
    
        Cd_thickness = (self.Cd[0] * Atten[0] + 
                        self.Cd[1] * Atten[1] +
                        self.Cd[2] * Atten[2] +
                        self.Cd[3] * Atten[3])* 1e-4  #cm

        return Cd_thickness
    
    def _trim_and_chop(self, T, I, E, DTC, chop2_phase):
        # trim off points that drop below the intensity threshold
        p = np.where(I >= self.Icut)
        T=T[p]
        I=I[p]
        E=E[p]

        #trim based on chop2 phase
        TCut = (chop2_phase * self.Chop2_cut_fn[0] + self.Chop2_cut_fn[1])/1000.0 + self.CutOffset
        p = np.where(T >= TCut)
        T=T[p]
        I=I[p]
        E=E[p]
        DTC=DTC[p]

        # trim off parts above the DTC threshold
        above = np.where(DTC > self.DTCcut)[0]
        if above.size == 0:
            p = np.where(DTC < self.DTCcut)[0]
        else:
            p0 = above[-1]
            p = np.where(DTC[p0+1:] < self.DTCcut)[0] + p0 +1

        T=T[p]
        I=I[p]
        E=E[p]
        DTC=DTC[p]

        return T, I, E, DTC