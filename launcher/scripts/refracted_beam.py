import sys
import os
import argparse
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.lines as mlines

import mantid
import mantid.simpleapi as api
mantid.kernel.config.setLogLevel(3)

import lmfit
from lmfit.models import GaussianModel, Model


SLD = dict(Si=2.07e-6, Quartz=4.18e-6)


class Refracted(object):

    def __init__(self, ws, material='Si', tof_bin=200, offset=0.01, pixel_size=0.00072):
        """
            Initial process of the raw data to organize it.
        """
        self.run_number = ws.getRun()['run_number'].value

        # Establish the useful TOF range
        self.tof_min = ws.getTofMin()+1100
        self.tof_max = ws.getTofMax()-1100
        print("TOF range used: %g %g" % (self.tof_min, self.tof_max))

        _ws = api.Rebin(InputWorkspace=ws, Params="%s,%s,%s" % (self.tof_min, tof_bin, self.tof_max))
        self.counts = _ws.extractY()
        
        # The x-axis (TOF in this case) in Mantid is stored with bin boundaries.
        self.tof_boundaries = _ws.extractX()[0]
        self.tof_center = (self.tof_boundaries[1:]+self.tof_boundaries[:-1])/2.0
        print("TOF bin boundaries: %s" % len(self.tof_boundaries))

        # Reshape the count array to map it to the detector geometry (256x304 pixels)
        y = np.reshape(self.counts, (256, 304, self.counts.shape[1]))

        # y pixel count vs TOF (integrate over the x pixels)
        self.p_vs_t = np.sum(y, axis=0)
        print("p vs t: %s" % str(self.p_vs_t.shape))

        # y pixel count total (total y counts for all TOF)
        self.y_counts = np.sum(self.p_vs_t, axis=1)

        # pixel number (center of the pixel)
        self.y_pixel = np.arange(self.y_counts.shape[0])+0.5

        # Collect useful meta-data
        self.sld = SLD[material]
        self.theta_sample = ws.getRun()['ths'].value[0]
        tthd = ws.getRun()['tthd'].value[0]
        self.wl_request = ws.getRun()['LambdaRequest'].value[0]

        _n=1-2.07e-6*self.wl_request**2/2/np.pi
        theta_c = self.wl_request*np.sqrt(self.sld/np.pi)
        self.qc = 4*np.sqrt(np.pi*self.sld)

        print("ths=%g; tthd=%g" % (self.theta_sample, tthd))
        print("n=%g" % _n)
        print("theta_c=%g" % theta_c)
        print("q_c=%g" % self.qc)

        # Computed values
        self.offset = offset
        #self.theta = self.theta_sample + self.offset
        self.pixel_size = pixel_size

        # Useful constants
        self.det_distance = 1.83
        source_sample_distance = 13.63
        source_detector_distance = source_sample_distance + self.det_distance

        h = 6.626e-34  # m^2 kg s^-1
        m = 1.675e-27  # kg
        self.constant = 1e-4 * m * source_detector_distance / h
        self.wl_boundaries = self.tof_boundaries / self.constant
        self.wl_center = self.tof_center / self.constant


    def fit_refracted_pixel(self):
        """ 
            Fit the refracted peak as a function of wavelength and obtain the center pixel and width.
            Using this we will be able to find the best offset to map the center pixel to theory
            as a function of wavelength
            
            TODO: auomatically determine the fitting range
        """
        refracted_pixel = []
        d_refracted_pixel = []

        # Compute a good starting point. When fitting the angle, the best guess is 
        #best_guess = -self.theta_sample/2.0
        best_guess = self.center_db + (self.center_r-self.center_db)/4.0
        
        # Determine a good fitting range
        y_pixel_min = int(np.ceil(self.center_db + 5))
        y_pixel_max = int(np.floor(self.center_r - 15))

        for i in range(self.p_vs_t.shape[1]):
            gauss = GaussianModel(prefix='g_')
            pars = gauss.make_params(amplitude=100, center=best_guess, sigma=5)

            fit = gauss.fit(self.p_vs_t[y_pixel_min:y_pixel_max,i], pars, method='leastsq',
                            x=self.y_pixel[y_pixel_min:y_pixel_max])

            a=fit.params['g_amplitude'].value
            c=fit.params['g_center'].value
            width=fit.params['g_sigma'].value

            refracted_pixel.append(c)
            d_refracted_pixel.append(width)

        # Refracted pixel position as a function of wavelength
        self.refracted_pixel = np.asarray(refracted_pixel)
        self.d_refracted_pixel = np.asarray(d_refracted_pixel)

        return refracted_pixel, d_refracted_pixel

    def pixel_prediction(self, wl, offset=0.01):
        """ Compute the expected pixel position given a wavelength and offset parameter """
        # First get the refracted angle
        _angle = self.compute_refracted_angle(wl, offset)

        # Compute the expected pixel position
        _theta_specular = self.theta_sample + offset
        return np.sin(np.pi/180*(_angle + _theta_specular))*self.det_distance/self.pixel_size+self.center_db

    def pixel_prediction_with_size(self, wl, offset=0.01, pixel_size=0.00072):
        """ Compute the expected pixel position given a wavelength and offset parameter """
        # First get the refracted angle
        _angle = self.compute_refracted_angle(wl, offset)

        # Compute the expected pixel position
        _theta_specular = self.theta_sample + offset
        return np.sin(np.pi/180*(_angle + _theta_specular))*self.det_distance/pixel_size+self.center_db

    def compute_angle_from_pixel(self, pixels, offset, pixel_size=None):
        """ Compute the scattering angle from the pixel coordinate """
        if pixel_size is None:
            pixel_size = self.pixel_size

        _theta_specular = self.theta_sample + offset
        return 180.0/np.pi * np.arcsin((pixel_size*(pixels-self.center_db)/self.det_distance)) - _theta_specular

    def compute_refracted_angle(self, wl, offset):
        """ Compute refracted angle as a function of wavelength
        """
        _theta_specular = self.theta_sample + offset
        n = 1 - self.sld * wl**2 / 2 / np.pi
        return -np.arccos(np.cos(np.pi/180*_theta_specular)/n) * 180/np.pi

    def residuals(self, pars, wl, data, err):
        parvals = pars.valuesdict()
        offset = parvals['offset']
        pixel_size = parvals['pixel_size']

        # First get the refracted angle
        _angle = self.compute_refracted_angle(wl, offset)

        # Compute the expected pixel position
        _theta_specular = self.theta_sample + offset
        _refr_pixel = np.sin(np.pi/180*(_angle + _theta_specular))*self.det_distance/pixel_size+self.center_db

        # Compute the reflected pixel position
        _theta_specular = self.theta_sample + offset
        _spec_pixel = np.sin(np.pi/180*(2*_theta_specular))*self.det_distance/pixel_size+self.center_db

        resid = (_refr_pixel-data)**2/err**2 + (_spec_pixel - self.center_r)**2/self.d_center_r**2
        return np.sqrt(resid)

    def fit_offset_and_pixel(self, margin=0.5):
        """
            Fit the position of the refracted beam to extract the calibration offset to ths
            and the pixel size
            
            - margin is a range close to the critical edge where computing the refracted angle
            might return a NaN.
        """
        # Compute wavelength for Qc
        wl_c = 4*np.pi*np.sin(np.pi/180*(self.theta_sample + self.offset))/self.qc - margin
        wl = self.tof_center / self.constant
        index_c = np.max(np.where(wl<wl_c))

        print("Wavelength cutoff = %g  [%s]" % (wl_c, index_c))

        fit_params = lmfit.Parameters()
        fit_params.add('offset', value=self.offset, min=-0.5, max=0.5)
        fit_params.add('pixel_size', value=self.pixel_size, min=0.00065, max=0.00075)

        engine = 'nelder'
        engine = 'leastsq'
        #engine = 'emcee'
        fit = lmfit.minimize(self.residuals, fit_params, args=(wl[:index_c],),
                             method=engine,
                             kws={'data': self.refracted_pixel[:index_c],
                                  'err': self.d_refracted_pixel[:index_c]})

        print(lmfit.fit_report(fit))

        chi2 = fit.chisqr

        _offset=fit.params['offset']
        _pixel_size=fit.params['pixel_size']
        # We may not always get uncertainties
        if hasattr(fit, 'covar'):
            _d_offset = np.sqrt(fit.covar[0][0])
            _d_pixel_size = np.sqrt(fit.covar[1][1])
        else:
            _d_offset = 0
            _d_pixel_size = 0
        print("Chi2 = %g" % chi2)
        print("Fitted offset:     %g +- %g" % (_offset, _d_offset))
        print("Fitted pixel size: %g +- %g" % (_pixel_size, _d_pixel_size))

        return _offset, _d_offset, _pixel_size, _d_pixel_size

    def plot_refraction_map(self, output_dir=None):
        """
            Plot angle vs wavelength, including the refracted beam
        """
        fig, ax = plt.subplots(figsize=(8,8))
        cmap = plt.get_cmap('jet')
        
        # Plot the count data
        pixel_boundaries = np.arange(self.p_vs_t.shape[0]+1)
        _theta_boundaries = self.compute_angle_from_pixel(pixel_boundaries, 
                                                          pixel_size=self.pixel_size,
                                                          offset=self.offset)
        im = ax.pcolormesh(self.wl_boundaries, _theta_boundaries, np.log(self.p_vs_t), cmap=cmap)

        # Plot the specular reflection angle, and the refracted angle as a function of wavelength
        # Reflection (constant)
        _theta_specular = self.theta_sample + self.offset
        r_angle =_theta_specular * np.ones(len(self.wl_boundaries))
        ax.scatter(self.wl_boundaries, r_angle, zorder=1)
        ax.scatter(self.wl_boundaries, -r_angle, zorder=1)
        refr_angle = self.compute_refracted_angle(self.wl_boundaries, self.offset)
        ax.scatter(self.wl_boundaries, refr_angle, zorder=1)

        plt.ylabel("Scattering angle")
        plt.xlabel("wavelength")
        ax.set_ylim([-1, 1])
        plt.title("Offset=%g, Pixel size=%g" % (self.offset, self.pixel_size))
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, "%s-detector-map.png" % self.run_number))
        plt.show()

    def find_initial_peaks(self, count_threshold=100):
        """ Find rough estimate for direct beam and reflected beam """
        peak_left = 0
        peak_right = 0

        for i in range(self.y_counts.shape[0]):
            if self.y_counts[i] > count_threshold:
                peak_left = i
                break

        for i in reversed(range(self.y_counts.shape[0])):
            if self.y_counts[i] > count_threshold:
                peak_right = i
                break

        # Depending on the orientation, set the appropriate peak locations
        #   center_db is the center of the direct beam peak
        #   center_r is the center of the reflected beam
        if self.theta_sample > 0:
            self.center_db = peak_left
            self.center_r = peak_right
        else:
            self.center_db = peak_right
            self.center_r = peak_left

        print("Peak estimates: %g %g" % (peak_left, peak_right))
        return peak_left, peak_right

    def improve_peaks(self, output_dir=None):
        """
            Fit the reflected and direct beam pixels
        """
        # Direct beam
        gauss = GaussianModel(prefix='g_')
        x_min=int(self.center_db-10)
        x_max=int(self.center_db+10)

        pars = gauss.make_params(amplitude=100, center=self.center_db, sigma=5)

        fit = gauss.fit(self.y_counts[x_min:x_max], pars, method='leastsq',
                        x=self.y_pixel[x_min:x_max])

        a=fit.params['g_amplitude']
        c=fit.params['g_center']
        width=fit.params['g_sigma']

        print("Center: %g\t Width: %g" % (c, width))

        # Update the direct beam center
        self.center_db = c

        # Reflected peak
        x_min_r=int(self.center_r-20)
        x_max_r=int(self.center_r+20)

        pars = gauss.make_params(amplitude=100, center=self.center_r, sigma=5)

        fit_r = gauss.fit(self.y_counts[x_min_r:x_max_r], pars, method='leastsq',
                          x=self.y_pixel[x_min_r:x_max_r])

        a=fit_r.params['g_amplitude']
        c=fit_r.params['g_center']
        width=fit_r.params['g_sigma']

        print("Center: %g\t Width: %g" % (c, width))

        if output_dir:
            fig, ax = plt.subplots(figsize=(8,5))
            plt.plot(self.y_pixel, self.y_counts)
            plt.plot(self.y_pixel[x_min:x_max], fit.best_fit, label='direct beam')
            plt.plot(self.y_pixel[x_min_r:x_max_r], fit_r.best_fit, label='reflected')
            plt.xlabel('Counts')
            plt.ylabel('pixel')
            plt.legend()
            plt.savefig(os.path.join(output_dir, "%s-peaks.png" % self.run_number))
            plt.show()

        # Update the reflected beam center
        self.center_r = c
        self.d_center_r = width

    def plot_pixel_vs_tof(self, output_dir=None):
        """
            Plot pixel vs TOF and overlay found peaks
        """
        fig, ax = plt.subplots(figsize=(8,8))
        cmap = plt.get_cmap('jet')
        # Colormesh plots axes using bin boundaries
        im = ax.pcolormesh(self.tof_boundaries, np.arange(self.p_vs_t.shape[0]+1),
                           np.log(self.p_vs_t), cmap=cmap)
        ax.set_title('pixel vs TOF')

        ax.scatter(self.tof_boundaries, np.ones(self.tof_boundaries.shape)*self.center_r, zorder=1, linewidth=2, s=0.2)
        ax.scatter(self.tof_boundaries, np.ones(self.tof_boundaries.shape)*self.center_db, zorder=1, linewidth=2, s=0.2)

        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, "%s-pixel-vs-tof.png" % self.run_number))

        plt.show()

    def compute_reflectivity(self, output_dir):
        print("Computing reflectivity for offset=%g and pixel size %g" % (self.offset, self.pixel_size))
        _theta_specular = self.theta_sample + self.offset

        db = np.sum(self.p_vs_t[125:195,:], axis=0)
        _theta = 180.0/np.pi * np.arcsin((self.pixel_size*(self.y_pixel-self.center_db)/self.det_distance)) - _theta_specular

        wl = self.tof_center/self.constant
        n2 = 1-self.sld*wl**2/np.pi

        # Generate a resolution function
        def _gauss(x, mu):
            _sigma = 0.0004/2.35
            return np.exp(-(x-mu)**2/(2*_sigma**2))/np.sqrt(2*np.pi*_sigma**2)

        sin_theta = np.sin(np.pi/180*_theta_specular)
        q = 4 * np.pi * sin_theta / wl

        alpha = np.sqrt(n2)/sin_theta*np.sqrt(1-(1-sin_theta**2)/n2+0j)
        _theory = (1-2*alpha.real+alpha**2)/(1+2*alpha.real+alpha**2)

        # Because we don't have a linear binning, we're doing the Q resolution brute force
        theory = []
        for i in range(len(q)):
            _pt_sum = 0
            _total = 0
            for j in range(1, len(q)):
                _frac = _gauss(q[i], q[j])
                _pt_sum += _frac*_theory[j]
                _total += _frac
            theory.append(_pt_sum/_total)

        sc = np.sum(self.p_vs_t[175:195,:], axis=0)
        ref = sc/db
        err = np.sqrt(db)/db

        np.savetxt(os.path.join(output_dir, "%s-reflectivity.txt" % self.run_number),
                   np.asarray([q, ref, err]).T)

        fig, ax = plt.subplots(figsize=(8,8))
        plt.errorbar(q, ref, yerr=err)
        plt.plot(q, theory)
        plt.savefig(os.path.join(output_dir, "%s-reflectivity.png" % self.run_number))
        plt.show()


def process_run(run_number, output_dir, material='Si'):
    """
        Process file and extract offset and pixel size
    """
    ws_sc = api.LoadEventNexus("REF_L_%s.nxs.h5" % run_number)

    refr = Refracted(ws_sc, material=material, offset=0.0, pixel_size=0.0007)

    # Find rough peak positions
    refr.find_initial_peaks()

    # Find final peak positions
    refr.improve_peaks(output_dir=output_dir)

    #refr.plot_pixel_vs_tof(output_dir=output_dir)

    # Fit the location of the refracted peak as a function of wavelength
    refracted_pixel, d_refracted_pixel = refr.fit_refracted_pixel()

    _offset, _d_offset, _pixel_size, _d_pixel_size = refr.fit_offset_and_pixel()

    # Store new values
    refr.offset = _offset
    refr.pixel_size = _pixel_size

    refr.plot_refraction_map(output_dir=output_dir)

    fitted_pixel_prediction = refr.pixel_prediction_with_size(refr.wl_center,
                                                              offset=_offset,
                                                              pixel_size=_pixel_size)

    fig, ax = plt.subplots(figsize=(8,5))
    plt.errorbar(refr.wl_center, refracted_pixel, yerr=d_refracted_pixel, label='data')
    plt.plot(refr.wl_center, fitted_pixel_prediction, label='fitted')
    ax.set_ylim([refr.center_db, refr.center_db+0.7*(refr.center_r-refr.center_db)])
    plt.legend()
    plt.savefig(os.path.join(output_dir, "%s-refracted-pixel.png" % refr.run_number))

    refr.compute_reflectivity(output_dir=output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument('run_number', type=str,
                        help='Run number to process')
    parser.add_argument('output_dir', type=str,
                        help='Output directory')
    parser.add_argument('--material', type=str, dest='material',
                        help='Si or Quartz', required=False, default='Si')


    # Parse arguments
    args = parser.parse_args()
    process_run(args.run_number, output_dir=args.output_dir, material=args.material)
