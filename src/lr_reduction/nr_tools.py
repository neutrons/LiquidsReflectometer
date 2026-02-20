"""
Subset of functions which may be useful outside of the NRReduction calculation class.
# TODO: Lots of this files still needs commenting.
"""

import numpy as np
from scipy.optimize import curve_fit
from pathlib import Path

def safe_divide(numerator, denominator):
    """Wrapper for safe division with zeros where denominator is zero.
    
    Args:
        numerator: Numerator array
        denominator: Denominator array
    
    Returns:
        Result array with zeros where denominator is zero
    """
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator != 0)

def divide_propagate_error(num, num_err, denom, denom_err):
    """Wrapper to proagate errors through division consistently: (num/denom).
    
    Args:
        num: Numerator array
        num_err: Numerator error array
        denom: Denominator array
        denom_err: Denominator error array
    
    Returns:
        (result, error): Divided array and propagated error
    """
    result = safe_divide(num, denom)
    
    # Error propagation: (num/denom) has error sqrt((num_err/num)^2 + (denom_err/denom)^2) * (num/denom)
    term1 = safe_divide(num_err, num)
    term2 = safe_divide(denom_err, denom)
    error = np.abs(result) * np.sqrt(term1**2 + term2**2)
    
    return result, error

def find_pixel_index(pixels, value, occurence=0):
    """Wrapper to find specific pixel index
    
    Args:
        pixels: Pixel array
        value: Value to find
        occurence: specific occurence or all (e.g. first, last). Use None for full array. default is first.
    
    Returns:
        Index depending on occurrence
    """
    if occurence is not None:
        return np.where(pixels == value)[0][occurence]
    else:
        return np.where(pixels == value)[0]

def log_qvector(qmin, qmax, dqbin):
    """
    Generate logarithmically spaced q vector of defined bin size.
    
    qmin: minimum q
    qmax: maximum q
    dqbin: bin size
    """
    n = int(np.floor(np.log(qmax/qmin) / np.log(1 + dqbin))) + 1
    return qmin * (1 + dqbin) ** np.arange(n + 1)

def weighted_mean(y1,y2,e1,e2, sigma_mask=3):       
    #find the weighted average
    v=y1/y2
    sigma_v = np.sqrt((e1 / y2)**2 + (e2 * y1 / y2**2)**2)
    w=1/sigma_v**2
    # initial weighted mean
    mean=np.sum(v*w)/np.sum(w)

    # identify outlier
    mask = np.abs(v - mean) < sigma_mask * sigma_v
    v = v[mask]
    w = w[mask]

    # Recompute
    mean = np.sum(v * w) / np.sum(w)
    sigma_mean = np.sqrt(1 / np.sum(w))

    return mean, sigma_mean 

# Define functions for peak_fitting:
def gaussian(x,a,x0,sig):
    return a*np.exp(-((x-x0)**2)/(2*sig**2))

def gaussian_constant(x, a, x0, sig, c):
    return gaussian(x, a, x0, sig) + c

def gaussian_slope(x, a, x0, sig, m, b):
    return gaussian(x, a, x0, sig) + m*x + b

def super_gaussian(x, a, x0, sig, ex):
    return a * np.exp(-((np.abs(x - x0)/sig)**ex))

def super_gaussian_slope(x, a, x0, sig, ex, m, b):
    return super_gaussian(x, a, x0, sig, ex) + m*x + b

def super_gaussian_constant(x, a, x0, sig, ex, c):
    return super_gaussian(x, a, x0, sig, ex) + c

def fit_peak(ypix, iY, peaktype="gauss", bkgtype="none"):
    """
    Fit a peak using a gaussian or super-gauss 
    Returns
    -------
    par : Fit parameters
    fit : y and fitted curve
    bkg : Background from fit 
    """
    # TODO: can add check on goodness of fit later.

    # ---- initial guess ----
    mean  = np.sum(iY * ypix) / np.sum(iY)
    var   = np.sum(iY * (ypix - mean)**2) / np.sum(iY)
    sigma = np.sqrt(var)
    A0 = np.max(iY)
    yv = np.linspace(ypix.min(), ypix.max(), 200)

    # --- Select model and initial parameters ---
    if peaktype == "gauss":
        if bkgtype == "none":
            model = gaussian
            p0 = [A0, mean, sigma]
        elif bkgtype == "constant":
            model = gaussian_constant
            p0 = [A0, mean, sigma, 0.0]
        elif bkgtype == "linear":
            model = gaussian_slope
            p0 = [A0, mean, sigma, 0.0, 0.0]
        else:
            raise ValueError("bkgtype must be 'none', 'constant', or 'linear'")
        par, cov = curve_fit(model, ypix, iY, p0=p0)

    elif peaktype == "supergauss":
        n0 = 2.0  # start at Gaussian
        if bkgtype == "none":
            model = super_gaussian
            p0 = [A0, mean, sigma, n0]
            bounds = ([A0/10, -np.inf, 1e-6, 2.0],
                      [A0*10, np.inf, np.inf, 10.0])
        elif bkgtype == "constant":
            model = super_gaussian_constant
            p0 = [A0, mean, sigma, n0, 0.0]
            bounds = ([A0/10, -np.inf, 1e-6, 2.0, -np.inf],
                      [A0*10, np.inf, np.inf, 10.0, np.inf])
        elif bkgtype == "linear":
            model = super_gaussian_slope
            p0 = [A0, mean, sigma, n0, 0.0, 0.0]
            bounds = ([A0/10, -np.inf, 1e-6, 2.0, -np.inf, -np.inf],
                      [A0*10, np.inf, np.inf, 10.0, np.inf, np.inf])
        else:
            raise ValueError("bkgtype must be 'none', 'constant', or 'linear'")
        par, cov = curve_fit(model, ypix, iY, p0=p0, bounds=bounds)

    else:
        raise ValueError("peaktype must be 'gauss' or 'supergauss'")

    # --- Evaluate fit and extract background ---
    f = model(yv, *par)
    if peaktype == "gauss":
        core = gaussian(yv, par[0], par[1], par[2])
    else:
        core = super_gaussian(yv, par[0], par[1], par[2], par[3])

    bkg = f - core if bkgtype != "none" else np.zeros_like(f)
    fit = np.column_stack((yv, f))

    return par, fit, bkg

def calc_beam_geometry_from_slits(si_H, s1_H, d_s1_si, d_si_sam, d_sam_det, radians=True):
    # TODO: fix the descrepancies so the radians tag isn't needed.
    # based on slits and instrument geometry, calculate beam profile
    L = d_si_sam + d_sam_det

    # define positive angles and up and negative angles as down
    # a1, a4 are the max and min angle passing through the top of the slit
    # a2, a3 are the max and min angles passing through the bottom of the slit 
    # Defined here as a=[a1, a2, a3, a4]
    
    # angular ranges at the incident slits
    a_plus  = np.arctan((s1_H + si_H)/(2*d_s1_si))
    a_minus = np.arctan((s1_H - si_H)/(2*d_s1_si))

    a = np.array([ a_plus, a_minus, -a_plus, -a_minus ])

    # propagate vertical position of the beam to the detector
    y = np.array([
        si_H/2 + L*np.tan(a[0]),
       -si_H/2 + L*np.tan(a[1]),
       -si_H/2 + L*np.tan(a[2]),
        si_H/2 + L*np.tan(a[3])
    ])

    # Toggle to have the calculation in degrees. Some point can alter other functions so this isn't needed!
    if not radians:
       a = np.degrees(a)

    #calculate the beam position (Y) as a function of beam angle at Si (X)
    #calculate the slope and offset for the lines that define the parallegram in these coordinates
    # 4 corners of the parralleagram are: (1) top right, (2) bot right, (3) bot left, (4) top left
    # note: intercept, b, is the same for all propagation distances
    #       m12 = m34, m23 = m41, b12 = -b34, b23 = -b41
    m = (y - np.roll(y, -1)) / (a - np.roll(a, -1))
    b = y - m*a

    return a, y, m, b

def calc_beam_on_detector(Ypix, CenPix, Si, S1, dS1Si, dSiSam, dSamDet, mmpix, DetRes, DetResFn):
    # based on slits and instrument geometry, calculate beam profile on the detector
    # read in the Ypix array to calculate over
    # adjust by the center pixel 
    
    Ymm = (Ypix-CenPix) * mmpix

    _, y, m, b = calc_beam_geometry_from_slits(Si, S1, dS1Si, dSiSam, dSamDet)

    #create the intensity distribution
    # calculate the spread of vertical position (height of the polyon) as a function of angle coming out of slits
    # this is the beam height as a function of angle, which scales as intensity as a function of angle
    # this is independent of the propagtion distance
    ss = np.maximum((Ymm - b[3]) / m[3],
                    (Ymm - b[2]) / m[2])

    fs = np.minimum((Ymm - b[0]) / m[0],
                    (Ymm - b[1]) / m[1])

    I = fs - ss
    
    # remove points below y3 and above y1, limits of the polygon
    I[(Ymm < y.min()) | (Ymm > y.max())] = 0
    I[I < 0] = 0
    
    # grid spacing
    dy = Ymm[1] - Ymm[0]
    
    # Convolve with the detector resolution, works better as tophat than gauss!
    if DetResFn == 'rectangular':
        width = DetRes * np.sqrt(12)
        n = max(1, int(np.round(width / dy)))
        kernel = np.ones(n) / n
        I = np.convolve(I, kernel, mode="same")
        
    # Gaussian kernel half-width (~±4σ)
    elif DetResFn == 'gaussian':
        half_width = int(np.ceil(4 * DetRes / dy))
        xk = np.arange(-half_width, half_width + 1) * dy
        kernel = np.exp(-0.5 * (xk / DetRes)**2)
        kernel /= kernel.sum()
        I = np.convolve(I, kernel, mode="same")
    
    elif DetResFn not in ["none", None]:
         raise ValueError("DetResFn must be 'rectangular', 'gaussian', or 'none'")
    
    # normalize
    if I.max() <= 0:
        raise ValueError("Maximum I has calculated negative")
    else:
        I = I/I.max()
        return I
    
def intersect(m1, b1, m2, b2):
    # find intersection of two lines
    x = (b2 - b1) / (m1 - m2)
    y = m1 * x + b1
    return x, y

# TODO: Link these resolution functions to the lr_reduction ones...
def dTheta_Sigma(Si, S1, dS1Si):
    
    # Trapezoidal half-width angles (deg)
    HW_bot = np.degrees(np.arctan((S1 + Si) / (2 * dS1Si)))     # full half-width
    HW_top = np.degrees(np.arctan((S1 - Si) / (2 * dS1Si)))     # flat-top half-width
    
    a=-HW_bot
    b=-HW_top
    c=HW_top
    d=HW_bot
    
    #analytical equation for the sigma of a trapezoidal function
    Sigma=np.sqrt(((d-a)**2+(c-b)**2+(d-a)*(c-b))/24)
     
    return Sigma

def dLambda_Sigma(Lambda):
    #exponential fit to dLAM
    L = 0.07564423
    A = 0.13093263
    k = 0.34918918
    
    dLAM = L - A*np.exp(-k*Lambda)

    return dLAM

# TODO: Link this into the gravity_correction in lr_reduction...
def gravity_correct(LAM, ThetaIn, dSamp, dSlit):
    
    dSamp=dSamp/1000    #dSamp is m from sample to incident slit
    dSlit=dSlit/1000    #dSlit is m between slits
    
    #calculation from the ILL paper. this works for inclined beams.
    g=9.8067                #m/s^2
    h=6.6260715e-34         #Js=kg m^2/s
    mn=1.67492749804e-27    #kg
    V=h/(mn*LAM*1e-10)
    K=g/(2*V**2)
       
    #define the sample position as x=0, y=0. increasing x is towards moderator
    xs=0
    
    #positions of slits
    x1=dSamp
    x2=(dSamp+dSlit)
    
    #height of slits determined by incident theta, y=0 is the sample height
    y1=x1*np.tan(ThetaIn*np.pi/180)
    y2=x2*np.tan(ThetaIn*np.pi/180)
    
    #this is the location of the top of the parabola
    x0=(y1-y2+K*(x1**2-x2**2))/(2*K*(x1-x2))
    y0=y2+K*(x2-x0)**2
    
    xs=x0-np.sqrt(y0/K)
    
    ThetaSam=np.arctan(2*K*(x0-xs))*180/np.pi                 #angle is arctan(dy/dx) at sample
    
    dTheta=ThetaSam-ThetaIn
    
    return dTheta


def load_db_file(db_dir, db_name, comment="#"):
    """
    Load header metadata (key=value in commented lines) and numeric data (3 columns)
    from a single text file. Here applied for loading DB.

    Parameters
    ----------
    db_dir : str or Path
        Directory containing the file.
    db_name : str
        File name.
    comment : str
        Comment character used in header lines.

    Returns
    -------
    lDB, iDB, eDB : np.ndarray
        Numeric arrays unpacked from the file (assumes 3 columns).
    meta : dict[str, str]
        Parsed header metadata as strings.
    """
    path = Path(db_dir) / db_name
    meta = {}

    with path.open("r", encoding="utf-8") as f:
        # Read header lines once; stop at the first non-comment line.
        while True:
            pos = f.tell()
            line = f.readline()
            if not line:
                # End of file (no data section present)
                break

            s = line.strip()
            if not s.startswith(comment):
                # Rewind to the first data line so np.loadtxt can read it
                f.seek(pos)
                break

            # Parse key=value pairs from header comment lines
            s = s.lstrip(comment).strip()
            if s and "=" in s:
                key, value = s.split("=", 1)
                meta[key.strip()] = value.strip()

        # TODO: Add error handling if doesn't have the 3 columns.
        # Load numeric data (expects 3 columns to unpack)
        col1, col2, col3 = np.loadtxt(f, unpack=True)

    return col1, col2, col3, meta

def get_lam_range(chopper_lam, chopper_speed, scaled_width=3.5):
    """
    Determine lambda range from chopper settings

    Parameters
    ----------
    chopper_lam
        mid-point of lambda from log value (LambdaRequest)
    chopper_speed
        speed of chopper
    scaled_width
        equivalent bandwidth to 60Hz to use (default 3.5A)
        # TODO: Check this range. It is wider that the lr_reduction allowed before.

    Returns
    -------
      list
        [min, max] wavelength range
    """

    # Cut the edges by using a width of 2.6 A
    wl_min = chopper_lam - (scaled_width / 2) * 60.0 / chopper_speed
    wl_max = chopper_lam + (scaled_width / 2) * 60.0 / chopper_speed

    return [wl_min, wl_max]