import warnings

import numpy as np

warnings.filterwarnings('ignore', module='numpy')
warnings.filterwarnings('ignore')

import mantid.simpleapi as api
from lmfit.models import GaussianModel


def process_data(workspace, summed=True, tof_step=200):
    """
        Process a Mantid workspace to extract counts vs pixel.

        Parameters
        ----------
        workspace : Mantid workspace
            The Mantid workspace to process.
        summed : bool, optional
            If True, the x pixels will be summed (default is True).
        tof_step : int, optional
            The TOF bin size (default is 200).

        Returns
        -------
        tuple
            A tuple containing:
            - tof : numpy.ndarray
                The time-of-flight values.
            - _x : numpy.ndarray
                The pixel indices.
            - _y : numpy.ndarray
                The summed counts for each pixel.
    """
    tof_min = workspace.getTofMin()
    tof_max = workspace.getTofMax()
    _ws = api.Rebin(InputWorkspace=workspace, Params="%s,%s,%s" % (tof_min, tof_step, tof_max))
    y=_ws.extractY()
    y = np.reshape(y, (256, 304, y.shape[1]))

    tof=_ws.extractX()[0]
    tof = (tof[:-1]+tof[1:])/2.0

    if summed:
        y = np.sum(y, axis=2)

    _y = np.sum(y, axis=0)
    _x = np.arange(304)
    return tof, _x, _y


def fit_signal_flat_bck(x, y, x_min=110, x_max=170, center=None, sigma=None,
                        background=None):
    """
    Fit a Gaussian peak.

    Parameters
    ----------
    x : list
        List of x values.
    y : list
        List of y values.
    x_min : int, optional
        Start index of the list of points, by default 110.
    x_max : int, optional
        End index of the list of points, by default 170.
    center : float, optional
        Estimated center position, by default None.
    sigma : float, optional
        If provided, the sigma will be fixed to the given value, by default None.
    background : float, optional
        If provided, the value will be subtracted from y, by default None.

    Returns
    -------
    c : float
        Fitted center position of the Gaussian peak.
    width : float
        Fitted width (sigma) of the Gaussian peak.
    fit : lmfit.model.ModelResult
        The result of the fit.
    """
    gauss = GaussianModel(prefix='g_')

    amplitude_guess = np.max(y[x_min:x_max])
    if background is None:
        background = np.min(y[x_min:x_max])

    _center = 140
    _sigma = 2
    if center is not None:
        _center = center
    if sigma is not None:
        _sigma = sigma

    pars = gauss.make_params(amplitude=amplitude_guess, center=_center, sigma=_sigma)

    if sigma is not None:
        pars['g_sigma'].vary=False
    pars['g_amplitude'].min=0
    pars['g_center'].min=_center-2
    pars['g_center'].max=_center+2

    weights=1/np.sqrt(y)
    weights[y<1]=1

    fit = gauss.fit(y[x_min:x_max]-background, pars, method='leastsq',
                    x=x[x_min:x_max],
                    weights=1/weights[x_min:x_max])

    fit.params['g_amplitude']
    c=fit.params['g_center']
    width=fit.params['g_sigma']

    return c, width, fit
