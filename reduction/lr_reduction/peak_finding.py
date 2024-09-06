import warnings

import numpy as np

warnings.filterwarnings('ignore', module='numpy')
warnings.filterwarnings('ignore')

import mantid.simpleapi as api
from lmfit.models import GaussianModel


def process_data(workspace, summed=True, tof_step=200):
    r"""
        Process a Mantid workspace to extract counts vs pixel.
        :param workspace: Mantid workspace
        :param summed: if True, the x pixels with be summed
        :param tof_step: TOF bin size
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
    r"""
        Fit a Gaussian peak.
        :param x: list of x values
        :param y: list of y values
        :param x_min: start index of the list of points
        :param x_max: end index of the list of points
        :param center: estimated center position
        :param sigma: if provided, the sigma will be fixed to the given value
        :param background: if provided, the value will be subtracted from y
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
