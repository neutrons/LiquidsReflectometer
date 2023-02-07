import numpy as np
import warnings
warnings.filterwarnings('ignore', module='numpy')
warnings.filterwarnings('ignore')

from lmfit.models import GaussianModel, LinearModel, ConstantModel, RectangleModel, QuadraticModel

import mantid
import mantid.simpleapi as api


def process_data(workspace, summed=True, tof_step=200):
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


def fit_signal_flat_bck(x, y, x_min=110, x_max=170, center=None, sigma=None):
    gauss = GaussianModel(prefix='g_')
    linear = LinearModel(prefix='l_')
    quadratic = QuadraticModel(prefix='q_')
    rectangular = RectangleModel(prefix='r_')

    amplitude_guess = np.max(y[x_min:x_max])

    _center = 140
    _sigma = 2
    if center is not None:
        _center = center
    if sigma is not None:
        _sigma = sigma
        
    pars = gauss.make_params(amplitude=amplitude_guess, center=_center, sigma=_sigma)
    pars.update(linear.make_params(a=0, b=0))

    #if center is not None:
    #    pars['g_center'].vary=False
    if sigma is not None:
        pars['g_sigma'].vary=False
    pars['g_amplitude'].min=0
    pars['g_center'].min=_center-2
    pars['g_center'].max=_center+2
    
    weights=1/np.sqrt(y)
    weights[y<1]=1
    
    model = gauss + linear
    fit = model.fit(y[x_min:x_max], pars, method='leastsq',
                    x=x[x_min:x_max], 
                    weights=1/weights[x_min:x_max])
    #print(fit.fit_report())

    a=fit.params['g_amplitude']
    c=fit.params['g_center']
    width=fit.params['g_sigma']
    #print("Gaussian: \t %5.4g +- %5.4g \t %3.3g +- %3.3g \t %3.3g +- %3.3g" % (a.value, a.stderr, c.value, c.stderr, width.value, width.stderr))
    return c, width, fit