"""
"""
import numpy as np
from scipy.special import erf, erfinv
from scipy.stats import norm


def _percentile_from_z_score(z_score):
    return 0.5*(1 + erf(z_score/np.sqrt(2)))


def _z_score_from_percentile(percentile):
    return np.sqrt(2)*erfinv(2*percentile-1)


def kernel(model):
    mean_sm = model.mean_stellar_mass_centrals(prim_haloprop=model.mock.galaxy_table['halo_mpeak'])
    mean_logsm = np.log10(mean_sm)
    dlogsm = (mean_logsm - model.threshold)
    dlogsm_zscore = dlogsm/model.param_dict[u'scatter_model_param1']
    dlogsm_percentile = _percentile_from_z_score(dlogsm_zscore)

    u = 1 - np.random.uniform(0, dlogsm_percentile)
    mc_logsm = norm.isf(1 - u, loc=mean_logsm,
                  scale=model.param_dict[u'scatter_model_param1'])
    return mc_logsm
