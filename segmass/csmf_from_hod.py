"""
"""
import numpy as np
from scipy.special import erf, erfinv
from scipy.stats import norm
from halotools.empirical_models import Leauthaud11Cens


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


class StellarMassLeauthaud11HOD(Leauthaud11Cens):

    def __init__(self, gal_type, prim_haloprop_key='halo_mpeak', **kwargs):
        Leauthaud11Cens.__init__(self, prim_haloprop_key=prim_haloprop_key, **kwargs)

        self.gal_type = gal_type

        # The _mock_generation_calling_sequence determines which methods
        # will be called during mock population, as well as in what order they will be called
        self._mock_generation_calling_sequence = ['assign_stellar_mass']
        self._galprop_dtypes_to_allocate = np.dtype([('stellar_mass', 'f8')])

    def assign_stellar_mass(self, **kwargs):
        if 'table' in kwargs.keys():
            table = kwargs['table']
            prim_haloprop = table[self.prim_haloprop_key]
        else:
            try:
                prim_haloprop = kwargs['prim_haloprop']
            except KeyError:
                raise KeyError("Input ``assign_stellar_mass`` accepts either \n"
                    "``table`` or ``prim_haloprop`` keyword arguments")

        mean_logsm = np.log10(self.mean_stellar_mass(prim_haloprop=prim_haloprop))
        scatter = self.param_dict[u'scatter_model_param1']
        dlogsm = (mean_logsm - self.threshold)
        dlogsm_zscore = dlogsm/scatter
        dlogsm_percentile = _percentile_from_z_score(dlogsm_zscore)
        u = 1 - np.random.uniform(0, dlogsm_percentile)
        mc_logsm = norm.isf(1 - u, loc=mean_logsm, scale=scatter)

        if 'table' in kwargs.keys():
            table['stellar_mass'][:] = 10**mc_logsm
        else:
            return 10**mc_logsm





