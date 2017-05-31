"""
"""
import numpy as np
from halotools.empirical_models import BiasedNFWPhaseSpace


def gal_conc_bias(mu, bias_low_mu, bias_high_mu, low_mu, high_mu):
    return np.interp(np.log10(mu), [np.log10(low_mu), np.log10(high_mu)],
                     [bias_low_mu, bias_high_mu])


def galaxy_concentration(mpeak, mhost, conc_host, bias_low_mu, bias_high_mu):
    mu = mpeak/mhost
    conc_bias = gal_conc_bias(mu, bias_low_mu, bias_high_mu)
    gal_conc = conc_bias*conc_host
    return gal_conc


class NFWMassSegregation(BiasedNFWPhaseSpace):
    """
    """
    def __init__(self, *args, **kwargs):
        """
        """
        BiasedNFWPhaseSpace.__init__(self, *args, **kwargs)

        self.param_dict['bias_low_mu'] = 2.5
        self.param_dict['bias_high_mu'] = 0.5

    def calculate_conc_gal_bias(self, seed=None, **kwargs):
        r""" Calculate the ratio of the galaxy concentration to the halo concentration,
        :math:`c_{\rm gal}/c_{\rm halo}`.

        Parameters
        ----------
        prim_haloprop : array, optional
            Array storing the mass-like variable, e.g., ``halo_mvir``.

            If ``prim_haloprop`` is not passed,
            then ``table`` keyword argument must be passed.

        table : object, optional
            `~astropy.table.Table` storing the halo catalog.

            If your NFW model is based on the virial definition,
            then ``halo_mvir`` must appear in the input table,
            and likewise for other halo mass definitions.

            If ``table`` is not passed,
            then ``prim_haloprop`` keyword argument must be passed.

        Returns
        -------
        conc_gal_bias : array_like
            Ratio of the galaxy concentration to the halo concentration,
            :math:`c_{\rm gal}/c_{\rm halo}`.
        """
        if 'table' in list(kwargs.keys()):
            table = kwargs['table']
            subhalo_mpeak = table['halo_mpeak']
            host_halomass = table['halo_mvir_host_halo']
        elif 'subhalo_segprop' in list(kwargs.keys()):
            subhalo_mpeak = np.atleast_1d(kwargs['subhalo_mpeak']).astype('f4')
            host_halomass = np.atleast_1d(kwargs['host_halomass']).astype('f4')
        else:
            msg = ("\nYou must pass either a ``table``, \n"
                "or ``subhalo_mpeak`` and ``host_halomass`` arguments, \n"
                "to the ``assign_conc_gal_bias`` function of the ``NFWMassSegregation`` class.\n")
            raise KeyError(msg)

        mu = subhalo_mpeak/host_halomass
        result = self.gal_bias_factor(mu)

        #  Now that conc_gal_bias has been calculated, write it to the mock galaxy_table
        if 'table' in list(kwargs.keys()):
            table['conc_gal_bias'][:] = result
            halo_conc = table['conc_NFWmodel']
            gal_conc = self._clipped_galaxy_concentration(halo_conc, result)
            table['conc_galaxy'][:] = gal_conc
        else:
            return result

    def gal_bias_factor(self, mu):
        """
        """
        return gal_conc_bias(mu, self.param_dict['bias_low_mu'],
                self.param_dict['bias_high_mu'], self.low_mu, self.high_mu)




