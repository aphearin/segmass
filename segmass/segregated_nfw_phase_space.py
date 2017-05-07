"""
"""
from __future__ import division, print_function, absolute_import

import numpy as np

from .biased_nfw_phase_space import BiasedNFWPhaseSpace


__author__ = ('Andrew Hearin', )
__all__ = ('SegregatedNFWPhaseSpace', )


missing_quiescent_key_msg = ("The `SegregatedNFWPhaseSpace` class "
    "can only be used to make mocks \nin concert"
    "with some other component model that is responsible for \nmodeling an"
    "``quiescent`` property of the ``galaxy_table``.\n")


class SegregatedNFWPhaseSpace(BiasedNFWPhaseSpace):
    """
    """
    def __init__(self, **kwargs):
        BiasedNFWPhaseSpace.__init__(self, **kwargs)

        self._initialize_conc_bias_param_dict()

    def _initialize_conc_bias_param_dict(self, **kwargs):
        r""" Set up the appropriate number of keys in the parameter dictionary
        and give the keys standardized names.
        """
        self.param_dict['conc_bias_1'] = 1
        self.param_dict['conc_bias_2'] = 1

        self._conc_bias_abscissa_1 = 9
        self._conc_bias_abscissa_2 = 11.5

    def calculate_conc_gal_bias(self, seed=None, **kwargs):
        """
        """
        if 'table' in list(kwargs.keys()):
            table = kwargs['table']
            try:
                stellar_mass = table['stellar_mass']
            except KeyError:
                raise KeyError("Mocks populated with the SegregatedNFWPhaseSpace require \n"
                    "some other component model for ``stellar_mass``")
        else:
            try:
                stellar_mass = np.atleast_1d(kwargs['stellar_mass'])
            except:
                msg = ("\nYou must pass either a ``table`` or ``stellar_mass`` argument \n"
                    "to the ``calculate_conc_gal_bias`` function of the ``SegregatedNFWPhaseSpace`` class.\n")
                raise KeyError(msg)

        try:
            logsm = np.log10(stellar_mass)
        except:
            raise ValueError("Input ``stellar_mass`` must not be zero-valued")

        xp = (self._conc_bias_abscissa_1, self._conc_bias_abscissa_2)
        yp = (self.param_dict['conc_bias_1'], self.param_dict['conc_bias_2'])
        conc_gal_bias = np.interp(logsm, xp, yp)

        if 'table' in list(kwargs.keys()):
            table['conc_gal_bias'][:] = conc_gal_bias
            halo_conc = table['conc_NFWmodel']
            gal_conc = self._clipped_galaxy_concentration(halo_conc, conc_gal_bias)
            table['conc_galaxy'][:] = gal_conc
        else:
            return conc_gal_bias
