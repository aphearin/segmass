"""
"""
import numpy as np
from halotools.empirical_models import randomly_resort

__all__ = ('segregated_luminosities', )


def segregated_luminosities(gals, sigma, **kwargs):
    """
    """
    segregated_luminosity = np.copy(gals['luminosity'].data)

    mhost_bins = kwargs.get('mhost_bins', np.logspace(11.5, 15, 15))

    for low_mhost, high_mhost in zip(mhost_bins[:-1], mhost_bins[1:]):
        mask = (gals['halo_mvir_host_halo'] >= low_mhost)
        mask *= (gals['halo_mvir_host_halo'] < high_mhost)
        mask *= gals['gal_type'] == 'satellites'

        num_sample = np.count_nonzero(mask)
        if num_sample > 10:
            sample = gals[mask]
            sample_luminosity = np.copy(sample['luminosity'].data)
            idx_sorted_mpeak = np.argsort(sample['halo_mpeak'])
            idx_sorted_luminosity = np.argsort(sample['luminosity'])
            sorted_luminosity = sample_luminosity
            sorted_luminosity[idx_sorted_mpeak] = sample_luminosity[idx_sorted_luminosity]
            new_sample_luminosity = randomly_resort(sorted_luminosity, sigma)
            segregated_luminosity[mask] = new_sample_luminosity

    return segregated_luminosity
