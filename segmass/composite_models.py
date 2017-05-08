import numpy as np

from halotools.empirical_models import ZuMandelbaum15Cens, ZuMandelbaum15Sats
from halotools.empirical_models import TrivialPhaseSpace
from halotools.empirical_models import SubhaloPhaseSpace
from halotools.empirical_models import HodModelFactory

from csmf_from_hod import StellarMassFromHOD
from segregated_nfw_phase_space import SegregatedNFWPhaseSpace


def composite_model(threshold, host_haloprop_bins=np.logspace(11, 15.15, 8)):

    cens_hod = ZuMandelbaum15Cens(threshold=threshold, prim_haloprop_key='halo_mpeak', redshift=0)
    sats_hod = ZuMandelbaum15Sats(threshold=threshold, prim_haloprop_key='halo_mpeak', redshift=0)

    csmf_centrals = StellarMassFromHOD('centrals', cens_hod.threshold, cens_hod)
    csmf_satellites = StellarMassFromHOD('satellites', cens_hod.threshold, cens_hod)
    csmf_satellites._suppress_repeated_param_warning = True

    seg_nfw = SegregatedNFWPhaseSpace(prim_haloprop_key='halo_mpeak')

    d = dict(halo_id=('halo_id', 'i8'), halo_mpeak=('halo_mpeak', 'f8'))
    satellites_mpeak = SubhaloPhaseSpace('satellites', host_haloprop_bins,
            inherited_subhalo_props_dict=d, binning_key='halo_mpeak_host_halo')

    cens_profile = TrivialPhaseSpace()

    model_dictionary = dict(centrals_occupation=cens_hod,
                           centrals_profile=cens_profile,
                           centrals_stellar_mass=csmf_centrals,
                           satellites_occupation=sats_hod,
                           satellites_mpeak=satellites_mpeak,
                           satellites_stellar_mass=csmf_satellites,
                           satellites_profile=seg_nfw)

    return HodModelFactory(**model_dictionary)
