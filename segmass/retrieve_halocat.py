"""
"""
from halotools.sim_manager import CachedHaloCatalog, UserSuppliedHaloCatalog
from halotools.utils import crossmatch


__all__ = ('retrieve_halocat', )


def retrieve_halocat(subhalo_mpeak_cut, halo_mpeak_host_halo_cut,
        simname='bolplanck', redshift=0, **kwargs):

    halocat = CachedHaloCatalog(simname=simname, redshift=redshift, **kwargs)
    halocat.halo_table['halo_mpeak_host_halo'] = halocat.halo_table['halo_mpeak']

    idxA, idxB = crossmatch(halocat.halo_table['halo_hostid'], halocat.halo_table['halo_id'])
    halocat.halo_table['halo_mpeak_host_halo'][idxA] = halocat.halo_table['halo_mpeak'][idxB]

    mask = halocat.halo_table['halo_mpeak'] > subhalo_mpeak_cut
    mask *= halocat.halo_table['halo_mpeak_host_halo'] > halo_mpeak_host_halo_cut

    halo_catalog_columns = {key: halocat.halo_table[key][mask] for key in halocat.halo_table.keys()}

    metadata = dict(Lbox=halocat.Lbox,
                   particle_mass=halocat.particle_mass,
                   redshift=halocat.redshift)
    halo_catalog_columns.update(metadata)

    return UserSuppliedHaloCatalog(**halo_catalog_columns)
