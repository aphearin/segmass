"""
"""
from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import PrebuiltSubhaloModelFactory


def load_behroozi10_mock(simname='bolplanck', redshift=0, logsm_threshold=10):
    halocat = CachedHaloCatalog(simname=simname, redshift=redshift)
    model = PrebuiltSubhaloModelFactory('behroozi10', redshift=redshift)
    model.populate_mock(halocat)

    mask = model.mock.galaxy_table['stellar_mass'] > 10**logsm_threshold
    return model.mock.galaxy_table[mask]
