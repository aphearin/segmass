"""
"""
import numpy as np
from halotools.empirical_models.phase_space_models.analytic_models.satellites.nfw.kernels.mc_generate_nfw_radial_positions import mc_generate_nfw_radial_positions


def digitize_concentration(conc, conc_bin_edges):
    """
    """
    num_conc_bin_edges = len(conc_bin_edges)
    conc_bin_midpoints = 0.5*(conc_bin_edges[:-1] + conc_bin_edges[1:])

    iconc = np.digitize(conc, conc_bin_edges) - 1
    iconc = np.where(iconc >= num_conc_bin_edges-2, num_conc_bin_edges-2, iconc)
    iconc = np.where(iconc < 0, 0, iconc)
    return iconc, conc_bin_midpoints


def assign_r_by_rvir(galaxy_concentration, conc_bin_edges):
    """

    Examples
    --------
    >>> num_gals = int(1e3)
    >>> cmin, cmax = 1.5, 30.
    >>> galaxy_concentration = np.random.uniform(cmin, cmax, num_gals)
    >>> conc_bin_edges = np.linspace(2, 25, 25)
    >>> r_by_rvir = assign_r_by_rvir(galaxy_concentration, conc_bin_edges)
    >>> assert np.all(r_by_rvir >= 0)
    >>> assert np.all(r_by_rvir <= 1)
    """
    iconc, conc_bin_midpoints = digitize_concentration(galaxy_concentration, conc_bin_edges)

    r_by_rvir = np.zeros_like(galaxy_concentration)
    num_bin_midpoints = len(conc_bin_midpoints)

    for i in range(num_bin_midpoints):
        approx_conc = conc_bin_midpoints[i]
        iconc_mask = iconc == i
        num_gals_iconc = np.count_nonzero(iconc_mask)

        if num_gals_iconc > 0:
            r_by_rvir[iconc_mask] = mc_generate_nfw_radial_positions(
                num_pts=num_gals_iconc, conc=approx_conc, halo_radius=1.)

    return r_by_rvir
