"""
"""
import numpy as np
from halotools.empirical_models.phase_space_models.analytic_models.satellites.nfw.kernels.mc_generate_nfw_radial_positions import mc_generate_nfw_radial_positions
from halotools.empirical_models.phase_space_models.analytic_models.satellites.nfw.kernels.biased_isotropic_velocity import dimensionless_radial_velocity_dispersion


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


def assign_posvel(halo_conc, gal_conc, conc_bin_edges):
    """
    Examples
    --------
    >>> num_gals = int(1e4)
    >>> cmin, cmax = 2, 25
    >>> halo_conc = np.random.uniform(cmin, cmax, num_gals)
    >>> gal_conc = np.random.uniform(0.5, 2, num_gals)*halo_conc
    >>> conc_bin_edges = np.linspace(cmin, cmax, 5)
    >>> r_by_rvir, vdisp = assign_posvel(halo_conc, gal_conc, conc_bin_edges)
    """

    r_by_rvir = np.zeros_like(halo_conc)
    vdisp = np.zeros_like(halo_conc)

    idx_gal_conc, conc_bin_midpoints = digitize_concentration(gal_conc, conc_bin_edges)
    idx_halo_conc, conc_bin_midpoints = digitize_concentration(halo_conc, conc_bin_edges)
    num_bin_midpoints = len(conc_bin_midpoints)

    #  Loop over galaxy concentrations
    for i in range(num_bin_midpoints):
        approx_gal_conc = conc_bin_midpoints[i]

        #  Identify all objects with gal_conc in bin i
        idx_gal_conc_mask = idx_gal_conc == i

        num_gals_ibin = np.count_nonzero(idx_gal_conc_mask)

        if num_gals_ibin > 0:
            #  For all objects with gal_conc in bin i, calculate r/Rvir
            r_by_rvir_ibin = mc_generate_nfw_radial_positions(
                    num_pts=num_gals_ibin, conc=approx_gal_conc, halo_radius=1.)
            r_by_rvir[idx_gal_conc_mask] = r_by_rvir_ibin

            vdisp_ibin = np.zeros(num_gals_ibin)
            #  Loop over halo concentrations
            for j in range(num_bin_midpoints):
                approx_halo_conc = conc_bin_midpoints[j]

                #  Identify all objects with gal_conc in bin i
                idx_halo_conc_mask = idx_halo_conc[idx_gal_conc_mask] == j

                num_gals_ijbin = np.count_nonzero(idx_halo_conc_mask)
                if num_gals_ijbin > 0:
                    r_by_rvir_ijbin = r_by_rvir_ibin[idx_halo_conc_mask]
                    vdisp_ijbin = dimensionless_radial_velocity_dispersion(
                        r_by_rvir_ijbin, approx_halo_conc, approx_gal_conc)
                    vdisp_ibin[idx_halo_conc_mask] = vdisp_ijbin

            vdisp[idx_gal_conc_mask] = vdisp_ibin

    return r_by_rvir, vdisp

