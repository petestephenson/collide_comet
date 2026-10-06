"""Derived physical quantities.  Ports of the ``asymPlots.calc*`` functions.

Conventions
-----------
* Spatial arrays are indexed ``[ix, iy, iz]``; the energy-resolved electron
  distribution is ``plas_dens_E[ix, iy, iz, iE, i_gen]``.
* ``gens`` selects test-particle generations with **0-based** indices:
  MATLAB ``inc_gens = [1, 2]`` is ``gens=(0, 1)``.
* Units follow the MATLAB code: lengths in km, densities in cm^-3, energies
  in eV, outgassing rates in s^-1.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
from scipy import ndimage

from . import constants as c
from .io import CrossSections, Grids, Outgassing

__all__ = [
    "interp1",
    "local_outgassing_from_cos",
    "local_outgassing",
    "neutral_density",
    "exobase_radius",
    "box_smooth",
    "plasma_densities",
    "smooth_plas_dens_E",
    "eii_frequency",
    "eii_production",
    "cold_electron_production",
    "temperature",
    "combine_run_velocities",
    "pressure",
    "ambipolar_field",
]


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def interp1(x: np.ndarray, y: np.ndarray, xq: np.ndarray) -> np.ndarray:
    """Linear interpolation returning NaN outside ``x`` (MATLAB ``interp1``)."""
    return np.interp(xq, x, y, left=np.nan, right=np.nan)


_trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy 2 renamed trapz


def _select_gens(plas_dens_E: np.ndarray, gens: Sequence[int]) -> np.ndarray:
    gens = np.atleast_1d(np.asarray(gens, dtype=int))
    if gens.max() >= plas_dens_E.shape[-1]:
        raise IndexError(
            f"gens={tuple(gens)} requested but the array only has "
            f"{plas_dens_E.shape[-1]} generation(s)"
        )
    return plas_dens_E[..., gens]


# ---------------------------------------------------------------------------
# neutral coma
# ---------------------------------------------------------------------------
def local_outgassing_from_cos(cos_to_max, q_uniform: float, q_dayside: float) -> np.ndarray:
    """Local outgassing rate for a given cosine of the angle to ``q_vec``.

    ``Q = q_uniform + 4 q_dayside max(cos, 0)``.
    Port of ``calc_local_q_from_angle``.
    """
    cos_to_max = np.asarray(cos_to_max, dtype=float)
    return q_uniform + 4.0 * q_dayside * np.clip(cos_to_max, 0.0, None)


def local_outgassing(x, y, z, outgassing: Outgassing) -> Tuple[np.ndarray, np.ndarray]:
    """Local outgassing rate at positions ``(x, y, z)``.

    Port of ``calc_local_q_from_pos``.  Returns ``(Q, cos_to_max)``.
    As in the MATLAB code ``q_vec`` is used as stored (it is not normalised).
    """
    x, y, z = np.broadcast_arrays(*(np.asarray(a, dtype=float) for a in (x, y, z)))
    r = np.sqrt(x**2 + y**2 + z**2)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_to_max = (
            x * outgassing.q_vec[0] + y * outgassing.q_vec[1] + z * outgassing.q_vec[2]
        ) / r
    q = local_outgassing_from_cos(cos_to_max, outgassing.q_uniform, outgassing.q_dayside)
    return q, cos_to_max


def neutral_density(x, y, z, outgassing: Outgassing, u_gas: float = c.U_GAS) -> np.ndarray:
    """Water density of the (asymmetric) Haser-like coma [cm^-3].

    Port of ``asymPlots.calcNeutralDens``:
    ``n = Q_local / (4 pi u_gas r^2 1e15)`` with ``r`` in km and ``u_gas`` in
    km/s.  ``x, y, z`` broadcast against each other.
    """
    x, y, z = (np.asarray(a, dtype=float) for a in (x, y, z))
    r2 = x**2 + y**2 + z**2
    q, _ = local_outgassing(x, y, z, outgassing)
    with np.errstate(divide="ignore", invalid="ignore"):
        return q / (4.0 * np.pi * u_gas * r2 * 1e15)


def neutral_density_on_grid(grids: Grids, outgassing: Outgassing, u_gas: float = c.U_GAS):
    """Neutral density at the cell centres of ``grids``, shape ``(nx, ny, nz)``."""
    xc, yc, zc = grids.centres
    return neutral_density(
        xc[:, None, None], yc[None, :, None], zc[None, None, :], outgassing, u_gas
    )


def exobase_radius(q_local, sigma: float = c.SIGMA_EXOBASE, u_gas: float = c.U_GAS):
    """Electron exobase distance [km] for a local outgassing rate [s^-1].

    ``r_exo = 100 Q sigma 1e-10 / (4 pi u_gas)`` -- the expression used in the
    MATLAB plotting functions (``sigma`` in cm^2, ``u_gas`` in km/s).
    """
    return 100.0 * np.asarray(q_local, dtype=float) * sigma * 1e-10 / (4.0 * np.pi * u_gas)


# ---------------------------------------------------------------------------
# electron densities
# ---------------------------------------------------------------------------
def box_smooth(arr: np.ndarray) -> np.ndarray:
    """3x3x3 box average over the first three axes with edge replication.

    Equivalent to the explicit triple loop in ``calculatePlasmaDensities`` /
    ``smooth_plas_dens_E``.
    """
    size = (3, 3, 3) + (1,) * (arr.ndim - 3)
    return ndimage.uniform_filter(np.asarray(arr, dtype=float), size=size, mode="nearest")


def energy_band_indices(energy_edges: np.ndarray, min_energy: float, max_energy: float) -> slice:
    """Energy bins covered by ``[min_energy, max_energy]`` as a slice.

    Reproduces the MATLAB selection
    ``find(E > minEng, 1) - 1 : find(E <= maxEng, 1, 'last') - 1``.
    """
    above = np.flatnonzero(energy_edges > min_energy)
    below = np.flatnonzero(energy_edges <= max_energy)
    if above.size == 0 or below.size == 0:
        return slice(0, 0)
    first = max(above[0] - 1, 0)
    last = below[-1] - 1  # inclusive
    return slice(first, last + 1)


def plasma_densities(
    grids: Grids,
    plas_dens_E: np.ndarray,
    gens: Sequence[int] = c.ALL_GENS,
    min_energy: float = 0.0,
    max_energy: float = 1.0,
    smooth: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """Total and energy-band electron densities [cm^-3].

    Port of ``asymPlots.calculatePlasmaDensities``.

    Returns
    -------
    plas_dens : (nx, ny, nz) density integrated over all energies
    plas_dens_band : (nx, ny, nz) density between ``min_energy`` and
        ``max_energy`` (the "cold" population for the default 0-1 eV).

    Both are 3x3x3 box-smoothed unless ``smooth=False``.
    """
    dE = grids.dE
    f = _select_gens(plas_dens_E, gens)
    band = energy_band_indices(grids.energy, min_energy, max_energy)

    dens = np.einsum("xyzeg,e->xyz", f, dE)
    dens_band = np.einsum("xyzeg,e->xyz", f[:, :, :, band, :], dE[band])
    if smooth:
        dens, dens_band = box_smooth(dens), box_smooth(dens_band)
    return dens, dens_band


def smooth_plas_dens_E(plas_dens_E: np.ndarray) -> np.ndarray:
    """Spatially box-smooth an energy distribution.  Port of ``smooth_plas_dens_E``."""
    return box_smooth(plas_dens_E)


# ---------------------------------------------------------------------------
# electron-impact rates
# ---------------------------------------------------------------------------
def eii_frequency(
    plas_dens_E: np.ndarray,
    energy_edges: np.ndarray,
    xs: CrossSections,
    gens: Sequence[int] = c.ALL_GENS,
    ix_ioni: int = c.IX_IONI,
) -> np.ndarray:
    """Electron-impact ionisation frequency [s^-1], shape ``(nx, ny, nz)``.

    Port of ``asymPlots.calcEIoniFreq``:
    ``nu = int sigma_ioni(E) f(E) v(E) dE`` (trapezoidal rule on the energy
    bin mid-points).
    """
    e_mid = 0.5 * (energy_edges[:-1] + energy_edges[1:])
    sigma = interp1(xs.energy, xs.xs[:, ix_ioni:].sum(axis=1), e_mid)
    sigma = np.where(np.isnan(sigma), 0.0, sigma)
    speed = np.sqrt(e_mid * 2.0 * c.QE / c.ME) * 1e2  # cm/s
    flux = _select_gens(plas_dens_E, gens).sum(axis=-1) * speed
    return _trapz(sigma * flux, e_mid, axis=3)


def eii_production(
    plas_dens_E: np.ndarray,
    grids: Grids,
    outgassing: Outgassing,
    xs: CrossSections,
    gens: Sequence[int] = c.ALL_GENS,
) -> np.ndarray:
    """Electron-impact ionisation production rate [cm^-3 s^-1].

    Port of ``asymPlots.calculateEIIProduction`` (``n_H2O * nu_ioni``).
    """
    n_neutral = neutral_density_on_grid(grids, outgassing)
    return n_neutral * eii_frequency(plas_dens_E, grids.energy, xs, gens)


def cold_electron_cross_section(
    energy_edges: np.ndarray,
    xs: CrossSections,
    e_cold: float = 1.0,
    ix_ioni: int = c.IX_IONI,
) -> Tuple[np.ndarray, np.ndarray]:
    """Effective cross section for producing an electron below ``e_cold``.

    For each inelastic process an electron ends up cold if its energy lies in
    ``[max(E_th, E_cold), E_th + E_cold]`` (``+ E_cold`` again for ionisation,
    where the energy is shared).  The partially covered first and last energy
    bins are weighted by the covered fraction.

    Returns ``(xs_total, energy_nodes)``: the summed cross section per energy
    bin [cm^2] and the energies at which the last process was evaluated (see
    ``cold_electron_production`` for why the latter is returned).
    """
    eg = np.asarray(energy_edges, dtype=float)
    n_eng = eg.size - 1
    e_mid = 0.5 * (eg[1:] + eg[:-1])
    xs_total = np.zeros(n_eng)
    nodes = e_mid.copy()

    for i_xs, e_th in enumerate(xs.e_th):
        if not e_th > 0:
            continue
        lim_lo = max(e_th, e_cold)
        lim_hi = e_cold + e_th
        if i_xs >= ix_ioni:
            lim_hi += e_cold

        above_lo = np.flatnonzero(eg > lim_lo)
        above_hi = np.flatnonzero(eg > lim_hi)
        if above_lo.size == 0 or above_hi.size == 0:
            raise ValueError(
                f"Energy grid (max {eg[-1]:g} eV) does not cover the cold-electron "
                f"window [{lim_lo:g}, {lim_hi:g}] eV of process {i_xs}."
            )
        j_lo = above_lo[0] - 1  # bin containing the lower limit
        j_hi = above_hi[0] - 1  # bin containing the upper limit

        eng_low = 0.5 * (lim_lo + eg[j_lo + 1])
        eng_upp = 0.5 * (lim_hi + eg[j_hi])
        low_factor = (eg[j_lo + 1] - lim_lo) / (eg[j_lo + 1] - eg[j_lo])
        upp_factor = (lim_hi - eg[j_hi]) / (eg[j_hi + 1] - eg[j_hi])

        nodes = e_mid.copy()
        nodes[j_lo] = eng_low
        nodes[j_hi] = eng_upp  # wins if both limits share a bin (as in MATLAB)
        xs_interp = interp1(xs.energy, xs.xs[:, i_xs], nodes)
        if j_hi == j_lo:
            # MATLAB: x([i i]) = x([i i]) .* [lo up]  ->  last assignment wins
            xs_interp[j_hi] *= upp_factor
        else:
            xs_interp[j_lo] *= low_factor
            xs_interp[j_hi] *= upp_factor
        xs_interp[:j_lo] = 0.0
        xs_interp[j_hi + 1 :] = 0.0
        xs_total += xs_interp
    return xs_total, nodes


def cold_electron_production(
    plas_dens_E: np.ndarray,
    grids: Grids,
    outgassing: Outgassing,
    xs: CrossSections,
    e_cold: float = 1.0,
    gens: Sequence[int] = c.ALL_GENS,
    legacy_speed_grid: bool = True,
) -> np.ndarray:
    """Production rate of cold (< ``e_cold`` eV) electrons [cm^-3 s^-1].

    Port of ``asymPlots.calculateColdeProduction``.

    Parameters
    ----------
    legacy_speed_grid : bool
        The MATLAB function evaluates the electron speed on the energy nodes
        left over from the *last* collision process in its loop (bin
        mid-points with two entries shifted to that process's window edges)
        rather than on the plain bin mid-points.  ``True`` (default) keeps
        that behaviour so results match MATLAB exactly; ``False`` uses the
        bin mid-points.  The difference is confined to two energy bins.
    """
    xs_total, nodes = cold_electron_cross_section(grids.energy, xs, e_cold)
    e_for_speed = nodes if legacy_speed_grid else grids.energy_mid
    speed = np.sqrt(2.0 * e_for_speed * abs(c.QE) / c.ME) * 1e2  # cm/s
    weight = xs_total * grids.dE * speed
    per_molecule = np.einsum("xyzeg,e->xyz", _select_gens(plas_dens_E, gens), weight)
    return per_molecule * neutral_density_on_grid(grids, outgassing)


# ---------------------------------------------------------------------------
# temperatures, pressure, fields
# ---------------------------------------------------------------------------
def temperature(
    square_vel: np.ndarray, bulk_vel: np.ndarray, me: float, component: Optional[int] = None
) -> np.ndarray:
    """Electron temperature [eV] from velocity moments.

    Port of ``asymPlots.calculateTemperature``.

    Parameters
    ----------
    square_vel, bulk_vel : (nx, ny, nz, 3)
        ``<v_i^2>`` and ``<v_i>`` for one population.
    me : float
        Electron mass used by the run [kg].
    component : {None, 0, 1, 2}
        ``None`` gives the scalar temperature ``Tr(T)/3`` (MATLAB
        ``i_dim = 0``); 0/1/2 give the x/y/z component.
    """
    if component is None:
        sq = np.nansum(square_vel, axis=3)
        bulk_sq = np.nansum(bulk_vel**2, axis=3)
        return me * (sq - bulk_sq) / c.QE / 3.0
    return me * (square_vel[..., component] - bulk_vel[..., component] ** 2) / c.QE


def combine_run_velocities(
    pe_bulk_vel: np.ndarray,
    pe_square_vel: np.ndarray,
    plas_dens_E_pe: np.ndarray,
    sw_bulk_vel: np.ndarray,
    sw_square_vel: np.ndarray,
    plas_dens_E_sw: np.ndarray,
    dE: np.ndarray,
    gens: Sequence[int] = c.ALL_GENS,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Density-weighted velocity moments of the combined PE + SW population.

    Port of ``asymPlots.combine_run_velocities``.

    Returns
    -------
    tot_plas_dens : (nx, ny, nz, 1)
    ave_bulk_vel, ave_square_vel : (nx, ny, nz, 3)
    """
    gens = np.atleast_1d(np.asarray(gens, dtype=int))

    def moments(bulk, square, f):
        dens_gen = np.einsum("xyzeg,e->xyzg", f, dE)[..., gens]  # (x,y,z,g)
        w = dens_gen[:, :, :, np.newaxis, :]
        sq = np.nansum(square[..., gens] * w, axis=-1)
        bk = np.nansum(bulk[..., gens] * w, axis=-1)
        return np.nansum(dens_gen, axis=-1), bk, sq

    dens_pe, bulk_pe, sq_pe = moments(pe_bulk_vel, pe_square_vel, plas_dens_E_pe)
    dens_sw, bulk_sw, sq_sw = moments(sw_bulk_vel, sw_square_vel, plas_dens_E_sw)

    tot_dens = (dens_pe + dens_sw)[..., np.newaxis]
    with np.errstate(divide="ignore", invalid="ignore"):
        ave_square = (sq_pe + sq_sw) / tot_dens
        ave_bulk = (bulk_pe + bulk_sw) / tot_dens
    return tot_dens, ave_bulk, ave_square


def pressure(plas_dens: np.ndarray, temperature_eV: np.ndarray) -> np.ndarray:
    """Electron pressure ``n k T`` [Pa] from density [cm^-3] and T [eV].

    Note: MATLAB ``calculatePressure`` additionally multiplies by ``k_B`` (and
    ``calcAmbipolarfield`` divides it out again).  This returns the pressure
    in pascals; :func:`ambipolar_field` is adjusted so the field is unchanged.
    """
    return np.asarray(plas_dens) * 1e6 * np.asarray(temperature_eV) * c.QE_PRESSURE


def ambipolar_field(pressure_pa: np.ndarray, plas_dens: np.ndarray, grids: Grids) -> np.ndarray:
    """Ambipolar electric field ``-grad(P_e) / (n e)`` [V/m].

    Port of ``asymPlots.calcAmbipolarfield``.  ``pressure_pa`` in Pa (see
    :func:`pressure`), ``plas_dens`` in cm^-3; gradients are taken on the cell
    centres (km, converted to m).

    Returns an array of shape ``(3, nx, ny, nz)``.
    """
    xc, yc, zc = grids.centres
    grads = np.gradient(np.asarray(pressure_pa, dtype=float), xc, yc, zc)
    with np.errstate(divide="ignore", invalid="ignore"):
        field = np.stack([-g / (1e6 * plas_dens * c.QE_PRESSURE) for g in grads])
    return field / 1e3
