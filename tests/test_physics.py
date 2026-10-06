"""Check the vectorised ports against literal transliterations of the MATLAB code."""

import numpy as np
import pytest

import collide_comet as cc
from collide_comet import physics
from collide_comet.io import Grids, Outgassing

QE, ME = 1.6e-19, 9.11e-31


@pytest.fixture(scope="module")
def setup():
    rng = np.random.default_rng(3)
    grids = Grids(
        x=np.linspace(-60, 60, 7), y=np.linspace(-50, 50, 6), z=np.linspace(-40, 40, 5),
        energy=np.r_[np.arange(0, 10, 0.5), np.arange(10, 60, 5), np.arange(60, 201, 20)],
    )
    f = rng.random(grids.shape + (grids.n_energy, 2)) + 0.01
    og = Outgassing(1e27, 0.75e27, 0.25e27, np.array([0.0, 1.0, 0.0]))
    return grids, f, og, cc.load_cross_sections()


def matlab_plasma_densities(grids, f, inc_gens, min_eng, max_eng):
    """calculatePlasmaDensities.m, written with MATLAB's 1-based indices."""
    E = grids.energy
    dE = E[1:] - E[:-1]
    ix_min = np.flatnonzero(E > min_eng)[0] + 1 - 1  # find(...,'first') - 1
    ix_max = np.flatnonzero(E <= max_eng)[-1] + 1 - 1  # find(...,'last') - 1
    inc = np.arange(ix_min, ix_max + 1) - 1  # to 0-based
    g = np.array(inc_gens) - 1
    dens = (f[..., g] * dE[None, None, None, :, None]).sum(axis=(3, 4))
    cold = (f[:, :, :, inc][..., g] * dE[inc][None, None, None, :, None]).sum(axis=(3, 4))
    out = []
    for arr in (dens, cold):
        nx, ny, nz = arr.shape
        pad = arr[np.r_[0, 0:nx, nx - 1]][:, np.r_[0, 0:ny, ny - 1]][:, :, np.r_[0, 0:nz, nz - 1]]
        acc = np.zeros_like(arr)
        for ix in range(3):
            for iy in range(3):
                for iz in range(3):
                    acc += pad[ix : ix + nx, iy : iy + ny, iz : iz + nz]
        out.append(acc / 27)
    return out


def matlab_neutral_density(grids, og):
    xc, yc, zc = grids.centres
    n = np.zeros(grids.shape)
    for i, x in enumerate(xc):
        for j, y in enumerate(yc):
            for k, z in enumerate(zc):
                r = np.sqrt(x * x + y * y + z * z)
                cos = (x * og.q_vec[0] + y * og.q_vec[1] + z * og.q_vec[2]) / r
                q = og.q_uniform + 4 * og.q_dayside * max(cos, 0.0)
                n[i, j, k] = q / (4 * np.pi * 1.0 * r**2 * 1e15)
    return n


def matlab_cold_production(grids, f, og, xs, e_cold, inc_gens):
    """calculateColdeProduction.m, statement by statement (1-based indices)."""
    Eg = grids.energy
    n_eng = Eg.size - 1
    g1 = lambda i: Eg[i - 1]  # noqa: E731  MATLAB grids{4}(i)
    E_mid = 0.5 * (Eg[1:] + Eg[:-1])
    dE = Eg[1:] - Eg[:-1]
    xs_total = np.zeros(n_eng)
    for i_xs in range(1, xs.e_th.size + 1):
        e_th = xs.e_th[i_xs - 1]
        if e_th > 0:
            lims = [max(e_th, e_cold), e_cold + e_th]
            if i_xs >= 7:
                lims[1] += e_cold
            ix_th = (np.flatnonzero(Eg > lims[0])[0] + 1) - 1
            ix_cold = (np.flatnonzero(Eg > lims[1])[0] + 1) - 1
            eng_low = 0.5 * (lims[0] + g1(ix_th + 1))
            eng_upp = 0.5 * (lims[1] + g1(ix_cold))
            low_factor = (g1(ix_th + 1) - lims[0]) / (g1(ix_th + 1) - g1(ix_th))
            upp_factor = (lims[1] - g1(ix_cold)) / (g1(ix_cold + 1) - g1(ix_cold))
            xs_eng2 = E_mid.copy()
            xs_eng2[ix_th - 1] = eng_low
            xs_eng2[ix_cold - 1] = eng_upp
            interp = np.interp(xs_eng2, xs.energy, xs.xs[:, i_xs - 1], left=np.nan, right=np.nan)
            vals = interp[[ix_th - 1, ix_cold - 1]] * np.array([low_factor, upp_factor])
            interp[ix_th - 1] = vals[0]
            interp[ix_cold - 1] = vals[1]
            interp[: ix_th - 1] = 0
            interp[ix_cold:] = 0
            xs_total += interp
    v_eng = np.sqrt(2 * xs_eng2 * QE / ME)
    temp = xs_total * dE * v_eng * 1e2
    g = np.array(inc_gens) - 1
    per_mol = (f[..., g] * temp[None, None, None, :, None]).sum(axis=(3, 4))
    return per_mol * matlab_neutral_density(grids, og)


def test_plasma_densities(setup):
    grids, f, _, _ = setup
    for gens, lo, hi in [((0, 1), 0, 1), ((0,), 3, 350), ((1,), 0, 1), ((0, 1), 0.5, 12)]:
        dens, band = physics.plasma_densities(grids, f, gens, lo, hi)
        ref_dens, ref_band = matlab_plasma_densities(grids, f, [g + 1 for g in gens], lo, hi)
        np.testing.assert_allclose(dens, ref_dens, rtol=1e-12)
        np.testing.assert_allclose(band, ref_band, rtol=1e-12)


def test_neutral_density(setup):
    grids, _, og, _ = setup
    np.testing.assert_allclose(
        physics.neutral_density_on_grid(grids, og), matlab_neutral_density(grids, og), rtol=1e-12
    )
    # symmetric coma: Haser density Q / (4 pi u r^2)
    n = physics.neutral_density(100.0, 0.0, 0.0, og.symmetric())
    assert n == pytest.approx(1e27 / (4 * np.pi * 1.0 * 100.0**2 * 1e15))


def test_eii_frequency(setup):
    grids, f, _, xs = setup
    E_mid = grids.energy_mid
    sigma = np.interp(E_mid, xs.energy, xs.xs[:, 6:].sum(axis=1), left=np.nan, right=np.nan)
    sigma[np.isnan(sigma)] = 0
    flux = f.sum(axis=-1) * np.sqrt(E_mid * 2 * QE / ME) * 1e2
    integrand = sigma * flux
    expected = (0.5 * (integrand[..., 1:] + integrand[..., :-1]) * np.diff(E_mid)).sum(axis=-1)
    np.testing.assert_allclose(physics.eii_frequency(f, grids.energy, xs), expected, rtol=1e-12)


@pytest.mark.parametrize("e_cold", [1.0, 0.3, 2.5])
def test_cold_electron_production(setup, e_cold):
    grids, f, og, xs = setup
    got = physics.cold_electron_production(f, grids, og, xs, e_cold, (0, 1))
    np.testing.assert_allclose(got, matlab_cold_production(grids, f, og, xs, e_cold, [1, 2]), rtol=1e-12)
    assert np.all(got > 0)
    # the documented deviation only touches two energy bins
    clean = physics.cold_electron_production(f, grids, og, xs, e_cold, (0, 1), legacy_speed_grid=False)
    np.testing.assert_allclose(clean, got, rtol=0.2)


def test_temperature_of_maxwellian():
    T_eV, drift = 7.5, np.array([3e4, -1e4, 2e4])
    bulk = np.tile(drift, (2, 2, 2, 1))
    square = bulk**2 + QE * T_eV / ME
    np.testing.assert_allclose(physics.temperature(square, bulk, ME), T_eV)
    np.testing.assert_allclose(physics.temperature(square, bulk, ME, component=1), T_eV)


def test_combine_run_velocities_single_population(setup):
    grids, f, _, _ = setup
    rng = np.random.default_rng(0)
    bulk = rng.normal(size=grids.shape + (3, 2))
    square = rng.random(grids.shape + (3, 2)) + 2
    dens, ave_bulk, ave_square = physics.combine_run_velocities(
        bulk, square, f, bulk, square, np.zeros_like(f), grids.dE, (0,)
    )
    np.testing.assert_allclose(ave_bulk, bulk[..., 0])
    np.testing.assert_allclose(ave_square, square[..., 0])
    np.testing.assert_allclose(dens[..., 0], (f[..., 0] * grids.dE).sum(axis=-1))


def test_ambipolar_field_of_isothermal_gradient():
    grids = Grids(*(np.linspace(-100, 100, 21),) * 3, energy=np.array([0.0, 1.0]))
    X, _, _ = grids.mesh()
    L, T = 50.0, 5.0  # km, eV
    n = 100.0 * np.exp(-X / L)
    E = physics.ambipolar_field(physics.pressure(n, T), n, grids)
    # E = -(kT/e) dln(n)/dx = T[V] / L
    np.testing.assert_allclose(E[0, 2:-2], T / (L * 1e3), rtol=2e-2)
    np.testing.assert_allclose(E[1:], 0, atol=1e-12)
