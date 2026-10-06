import numpy as np
import pytest
from scipy import stats

import collide_comet as cc
from collide_comet import analysis as an
from collide_comet.io import Grids


def test_bin_index_matches_histcounts_convention():
    edges = np.array([0.0, 1.0, 2.0, 3.0])
    vals = np.array([-0.1, 0.0, 0.99, 1.0, 2.5, 3.0, 3.01, np.nan])
    np.testing.assert_array_equal(an.bin_index(vals, edges), [-1, 0, 0, 1, 2, 2, -1, -1])


def test_mc_points_follow_matlab_rng_stream():
    # MATLAB: rng(1); rand(1,5) -> 0.4170 0.7203 0.0001 0.3023 0.1468
    u = np.random.RandomState(1).rand(5)
    np.testing.assert_allclose(u, [0.4170, 0.7203, 0.0001, 0.3023, 0.1468], atol=5e-5)
    x, y, z, r = an.shell_mc_points(20000, 45.0, 55.0, seed=1)
    np.testing.assert_allclose(np.sqrt(x**2 + y**2 + z**2), r)
    assert r.min() >= 45 and r.max() <= 55
    np.testing.assert_allclose([x.mean(), y.mean(), z.mean()], 0, atol=1.0)
    x2, *_ = an.shell_mc_points(20000, 45.0, 55.0, seed=1)
    np.testing.assert_array_equal(x, x2)


@pytest.mark.parametrize("axis", [(1, 0, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1), (1, 1, 0)])
def test_cone_points_stay_in_cone(axis):
    x, y, z, r = an.cone_mc_points(5000, axis, 10.0, 300.0, seed=2)
    angle = an.angle_to_direction(x, y, z, axis)
    assert angle.max() <= 10.0 + 1e-9
    assert r.max() <= 300.0


def test_weighted_quantiles_and_histogram():
    vals = np.array([1.0, 2.0, 3.0, 4.0])
    w = np.array([1.0, 1.0, 1.0, 5.0])
    np.testing.assert_allclose(an.weighted_quantiles(vals, w, [0.25, 0.5, 0.99]), [2, 4, 4])
    hist, centres = an.weighted_histogram(vals, w, [0.5, 2.5, 4.5], "probability")
    np.testing.assert_allclose(hist, [0.25, 0.75])
    hist, _ = an.weighted_histogram(vals, w, [0.5, 2.5, 4.5], "pdf")
    np.testing.assert_allclose(hist, [0.125, 0.375])


def test_radial_profile_of_power_law():
    grids = Grids(*(np.linspace(-200, 200, 81),) * 3, energy=np.array([0.0, 1.0]))
    r = grids.radius()
    prof = an.radial_profiles(grids, 1000.0 / r, 100.0 / r, np.arange(20.0, 200.0, 10.0))
    np.testing.assert_allclose(prof["mean_total"] * prof["r_centres"], 1000.0, rtol=0.03)
    np.testing.assert_allclose(prof["frac_cold"], 0.1, rtol=1e-12)
    shell = 4 / 3 * np.pi * (prof["r_edges"][1:] ** 3 - prof["r_edges"][:-1] ** 3)
    np.testing.assert_allclose(prof["vol_per_bin"], shell, rtol=0.06)


def test_r_angle_map_against_scipy():
    rng = np.random.default_rng(0)
    r, ang = rng.uniform(0, 300, 20000), rng.uniform(0, 180, 20000)
    vals = rng.lognormal(size=20000)
    r_edges, a_edges = np.linspace(0, 300, 11), np.linspace(0, 180, 7)
    for stat in ("median", "mean"):
        got, n = an.mc_r_angle_map(r, ang, vals, r_edges, a_edges, stat)
        want = stats.binned_statistic_2d(r, ang, vals, stat, bins=[r_edges, a_edges]).statistic
        np.testing.assert_allclose(got, want)
    assert n.sum() == 20000
    geo, _ = an.mc_r_angle_map(r, ang, vals, r_edges, a_edges, "geomean")
    assert np.all(geo < an.mc_r_angle_map(r, ang, vals, r_edges, a_edges, "mean")[0])


def test_angle_profiles_recover_cosine_dependence():
    x, y, z, _ = an.shell_mc_points(200000, 49, 51, seed=1)
    angle = an.angle_to_direction(x, y, z, (0, 0, 1))
    total = 10.0 * (2.0 + np.cos(np.radians(angle)))
    prof = an.mc_angle_profiles(angle, total, 0.2 * total, np.linspace(0, 180, 19))
    expected = 10.0 * (2.0 + np.cos(np.radians(prof["angle_centres"])))
    np.testing.assert_allclose(prof["mean_total"], expected, rtol=0.01)
    np.testing.assert_allclose(prof["frac_cold"], 0.2, rtol=1e-9)
    assert np.all(prof["err_low_total"] >= 0)


def test_shell_angle_profile_uniform_field():
    grids = Grids(*(np.linspace(-100, 100, 41),) * 3, energy=np.array([0.0, 1.0]))
    prof, n = an.shell_angle_profile(grids, np.full(grids.shape, 3.0), (1, 0, 0), 40, 60, np.linspace(0, 180, 10))
    np.testing.assert_allclose(prof, 3.0)
    assert n.sum() > 0


def matlab_median_depth(rad_edges, depth_grid, dist, i_dim, i_tau, level):
    """compute_medianDepth_grid.m with explicit loops (1-based arguments)."""
    slab = dist[1:, i_tau - 1, :, :, i_dim - 1]
    cum = np.cumsum(slab, axis=0)
    n_r, n_t = slab.shape[1:]
    out = np.full((n_r, n_t), np.nan)
    for it in range(n_t):
        for ir in range(n_r):
            denom = cum[-1, ir, it]
            if denom == 0:
                continue
            ix = np.flatnonzero(cum[:, ir, it] / denom >= level)
            if ix.size:
                out[ir, it] = depth_grid[(ix[0] + 1) + 1 - 1]  # elDepthGrid(ix + 1)
    return out


def test_median_depth_grid(cases):
    rad, grid, dist = cases["+Z"].electron_depth("total")
    dist = dist.copy()
    dist[:, :, 3, 4, :] = 0  # an empty bin must give NaN
    for i_dim in (0, 2):
        got, r_mids, theta = an.median_depth_grid(rad, grid, dist, i_dim, i_tau=2)
        np.testing.assert_array_equal(got, matlab_median_depth(rad, grid, dist, i_dim + 1, 3, 0.5))
        assert np.isnan(got[3, 4])
    assert theta[0] > theta[-1]  # runs from ~pi to ~0
    # '+Z' outgasses towards theta = 0 of the z axis -> deeper there
    assert np.nanmean(got[:, -3:]) > np.nanmean(got[:, :3])


def test_exobase_from_depth(cases):
    rad, grid, dist = cases["Ref"].electron_depth("total")
    out = an.median_depth_and_exobase(rad, grid, dist, 0)
    assert out["theta_mids"][0] == pytest.approx(np.pi) and out["theta_mids"][-1] == 0
    env = out["exo"]["envelope"]
    assert np.nanstd(env["r"]) < 0.2 * np.nanmean(env["r"])  # symmetric case: ~circular
    assert out["exo"]["x_smooth"][0] == out["exo"]["x_smooth"][-1]  # closed curve


def test_temperature_vs_radius_fit(cases):
    case = cases["Ref"]
    data = an.temperature_vs_radius(case.temperature(), case.densities()[0], case.grids, n_bins=20)
    assert np.isfinite(data["T_weighted"]).sum() > 10
    assert data["fit_function"] is None or np.all(np.isfinite(data["fit_function"](data["r_mids"])))


def test_movmean():
    x = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    np.testing.assert_allclose(an.movmean(x, 3, omitnan=True), [1.5, 1.5, 3.0, 4.5, 4.5])
    np.testing.assert_allclose(an.movmean(np.arange(5.0), 3), [0.5, 1, 2, 3, 3.5])
