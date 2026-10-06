import numpy as np
import pytest

import collide_comet as cc
from collide_comet import io


def test_bundled_cross_sections_match_source_table():
    xs = cc.load_cross_sections()
    assert xs.xs.shape == (103, 12)
    assert xs.energy[0] == pytest.approx(0.001) and xs.energy[-1] == pytest.approx(1000.0)
    # thresholds from XS/H2O_coll_E_th_table.csv in the repository
    np.testing.assert_allclose(
        xs.e_th, [0, 0.004, 0.198, 0.453, 10, 7, 12.62, 18.1, 19, 20.7, 16.95, 54.1]
    )
    # first rows of the elastic column in XS/H2O_coll_xs_table.csv
    np.testing.assert_allclose(xs.xs[:2, 0], [4.82e-12, 2.36e-12])


@pytest.mark.parametrize("name,n_process", [("Adv_scat_Coll_XS.nc", 12), ("Combined_XS.nc", 4)])
def test_other_bundled_tables(name, n_process):
    xs = cc.load_cross_sections(name)
    assert xs.xs.shape == (103, n_process)
    if name.startswith("Adv"):
        assert xs.scat_cdf.shape == (xs.scat_energy.size, xs.scat_theta.size)


def test_run_file_layout(cases):
    case = cases["-X"]
    grids, f = cc.load_run_data(case.pe_path)
    assert f.shape == grids.shape + (grids.n_energy, 2)
    bulk, square = cc.load_velocity_arrays(case.pe_path)
    assert bulk.shape == square.shape == grids.shape + (3, 2)
    rad, depth, dist = cc.load_electron_depth_dist(case.pe_path)
    assert dist.shape == (depth.size, 3, rad.size - 1, 18, 3)
    np.testing.assert_allclose(depth, np.logspace(-5, 4, depth.size))
    og = cc.load_outgassing(case.pe_path)
    np.testing.assert_allclose(og.q_vec, [-1, 0, 0])
    assert og.q_uniform + og.q_dayside == pytest.approx(og.q_spherical)
    assert io.load_electron_mass(case.pe_path) == pytest.approx(9.11e-31)
    dt_edges, dt_count = io.read_timestep_grid(case.pe_path)
    assert dt_edges.size == dt_count.size + 1


def test_first_axis_is_x(cases):
    """The '-X' case outgasses towards -x, so density must be higher there."""
    n_e, _ = cases["-X"].densities("pe", smooth=False)
    nx = n_e.shape[0]
    assert n_e[: nx // 2].sum() > 1.5 * n_e[nx // 2 :].sum()
    np.testing.assert_allclose(n_e[:, : nx // 2].sum(), n_e[:, nx // 2 :].sum(), rtol=0.05)
    # ... and the '+Z' case along the third axis
    n_e, _ = cases["+Z"].densities("pe", smooth=False)
    assert n_e[:, :, nx // 2 :].sum() > 1.5 * n_e[:, :, : nx // 2].sum()


def test_missing_variable_message(cases):
    with io._open(cases[0].pe_path) as ds:
        with pytest.raises(KeyError, match="Available variables"):
            io._get(ds, "Not A Variable")


def test_describe(cases, capsys):
    info = cc.describe(cases[0].pe_path)
    assert info["Plasma Density"][:3] == (12, 12, 12)
    assert "Plasma Density" in capsys.readouterr().out
