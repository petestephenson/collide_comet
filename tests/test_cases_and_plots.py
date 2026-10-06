import numpy as np
import pytest

import collide_comet as cc
from collide_comet import plotting as ccp


def test_case_set_basics(cases):
    assert cases.names == ["-X", "+Y", "+Z", "-Z", "Ref"]
    assert cases.reference.name == "Ref" and len(cases.asymmetric) == 4
    assert cases.missing_files() == []
    assert cases[0] is cases["-X"]
    with pytest.raises(KeyError, match="available"):
        cases["+X"]
    # the reference case drops the dayside term, the others keep it
    assert cases["Ref"].outgassing.q_dayside == 0
    assert cases["-X"].outgassing.q_dayside > 0


def test_large_arrays_are_released(cases):
    for case in cases:
        case.densities()
    loaded = [case.name for case in cases if case._raw is not None]
    assert len(loaded) <= cases.max_loaded
    # derived quantities survive and are cached
    assert cases["-X"].densities()[0] is cases["-X"].densities()[0]


def test_total_is_pe_plus_sw(cases):
    case = cases["+Y"]
    total = case.densities("total", smooth=False)[0]
    parts = case.densities("pe", smooth=False)[0] + case.densities("sw", smooth=False)[0]
    np.testing.assert_allclose(total, parts, rtol=1e-6)
    cometary = case.densities("total", cometary_only=True, smooth=False)[0]
    assert np.all(cometary < total)


@pytest.mark.parametrize("quantity", list(cc.QUANTITIES))
def test_every_quantity_plots(cases, quantity):
    fig, axes, cbar = ccp.grouped_slice_plot(cases, quantity, view="XZ")
    assert set(axes) == set(cases.names)
    fig, axes, cbar = ccp.grouped_ratio_plot(cases, quantity, view="XY")
    assert set(axes) == {"-X", "+Y", "+Z", "-Z"}


@pytest.mark.parametrize("view", ["YZ", "XZ", "XY"])
def test_slice_orientation(cases, view):
    """A field equal to the x coordinate must increase along the X axis."""
    grids = cases[0].grids
    X, Y, Z = grids.mesh()
    for field, letter in ((X, "X"), (Y, "Y"), (Z, "Z")):
        fig, ax, mesh, _ = ccp.plot_slice(field, grids, view=view, plot_lim=None)
        arr = mesh.get_array().reshape(len(mesh.get_coordinates()) - 1, -1)
        if ax.get_xlabel().startswith(letter):
            assert np.all(np.diff(arr, axis=1) > 0) and np.allclose(np.diff(arr, axis=0), 0)
        elif ax.get_ylabel().startswith(letter):
            assert np.all(np.diff(arr, axis=0) > 0) and np.allclose(np.diff(arr, axis=1), 0)
        else:  # the cut direction: constant, and it is the cell just above 0
            assert np.allclose(arr, arr.flat[0]) and 0 < arr.flat[0] < np.diff(grids.x)[0]


def test_single_case_plots(cases):
    case = cases["-X"]
    ccp.plot_quantity(case, "cold_density", view="XY")
    f = case.plasma_density_E("sw")[..., 0]
    ccp.plot_edf(f, case.grids, axis="X")
    ccp.plot_edf_ratio(cases["Ref"].plasma_density_E("sw")[..., 0], f, case.grids, axis=2)
    ccp.temperature_vs_radius_plot(case)
    ccp.timestep_histogram(case.pe_path)
    ccp.axis_cone_profiles_plot(case, n_mc=5000)
    ccp.axis_cone_profiles_plot(case, n_mc=5000, scale_by_local_exobase=True)
    ccp.plot_el_depth_2d(*case.electron_depth(), i_dim="X")


def test_grouped_figures(cases):
    ccp.grouped_edf_plot(cases, axis="X")
    ccp.grouped_edf_ratio_plot(cases, axis="Z")
    _, _, profiles = ccp.radial_profiles_plot(cases)
    assert set(profiles) == set(cases.names)
    _, _, data = ccp.angle_shell_profiles_plot(cases, r0=60, dr=40, n_mc=20000, n_angle_bins=9)
    assert "Ref relative to -X" in data["profiles"] and "+Z" in data["profiles"]
    for quantity in ("total", "cold", "fraction"):
        _, axes, data = ccp.r_angle_maps_plot(cases, quantity, n_mc=20000, n_rad_bins=10, n_angle_bins=6,
                                              show_boundary=True)
        assert len(data["panels"]) == 7  # 4 asymmetric + 3 reference panels
    _, axes, summary = ccp.cold_production_summary_plot(cases, r0=80, dr=60, n_angle_bins=6)
    assert "summary" in axes and len(summary["profiles"]) == 8
    ccp.cold_production_histogram(cases, r_range=(0, 150), yscale="log")
    ccp.cold_production_histogram(cases, log=False, normalization="pdf", percentiles=None)
    ccp.grouped_el_depth_2d_plot(cases, i_dim=2)
    _, axes, data = ccp.el_depth_angle_maps_plot(cases, r0=[50, 100], dr=20)
    assert data["cases"]["-X"]["sign"] == -1 and data["cases"]["+Z"]["axis"] == 2
    _, _, data = ccp.el_depth_distribution_plot(cases)
    assert all(d["q1"] <= d["median"] <= d["q3"] for d in data.values())


def test_asymmetric_case_peaks_towards_outgassing(cases):
    _, _, data = ccp.angle_shell_profiles_plot(cases, r0=60, dr=40, n_mc=40000, n_angle_bins=6)
    prof = data["profiles"]["-X"]["mean_cold"]
    assert prof[0] > prof[-1]  # more cold electrons at small angle to q_vec


def test_cross_section_and_trajectory_plots():
    xs = cc.load_cross_sections()
    ccp.plot_cross_sections(xs)
    ccp.plot_loss_rates(xs)
    ccp.plot_cross_sections(cc.load_cross_sections("Combined_XS.nc"))
    t = np.linspace(0, 20, 400)
    X = np.column_stack([np.cos(t), np.sin(t), 0.1 * t]) * 50
    V = np.gradient(X, t, axis=0) * 1e5
    ccp.plot_trajectory_3d(X, V, comet_radius=2.0)


def test_save_figure(cases, tmp_path, monkeypatch):
    fig, *_ = ccp.plot_quantity(cases["Ref"], "total_density")
    (tmp_path / "pyproject.toml").write_text("")
    work = tmp_path / "notebooks"
    work.mkdir()
    monkeypatch.chdir(work)
    monkeypatch.delenv("COLLIDE_COMET_FIGURES", raising=False)
    cc.set_figure_dir(None)
    written = cc.save_figure(fig, "density_r0.5", subfolder="demo", formats=["png", "pdf"])
    assert [p.relative_to(tmp_path).as_posix() for p in written] == [
        "Figures/demo/density_r0.5.png", "Figures/demo/density_r0.5.pdf"]
    assert all(p.stat().st_size > 0 for p in written)
    assert cc.save_figure(fig, "explicit.svg")[0].name == "explicit.svg"
    cc.set_figure_dir(tmp_path / "elsewhere")
    assert cc.save_figure(fig, "a")[0].parent == tmp_path / "elsewhere"
    cc.set_figure_dir(None)


def test_cold_production_histogram_options(cases, capsys):
    """Log y axis and volume-weighted percentile lines (MATLAB opts.yScale / opts.percentiles)."""
    import numpy as np

    from collide_comet import analysis as an

    _, ax, data = ccp.cold_production_histogram(
        cases, yscale="log", y_min_log=1e-5, percentiles=(50, 75, 90),
        percentile_linestyles=("-", "--"), print_percentiles=True, xlim=(-8, 2),
    )
    assert ax.get_yscale() == "log" and ax.get_ylim()[0] >= 1e-5
    assert ax.get_xlim() == (-8, 2)
    assert data["percentiles"] == [50, 75, 90]
    assert "weighted percentiles" in capsys.readouterr().out
    for entry in data["cases"].values():
        p50, p75, p90 = entry["percentiles"]
        assert p50 <= p75 <= p90
        assert not np.any(entry["hist"] <= 0)  # empty bins are NaN on a log axis
        # definition used in MATLAB weightedPercentiles: first value whose cumulative weight >= p
        order = np.argsort(entry["values"], kind="stable")
        cdf = np.cumsum(entry["volumes"][order]) / entry["volumes"].sum()
        assert p75 == entry["values"][order][np.flatnonzero(cdf >= 0.75)[0]]
    # one vertical line per case and percentile, plus the step curves
    assert sum(len(line.get_xdata()) == 2 and line.get_xdata()[0] == line.get_xdata()[1]
               for line in ax.lines) == 3 * len(data["cases"])
    assert list(an.weighted_quantiles([1, 2, 3, 4], [1, 1, 1, 1], [0.5, 0.75, 0.9])) == [2, 3, 4]


def test_r_angle_maps_exobase_and_layout(cases):
    """Exobase overlay and reference-first layout of the radius-angle maps."""
    import numpy as np

    from collide_comet import constants as c

    kw = dict(n_mc=20000, n_rad_bins=10, n_angle_bins=6)
    _, axes, data = ccp.r_angle_maps_plot(cases, "cold", **kw)
    names = [p["case"] for p in data["panels"]]
    assert names == ["Ref", "Ref", "Ref", "-X", "+Y", "+Z", "-Z"]  # reference in the top-left
    assert axes[0].get_title() == "Ref: +X axis"
    scale = 100 * c.SIGMA_EXOBASE * 1e-10 / (4 * np.pi * c.U_GAS)
    ref = data["panels"][0]
    assert list(ref["exobase_angle"]) == [0, 180]
    np.testing.assert_allclose(ref["exobase_radius"], scale * cases.reference.outgassing_raw.q_spherical)
    asym, og = data["panels"][3], cases["-X"].outgassing_raw
    expected = scale * (og.q_uniform + 4 * og.q_dayside * np.clip(np.cos(np.radians(data["angle_centres"])), 0, None))
    np.testing.assert_allclose(asym["exobase_radius"], expected)
    assert asym["exobase_radius"][0] > asym["exobase_radius"][-1]  # exobase furthest out towards q_vec
    assert all(len(ax.lines) == 1 for ax in axes)  # one exobase curve per panel

    # the options: natural order, no overlay, one reference panel, cell binning
    _, axes, data = ccp.r_angle_maps_plot(cases, "cold", reference_first=False, show_exobase=False,
                                          reference_axes=False, **kw)
    assert [p["case"] for p in data["panels"]] == ["-X", "+Y", "+Z", "-Z", "Ref"]
    assert data["panels"][0]["exobase_radius"] is None and not axes[0].lines
    _, _, cells = ccp.r_angle_maps_plot(cases, "fraction", method="cells", n_rad_bins=6, n_angle_bins=4)
    n_cells = np.prod(cases["-X"].grids.shape)
    assert 0 < cells["panels"][3]["n_per_bin"].sum() <= n_cells


def test_el_depth_summary_plot(cases):
    """Port of make_eDepth_QVec_summaryPlot: shell profile per case + one radius-angle map."""
    import numpy as np

    _, (ax1, ax2), data = ccp.el_depth_summary_plot(cases, r0=50, dr=40)
    assert data["case_2d"] == "Ref" and len(ax1.lines) == len(cases)
    assert data["cases"]["-X"]["axis"] == 0 and data["cases"]["-X"]["sign"] == -1
    assert data["cases"]["Ref"]["axis"] == 2 and data["cases"]["Ref"]["sign"] == 1
    entry = data["cases"]["+Z"]
    assert entry["shell_profile"].shape == entry["theta_deg"].shape
    assert data["plot_values"].shape == data["cases"]["Ref"]["median_depth"].shape

    _, _, rebinned = ccp.el_depth_summary_plot(cases, r0=50, dr=40, n_angle_bins=4, case_2d="-X", log=False)
    prof = rebinned["cases"]["+Z"]
    assert prof["shell_profile"].shape == (4,) and list(prof["shell_angle_centres"]) == [22.5, 67.5, 112.5, 157.5]
    native = entry["shell_profile"][(entry["theta_deg"] < 45) & (entry["shell_profile"] > 0)]
    if native.size:
        assert prof["shell_profile"][0] == np.median(native)


def test_profile_error_styles_and_map_options(cases):
    kw = dict(r0=60, dr=40, n_mc=20000, n_angle_bins=6)
    _, axes, _ = ccp.angle_shell_profiles_plot(cases, error_style="bars", **kw)
    assert axes[0].containers and not axes[0].collections == []  # error-bar containers drawn
    _, axes, _ = ccp.angle_shell_profiles_plot(cases, error_style="none", **kw)
    assert not axes[0].containers and not axes[0].collections
    _, axes, data = ccp.r_angle_maps_plot(cases, "cold", method="cells", r_max=None, n_rad_bins=6,
                                          n_angle_bins=4, show_boundary=True, boundary_linestyle="w:",
                                          show_exobase=False)
    assert data["r_edges"][-1] == max(case.grids.radius().max() for case in cases)
    assert axes[0].lines[0].get_linestyle() == ":"
