"""Electron (collisional) depth figures.

Ports of ``asymPlots.make_2dElDepthPlot``, ``make_grouped_ElDepth2d_Plot``,
``make_eDepth_QVec_2d_AllCases_Plot``, ``make_eDepth_QVec_summaryPlot`` and
``plot_eDepth_dist``.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .. import analysis as an
from ..cases import CaseData, CaseSet
from ..colormaps import get_cmap
from ..physics import exobase_radius, local_outgassing_from_cos
from .profiles import _color, _style
from .slices import _AXIS_NAMES, _finish_group, _grouped_axes

__all__ = [
    "plot_el_depth_2d",
    "grouped_el_depth_2d_plot",
    "el_depth_angle_maps_plot",
    "el_depth_summary_plot",
    "el_depth_distribution_plot",
]

_SHELL_LINESTYLES = ("-", "--", ":", "-.")


def _dim_index(i_dim) -> int:
    if isinstance(i_dim, str):
        return _AXIS_NAMES.index(i_dim.upper())
    return int(i_dim)


def plot_el_depth_2d(
    rad_edges: np.ndarray,
    depth_grid: np.ndarray,
    depth_dist: np.ndarray,
    i_dim=0,
    i_tau: int = 2,
    median_level: float = 0.5,
    exo_thresh_log10: float = 0.0,
    clim: Sequence[float] = (-2.0, 2.0),
    cmap=None,
    r_max: Optional[float] = None,
    show_angle_rays: bool = True,
    show_hull: bool = True,
    ax: Optional[plt.Axes] = None,
    colorbar: bool = True,
    fontsize: float = 11,
    **hull_kwargs,
):
    """Half-disc map of the median electron depth about one axis.

    Port of ``asymPlots.make_2dElDepthPlot``.  The horizontal axis is the
    distance along the polar axis ``i_dim`` (0/1/2 or "X"/"Y"/"Z"), the
    vertical one the distance from it.  The black curve is the smoothed
    convex hull of the region where ``log10(median depth) >
    exo_thresh_log10``, i.e. an estimate of the electron exobase.

    Returns ``fig, ax, cbar, data``.
    """
    i_dim = _dim_index(i_dim)
    data = an.median_depth_and_exobase(
        rad_edges, depth_grid, depth_dist, i_dim, i_tau, median_level, exo_thresh_log10,
        pad_theta=False, **hull_kwargs,
    )
    rad_edges = np.asarray(rad_edges, dtype=float)
    r_max = float(rad_edges.max()) if r_max is None else r_max
    tb = an.theta_bounds(depth_dist.shape[3])

    if ax is None:
        fig, ax = plt.subplots(figsize=(8.0, 4.6))
    else:
        fig = ax.figure
    xv = rad_edges[:, None] * np.cos(tb)[None, :]
    yv = rad_edges[:, None] * np.sin(tb)[None, :]
    mesh = ax.pcolormesh(xv, yv, np.ma.masked_invalid(data["log10_median_depth"]),
                         cmap=get_cmap(cmap), shading="flat", rasterized=True)
    mesh.set_clim(*clim)
    ax.set_aspect("equal")
    ax.set_xlim(-r_max, r_max)
    ax.set_ylim(0, r_max)
    ax.set_xlabel(f"{_AXIS_NAMES[i_dim].lower()} [km]", fontsize=fontsize)
    ax.yaxis.set_visible(False)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=fontsize)
    ax.patch.set_alpha(0.0)

    exo = data["exo"]
    if show_hull and "x_smooth" in exo:
        ax.plot(exo["x_smooth"], exo["y_smooth"], "k-", linewidth=1.7)
    if show_angle_rays:
        for ang_deg in range(0, 181, 30):
            ang = np.radians(ang_deg)
            xe, ye = r_max * np.cos(ang), r_max * np.sin(ang)
            ax.plot([0, xe], [0, ye], "--", color=(0.2, 0.2, 0.2), linewidth=0.7)
            ax.text(1.06 * xe, 1.06 * ye, f"{ang_deg}\N{DEGREE SIGN}", ha="center", va="center",
                    fontsize=fontsize - 2, clip_on=False)

    cbar = None
    if colorbar:
        cbar = fig.colorbar(mesh, ax=ax, shrink=0.85)
        cbar.set_label(r"$\log_{10}$(median electron depth)", fontsize=fontsize)
    data["mesh"] = mesh
    return fig, ax, cbar, data


def grouped_el_depth_2d_plot(
    cases: CaseSet,
    i_dim=0,
    source: str = "total",
    clim: Sequence[float] = (-2.0, 0.5),
    cmap=None,
    titles: bool = True,
    figsize: Tuple[float, float] = (10.0, 6.0),
    fontsize: float = 11,
    **kwargs,
):
    """Median electron-depth half-discs for every case.

    Port of ``make_grouped_ElDepth2d_Plot``.  Returns ``fig, axes, cbar, data``.
    """
    fig = plt.figure(figsize=figsize)
    case_list = list(cases)
    axes, bottom_row = _grouped_axes(fig, case_list, ratio=False)
    data, mesh = {}, None
    for case in case_list:
        rad_edges, depth_grid, dist = case.electron_depth(source)
        _, _, _, out = plot_el_depth_2d(rad_edges, depth_grid, dist, i_dim, clim=clim, cmap=cmap,
                                        ax=axes[case.name], colorbar=False, fontsize=fontsize, **kwargs)
        data[case.name] = out
        mesh = out["mesh"]
    cbar = _finish_group(fig, axes, bottom_row, mesh, r"$\log_{10}$(median electron depth)",
                         False, fontsize)
    for name, ax in axes.items():
        ax.set_xlabel(f"{_AXIS_NAMES[_dim_index(i_dim)].lower()} [km]", fontsize=fontsize)
        if titles:
            ax.set_title(name, fontsize=fontsize, pad=14)
    return fig, axes, cbar, data


def _case_axis(case: CaseData, index: int) -> Tuple[int, int]:
    """Polar axis and sign of a case's outgassing direction.

    Port of ``get_case_asymmetry_dim``: taken from the case name (``"-X"`` ->
    axis 0, sign -1); the reference case is analysed about -X.
    """
    if case.is_reference:
        return 0, -1
    try:
        direction = an.axis_direction(case.name)
    except ValueError:
        q_hat = case.outgassing_raw.q_hat
        if not np.all(np.isfinite(q_hat)):
            return 0, -1
        direction = q_hat
    axis = int(np.argmax(np.abs(direction)))
    return axis, (1 if direction[axis] >= 0 else -1)


def _oriented_median_depth(case, axis, sign, source, i_tau, median_level):
    rad_edges, depth_grid, dist = case.electron_depth(source)
    median, r_mids, theta_mids = an.median_depth_grid(rad_edges, depth_grid, dist, axis, i_tau, median_level)
    theta_deg = np.degrees(theta_mids)
    if sign < 0:
        theta_deg = 180.0 - theta_deg
    order = np.argsort(theta_deg)
    return median[:, order], rad_edges, r_mids, theta_deg[order]


def el_depth_angle_maps_plot(
    cases: CaseSet,
    r0: Sequence[float] = 50.0,
    dr: Sequence[float] = 10.0,
    source: str = "total",
    i_tau: int = 2,
    median_level: float = 0.5,
    clim: Sequence[float] = (-2.0, 0.5),
    log: bool = True,
    cmap=None,
    r_max: Optional[float] = None,
    depth_thresh: float = 1.0,
    show_depth_boundary: bool = True,
    show_exobase: bool = True,
    figsize: Tuple[float, float] = (14.0, 8.0),
    fontsize: float = 11,
):
    """Median electron depth versus radius and angle to the outgassing axis.

    Port of ``make_eDepth_QVec_2d_AllCases_Plot``: one map per case (angle
    measured from the direction of maximum outgassing) and a summary panel
    with the median depth in the shell(s) ``r0 +/- dr/2`` against angle, each
    asymmetric case compared with the reference case about the same axis.

    Parameters
    ----------
    r0, dr : float or sequences
        Shell centre(s) and width(s) [km] for the summary panel.
    depth_thresh : float
        The solid line marks the outermost radius where the median depth
        exceeds this value (1 = collisional).
    show_exobase : bool
        Dashed line: analytic exobase distance for the local outgassing rate.

    Returns ``fig, axes, data``; ``axes`` is a dict keyed by case name plus
    ``"summary"``.
    """
    r0 = np.atleast_1d(np.asarray(r0, dtype=float))
    dr = np.atleast_1d(np.asarray(dr, dtype=float))
    r0, dr = np.broadcast_arrays(r0, dr)
    case_list = list(cases)
    asym = [case for case in case_list if not case.is_reference]
    refs = [case for case in case_list if case.is_reference]
    if len(asym) > 4:
        raise ValueError("The fixed 2 x 3 layout supports at most four asymmetric cases")

    fig, grid = plt.subplots(2, 3, figsize=figsize, constrained_layout=True)
    tiles = {**{case.name: ax for case, ax in zip(asym, (grid[0, 1], grid[0, 2], grid[1, 1], grid[1, 2]))}}
    if refs:
        tiles[refs[0].name] = grid[0, 0]
    ax_sum = grid[1, 0]
    used = set(tiles.values()) | {ax_sum}
    for ax in grid.ravel():
        if ax not in used:
            ax.set_visible(False)

    ax_sum.set_yscale("log")
    _style(ax_sum, fontsize)
    data = dict(r0=r0, dr=dr, cases={})
    ref = refs[0] if refs else None
    meshes = []
    for i, case in enumerate(case_list):
        axis, sign = _case_axis(case, i)
        median, rad_edges, r_mids, theta_deg = _oriented_median_depth(case, axis, sign, source, i_tau, median_level)
        entry = dict(axis=axis, sign=sign, r_mids=r_mids, theta_deg=theta_deg, median_depth=median,
                     shell_profiles=[])

        if not case.is_reference:
            ref_median = None
            if ref is not None:
                ref_median, _, ref_r, ref_theta = _oriented_median_depth(ref, axis, sign, source, i_tau, median_level)
            for k, (r_c, d_r) in enumerate(zip(r0, dr)):
                suffix = "" if r0.size == 1 else f", r = {r_c:g} km"
                profile = an.shell_median_profile(median, r_mids, r_c, d_r)
                entry["shell_profiles"].append(profile)
                good = np.isfinite(profile) & (profile > 0)
                ax_sum.plot(theta_deg[good], profile[good], color=_color(i),
                            linestyle=_SHELL_LINESTYLES[k % 4], marker="o", markersize=3.5,
                            linewidth=1.5, label=case.name + suffix)
                if ref_median is not None:
                    ref_profile = an.shell_median_profile(ref_median, ref_r, r_c, d_r)
                    good = np.isfinite(ref_profile) & (ref_profile > 0)
                    ax_sum.plot(ref_theta[good], ref_profile[good], color=_color(i), alpha=0.5,
                                linestyle=_SHELL_LINESTYLES[(k + 1) % 4], linewidth=1.2,
                                label=f"{ref.name} {case.name}{suffix}")

        ax = tiles[case.name]
        _style(ax, fontsize)
        plot_vals = median.copy()
        if log:
            with np.errstate(invalid="ignore", divide="ignore"):
                plot_vals = np.log10(np.where(plot_vals > 0, plot_vals, np.nan))
        # theta bins are uniform in cos(theta); rebuild their edges in degrees
        tb = np.degrees(an.theta_bounds(theta_deg.size))[::-1]
        mesh = ax.pcolormesh(tb, rad_edges, np.ma.masked_invalid(plot_vals), cmap=get_cmap(cmap),
                             shading="flat", rasterized=True)
        meshes.append(mesh)
        entry["plot_values"] = plot_vals

        if show_exobase:
            og = case.outgassing_raw
            if case.is_reference:
                exo_theta = np.array([0.0, 180.0])
                exo_r = exobase_radius(og.q_spherical) * np.ones(2)
            else:
                exo_theta = theta_deg
                q_local = local_outgassing_from_cos(np.cos(np.radians(theta_deg)), og.q_uniform, og.q_dayside)
                exo_r = exobase_radius(q_local)
            ax.plot(exo_theta, exo_r, "k--", linewidth=1.0, label="analytic exobase")
            entry.update(exobase_theta=exo_theta, exobase_radius=exo_r)
        if show_depth_boundary:
            r_boundary = an.largest_radius_above_threshold(median, r_mids, depth_thresh)
            ok = np.isfinite(r_boundary)
            ax.plot(theta_deg[ok], r_boundary[ok], "k-", linewidth=1.0,
                    label=rf"median depth $>$ {depth_thresh:g}")
            entry["r_boundary"] = r_boundary

        ax.set_xlim(0, 180)
        ax.set_xticks(np.arange(0, 181, 30))
        if r_max is not None:
            ax.set_ylim(0, r_max)
        ax.set_title(case.name, fontsize=fontsize)
        ax.set_xlabel("Angle to max outgassing [deg]", fontsize=fontsize)
        ax.set_ylabel("Radius [km]", fontsize=fontsize)
        data["cases"][case.name] = entry

    if log:
        shared = tuple(clim)
    else:
        finite = np.concatenate([e["plot_values"][np.isfinite(e["plot_values"])] for e in data["cases"].values()])
        shared = (finite.min(), finite.max()) if finite.size else (0.0, 1.0)
    for mesh in meshes:
        mesh.set_clim(*shared)

    ax_sum.set_xlabel("Angle to max outgassing [deg]", fontsize=fontsize)
    ax_sum.set_ylabel("Median electron depth", fontsize=fontsize)
    ax_sum.set_xlim(0, 180)
    ax_sum.set_xticks(np.arange(0, 181, 30))
    if ax_sum.get_legend_handles_labels()[0]:
        ax_sum.legend(ncol=2, fontsize=fontsize - 3, loc="best")
    cbar = fig.colorbar(meshes[0], ax=list(tiles.values()), location="right", shrink=0.9)
    cbar.set_label(r"$\log_{10}$(median electron depth)" if log else "Median electron depth",
                   fontsize=fontsize)
    data["colorbar"] = cbar
    return fig, {**tiles, "summary": ax_sum}, data


def _rebin_profile_median(theta_deg: np.ndarray, values: np.ndarray, angle_edges: np.ndarray) -> np.ndarray:
    """Median of the positive ``values`` in each angle bin (``rebin_angle_profile_median``)."""
    out = np.full(angle_edges.size - 1, np.nan)
    idx = an.bin_index(theta_deg, angle_edges)
    good = (idx >= 0) & np.isfinite(values) & (values > 0)
    for k in np.unique(idx[good]):
        out[k] = np.median(values[good & (idx == k)])
    return out


def el_depth_summary_plot(
    cases: CaseSet,
    r0: float = 50.0,
    dr: float = 10.0,
    n_angle_bins: Optional[int] = None,
    case_2d: Optional[str] = None,
    source: str = "total",
    i_tau: int = 2,
    median_level: float = 0.5,
    clim: Sequence[float] = (-2.0, 0.5),
    log: bool = True,
    cmap=None,
    figsize: Tuple[float, float] = (12.5, 5.2),
    fontsize: float = 12,
):
    """Two-panel summary of the electron depth.

    Port of ``make_eDepth_QVec_summaryPlot``.

    * Left: median electron depth in the shell ``r0 +/- dr/2`` against the
      angle to the direction of maximum outgassing, one line per case.
    * Right: rectangular map of the median electron depth for one case, with
      the angle to maximum outgassing on x and the radius on y.

    Each asymmetric case is analysed about its own outgassing axis.  As in the
    MATLAB function, the reference case is analysed about +Z.

    Parameters
    ----------
    n_angle_bins : int, optional
        Re-bin the shell profile into this many equal angle bins (median of
        the native bins in each).  ``None`` keeps the native bins, which are
        uniform in ``cos(angle)``.
    case_2d : str, optional
        Name of the case shown in the map (default: the reference case, or
        the last case if there is none).
    i_tau : int
        0-based index of the cross section used for the depth (MATLAB
        ``iTau = 3`` is ``i_tau=2``).

    Returns ``fig, (ax_profile, ax_map), data``.
    """
    case_list = list(cases)
    if case_2d is None:
        case_2d = next((case.name for case in case_list if case.is_reference), case_list[-1].name)
    if case_2d not in [case.name for case in case_list]:
        raise KeyError(f"case_2d={case_2d!r} is not one of {[case.name for case in case_list]}")

    def axis_of(case, index):
        return (2, 1) if case.is_reference else _case_axis(case, index)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, constrained_layout=True)
    _style(ax1, fontsize)
    _style(ax2, fontsize)
    ax1.set_yscale("log")
    ax1.grid(True, alpha=0.3)

    data = dict(r0=r0, dr=dr, cases={}, case_2d=case_2d)
    angle_edges = None if n_angle_bins is None else np.linspace(0.0, 180.0, n_angle_bins + 1)
    for i, case in enumerate(case_list):
        axis, sign = axis_of(case, i)
        median, rad_edges, r_mids, theta_deg = _oriented_median_depth(case, axis, sign, source, i_tau, median_level)
        profile = an.shell_median_profile(median, r_mids, r0, dr)
        if angle_edges is None:
            angle_centres = theta_deg
        else:
            angle_centres = 0.5 * (angle_edges[:-1] + angle_edges[1:])
            profile = _rebin_profile_median(theta_deg, profile, angle_edges)
        good = np.isfinite(profile) & (profile > 0)
        ax1.plot(angle_centres[good], profile[good], "o-", color=_color(i), markersize=5,
                 linewidth=1.5, label=case.name)
        data["cases"][case.name] = dict(axis=axis, sign=sign, r_mids=r_mids, rad_edges=rad_edges,
                                        theta_deg=theta_deg, median_depth=median,
                                        shell_angle_centres=angle_centres, shell_profile=profile)

    ax1.set_xlabel("Angle to max outgassing [deg]", fontsize=fontsize)
    ax1.set_ylabel("Median electron depth", fontsize=fontsize)
    ax1.set_title(rf"Median electron depth at $r = {r0:g} \pm {0.5 * dr:g}$ km", fontsize=fontsize)
    ax1.set_xlim(0, 180)
    ax1.set_xticks(np.arange(0, 181, 30))
    if ax1.get_legend_handles_labels()[0]:
        ax1.legend(loc="best", fontsize=fontsize - 2)

    entry = data["cases"][case_2d]
    plot_vals = entry["median_depth"].copy()
    if log:
        with np.errstate(invalid="ignore", divide="ignore"):
            plot_vals = np.log10(np.where(plot_vals > 0, plot_vals, np.nan))
    # theta bins are uniform in cos(theta); rebuild their edges in degrees
    tb = np.degrees(an.theta_bounds(entry["theta_deg"].size))[::-1]
    mesh = ax2.pcolormesh(tb, entry["rad_edges"], np.ma.masked_invalid(plot_vals), cmap=get_cmap(cmap),
                          shading="flat", rasterized=True)
    if log:
        mesh.set_clim(*clim)
    ax2.set_xlim(0, 180)
    ax2.set_xticks(np.arange(0, 181, 30))
    ax2.set_xlabel("Angle to max outgassing [deg]", fontsize=fontsize)
    ax2.set_ylabel("Radius [km]", fontsize=fontsize)
    axis_name = "XYZ"[entry["axis"]]
    ax2.set_title(f"Radius-angle electron depth map, {case_2d} "
                  f"(about {'+' if entry['sign'] > 0 else '-'}{axis_name})", fontsize=fontsize)
    cbar = fig.colorbar(mesh, ax=ax2)
    cbar.set_label(r"$\log_{10}$(median electron depth)" if log else "Median electron depth",
                   fontsize=fontsize)
    fig.suptitle("Electron depth summary vs radius and angle", fontsize=fontsize)
    data.update(plot_values=plot_vals, colorbar=cbar)
    return fig, (ax1, ax2), data


def el_depth_distribution_plot(
    cases: CaseSet,
    i_xs=2,
    source: str = "pe",
    i_dim: int = 0,
    max_radius: float = 200.0,
    xlim: Tuple[float, float] = (1e-4, 1e3),
    figsize: Tuple[float, float] = (5.5, 4.4),
    fontsize: float = 10,
):
    """Distribution of electron depth for every case, with median and IQR.

    Port of ``asymPlots.plot_eDepth_dist``.

    Parameters
    ----------
    i_xs : int or sequence of int
        Index/indices along the collision-process axis of the depth
        distribution (0-based; MATLAB ``ixXs = 3`` is ``i_xs = 2``).
    source : {"pe", "sw", "total"}
        The MATLAB function uses the photoelectron run only (default).
    max_radius : float
        Only electrons produced inside this cometocentric distance [km].

    Returns ``fig, (ax_top, ax), data``.
    """
    i_xs = np.atleast_1d(np.asarray(i_xs, dtype=int))
    case_list = list(cases)
    n_cases = len(case_list)
    fig = plt.figure(figsize=figsize)
    ax_top = fig.add_axes([0.14, 0.73, 0.81, 0.20])
    ax = fig.add_axes([0.14, 0.12, 0.81, 0.55], sharex=ax_top)
    ax.set_xscale("log")

    data = {}
    for i, case in enumerate(case_list):
        rad_edges, depth_grid, dist = case.electron_depth(source)
        r_mids = 0.5 * (rad_edges[1:] + rad_edges[:-1])
        in_r = r_mids < max_radius
        w = dist[..., i_dim][:, i_xs][:, :, in_r, :].sum(axis=(1, 2, 3))
        tau = np.sqrt(depth_grid[1:] * depth_grid[:-1])
        n = min(tau.size, w.size)
        tau, w = tau[:n], w[:n].astype(float)
        if n >= 2:  # drop the first (underflow) bin
            tau, w = tau[1:], w[1:]
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
        if w.sum() > 0:
            w = w / w.sum()
        median, q1, q3 = an.weighted_quantiles(tau, w, [0.50, 0.25, 0.75])

        ax.plot(tau, w, color=_color(i), linewidth=1.5, label=case.name)
        row = i + 1
        ax_top.plot([q1, q3], [row, row], color=_color(i), linewidth=2)
        ax_top.plot(median, row, "o", color=_color(i), markerfacecolor="w", markersize=5)
        data[case.name] = dict(tau=tau, weight=w, median=median, q1=q1, q3=q3)

    for a in (ax, ax_top):
        a.axvline(1.0, color="k", linestyle="--", linewidth=1.5)
        a.tick_params(direction="out", labelsize=fontsize)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
    ax.set_xlim(*xlim)
    ax.set_ylabel("Electron distribution", fontsize=fontsize)
    ax.set_xlabel(r"$\tau_e$", fontsize=fontsize + 1)
    ax.legend(loc="best", fontsize=fontsize - 1)
    ax_top.set_ylabel(r"Median $\tau_e$ (IQR)", fontsize=fontsize)
    ax_top.set_ylim(0.5, n_cases + 0.5)
    ax_top.set_yticks(range(1, n_cases + 1))
    ax_top.set_yticklabels([case.name for case in case_list])
    plt.setp(ax_top.get_xticklabels(), visible=False)
    return fig, (ax_top, ax), data
