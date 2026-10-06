"""Line profiles, Monte-Carlo maps and histograms comparing the cases.

Ports of ``make_electrondens_radial_AllCases_Plot``,
``make_electrondens_angleToQVec_shell_MC_AllCases_Plot_v2``,
``make_electrondens_axisCone_MC_scaledExobase_Plot``,
``make_electrondens_rAngle_QVec_2d_MC_AllCases_Plot_v2``,
``make_grouped_cold_e_prod_Plot_v3``,
``make_cold_e_prod_volumeWeighted_Hist_AllCases`` and
``makeTemperature_vs_r_Plot_v2``.

Every function returns ``(fig, axes, data)`` where ``data`` holds the plotted
numbers (the MATLAB ``prof`` / ``out`` / ``summary`` structs).
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np

from .. import analysis as an
from .. import constants as c
from .. import io
from ..cases import CaseData, CaseSet
from ..colormaps import get_cmap
from ..physics import exobase_radius, local_outgassing_from_cos
from .slices import _CBAR_POS, _GROUP_BOTTOM, _GROUP_LEFT, _GROUP_REF, _GROUP_SIZE, plot_slice

__all__ = [
    "MATLAB_LINES",
    "radial_profiles_plot",
    "angle_shell_profiles_plot",
    "axis_cone_profiles_plot",
    "r_angle_maps_plot",
    "cold_production_summary_plot",
    "cold_production_histogram",
    "temperature_vs_radius_plot",
    "timestep_histogram",
]

#: MATLAB's default ``lines`` colour order.
MATLAB_LINES = (
    (0.0000, 0.4470, 0.7410),
    (0.8500, 0.3250, 0.0980),
    (0.9290, 0.6940, 0.1250),
    (0.4940, 0.1840, 0.5560),
    (0.4660, 0.6740, 0.1880),
    (0.3010, 0.7450, 0.9330),
    (0.6350, 0.0780, 0.1840),
)
_CASE_MARKERS = ("o", "*", "^", "v", "s", "d", "p", "h")
_CASE_LINESTYLES = ("-", "-", "-", "-", "--", ":", "-.")


def _color(i: int):
    return MATLAB_LINES[i % len(MATLAB_LINES)]


def _axis_marker(name: str) -> str:
    upper = str(name).upper()
    for letter, marker in (("X", "o"), ("Y", "*"), ("Z", "^")):
        if letter in upper:
            return marker
    return "s"


def _style(ax, fontsize=12):
    ax.tick_params(direction="out", labelsize=fontsize)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)


def _density_triptych(figsize=(14.5, 5.0), fontsize=12):
    fig, axes = plt.subplots(1, 3, figsize=figsize, constrained_layout=True)
    for ax in axes:
        _style(ax, fontsize)
    return fig, axes


def _label_triptych(axes, xlabel, fontsize=12, titles=None):
    titles = titles or ("Total electron density", "Cold electron density", "Cold electron fraction")
    ylabels = (r"Total $n_e$ [cm$^{-3}$]", r"Cold $n_e$ [cm$^{-3}$]", "Cold / total")
    for ax, title, ylabel in zip(axes, titles, ylabels):
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.set_title(title, fontsize=fontsize)


# ---------------------------------------------------------------------------
# radial profiles
# ---------------------------------------------------------------------------
def radial_profiles_plot(
    cases: CaseSet,
    r_edges: Optional[Sequence[float]] = None,
    cometary_only: bool = False,
    gens: Sequence[int] = c.ALL_GENS,
    min_energy: float = 0.0,
    max_energy: float = 1.0,
    r_max: Optional[float] = None,
    markers: bool = True,
    marker_step: int = 3,
    linewidth: float = 1.5,
    yscale_density: str = "log",
    yscale_fraction: str = "log",
    fraction_ylim: Optional[Tuple[float, float]] = (1e-3, 1.0),
    legend: bool = True,
    fontsize: float = 12,
):
    """Volume-weighted radial profiles of total and cold electron density.

    Port of ``make_electrondens_radial_AllCases_Plot``.

    Returns ``fig, axes, profiles`` with ``profiles[case_name]`` the dict
    from :func:`collide_comet.analysis.radial_profiles`.
    """
    r_edges = np.arange(0.0, 300.0 + 5.0, 5.0) if r_edges is None else np.asarray(r_edges, float)
    fig, axes = _density_triptych(fontsize=fontsize)
    axes[0].set_yscale(yscale_density)
    axes[1].set_yscale(yscale_density)
    axes[2].set_yscale(yscale_fraction)

    profiles = {}
    for i, case in enumerate(cases):
        total, cold = case.densities("total", gens, min_energy, max_energy, cometary_only)
        prof = an.radial_profiles(case.grids, total, cold, r_edges)
        profiles[case.name] = prof
        style = dict(
            color=_color(i), linestyle=_CASE_LINESTYLES[i % len(_CASE_LINESTYLES)],
            marker=_CASE_MARKERS[i % len(_CASE_MARKERS)] if markers else None,
            markersize=5, markevery=marker_step, linewidth=linewidth, label=case.name,
        )
        for ax, key in zip(axes, ("mean_total", "mean_cold", "frac_cold")):
            good = np.isfinite(prof[key]) & (prof[key] > 0)
            ax.plot(prof["r_centres"][good], prof[key][good], **style)

    _label_triptych(axes, "Radius [km]", fontsize)
    for ax in axes:
        ax.set_xlim(r_edges[0] if r_max is None else 0.0, r_edges[-1] if r_max is None else r_max)
    if fraction_ylim is not None:
        axes[2].set_ylim(*fraction_ylim)
    if legend:
        axes[2].legend(loc="best", fontsize=fontsize - 1)
    return fig, axes, profiles


# ---------------------------------------------------------------------------
# angular profiles in a shell (Monte Carlo)
# ---------------------------------------------------------------------------
def _fill_profile(ax, x, y, err_low, err_high, color, linestyle, marker, label, alpha, fill_alpha=0.1,
                  error_style="band"):
    """Draw a binned profile with its spread as a shaded band or as error bars."""
    if error_style not in ("band", "bars", "none"):
        raise ValueError('error_style must be "band", "bars" or "none"')
    y_low, y_high = y - err_low, y + err_high
    good = np.isfinite(x) & np.isfinite(y) & np.isfinite(y_low) & np.isfinite(y_high)
    if ax.get_yscale() == "log":
        good &= (y > 0) & (y_high > 0)
        y_low = np.maximum(np.where(np.isfinite(y_low), y_low, 1e-30), 1e-30)
    if not good.any():
        return
    if error_style == "bars":  # scatter points with error bars (MATLAB errorbar(..., 'o-'))
        ax.errorbar(x[good], y[good], yerr=[y[good] - y_low[good], y_high[good] - y[good]],
                    color=color, linestyle=linestyle, marker=marker, markersize=5, linewidth=1.5,
                    capsize=4, alpha=alpha, label=label)
        return
    if error_style == "band":
        ax.fill_between(x[good], y_low[good], y_high[good], color=color, alpha=fill_alpha, linewidth=0)
    ax.plot(x[good], y[good], color=color, linestyle=linestyle, marker=marker, markersize=5,
            markevery=3, linewidth=1.5, alpha=alpha, label=label)


def angle_shell_profiles_plot(
    cases: CaseSet,
    r0: float = 50.0,
    dr: float = 10.0,
    n_mc: int = 200_000,
    n_angle_bins: int = 36,
    cometary_only: bool = False,
    gens: Sequence[int] = c.ALL_GENS,
    min_energy: float = 0.0,
    max_energy: float = 1.0,
    seed: int = 1,
    fraction_ylim: Optional[Tuple[float, float]] = (1e-3, 0.2),
    error_style: str = "band",
    fontsize: float = 12,
):
    """Electron density versus angle from the outgassing axis in a shell.

    Port of ``make_electrondens_angleToQVec_shell_MC_AllCases_Plot_v2``.
    Monte-Carlo points are drawn uniformly in the shell ``r0 +/- dr/2`` and
    assigned the density of the nearest cell.  Asymmetric cases are binned by
    angle to their ``q_vec``; the reference case is binned by angle to each
    asymmetric case's axis (dashed lines), giving a like-for-like baseline.

    ``error_style`` sets how the spread of the points in each angle bin (one
    standard deviation, in log space for the densities) is shown: ``"band"``
    (shaded, as in the ``_v2`` MATLAB function), ``"bars"`` (points with error
    bars, as in ``make_electrondens_angleToQVec_shell_MC_AllCases_Plot``) or
    ``"none"``.

    Returns ``fig, axes, data`` with ``data["profiles"]`` keyed by label.
    """
    r_min, r_max = r0 - 0.5 * dr, r0 + 0.5 * dr
    if r_min < 0:
        raise ValueError(f"r0 - dr/2 must be >= 0 (got {r_min:g} km)")
    angle_edges = np.linspace(0.0, 180.0, n_angle_bins + 1)
    rng = np.random.RandomState(seed)

    fig, axes = _density_triptych(fontsize=fontsize)
    for ax in axes:
        ax.set_yscale("log")

    names = cases.names
    data = {"r0": r0, "dr": dr, "angle_edges": angle_edges, "profiles": {}}
    for i, case in enumerate(cases):
        total, cold = case.densities("total", gens, min_energy, max_energy, cometary_only)
        x, y, z, _ = an.shell_mc_points(n_mc, r_min, r_max, rng)
        inside, (total_s, cold_s) = an.sample_nearest(case.grids, [total, cold], x, y, z)
        x, y, z = x[inside], y[inside], z[inside]
        if x.size == 0:
            continue

        if case.is_reference:
            targets = []
            for j, other in enumerate(cases):
                if other.is_reference:
                    continue
                try:
                    direction = an.axis_direction(other.name)
                except ValueError:
                    direction = other.outgassing_raw.q_hat
                targets.append((direction, f"{case.name} relative to {other.name}", _color(j),
                                "--", _axis_marker(other.name), 0.55))
        else:
            q_hat = case.outgassing_raw.q_hat
            if not np.all(np.isfinite(q_hat)):
                continue
            targets = [(q_hat, case.name, _color(i), "-", _axis_marker(case.name), 1.0)]

        for direction, label, color, linestyle, marker, alpha in targets:
            angle = an.angle_to_direction(x, y, z, direction)
            prof = an.mc_angle_profiles(angle, total_s, cold_s, angle_edges)
            prof.update(case=case.name, direction=np.asarray(direction), n_inside=int(x.size))
            data["profiles"][label] = prof
            for ax, key in zip(axes, ("total", "cold", "frac")):
                mean = prof["mean_" + key] if key != "frac" else prof["frac_cold"]
                _fill_profile(ax, prof["angle_centres"], mean, prof["err_low_" + key],
                              prof["err_high_" + key], color, linestyle, marker, label, alpha,
                              error_style=error_style)

    _label_triptych(axes, "Angle to reference direction [deg]", fontsize)
    for ax in axes:
        ax.set_xlim(0, 180)
        ax.set_xticks(np.arange(0, 181, 30))
    if fraction_ylim is not None:
        axes[2].set_ylim(*fraction_ylim)
    axes[2].legend(loc="best", fontsize=fontsize - 2)
    fig.suptitle(rf"Shell $r = {r0:g} \pm {0.5 * dr:g}$ km, $N = {int(n_mc):d}$ per case",
                 fontsize=fontsize)
    return fig, axes, data


# ---------------------------------------------------------------------------
# profiles along the six semi-axes (Monte-Carlo cones)
# ---------------------------------------------------------------------------
def axis_cone_profiles_plot(
    case: CaseData,
    n_mc: int = 200_000,
    scale_by_local_exobase: bool = False,
    cone_angle_deg: float = 10.0,
    r_max: float = 300.0,
    n_bins: int = 50,
    cometary_only: bool = True,
    gens: Sequence[int] = c.ALL_GENS,
    min_energy: float = 0.0,
    max_energy: float = 1.0,
    seed: int = 1,
    fontsize: float = 12,
):
    """Density profiles in narrow cones about +-X, +-Y and +-Z for one case.

    Port of ``make_electrondens_axisCone_MC_scaledExobase_Plot``.  With
    ``scale_by_local_exobase`` the distance of every sample is divided by the
    exobase distance for the local outgassing rate.

    Returns ``fig, axes, data`` with ``data[axis_label]`` per cone.
    """
    labels = ("+X", "-X", "+Y", "-Y", "+Z", "-Z")
    total, cold = case.densities("total", gens, min_energy, max_energy, cometary_only)
    og = case.outgassing_raw
    q_hat = og.q_hat
    rng = np.random.RandomState(seed)

    samples = {}
    for label in labels:
        x, y, z, r = an.cone_mc_points(n_mc, an.axis_direction(label), cone_angle_deg, r_max, rng)
        inside, (total_s, cold_s) = an.sample_nearest(case.grids, [total, cold], x, y, z)
        x, y, z, r = x[inside], y[inside], z[inside], r[inside]
        with np.errstate(invalid="ignore", divide="ignore"):
            cos_q = (x * q_hat[0] + y * q_hat[1] + z * q_hat[2]) / np.where(r == 0, np.nan, r)
        q_local = local_outgassing_from_cos(cos_q, og.q_uniform, og.q_dayside)
        dist = r / exobase_radius(q_local) if scale_by_local_exobase else r
        keep = np.isfinite(dist) & (dist >= 0) & (q_local > 0)
        samples[label] = dict(r=r[keep], q_local=q_local[keep], dist=dist[keep],
                              total=total_s[keep], cold=cold_s[keep])

    all_dist = np.concatenate([s["dist"] for s in samples.values()])
    if all_dist.size == 0:
        raise ValueError("No valid Monte-Carlo samples inside the grid")
    profile_lim = float(all_dist.max()) if scale_by_local_exobase else r_max
    edges = np.linspace(0.0, profile_lim, n_bins + 1)
    xlabel = r"$r / r_\mathrm{exo,local}$" if scale_by_local_exobase else "r [km]"

    fig, axes = _density_triptych(fontsize=fontsize)
    axes[0].set_yscale("log")
    axes[1].set_yscale("log")
    data = {"edges": edges, "q_hat": q_hat, "cone_angle_deg": cone_angle_deg}
    for i, label in enumerate(labels):
        s = samples[label]
        prof = an.mc_distance_profiles(s["dist"], s["total"], s["cold"], edges)
        prof.update(n_valid=int(s["dist"].size), mean_exobase_km=float(np.mean(exobase_radius(s["q_local"])))
                    if s["q_local"].size else np.nan)
        data[label] = prof
        style = dict(color=_color(i // 2), linestyle="-" if i % 2 == 0 else "--", linewidth=1.6,
                     label=f"{label}, N={s['dist'].size}")
        for ax, key in zip(axes, ("mean_total", "mean_cold", "frac_cold")):
            ax.plot(prof["centres"], prof[key], **style)

    _label_triptych(
        axes, xlabel, fontsize,
        titles=tuple(f"{t}, MC cones, {cone_angle_deg:g}\N{DEGREE SIGN}"
                     for t in (r"Total $n_e$", r"Cold $n_e$", "Cold fraction")),
    )
    for ax in axes:
        ax.set_xlim(0, profile_lim)
        ax.grid(True, alpha=0.3)
    axes[2].set_ylim(0, 1)
    axes[2].legend(loc="best", fontsize=fontsize - 2)
    suffix = ", distance scaled by local exobase" if scale_by_local_exobase else ""
    fig.suptitle(f"{case.name}: MC electron density along semi-axes{suffix}", fontsize=fontsize)
    return fig, axes, data


# ---------------------------------------------------------------------------
# radius-angle maps (Monte Carlo)
# ---------------------------------------------------------------------------
def r_angle_maps_plot(
    cases: CaseSet,
    quantity: str = "cold",
    statistic: str = "median",
    n_mc: int = 500_000,
    n_rad_bins: int = 80,
    n_angle_bins: int = 72,
    r_min: float = 0.0,
    r_max: Optional[float] = 300.0,
    log: Optional[bool] = None,
    clim: Optional[Sequence[float]] = None,
    cmap=None,
    cometary_only: bool = False,
    gens: Sequence[int] = c.ALL_GENS,
    min_energy: float = 0.0,
    max_energy: float = 1.0,
    seed: int = 1,
    show_boundary: bool = False,
    boundary_thresh: Optional[float] = None,
    boundary_linestyle: str = "k-",
    boundary_linewidth: float = 2.0,
    method: str = "mc",
    reference_axes: bool = True,
    reference_first: bool = True,
    show_exobase: bool = True,
    exobase_linestyle: str = "k--",
    exobase_linewidth: float = 2.0,
    exobase_sigma: float = c.SIGMA_EXOBASE,
    exobase_u_gas: float = c.U_GAS,
    fontsize: float = 11,
):
    """Maps of electron density versus radius and angle to the outgassing axis.

    Port of ``make_electrondens_rAngle_QVec_2d_MC_AllCases_Plot_v2`` with the
    exobase overlay and reference-first layout of
    ``make_electrondens_rAngle_QVec_2d_MC_AllCases_Plot`` and the cell binning
    of ``make_electrondens_rAngle_QVec_2d_AllCases_Plot``.

    Asymmetric cases get one panel (angle to ``q_vec``).  The reference case
    gets three (angle to +X, +Y, +Z), or a single one with
    ``reference_axes=False``.

    Parameters
    ----------
    quantity : {"total", "cold", "fraction"}
    statistic : {"median", "mean", "geomean"}
    method : {"mc", "cells"}
        ``"mc"`` bins ``n_mc`` Monte-Carlo points sampled uniformly in volume
        (nearest cell value); ``"cells"`` bins the grid-cell centres directly
        (``n_mc`` and ``seed`` are then unused).
    r_max : float, optional
        Outer radius of the map [km].  ``None`` uses the largest cell-centre
        radius of the grids (the corner of the box).
    show_boundary, boundary_thresh, boundary_linestyle, boundary_linewidth
        Draw the outermost radius at which the mapped quantity exceeds
        ``boundary_thresh`` (default 0.1 for the fraction, 1 otherwise).
    reference_axes : bool
        ``True``: reference case relative to +X, +Y and +Z.  ``False``: one
        reference panel relative to its stored ``q_vec`` (+X if that is zero).
    reference_first : bool
        Put the reference panel(s) first, i.e. in the top-left of the figure;
        the other cases follow in their usual order.
    show_exobase : bool
        Overlay the estimated electron exobase distance
        ``r_exo = 100 Q sigma 1e-10 / (4 pi u_gas)``.  Asymmetric cases use
        the local outgassing ``Q = q_uniform + 4 q_dayside max(cos(angle), 0)``;
        the reference case uses ``q_spherical`` (a horizontal line).
    exobase_linestyle, exobase_linewidth : matplotlib format string and width.
    exobase_sigma, exobase_u_gas : cross section [cm^2] and gas speed [km/s].

    Returns ``fig, axes, data`` with ``data["panels"]`` a list of dicts in the
    order the panels are drawn.  Each holds the map and, when the overlay is
    on, ``exobase_angle`` and ``exobase_radius``.
    """
    if quantity not in ("total", "cold", "fraction"):
        raise ValueError('quantity must be "total", "cold" or "fraction"')
    if method not in ("mc", "cells"):
        raise ValueError('method must be "mc" or "cells"')
    if r_max is None:
        r_max = max(float(case.grids.radius().max()) for case in cases)
    if r_max <= r_min or r_min < 0:
        raise ValueError("need 0 <= r_min < r_max")
    log = (quantity != "fraction") if log is None else log
    label = {"total": r"Total $n_e$ [cm$^{-3}$]", "cold": r"Cold $n_e$ [cm$^{-3}$]",
             "fraction": r"Cold / total $n_e$"}[quantity]
    r_edges = np.linspace(r_min, r_max, n_rad_bins + 1)
    angle_edges = np.linspace(0.0, 180.0, n_angle_bins + 1)
    r_centres = 0.5 * (r_edges[:-1] + r_edges[1:])
    angle_centres = 0.5 * (angle_edges[:-1] + angle_edges[1:])
    rng = np.random.RandomState(seed)

    # -- compute every panel (in case order, so the random stream is unchanged) --
    panels = []
    for case in cases:
        total, cold = case.densities("total", gens, min_energy, max_energy, cometary_only)
        if method == "mc":
            x, y, z, r = an.shell_mc_points(n_mc, r_min, r_max, rng)
            inside, (total_s, cold_s) = an.sample_nearest(case.grids, [total, cold], x, y, z)
            x, y, z, r = x[inside], y[inside], z[inside], r[inside]
        else:
            x, y, z = (a.ravel() for a in case.grids.mesh("centres"))
            r = np.sqrt(x**2 + y**2 + z**2)
            total_s, cold_s = total.ravel(), cold.ravel()
        if r.size == 0:
            continue
        with np.errstate(invalid="ignore", divide="ignore"):
            vals = {"total": total_s, "cold": cold_s, "fraction": cold_s / total_s}[quantity]

        q_hat = case.outgassing_raw.q_hat
        if case.is_reference and reference_axes:
            targets = [((1, 0, 0), "+X axis"), ((0, 1, 0), "+Y axis"), ((0, 0, 1), "+Z axis")]
        elif np.all(np.isfinite(q_hat)):
            targets = [(q_hat, "QVec")]
        elif case.is_reference:
            targets = [((1, 0, 0), "+X axis")]
        else:
            continue

        for direction, angle_name in targets:
            angle = an.angle_to_direction(x, y, z, direction)
            map_vals, n_per_bin = an.mc_r_angle_map(r, angle, vals, r_edges, angle_edges, statistic, quantity)
            plot_vals = map_vals.copy()
            if log:
                with np.errstate(invalid="ignore", divide="ignore"):
                    plot_vals = np.log10(np.where(plot_vals > 0, plot_vals, np.nan))

            r_boundary, thresh = None, None
            if show_boundary:
                thresh = boundary_thresh if boundary_thresh is not None else (0.1 if quantity == "fraction" else 1.0)
                r_boundary = an.largest_radius_above_threshold(map_vals, r_centres, thresh)

            exo_angle, exo_radius = None, None
            if show_exobase:
                og = case.outgassing  # reference: q_dayside = 0, q_uniform = q_spherical
                if case.is_reference:
                    exo_angle = np.array([0.0, 180.0])
                    q_local = np.full(2, og.q_spherical)
                else:
                    exo_angle = angle_centres
                    q_local = local_outgassing_from_cos(
                        np.cos(np.radians(angle_centres)), og.q_uniform, og.q_dayside)
                exo_radius = exobase_radius(q_local, exobase_sigma, exobase_u_gas)

            panels.append(dict(case=case.name, is_reference=case.is_reference, angle_name=angle_name,
                               direction=np.asarray(direction, float),
                               map=map_vals, plot_values=plot_vals, n_per_bin=n_per_bin,
                               n_inside=int(r.size), r_boundary=r_boundary, boundary_thresh=thresh,
                               exobase_angle=exo_angle, exobase_radius=exo_radius))

    if not panels:
        raise ValueError("No case produced any samples inside its grid")
    if reference_first:  # stable: reference panel(s) first, others keep their order
        panels.sort(key=lambda p: not p["is_reference"])

    # -- draw ---------------------------------------------------------------
    n_panels = len(panels)
    n_rows = 2 if n_panels > 1 else 1
    n_cols = int(np.ceil(n_panels / n_rows))
    fig, grid = plt.subplots(n_rows, n_cols, figsize=(3.4 * n_cols + 0.8, 6.5 if n_rows == 2 else 3.6),
                             constrained_layout=True, squeeze=False)
    flat_axes = list(grid.ravel())
    used = flat_axes[:n_panels]
    for ax in flat_axes[n_panels:]:
        ax.set_visible(False)

    for ax, panel in zip(used, panels):
        _style(ax, fontsize)
        ax.pcolormesh(angle_edges, r_edges, np.ma.masked_invalid(panel["plot_values"]),
                      cmap=get_cmap(cmap), shading="flat", rasterized=True)
        ax.set_xlim(0, 180)
        ax.set_xticks(np.arange(0, 181, 30))
        ax.set_ylim(r_min, r_max)
        ax.set_ylabel("Radius [km]", fontsize=fontsize)
        if panel["angle_name"] == "QVec":
            ax.set_xlabel("Angle to max outgassing [deg]", fontsize=fontsize)
            ax.set_title(panel["case"], fontsize=fontsize)
        else:
            ax.set_xlabel(f"Angle to {panel['angle_name']} [deg]", fontsize=fontsize)
            ax.set_title(f"{panel['case']}: {panel['angle_name']}", fontsize=fontsize)
        if panel["r_boundary"] is not None:
            ok = np.isfinite(panel["r_boundary"])
            ax.plot(angle_centres[ok], panel["r_boundary"][ok], boundary_linestyle,
                    linewidth=boundary_linewidth)
        if panel["exobase_radius"] is not None:
            ok = np.isfinite(panel["exobase_angle"]) & np.isfinite(panel["exobase_radius"])
            ax.plot(panel["exobase_angle"][ok], panel["exobase_radius"][ok], exobase_linestyle,
                    linewidth=exobase_linewidth, label="Estimated exobase")

    if clim is None:  # shared limits so the single colour bar is valid for every panel
        finite = np.concatenate([p["plot_values"][np.isfinite(p["plot_values"])] for p in panels])
        clim = (finite.min(), finite.max()) if finite.size else (0.0, 1.0)
    for ax in used:
        ax.collections[0].set_clim(*clim)
    cbar = fig.colorbar(used[0].collections[0], ax=used, location="right", shrink=0.9)
    cbar.set_label(rf"$\log_{{10}}$({label})" if log else label, fontsize=fontsize)
    if show_exobase:
        used[0].legend(loc="upper right", fontsize=fontsize - 2, frameon=True)
    how = f"MC {quantity}" if method == "mc" else f"Cell-binned {quantity}"
    tail = f", N = {int(n_mc)} per case" if method == "mc" else ""
    fig.suptitle(f"{how} electron-density map vs radius and angle ({statistic}{tail})", fontsize=fontsize)
    data = dict(r_edges=r_edges, angle_edges=angle_edges, r_centres=r_centres,
                angle_centres=angle_centres, panels=panels, colorbar=cbar)
    return fig, used, data


# ---------------------------------------------------------------------------
# cold-electron production
# ---------------------------------------------------------------------------
def cold_production_summary_plot(
    cases: CaseSet,
    view: str = "XZ",
    r0: float = 50.0,
    dr: float = 10.0,
    n_angle_bins: int = 36,
    e_cold: float = 1.0,
    clim: Sequence[float] = (-3.0, 3.0),
    cmap=None,
    plot_lim: float = 300.0,
    figsize: Tuple[float, float] = (10.0, 6.0),
    fontsize: float = 11,
):
    """Cold-electron production maps plus an angular summary in a shell.

    Port of ``make_grouped_cold_e_prod_Plot_v3``: the five slice panels of the
    grouped figure and, bottom-left, the volume-weighted mean production in
    the shell ``r0 +/- dr/2`` versus angle to the outgassing direction.

    Returns ``fig, axes, summary``; ``axes`` is a dict keyed by case name
    plus ``"summary"``.
    """
    r_min, r_max = r0 - 0.5 * dr, r0 + 0.5 * dr
    if r_min < 0:
        raise ValueError(f"r0 - dr/2 must be >= 0 (got {r_min:g} km)")
    angle_edges = np.linspace(0.0, 180.0, n_angle_bins + 1)
    angle_centres = 0.5 * (angle_edges[:-1] + angle_edges[1:])

    fig = plt.figure(figsize=figsize)
    asym = cases.asymmetric
    if len(asym) > 4:
        raise ValueError("This layout supports at most four asymmetric cases")
    axes = {}
    for i, case in enumerate(asym):
        axes[case.name] = fig.add_axes([_GROUP_LEFT[i], _GROUP_BOTTOM[i], *_GROUP_SIZE])
    ref = cases.reference if cases.reference_name else None
    if ref is not None:
        axes[ref.name] = fig.add_axes([*_GROUP_REF, *_GROUP_SIZE])
    ax_sum = fig.add_axes([_GROUP_REF[0], 0.10, *_GROUP_SIZE])
    axes["summary"] = ax_sum
    ax_sum.set_yscale("log")
    _style(ax_sum, fontsize)

    summary = dict(r0=r0, dr=dr, angle_edges=angle_edges, angle_centres=angle_centres, profiles={})
    mesh = None
    for i, case in enumerate(cases):
        prod = case.cold_electron_production(e_cold)
        ax = axes[case.name]
        _, _, mesh, _ = plot_slice(prod, case.grids, view=view, log=True, clim=clim, cmap=cmap,
                                   plot_lim=plot_lim, ax=ax, colorbar=False, fontsize=fontsize)
        ax.set_title(case.name, fontsize=fontsize)
        if case.is_reference or asym.index(case) < 2:  # top row
            ax.set_xlabel("")
        if not case.is_reference and asym.index(case) != 2:
            ax.set_ylabel("")

        if case.is_reference:
            targets = []
            for j, other in enumerate(cases):
                if other.is_reference:
                    continue
                try:
                    direction = an.axis_direction(other.name)
                except ValueError:
                    direction = other.outgassing_raw.q_hat
                targets.append((direction, f"{case.name}, {other.name}", _color(j), "--",
                                _CASE_MARKERS[j % len(_CASE_MARKERS)]))
        else:
            q_hat = case.outgassing_raw.q_hat
            if not np.all(np.isfinite(q_hat)):
                continue
            targets = [(q_hat, case.name, _color(i), "-", _CASE_MARKERS[i % len(_CASE_MARKERS)])]

        for direction, label, color, linestyle, marker in targets:
            profile, n_cells = an.shell_angle_profile(case.grids, prod, direction, r_min, r_max, angle_edges)
            summary["profiles"][label] = dict(profile=profile, n_cells=n_cells, direction=np.asarray(direction))
            good = np.isfinite(profile) & (profile > 0)
            ax_sum.plot(angle_centres[good], profile[good], color=color, linestyle=linestyle,
                        marker=marker, markersize=4, linewidth=1.3, label=label)

    ax_sum.set_xlabel("Angle to outgassing direction [deg]", fontsize=fontsize)
    ax_sum.set_ylabel(r"Mean $P_e$ [cm$^{-3}$ s$^{-1}$]", fontsize=fontsize)
    ax_sum.set_title(rf"$r = {r0:g} \pm {0.5 * dr:g}$ km", fontsize=fontsize)
    ax_sum.set_xlim(0, 180)
    ax_sum.set_xticks(np.arange(0, 181, 30))
    ax_sum.legend(ncol=2, fontsize=fontsize - 4, loc="best")
    cbar = fig.colorbar(mesh, cax=fig.add_axes(_CBAR_POS))
    cbar.set_label(r"$\log_{10}$($P_e$ [cm$^{-3}$ s$^{-1}$])", fontsize=fontsize)
    summary["colorbar"] = cbar
    return fig, axes, summary


def cold_production_histogram(
    cases: CaseSet,
    e_cold: float = 1.0,
    gens: Sequence[int] = c.ALL_GENS,
    n_bins: int = 60,
    bin_edges: Optional[Sequence[float]] = None,
    normalization: str = "probability",
    r_range: Optional[Tuple[float, float]] = None,
    cometary_only: bool = False,
    log: bool = True,
    yscale: str = "linear",
    y_min_log: float = 1e-6,
    percentiles: Optional[Sequence[float]] = (50, 75, 90),
    percentile_linestyles: Sequence[str] = ("-", "--", ":"),
    print_percentiles: bool = False,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    linewidth: float = 1.8,
    title: str = "Volume-weighted cold electron production",
    fontsize: float = 12,
):
    """Volume-weighted distribution of the cold-electron production rate.

    Port of ``make_cold_e_prod_volumeWeighted_Hist_AllCases``.

    Parameters
    ----------
    normalization : {"probability", "count", "pdf"}
    r_range : (r_min, r_max), optional
        Restrict to cells with centres in this range of cometocentric distance.
    percentiles : sequence, optional
        Volume-weighted percentiles marked by vertical lines (``None`` = off).
        MATLAB ``opts.showPercentiles`` / ``opts.percentiles``.
    percentile_linestyles : sequence of str
        Line style of each percentile, cycled (median solid, 75th dashed,
        90th dotted by default).  The lines take the colour of their case.
    print_percentiles : bool
        Print the percentile values of every case (and the linear production
        rate when ``log=True``).
    yscale : {"linear", "log"}
        Scale of the y axis.  With ``"log"`` empty bins are not drawn and the
        lower limit is kept above ``y_min_log``.
    xlim, ylim : (low, high), optional
        Axis limits, applied last.

    Returns ``fig, ax, data`` with ``data["cases"][name]`` per case; the
    percentile values are in ``data["cases"][name]["percentiles"]`` (in
    ``log10`` units when ``log=True``) and ``data["percentiles"]`` lists the
    percentages they belong to.
    """
    per_case = {}
    for case in cases:
        prod = case.cold_electron_production(e_cold, gens=gens, cometary_only=cometary_only)
        vol = case.grids.cell_volume
        valid = np.isfinite(prod) & (prod > 0) & (vol > 0)
        if r_range is not None:
            r = case.grids.radius()
            valid &= (r >= r_range[0]) & (r < r_range[1])
        prod_vals, vol_vals = prod[valid], vol[valid]
        plot_vals = np.log10(prod_vals) if log else prod_vals
        entry = dict(values=plot_vals, volumes=vol_vals, total_volume=float(vol_vals.sum()),
                     n_cells=int(prod_vals.size))
        if percentiles is not None:
            entry["percentiles"] = an.weighted_quantiles(plot_vals, vol_vals, np.asarray(percentiles) / 100.0)
            if print_percentiles:
                print(f"{case.name} weighted percentiles:")
                for pct, value in zip(percentiles, entry["percentiles"]):
                    if log:
                        print(f"  P{pct:02g} = {value:.3f} log10, {10.0 ** value:.3e} cm^-3 s^-1")
                    else:
                        print(f"  P{pct:02g} = {value:.3e} cm^-3 s^-1")
        per_case[case.name] = entry

    if bin_edges is None:
        everything = np.concatenate([e["values"] for e in per_case.values()])
        everything = everything[np.isfinite(everything)]
        if everything.size == 0:
            raise ValueError("No valid cold electron production values found")
        lo, hi = np.floor(everything.min()), np.ceil(everything.max())
        if lo == hi:
            lo, hi = lo - 0.5, hi + 0.5
        bin_edges = np.linspace(lo, hi, n_bins + 1)
    bin_edges = np.asarray(bin_edges, dtype=float)

    fig, ax = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)
    _style(ax, fontsize)
    pct_styles = tuple(percentile_linestyles)
    for i, (name, entry) in enumerate(per_case.items()):
        hist, centres = an.weighted_histogram(entry["values"], entry["volumes"], bin_edges, normalization)
        if yscale == "log":
            hist = np.where(hist > 0, hist, np.nan)
        ax.step(centres, hist, where="post", color=_color(i), linewidth=linewidth, label=name)
        entry.update(hist=hist, bin_centres=centres)
        if percentiles is not None:
            for k, value in enumerate(entry["percentiles"]):
                if np.isfinite(value):
                    ax.axvline(value, color=_color(i), linestyle=pct_styles[k % len(pct_styles)], linewidth=1.0, alpha=0.7)

    ax.set_yscale(yscale)
    if yscale == "log":
        low, high = ax.get_ylim()
        ax.set_ylim(max(low, y_min_log), high)
    if log:
        ax.set_xlabel(r"$\log_{10}$($P_e$ [cm$^{-3}$ s$^{-1}$])", fontsize=fontsize)
    else:
        ax.set_xlabel(r"$P_e$ [cm$^{-3}$ s$^{-1}$]", fontsize=fontsize)
        ax.set_xscale("log")
    ax.set_ylabel({"probability": "Volume fraction", "count": r"Cell volume per bin [km$^3$]",
                   "pdf": "Volume-weighted PDF"}[normalization], fontsize=fontsize)
    full_title = title
    if r_range is not None:
        full_title += "\n" + rf"${r_range[0]:g} \leq r < {r_range[1]:g}$ km"
    ax.set_title(full_title, fontsize=fontsize)
    handles, labels = ax.get_legend_handles_labels()
    if percentiles is not None:
        for k, pct in enumerate(percentiles):
            handles.append(plt.Line2D([], [], color="k", linestyle=pct_styles[k % len(pct_styles)], linewidth=1.0))
            labels.append(f"{pct:g}th percentile")
    ax.legend(handles, labels, fontsize=fontsize - 2, loc="best")
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    return fig, ax, dict(
        bin_edges=bin_edges,
        cases=per_case,
        percentiles=None if percentiles is None else list(percentiles),
    )


# ---------------------------------------------------------------------------
# temperature vs radius, time steps
# ---------------------------------------------------------------------------
def temperature_vs_radius_plot(
    case: CaseData,
    mode: str = "combined",
    n_bins: int = 100,
    fit: bool = True,
    ax: Optional[plt.Axes] = None,
    max_points: Optional[int] = 200_000,
    fontsize: float = 12,
):
    """Electron temperature of every cell against cometocentric distance.

    Port of ``makeTemperature_vs_r_Plot_v2``: scatter coloured by X, the
    density-weighted and unweighted radial means, the weighted 1-sigma
    spread and a ``rat22`` fit to the weighted mean.

    ``max_points`` thins the scatter (not the statistics) for large grids.

    Returns ``fig, ax, data`` (see
    :func:`collide_comet.analysis.temperature_vs_radius`).
    """
    T = case.temperature(mode)
    dens = case.densities("total")[0]
    data = an.temperature_vs_radius(T, dens, case.grids, n_bins, fit)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
    else:
        fig = ax.figure
    _style(ax, fontsize)
    step = 1 if not max_points else max(1, int(np.ceil(data["r"].size / max_points)))
    sc = ax.scatter(data["r"][::step], data["T"][::step], s=5, c=data["x"][::step], cmap="viridis",
                    rasterized=True, linewidths=0)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("X [km]", fontsize=fontsize)

    mids = data["r_mids"]
    ax.plot(mids, data["T_weighted"], "k-", linewidth=1.6, label=r"$\langle T_e \rangle_w$ (density-weighted)")
    ax.plot(mids, data["T_unweighted"], "k--", linewidth=1.4, label=r"$\langle T_e \rangle$ (unweighted)")
    ok = np.isfinite(data["T_weighted"]) & np.isfinite(data["T_std_weighted"])
    if ok.any():
        ax.errorbar(mids[ok], data["T_weighted"][ok], yerr=data["T_std_weighted"][ok], linestyle="none",
                    color=(0, 0, 0, 0.25), capsize=0, label=r"1$\sigma$ (weighted)")
    if data["fit_function"] is not None:
        r_fit = np.linspace(data["r_edges"][0], data["r_edges"][-1], 200)
        ax.plot(r_fit, data["fit_function"](r_fit), "b-", linewidth=1.6, label="Rational fit (rat22)")
    ax.set_xlabel("r [km]", fontsize=fontsize)
    ax.set_ylabel(r"$T_e$ [eV]", fontsize=fontsize)
    ax.legend(loc="upper left", fontsize=fontsize - 2)
    data["colorbar"] = cbar
    return fig, ax, data


def timestep_histogram(path, ax: Optional[plt.Axes] = None, fontsize: float = 12):
    """Histogram of the particle time steps stored in a run file.

    Port of the ``dT`` histogram at the top of ``asymmetricPlots.m``.
    Returns ``fig, ax, (dt_edges, dt_count)``.
    """
    dt_edges, dt_count = io.read_timestep_grid(path)
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
    else:
        fig = ax.figure
    n = min(dt_count.size, dt_edges.size - 1)
    ax.stairs(dt_count[:n], dt_edges[: n + 1], fill=True, alpha=0.7)
    ax.set_xscale("log")
    ax.set_xlabel("Time step [s]", fontsize=fontsize)
    ax.set_ylabel("Count", fontsize=fontsize)
    _style(ax, fontsize)
    return fig, ax, (dt_edges, dt_count)
