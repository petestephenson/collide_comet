"""Profiles, binned statistics, Monte-Carlo sampling and electron-depth analysis.

These are the numerical parts of the larger MATLAB plotting functions
(``make_electrondens_*``, ``make_cold_e_prod_*``, ``compute_2dElDepth*``,
``makeTemperature_vs_r_Plot_v2`` ...), separated from the drawing code so they
can be used and tested on their own.  All 3D arrays are indexed
``[ix, iy, iz]``.
"""

from __future__ import annotations

import warnings
from typing import Optional, Sequence, Tuple

import numpy as np
from scipy import optimize, signal, spatial

from .io import Grids

__all__ = [
    "bin_index",
    "movmean",
    "weighted_quantiles",
    "weighted_histogram",
    "radial_profiles",
    "shell_mc_points",
    "cone_mc_points",
    "nearest_index",
    "sample_nearest",
    "angle_to_direction",
    "axis_direction",
    "mc_angle_profiles",
    "mc_distance_profiles",
    "mc_r_angle_map",
    "shell_angle_profile",
    "largest_radius_above_threshold",
    "temperature_vs_radius",
    "rat22",
    "median_depth_grid",
    "shell_median_profile",
    "exobase_hull",
    "exobase_envelope",
    "median_depth_and_exobase",
]


# ---------------------------------------------------------------------------
# generic binning helpers
# ---------------------------------------------------------------------------
def bin_index(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """0-based bin index of ``values`` in ``edges``; -1 if outside.

    Bins are ``[e_i, e_{i+1})`` with the last bin closed on the right, which
    is the behaviour of MATLAB ``histcounts`` / ``discretize``.
    """
    values = np.asarray(values, dtype=float)
    edges = np.asarray(edges, dtype=float)
    idx = np.searchsorted(edges, values, side="right") - 1
    idx[values == edges[-1]] = edges.size - 2
    bad = (idx < 0) | (idx >= edges.size - 1) | ~np.isfinite(values)
    idx[bad] = -1
    return idx


def _bin_sum(idx: np.ndarray, weights: Optional[np.ndarray], n_bins: int) -> np.ndarray:
    return np.bincount(idx, weights=weights, minlength=n_bins).astype(float)


def _bin_mean_std(idx: np.ndarray, vals: np.ndarray, n_bins: int):
    """Per-bin count, mean and sample standard deviation (NaN if empty)."""
    n = _bin_sum(idx, None, n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = _bin_sum(idx, vals, n_bins) / n
        # sample std (N-1), as MATLAB std; 0 for a single sample
        resid = vals - mean[idx]
        var = _bin_sum(idx, resid**2, n_bins) / (n - 1)
    var[n == 1] = 0.0
    std = np.sqrt(var)
    mean[n == 0] = np.nan
    std[n == 0] = np.nan
    return n, mean, std


def movmean(x: np.ndarray, span: int, omitnan: bool = False) -> np.ndarray:
    """Centred moving mean with a window that shrinks at the ends.

    Mirrors MATLAB ``movmean(x, span)`` for odd ``span`` (optionally
    ``'omitnan'``).
    """
    x = np.asarray(x, dtype=float)
    half = int(span) // 2
    out = np.full(x.shape, np.nan)
    for i in range(x.size):
        window = x[max(0, i - half) : i + half + 1]
        if omitnan:
            window = window[np.isfinite(window)]
        if window.size:
            out[i] = window.mean()
    return out


def weighted_quantiles(values, weights, probs: Sequence[float]) -> np.ndarray:
    """Weighted quantiles: smallest value whose cumulative weight >= p.

    Port of ``local_weighted_quantiles`` / ``weightedPercentiles``.  ``probs``
    are fractions in [0, 1].  Returns NaN where undefined.
    """
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    probs = np.atleast_1d(np.asarray(probs, dtype=float))
    out = np.full(probs.shape, np.nan)
    good = np.isfinite(values) & np.isfinite(weights) & (weights >= 0)
    values, weights = values[good], weights[good]
    if values.size == 0 or weights.sum() <= 0:
        return out
    order = np.argsort(values, kind="stable")
    values, weights = values[order], weights[order]
    cdf = np.cumsum(weights) / weights.sum()
    for i, p in enumerate(probs):
        j = np.flatnonzero(cdf >= p)
        out[i] = values[j[0]] if j.size else values[-1]
    return out


def weighted_histogram(values, weights, edges, normalization: str = "probability"):
    """Histogram of ``values`` weighted by ``weights``.

    Port of ``weightedHistCountsManual``.  ``normalization`` is ``"count"``
    (sum of weights per bin), ``"probability"`` or ``"pdf"``.
    Returns ``(hist, bin_centres)``.
    """
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    edges = np.asarray(edges, dtype=float)
    n_bins = edges.size - 1
    centres = 0.5 * (edges[:-1] + edges[1:])

    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    idx = bin_index(values[good], edges)
    w = weights[good][idx >= 0]
    idx = idx[idx >= 0]
    if idx.size == 0:
        return np.full(n_bins, np.nan), centres

    hist = _bin_sum(idx, w, n_bins)
    total = hist.sum()
    if normalization == "probability":
        hist = hist / total if total > 0 else np.full(n_bins, np.nan)
    elif normalization == "pdf":
        hist = hist / total / np.diff(edges) if total > 0 else np.full(n_bins, np.nan)
    elif normalization != "count":
        raise ValueError("normalization must be 'count', 'probability' or 'pdf'")
    return hist, centres


# ---------------------------------------------------------------------------
# cell-based (volume-weighted) profiles
# ---------------------------------------------------------------------------
def radial_profiles(
    grids: Grids, dens_total: np.ndarray, dens_cold: np.ndarray, r_edges: np.ndarray
) -> dict:
    """Volume-weighted radial profiles of total and cold electron density.

    Port of ``radialVolumeWeightedProfiles_local``.

    Returns a dict with ``r_centres``, ``mean_total``, ``mean_cold``,
    ``frac_cold``, ``vol_per_bin`` and ``n_cells_per_bin``.
    """
    r_edges = np.asarray(r_edges, dtype=float)
    n_bins = r_edges.size - 1
    r = grids.radius().ravel()
    vol = grids.cell_volume.ravel()
    tot = np.asarray(dens_total, dtype=float).ravel()
    cold = np.asarray(dens_cold, dtype=float).ravel()

    good = (
        np.isfinite(r) & np.isfinite(vol) & np.isfinite(tot) & np.isfinite(cold)
        & (vol > 0) & (tot > 0) & (cold >= 0)
    )
    idx = bin_index(r[good], r_edges)
    in_bin = idx >= 0
    idx, vol, tot, cold = idx[in_bin], vol[good][in_bin], tot[good][in_bin], cold[good][in_bin]

    n_cells = _bin_sum(idx, None, n_bins)
    vol_bin = _bin_sum(idx, vol, n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_total = _bin_sum(idx, tot * vol, n_bins) / vol_bin
        mean_cold = _bin_sum(idx, cold * vol, n_bins) / vol_bin
        frac_cold = mean_cold / mean_total

    empty = n_cells <= 0
    vol_bin[empty] = np.nan
    mean_total[empty | ~(mean_total > 0)] = np.nan
    mean_cold[empty | ~(mean_cold > 0)] = np.nan
    frac_cold[empty | ~(frac_cold > 0)] = np.nan
    return {
        "r_edges": r_edges,
        "r_centres": 0.5 * (r_edges[:-1] + r_edges[1:]),
        "mean_total": mean_total,
        "mean_cold": mean_cold,
        "frac_cold": frac_cold,
        "vol_per_bin": vol_bin,
        "n_cells_per_bin": n_cells,
    }


def angle_to_direction(x, y, z, direction) -> np.ndarray:
    """Angle [deg] between position vectors and ``direction``."""
    d = np.asarray(direction, dtype=float).ravel()
    d = d / np.linalg.norm(d)
    x, y, z = (np.asarray(a, dtype=float) for a in (x, y, z))
    r = np.sqrt(x**2 + y**2 + z**2)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_angle = (x * d[0] + y * d[1] + z * d[2]) / np.where(r == 0, np.nan, r)
    return np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))


def axis_direction(name: str) -> np.ndarray:
    """Unit vector for a case/axis label such as ``"-X"`` or ``"+z"``.

    Port of ``caseAxisDirectionFromName``.
    """
    key = str(name).strip().upper()
    table = {
        "+X": (1, 0, 0), "-X": (-1, 0, 0),
        "+Y": (0, 1, 0), "-Y": (0, -1, 0),
        "+Z": (0, 0, 1), "-Z": (0, 0, -1),
    }
    if key not in table:
        raise ValueError(f"Cannot interpret {name!r} as an axis direction (+X, -Y, ...)")
    return np.array(table[key], dtype=float)


def shell_angle_profile(
    grids: Grids,
    values: np.ndarray,
    direction,
    r_min: float,
    r_max: float,
    angle_edges: np.ndarray,
    where: str = "lower",
) -> Tuple[np.ndarray, np.ndarray]:
    """Volume-weighted mean of ``values`` in a radial shell vs angle.

    Port of ``coldProdShellAngleProfile``.  Only cells with ``values > 0``
    contribute.  ``where="lower"`` positions cells at their lower boundaries,
    as the MATLAB function does; use ``"centres"`` for cell centres.

    Returns ``(profile, n_cells)`` with one entry per angle bin.
    """
    angle_edges = np.asarray(angle_edges, dtype=float)
    n_bins = angle_edges.size - 1
    X, Y, Z = grids.mesh(where)
    r = np.sqrt(X**2 + Y**2 + Z**2)
    angle = angle_to_direction(X, Y, Z, direction)
    vol = grids.cell_volume
    values = np.asarray(values, dtype=float)

    good = (
        (r >= r_min) & (r < r_max) & np.isfinite(angle) & np.isfinite(values)
        & (values > 0) & (vol > 0)
    )
    idx = bin_index(angle[good], angle_edges)
    in_bin = idx >= 0
    idx, v, w = idx[in_bin], values[good][in_bin], vol[good][in_bin]

    n_cells = _bin_sum(idx, None, n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        profile = _bin_sum(idx, v * w, n_bins) / _bin_sum(idx, w, n_bins)
    profile[(n_cells <= 0) | ~(profile > 0)] = np.nan
    return profile, n_cells


def largest_radius_above_threshold(map_vals: np.ndarray, r_centres: np.ndarray, thresh: float):
    """Outermost radius at which each column of ``map_vals`` exceeds ``thresh``.

    ``map_vals`` has shape ``(n_r, n_angle)``.  Returns NaN for columns that
    never exceed the threshold.
    """
    map_vals = np.asarray(map_vals, dtype=float)
    r_centres = np.asarray(r_centres, dtype=float)
    above = np.isfinite(map_vals) & (map_vals > thresh)
    r = np.where(above, r_centres[:, None], -np.inf).max(axis=0)
    return np.where(np.isfinite(r), r, np.nan)


# ---------------------------------------------------------------------------
# Monte-Carlo sampling
# ---------------------------------------------------------------------------
def _rng(seed):
    """Random generator matching MATLAB's ``rng(seed); rand(n, 1)`` stream."""
    if isinstance(seed, np.random.RandomState):
        return seed
    return np.random.RandomState(seed)


def shell_mc_points(n_mc: int, r_min: float, r_max: float, seed=1):
    """Points uniformly distributed in volume in a spherical shell.

    Port of ``generateShellMCPoints``.  With an integer ``seed`` the points
    are identical to those MATLAB generates after ``rng(seed)`` (both use the
    Mersenne Twister with the same double-precision conversion).

    Returns ``(x, y, z, r)``.
    """
    rng = _rng(seed)
    n_mc = int(n_mc)
    u_r, u_mu, u_phi = rng.rand(n_mc), rng.rand(n_mc), rng.rand(n_mc)
    r = (r_min**3 + u_r * (r_max**3 - r_min**3)) ** (1.0 / 3.0)
    mu = -1.0 + 2.0 * u_mu
    phi = 2.0 * np.pi * u_phi
    sin_theta = np.sqrt(1.0 - mu**2)
    return r * sin_theta * np.cos(phi), r * sin_theta * np.sin(phi), r * mu, r


def _rotation_from_z(target: np.ndarray) -> np.ndarray:
    """Rotation matrix taking +z onto ``target`` (Rodrigues formula)."""
    target = np.asarray(target, dtype=float).ravel()
    target = target / np.linalg.norm(target)
    cos = target[2]
    if cos > 1 - 1e-12:
        return np.eye(3)
    if cos < -1 + 1e-12:
        return np.diag([1.0, -1.0, -1.0])
    v = np.cross([0.0, 0.0, 1.0], target)
    s = np.linalg.norm(v)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - cos) / s**2)


def cone_mc_points(n_mc: int, axis, cone_angle_deg: float, r_max: float, seed=1):
    """Points uniformly distributed in volume in a cone about ``axis``.

    Port of ``generateConeMCPoints``.  Returns ``(x, y, z, r)``.
    """
    rng = _rng(seed)
    n_mc = int(n_mc)
    cos_cone = np.cos(np.radians(cone_angle_deg))
    u_r, u_mu, u_phi = rng.rand(n_mc), rng.rand(n_mc), rng.rand(n_mc)
    r = r_max * u_r ** (1.0 / 3.0)
    theta = np.arccos(cos_cone + (1 - cos_cone) * u_mu)
    phi = 2.0 * np.pi * u_phi
    local = np.column_stack(
        [r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)]
    )
    pts = local @ _rotation_from_z(axis).T
    return pts[:, 0], pts[:, 1], pts[:, 2], r


def nearest_index(centres: np.ndarray, query: np.ndarray) -> np.ndarray:
    """Index of the nearest grid centre for every query value (clamped)."""
    centres = np.asarray(centres, dtype=float)
    query = np.asarray(query, dtype=float)
    if centres.size < 2:
        return np.zeros(query.shape, dtype=int)
    boundaries = 0.5 * (centres[1:] + centres[:-1])
    return np.clip(np.searchsorted(boundaries, query, side="right"), 0, centres.size - 1)


def sample_nearest(grids: Grids, arrays: Sequence[np.ndarray], x, y, z):
    """Nearest-cell values of ``arrays`` at the points ``(x, y, z)``.

    Points outside the grid are dropped.  Returns ``(inside, samples)`` where
    ``inside`` is the boolean mask of retained points and ``samples`` is a
    list with one 1D array per input array.
    """
    x, y, z = (np.asarray(a, dtype=float) for a in (x, y, z))
    inside = (
        (x >= grids.x[0]) & (x <= grids.x[-1])
        & (y >= grids.y[0]) & (y <= grids.y[-1])
        & (z >= grids.z[0]) & (z <= grids.z[-1])
    )
    xc, yc, zc = grids.centres
    ix = nearest_index(xc, x[inside])
    iy = nearest_index(yc, y[inside])
    iz = nearest_index(zc, z[inside])
    return inside, [np.asarray(a)[ix, iy, iz] for a in arrays]


def _binned_log_stats(idx: np.ndarray, vals: np.ndarray, n_bins: int):
    """Geometric mean and asymmetric 1-sigma (in log10) errors per bin."""
    good = np.isfinite(vals) & (vals > 0)
    _, log_mean, log_std = _bin_mean_std(idx[good], np.log10(vals[good]), n_bins)
    mid = 10.0**log_mean
    err_low = mid - 10.0 ** (log_mean - log_std)
    err_high = 10.0 ** (log_mean + log_std) - mid
    bad = ~np.isfinite(mid) | ~np.isfinite(err_low) | ~np.isfinite(err_high) | (mid <= 0)
    for a in (mid, err_low, err_high):
        a[bad] = np.nan
    return mid, err_low, err_high


def mc_angle_profiles(angle_deg, total, cold, angle_edges) -> dict:
    """Angular profiles of total / cold density and cold fraction.

    Port of ``calcMCVolumeAngleProfiles``: geometric means with log-space
    standard deviations for the densities, and for the fraction the ratio of
    the geometric means with a point-wise standard deviation.
    """
    angle_edges = np.asarray(angle_edges, dtype=float)
    n_bins = angle_edges.size - 1
    angle_deg, total, cold = (np.asarray(a, dtype=float).ravel() for a in (angle_deg, total, cold))
    valid = np.isfinite(angle_deg) & np.isfinite(total) & np.isfinite(cold) & (total > 0) & (cold >= 0)
    idx = bin_index(angle_deg[valid], angle_edges)
    in_bin = idx >= 0
    idx, total, cold = idx[in_bin], total[valid][in_bin], cold[valid][in_bin]
    frac_pts = cold / total

    n_per_bin = _bin_sum(idx, None, n_bins)
    mean_total, lo_total, hi_total = _binned_log_stats(idx, total, n_bins)
    mean_cold, lo_cold, hi_cold = _binned_log_stats(idx, cold, n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = mean_cold / mean_total
    _, f_mean, f_std = _bin_mean_std(idx, frac_pts, n_bins)
    lo_frac = frac - np.maximum(f_mean - f_std, 0.0)
    hi_frac = np.minimum(f_mean + f_std, 1.0) - frac
    lo_frac[~(lo_frac >= 0)] = np.nan
    hi_frac[~(hi_frac >= 0)] = np.nan

    few = n_per_bin <= 1
    for mean, lo, hi in ((mean_total, lo_total, hi_total), (mean_cold, lo_cold, hi_cold)):
        bad = few | ~(mean > 0)
        mean[bad], lo[bad], hi[bad] = np.nan, np.nan, np.nan
    bad = few | ~np.isfinite(frac) | (frac < 0) | (frac > 1)
    frac[bad], lo_frac[bad], hi_frac[bad] = np.nan, np.nan, np.nan

    return {
        "angle_edges": angle_edges,
        "angle_centres": 0.5 * (angle_edges[:-1] + angle_edges[1:]),
        "mean_total": mean_total, "err_low_total": lo_total, "err_high_total": hi_total,
        "mean_cold": mean_cold, "err_low_cold": lo_cold, "err_high_cold": hi_cold,
        "frac_cold": frac, "err_low_frac": lo_frac, "err_high_frac": hi_frac,
        "n_per_bin": n_per_bin,
    }


def mc_distance_profiles(distance, total, cold, edges) -> dict:
    """Arithmetic-mean profiles of sampled densities vs distance.

    Port of ``calcMCProfiles`` (used by the axis-cone plots).
    """
    edges = np.asarray(edges, dtype=float)
    n_bins = edges.size - 1
    distance, total, cold = (np.asarray(a, dtype=float).ravel() for a in (distance, total, cold))
    valid = (
        np.isfinite(distance) & np.isfinite(total) & np.isfinite(cold)
        & (total > 0) & (cold >= 0) & (distance >= 0) & (distance <= edges[-1])
    )
    idx = bin_index(distance[valid], edges)
    in_bin = idx >= 0
    idx, total, cold = idx[in_bin], total[valid][in_bin], cold[valid][in_bin]
    n = _bin_sum(idx, None, n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_total = _bin_sum(idx, total, n_bins) / n
        mean_cold = _bin_sum(idx, cold, n_bins) / n
        frac = mean_cold / mean_total
    mean_total[(n <= 0) | ~(mean_total > 0)] = np.nan
    mean_cold[(n <= 0) | ~(mean_cold >= 0)] = np.nan
    frac[(n <= 0) | ~np.isfinite(frac)] = np.nan
    return {
        "edges": edges,
        "centres": 0.5 * (edges[:-1] + edges[1:]),
        "mean_total": mean_total,
        "mean_cold": mean_cold,
        "frac_cold": frac,
        "n_per_bin": n,
    }


def mc_r_angle_map(
    r, angle_deg, values, r_edges, angle_edges, statistic: str = "median", quantity: str = "total"
) -> Tuple[np.ndarray, np.ndarray]:
    """2D map of a sampled quantity versus radius and angle.

    Port of ``bin_MC_r_angle_values``.

    Parameters
    ----------
    statistic : {"median", "mean", "geomean"}
    quantity : {"total", "cold", "fraction"}
        Controls which samples are valid (``> 0`` for densities, within
        ``[0, 1]`` for the cold fraction).

    Returns ``(map_vals, n_per_bin)``, both of shape ``(n_r, n_angle)``.
    """
    r, angle_deg, values = (np.asarray(a, dtype=float).ravel() for a in (r, angle_deg, values))
    valid = np.isfinite(r) & np.isfinite(angle_deg) & np.isfinite(values)
    if quantity in ("total", "cold"):
        valid &= values > 0
    elif quantity == "fraction":
        valid &= (values >= 0) & (values <= 1)
    else:
        raise ValueError("quantity must be 'total', 'cold' or 'fraction'")
    r, angle_deg, values = r[valid], angle_deg[valid], values[valid]

    r_edges = np.asarray(r_edges, dtype=float)
    angle_edges = np.asarray(angle_edges, dtype=float)
    shape = (r_edges.size - 1, angle_edges.size - 1)
    ir, ia = bin_index(r, r_edges), bin_index(angle_deg, angle_edges)
    good = (ir >= 0) & (ia >= 0)
    flat = np.ravel_multi_index((ir[good], ia[good]), shape)
    values = values[good]
    n_flat = shape[0] * shape[1]
    n_per_bin = _bin_sum(flat, None, n_flat)

    if statistic == "mean":
        with np.errstate(invalid="ignore", divide="ignore"):
            out = _bin_sum(flat, values, n_flat) / n_per_bin
    elif statistic == "geomean":
        with np.errstate(invalid="ignore", divide="ignore"):
            out = 10.0 ** (_bin_sum(flat, np.log10(values), n_flat) / n_per_bin)
    elif statistic == "median":
        out = np.full(n_flat, np.nan)
        if flat.size:
            order = np.lexsort((values, flat))
            f_sorted, v_sorted = flat[order], values[order]
            starts = np.flatnonzero(np.r_[True, f_sorted[1:] != f_sorted[:-1]])
            counts = np.diff(np.r_[starts, f_sorted.size])
            lo = starts + (counts - 1) // 2
            hi = starts + counts // 2
            out[f_sorted[starts]] = 0.5 * (v_sorted[lo] + v_sorted[hi])
    else:
        raise ValueError("statistic must be 'median', 'mean' or 'geomean'")
    out[n_per_bin == 0] = np.nan
    return out.reshape(shape), n_per_bin.reshape(shape)


# ---------------------------------------------------------------------------
# temperature versus radius
# ---------------------------------------------------------------------------
def rat22(x, p1, p2, p3, q1, q2):
    """MATLAB ``rat22`` model: ``(p1 x^2 + p2 x + p3) / (x^2 + q1 x + q2)``."""
    return (p1 * x**2 + p2 * x + p3) / (x**2 + q1 * x + q2)


def temperature_vs_radius(
    temperature: np.ndarray, plas_dens: np.ndarray, grids: Grids, n_bins: int = 100, fit: bool = True
) -> dict:
    """Density-weighted radial temperature profile (with optional rat22 fit).

    Port of the numerical part of ``makeTemperature_vs_r_Plot_v2``.  Cells are
    positioned at their lower boundaries, as in MATLAB.
    """
    X, Y, Z = grids.mesh("lower")
    r = np.sqrt(X**2 + Y**2 + Z**2).ravel()
    T = np.asarray(temperature, dtype=float).ravel()
    W = np.asarray(plas_dens, dtype=float).ravel()
    Xf = X.ravel()
    good = np.isfinite(r) & np.isfinite(T) & np.isfinite(W) & (W >= 0)
    r, T, W, Xf = r[good], T[good], W[good], Xf[good]

    edges = np.linspace(r.min(), r.max(), n_bins + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])
    idx = bin_index(r, edges)
    n = _bin_sum(idx, None, n_bins)
    den = _bin_sum(idx, W, n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_w = _bin_sum(idx, W * T, n_bins) / den
        var_w = _bin_sum(idx, W * (T - mean_w[idx]) ** 2, n_bins) / den
        mean_u = _bin_sum(idx, T, n_bins) / n
        var_u = _bin_sum(idx, (T - mean_u[idx]) ** 2, n_bins) / n
    for a in (mean_w, var_w, mean_u, var_u):
        a[n == 0] = np.nan

    out = {
        "r": r, "T": T, "x": Xf,
        "r_edges": edges, "r_mids": mids,
        "T_weighted": mean_w, "T_unweighted": mean_u,
        "T_var_weighted": var_w, "T_std_weighted": np.sqrt(var_w),
        "T_var_unweighted": var_u,
        "fit_params": None, "fit_function": None,
    }
    ok = np.isfinite(mids) & np.isfinite(mean_w)
    if fit and ok.sum() >= 5:
        p0 = np.array([0.117417650855806, 0.296675873218327, 0.318778301925882,
                       0.424166759713807, 0.507858284661118])
        lower = np.array([-np.inf, -np.inf, 1e3, -np.inf, -np.inf])
        try:
            popt, _ = optimize.curve_fit(
                rat22, mids[ok], mean_w[ok], p0=np.maximum(p0, lower),
                bounds=(lower, np.inf), maxfev=20000,
            )
            out["fit_params"] = popt
            out["fit_function"] = lambda x, _p=popt: rat22(np.asarray(x, dtype=float), *_p)
        except (RuntimeError, ValueError):
            pass
    return out


# ---------------------------------------------------------------------------
# electron depth
# ---------------------------------------------------------------------------
def theta_bounds(n_theta: int) -> np.ndarray:
    """Polar-angle bin boundaries, uniform in cos(theta): pi ... 0 [rad]."""
    return np.arccos(np.linspace(-1.0, 1.0, n_theta + 1))


def median_depth_grid(
    rad_edges: np.ndarray,
    depth_grid: np.ndarray,
    depth_dist: np.ndarray,
    i_dim: int,
    i_tau: int = 2,
    median_level: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Median electron depth on the native (radius, polar angle) grid.

    Port of ``compute_medianDepth_grid``.

    Parameters
    ----------
    depth_dist : (n_depth, n_xs, n_r, n_theta, 3)
    i_dim : {0, 1, 2}
        Polar axis (x, y, z) -- MATLAB ``iDim`` minus one.
    i_tau : int
        Index along the second axis of ``depth_dist`` (MATLAB ``iTau = 3``
        is ``i_tau = 2``).
    median_level : float
        Quantile of the depth distribution to return (0.5 = median).

    Returns
    -------
    median_depth : (n_r, n_theta)
    r_mids : (n_r,) [km]
    theta_mids : (n_theta,) [rad], running from ~pi down to ~0
    """
    rad_edges = np.asarray(rad_edges, dtype=float)
    r_mids = 0.5 * (rad_edges[:-1] + rad_edges[1:])
    n_theta = depth_dist.shape[3]
    tb = theta_bounds(n_theta)
    theta_mids = 0.5 * (tb[:-1] + tb[1:])

    slab = np.asarray(depth_dist[1:, i_tau, :, :, i_dim], dtype=float)  # (n_depth-1, n_r, n_theta)
    cum = np.cumsum(slab, axis=0)
    denom = cum[-1]
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = cum / np.where(denom == 0, np.nan, denom)
    reached = frac >= median_level
    first = np.argmax(reached, axis=0)
    has = reached.any(axis=0)
    depth_grid = np.asarray(depth_grid, dtype=float)
    lookup = np.clip(first + 1, 0, depth_grid.size - 1)
    median_depth = np.where(has, depth_grid[lookup], np.nan)
    return median_depth, r_mids, theta_mids


def shell_median_profile(median_depth: np.ndarray, r_mids: np.ndarray, r0: float, dr: float):
    """Median over the radial bins of a shell, per polar angle.

    Port of ``median_shell_profile_from_native_theta``.
    """
    mask = (r_mids >= r0 - 0.5 * dr) & (r_mids < r0 + 0.5 * dr)
    if not mask.any():
        return np.full(median_depth.shape[1], np.nan)
    vals = np.where(median_depth[mask] > 0, median_depth[mask], np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN columns
        return np.nanmedian(vals, axis=0)


def exobase_hull(
    x: np.ndarray,
    y: np.ndarray,
    in_exo: np.ndarray,
    smooth: bool = True,
    smooth_method: str = "sgolay",
    smooth_span: int = 14,
    sgolay_order: int = 2,
    resample_n: int = 50,
) -> Optional[dict]:
    """Convex hull of the collisional region, optionally smoothed.

    Port of ``compute_smoothed_exobase_hull_xy`` /
    ``compute_raw_exobase_hull_xy``.  Returns ``None`` if fewer than three
    non-collinear points are flagged, otherwise a dict with closed curves
    ``x_raw, y_raw`` (and ``x_smooth, y_smooth`` when ``smooth``).
    """
    xs, ys = np.asarray(x)[in_exo], np.asarray(y)[in_exo]
    good = np.isfinite(xs) & np.isfinite(ys)
    pts = np.column_stack([xs[good], ys[good]])
    if pts.shape[0] < 3 or np.linalg.matrix_rank(pts - pts.mean(axis=0)) < 2:
        return None
    try:
        hull = spatial.ConvexHull(pts)
    except spatial.QhullError:
        return None
    k = np.r_[hull.vertices, hull.vertices[0]]  # closed, counter-clockwise
    xr, yr = pts[k, 0], pts[k, 1]
    out = {"x_raw": xr, "y_raw": yr}
    if not smooth:
        return out

    s = np.r_[0.0, np.cumsum(np.hypot(np.diff(xr), np.diff(yr)))]
    length = s[-1]
    if not np.isfinite(length) or length == 0:
        return None
    n = max(20, int(round(resample_n)))
    s_u = np.linspace(0.0, length, n + 1)[:-1]
    x_u, y_u = np.interp(s_u, s, xr), np.interp(s_u, s, yr)

    span = max(3, int(round(smooth_span)))
    span += span % 2 == 0
    half = span // 2
    x_w = np.r_[x_u[-half:], x_u, x_u[:half]]
    y_w = np.r_[y_u[-half:], y_u, y_u[:half]]
    if smooth_method.lower() == "sgolay":
        order = min(int(round(sgolay_order)), span - 2)
        x_s = signal.savgol_filter(x_w, span, order)
        y_s = signal.savgol_filter(y_w, span, order)
    else:
        x_s, y_s = movmean(x_w, span), movmean(y_w, span)
    x_s, y_s = x_s[half : half + n], y_s[half : half + n]
    out.update(
        x_resamp=np.r_[x_u, x_u[0]], y_resamp=np.r_[y_u, y_u[0]],
        x_smooth=np.r_[x_s, x_s[0]], y_smooth=np.r_[y_s, y_s[0]],
    )
    return out


def exobase_envelope(r_mids, theta_mids, in_exo, smooth: bool = True, smooth_span: int = 14) -> dict:
    """Outermost collisional radius for each polar angle.

    Port of ``compute_exobase_envelope_rt``.
    """
    r_mids = np.asarray(r_mids, dtype=float)
    theta_mids = np.asarray(theta_mids, dtype=float)
    in_exo = np.asarray(in_exo, dtype=bool)
    last = in_exo.shape[0] - 1 - np.argmax(in_exo[::-1], axis=0)
    r_exo = np.where(in_exo.any(axis=0), r_mids[last], np.nan)
    if smooth:
        span = max(3, int(round(smooth_span)))
        span += span % 2 == 0
        r_s = movmean(r_exo, span, omitnan=True)
    else:
        r_s = r_exo
    return {
        "theta": theta_mids, "r": r_exo, "r_smooth": r_s,
        "x": r_s * np.cos(theta_mids), "y": r_s * np.sin(theta_mids),
    }


def median_depth_and_exobase(
    rad_edges,
    depth_grid,
    depth_dist,
    i_dim: int,
    i_tau: int = 2,
    median_level: float = 0.5,
    exo_thresh_log10: float = 0.0,
    pad_theta: bool = True,
    smooth_hull: bool = True,
    **hull_kwargs,
) -> dict:
    """Median depth map plus exobase hull and envelope.

    Port of ``asymPlots.compute_2dElDepthMedianAndExobase``.  The collisional
    region is where ``log10(median depth) > exo_thresh_log10``.  With
    ``pad_theta`` the angular grid is extended to exactly 0 and pi by
    duplicating the end columns.
    """
    median, r_mids, theta_mids = median_depth_grid(
        rad_edges, depth_grid, depth_dist, i_dim, i_tau, median_level
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        log_median = np.log10(median)
    in_exo = log_median > exo_thresh_log10

    out = {
        "median_depth_unpad": median, "log10_median_depth_unpad": log_median,
        "in_exo_unpad": in_exo, "theta_mids_unpad": theta_mids,
    }
    if pad_theta:
        theta_mids = np.r_[np.pi, theta_mids, 0.0]
        median = np.column_stack([median[:, 0], median, median[:, -1]])
        log_median = np.column_stack([log_median[:, 0], log_median, log_median[:, -1]])
        in_exo = np.column_stack([in_exo[:, 0], in_exo, in_exo[:, -1]])

    x = r_mids[:, None] * np.cos(theta_mids)[None, :]
    y = r_mids[:, None] * np.sin(theta_mids)[None, :]
    exo = exobase_hull(x, y, in_exo, smooth=smooth_hull, **hull_kwargs)
    exo = {} if exo is None else exo
    exo["envelope"] = exobase_envelope(
        r_mids, theta_mids, in_exo, smooth_span=hull_kwargs.get("smooth_span", 14)
    )
    out.update(
        median_depth=median, log10_median_depth=log_median, in_exo=in_exo,
        r_mids=r_mids, theta_mids=theta_mids, x=x, y=y, exo=exo,
    )
    return out
