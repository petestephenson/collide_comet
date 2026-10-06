"""Electron energy distribution along a line through the comet.

Ports of ``asymPlots.make_edf_plot``, ``make_edfRatioplot`` and the grouped
``make_SW_edf_plot`` / ``make_SW_edfRatioplot``.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .. import constants as c
from ..cases import CaseSet
from ..io import Grids
from ..physics import smooth_plas_dens_E
from .slices import _AXIS_NAMES, _finish_group, _grouped_axes, slice_index

__all__ = ["edf_line", "plot_edf", "plot_edf_ratio", "grouped_edf_plot", "grouped_edf_ratio_plot"]

_EDF_LABEL = r"$\log_{10}$($f(E)$ [cm$^{-3}$ eV$^{-1}$])"
_EDF_RATIO_LABEL = r"$\log_{10}$($f_\mathrm{asym}(E) / f_\mathrm{sym}(E)$)"


def _axis_index(axis) -> int:
    if isinstance(axis, str):
        return _AXIS_NAMES.index(axis.upper())
    return int(axis)


def edf_line(
    plas_dens_E: np.ndarray, grids: Grids, axis=0, position: Sequence[float] = (0.0, 0.0, 0.0),
    smooth: bool = True,
) -> np.ndarray:
    """Energy distribution along a line parallel to ``axis``.

    Parameters
    ----------
    plas_dens_E : (nx, ny, nz, nE) -- one generation (or a sum over several)
    axis : 0/1/2 or "X"/"Y"/"Z"
    position : (x, y, z) [km] the line passes through (the entry for
        ``axis`` is ignored).

    Returns an array of shape ``(n_along, nE)``.
    """
    axis = _axis_index(axis)
    if plas_dens_E.ndim != 4:
        raise ValueError("plas_dens_E must have shape (nx, ny, nz, nE); select or sum generations first")
    f = smooth_plas_dens_E(plas_dens_E) if smooth else np.asarray(plas_dens_E, dtype=float)
    index = [slice_index(e, p) for e, p in zip(grids.edges, position)]
    index[axis] = slice(None)
    return f[tuple(index)]


def _draw_edf(ax, grids, axis, log_values, clim, cmap, fontsize):
    mesh = ax.pcolormesh(
        grids.edges[axis], grids.energy, np.ma.masked_invalid(log_values).T, cmap=cmap,
        shading="flat", rasterized=True,
    )
    if clim is not None:
        mesh.set_clim(*clim)
    ax.set_xlabel(f"{_AXIS_NAMES[axis]} [km]", fontsize=fontsize)
    ax.set_ylabel("Energy [eV]", fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    return mesh


def plot_edf(
    plas_dens_E: np.ndarray,
    grids: Grids,
    axis=0,
    position: Sequence[float] = (0.0, 0.0, 0.0),
    clim: Optional[Sequence[float]] = (-4.0, 1.0),
    cmap="jet",
    ax: Optional[plt.Axes] = None,
    colorbar: bool = True,
    fontsize: float = 12,
):
    """Energy spectrogram ``log10 f(E)`` along a line.  Port of ``make_edf_plot``.

    Returns ``fig, ax, mesh, cbar``.
    """
    axis = _axis_index(axis)
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.0, 4.5))
    else:
        fig = ax.figure
    with np.errstate(divide="ignore", invalid="ignore"):
        values = np.log10(edf_line(plas_dens_E, grids, axis, position))
    mesh = _draw_edf(ax, grids, axis, values, clim, cmap, fontsize)
    cbar = None
    if colorbar:
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label(_EDF_LABEL, fontsize=fontsize)
    return fig, ax, mesh, cbar


def plot_edf_ratio(
    plas_dens_E_ref: np.ndarray,
    plas_dens_E_asym: np.ndarray,
    grids: Grids,
    axis=0,
    position: Sequence[float] = (0.0, 0.0, 0.0),
    clim: Optional[Sequence[float]] = (-1.0, 1.0),
    cmap="jet",
    ax: Optional[plt.Axes] = None,
    colorbar: bool = True,
    fontsize: float = 12,
):
    """``log10`` ratio of two energy distributions along a line.

    Port of ``make_edfRatioplot``.  Returns ``fig, ax, mesh, cbar``.
    """
    axis = _axis_index(axis)
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.0, 4.5))
    else:
        fig = ax.figure
    with np.errstate(divide="ignore", invalid="ignore"):
        values = np.log10(
            edf_line(plas_dens_E_asym, grids, axis, position)
            / edf_line(plas_dens_E_ref, grids, axis, position)
        )
    mesh = _draw_edf(ax, grids, axis, values, clim, cmap, fontsize)
    cbar = None
    if colorbar:
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label(_EDF_RATIO_LABEL, fontsize=fontsize)
    return fig, ax, mesh, cbar


def _edf_source(case, source: str, gen: Optional[int]) -> np.ndarray:
    f = case.plasma_density_E(source)
    return f.sum(axis=-1) if gen is None else f[..., gen]


def grouped_edf_plot(
    cases: CaseSet,
    axis=0,
    source: str = "sw",
    gen: Optional[int] = c.GEN_PRIMARY,
    clim: Optional[Sequence[float]] = (-4.0, 1.0),
    cmap="jet",
    plot_lim: float = 300.0,
    energy_lim: Tuple[float, float] = (0.0, 350.0),
    titles: bool = True,
    figsize: Tuple[float, float] = (10.0, 6.0),
    fontsize: float = 12,
):
    """Energy spectrograms for every case.  Port of ``make_SW_edf_plot``.

    The defaults (``source="sw"``, primary generation) reproduce the MATLAB
    figure of the solar-wind electron distribution; ``gen=None`` sums the
    generations and ``source`` may be ``"pe"`` or ``"total"``.

    Returns ``fig, axes, cbar``.
    """
    axis = _axis_index(axis)
    fig = plt.figure(figsize=figsize)
    case_list = list(cases)
    axes, bottom_row = _grouped_axes(fig, case_list, ratio=False)
    mesh = None
    for case in case_list:
        ax = axes[case.name]
        _, _, mesh, _ = plot_edf(
            _edf_source(case, source, gen), case.grids, axis, clim=clim, cmap=cmap, ax=ax,
            colorbar=False, fontsize=fontsize,
        )
        ax.set_xlim(-plot_lim, plot_lim)
        ax.set_ylim(*energy_lim)
    cbar = _finish_group(fig, axes, bottom_row, mesh, _EDF_LABEL, titles, fontsize)
    return fig, axes, cbar


def grouped_edf_ratio_plot(
    cases: CaseSet,
    axis=0,
    source: str = "sw",
    gen: Optional[int] = c.GEN_PRIMARY,
    clim: Optional[Sequence[float]] = (-1.0, 1.0),
    cmap="jet",
    plot_lim: Optional[float] = None,
    energy_lim: Optional[Tuple[float, float]] = None,
    titles: bool = True,
    figsize: Tuple[float, float] = (10.0, 6.0),
    fontsize: float = 12,
):
    """Asymmetric / reference energy-distribution ratios.

    Port of ``make_SW_edfRatioplot``.  Returns ``fig, axes, cbar``.
    """
    axis = _axis_index(axis)
    ref = _edf_source(cases.reference, source, gen)
    fig = plt.figure(figsize=figsize)
    asym = cases.asymmetric
    axes, bottom_row = _grouped_axes(fig, asym, ratio=True)
    mesh = None
    for case in asym:
        ax = axes[case.name]
        _, _, mesh, _ = plot_edf_ratio(
            ref, _edf_source(case, source, gen), case.grids, axis, clim=clim, cmap=cmap, ax=ax,
            colorbar=False, fontsize=fontsize,
        )
        if plot_lim is not None:
            ax.set_xlim(-plot_lim, plot_lim)
        if energy_lim is not None:
            ax.set_ylim(*energy_lim)
    cbar = _finish_group(fig, axes, bottom_row, mesh, _EDF_RATIO_LABEL, titles, fontsize)
    return fig, axes, cbar
