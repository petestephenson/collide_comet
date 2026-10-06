"""2D cuts through 3D quantities, single and grouped by case.

Ports of ``asymPlots.makePlot`` and its many thin wrappers
(``makePlasDensPlot``, ``makeTemperaturePlot``, ``makeEIIFreqPlot`` ...), and
of the grouped five-panel figures (``make_totalelectrons_Plot``,
``make_grouped_EIIFreqPlot``, ``make_*_Ratio_Plot`` ...).

MATLAB drew a 3D ``slice`` and looked at it along an axis; here the plane is
drawn directly with ``pcolormesh`` on the true cell boundaries.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np

from ..cases import CaseData, CaseSet
from ..colormaps import get_cmap
from ..io import Grids
from ..quantities import Quantity, compute, get_quantity

__all__ = [
    "VIEWS",
    "slice_index",
    "plot_slice",
    "plot_quantity",
    "grouped_slice_plot",
    "grouped_ratio_plot",
]

#: view name -> (normal axis, horizontal axis, vertical axis).  The order
#: matches MATLAB ``vangleList = {[1,0,0], [0,-1,0], [0,0,1]}``.
VIEWS: Dict[str, Tuple[int, int, int]] = {"YZ": (0, 1, 2), "XZ": (1, 0, 2), "XY": (2, 0, 1)}
_AXIS_NAMES = ("X", "Y", "Z")

# Panel positions of the MATLAB figures (left, bottom) in figure fractions.
_GROUP_LEFT = (0.35, 0.65, 0.35, 0.65)
_GROUP_BOTTOM = (0.57, 0.57, 0.10, 0.10)
_GROUP_REF = (0.07, 0.57)
_GROUP_SIZE = (0.23, 0.35)
_RATIO_LEFT = (0.10, 0.50, 0.10, 0.50)
_RATIO_BOTTOM = (0.60, 0.60, 0.15, 0.15)
_RATIO_SIZE = (0.36, 0.35)
_CBAR_POS = (0.90, 0.10, 0.01, 0.80)


def _view(view) -> Tuple[str, Tuple[int, int, int]]:
    if isinstance(view, (int, np.integer)):
        name = list(VIEWS)[int(view)]
    else:
        name = str(view).upper()
        if name not in VIEWS:
            name = name[::-1]
    if name not in VIEWS:
        raise ValueError(f"view must be one of {list(VIEWS)} (got {view!r})")
    return name, VIEWS[name]


def slice_index(edges: np.ndarray, position: float) -> int:
    """Index of the cell containing ``position`` (clamped to the grid)."""
    idx = int(np.searchsorted(edges, position, side="right") - 1)
    return min(max(idx, 0), edges.size - 2)


def _to_plot_values(values: np.ndarray, log: bool) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if not log:
        return np.ma.masked_invalid(values)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.ma.masked_invalid(np.log10(np.abs(values)))


def _log_label(label: str, log: bool) -> str:
    return rf"$\log_{{10}}$({label})" if log and label else label


def plot_slice(
    values: np.ndarray,
    grids: Grids,
    view: str = "XZ",
    slice_at: float = 0.0,
    log: bool = False,
    clim: Optional[Sequence[float]] = None,
    cmap=None,
    label: str = "",
    plot_lim: Optional[float] = 300.0,
    ax: Optional[plt.Axes] = None,
    colorbar: bool = True,
    fontsize: float = 12,
):
    """Draw a planar cut through a 3D array.

    Port of ``asymPlots.makePlot``.

    Parameters
    ----------
    values : (nx, ny, nz) array
    grids : Grids
    view : {"YZ", "XZ", "XY"}
        Plane to show; the cut is taken at ``slice_at`` [km] along the third
        axis.  (MATLAB ``iAngle`` 1, 2, 3.)
    log : bool
        Plot ``log10(|values|)``.
    clim : (low, high), optional
        Colour limits (of the log if ``log``).
    cmap : colour map or name; default the MATLAB ``redblue`` map.
    label : colour-bar label; wrapped in ``log10(...)`` when ``log``.
    plot_lim : float or None
        Axis limits ``[-plot_lim, plot_lim]`` km; ``None`` for the full grid.

    Returns
    -------
    fig, ax, mesh, cbar
    """
    name, (normal, h, v) = _view(view)
    values = np.asarray(values)
    if values.shape != grids.shape:
        raise ValueError(f"values has shape {values.shape}; grid has {grids.shape} cells")

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.6, 4.6))
    else:
        fig = ax.figure

    edges = grids.edges
    plane = np.take(values, slice_index(edges[normal], slice_at), axis=normal)
    # after removing `normal`, the remaining axes are in increasing order;
    # pcolormesh wants C[vertical, horizontal]
    if h < v:
        plane = plane.T
    mesh = ax.pcolormesh(
        edges[h], edges[v], _to_plot_values(plane, log), cmap=get_cmap(cmap), shading="flat",
        rasterized=True,
    )
    if clim is not None:
        mesh.set_clim(*clim)

    ax.set_xlabel(f"{_AXIS_NAMES[h]} [km]", fontsize=fontsize)
    ax.set_ylabel(f"{_AXIS_NAMES[v]} [km]", fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax.set_aspect("equal")
    if plot_lim is not None:
        ax.set_xlim(-plot_lim, plot_lim)
        ax.set_ylim(-plot_lim, plot_lim)

    cbar = None
    if colorbar:
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label(_log_label(label, log), fontsize=fontsize)
    return fig, ax, mesh, cbar


def plot_quantity(
    case: CaseData,
    quantity: Union[str, Quantity] = "total_density",
    view: str = "XZ",
    log: Optional[bool] = None,
    clim: Optional[Sequence[float]] = None,
    **kwargs,
):
    """Slice plot of a registered quantity for one case.

    ``quantity`` is a name from :data:`collide_comet.quantities.QUANTITIES`
    (e.g. ``"total_density"``, ``"cold_density"``, ``"eii_frequency"``,
    ``"temperature"``).  Label, log scaling and colour limits default to the
    values used by the corresponding MATLAB function.  Remaining keyword
    arguments go to :func:`plot_slice`.
    """
    q = get_quantity(quantity)
    log = q.log if log is None else log
    if clim is None and log == q.log:
        clim = q.clim
    kwargs.setdefault("label", q.label)
    return plot_slice(compute(case, q), case.grids, view=view, log=log, clim=clim, **kwargs)


# ---------------------------------------------------------------------------
# grouped figures
# ---------------------------------------------------------------------------
def _grouped_axes(fig, cases: Sequence[CaseData], ratio: bool):
    """Axes laid out as in the MATLAB grouped figures, one per case.

    Returns the axes (dict by case name) and the names on the bottom row.
    Axis labels are kept on the bottom row and the first panel of each row
    only, so they do not run into the neighbouring panels.
    """
    asym = [case for case in cases if not case.is_reference]
    refs = [case for case in cases if case.is_reference]
    axes = {}
    if len(asym) <= 4 and len(refs) <= 1:
        lefts, bottoms, size = (
            (_RATIO_LEFT, _RATIO_BOTTOM, _RATIO_SIZE) if ratio
            else (_GROUP_LEFT, _GROUP_BOTTOM, _GROUP_SIZE)
        )
        for i, case in enumerate(asym):
            axes[case.name] = fig.add_axes([lefts[i], bottoms[i], *size])
        for case in refs:
            axes[case.name] = fig.add_axes([*_GROUP_REF, *_GROUP_SIZE])
        bottom_row = {case.name for i, case in enumerate(asym) if i >= 2} or {c.name for c in asym}
        bottom_row |= {case.name for case in refs}
        if ratio or not refs:
            left_col = {case.name for i, case in enumerate(asym) if i % 2 == 0}
        else:  # the reference panel starts the top row
            left_col = {case.name for i, case in enumerate(asym) if i == 2} | {c.name for c in refs}
    else:  # generic fallback for other numbers of cases
        n = len(cases)
        ncols = int(np.ceil(np.sqrt(n)))
        nrows = int(np.ceil(n / ncols))
        gs = fig.add_gridspec(nrows, ncols, left=0.07, right=0.86, bottom=0.1, top=0.92,
                              wspace=0.35, hspace=0.35)
        for i, case in enumerate(cases):
            axes[case.name] = fig.add_subplot(gs[i // ncols, i % ncols])
        bottom_row = {case.name for i, case in enumerate(cases) if i + ncols >= n}
        left_col = {case.name for i, case in enumerate(cases) if i % ncols == 0}
    for name, ax in axes.items():
        ax._collide_label_x = name in bottom_row
        ax._collide_label_y = name in left_col
    return axes, bottom_row


def _finish_group(fig, axes, bottom_row, mesh, label, titles, fontsize):
    for name, ax in axes.items():
        if name not in bottom_row:
            ax.set_xlabel("")
        if not getattr(ax, "_collide_label_y", True):
            ax.set_ylabel("")
        if titles:
            ax.set_title(name, fontsize=fontsize)
    cax = fig.add_axes(_CBAR_POS)
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.set_label(label, fontsize=fontsize)
    return cbar


def grouped_slice_plot(
    cases: CaseSet,
    quantity: Union[str, Quantity] = "total_density",
    view: str = "XZ",
    log: Optional[bool] = None,
    clim: Optional[Sequence[float]] = None,
    cmap=None,
    plot_lim: Optional[float] = 300.0,
    slice_at: float = 0.0,
    titles: bool = True,
    figsize: Tuple[float, float] = (10.0, 6.0),
    fontsize: float = 12,
):
    """One panel per case with a shared colour bar.

    Replaces the ``make_*_Plot`` / ``make_grouped_*Plot`` family: the
    reference case sits top-left and the asymmetric cases form a 2 x 2 block
    on the right, as in the MATLAB figures.

    Parameters
    ----------
    cases : CaseSet
    quantity : str or Quantity
        See :func:`collide_comet.list_quantities`.
    view : {"YZ", "XZ", "XY"}
    log, clim : override the quantity's defaults.

    Returns
    -------
    fig, axes, cbar
        ``axes`` is a dict keyed by case name.
    """
    q = get_quantity(quantity)
    use_log = q.log if log is None else log
    if clim is None and use_log == q.log:
        clim = q.clim

    fig = plt.figure(figsize=figsize)
    case_list = list(cases)
    axes, bottom_row = _grouped_axes(fig, case_list, ratio=False)
    mesh = None
    for case in case_list:
        _, _, mesh, _ = plot_slice(
            compute(case, q), case.grids, view=view, slice_at=slice_at, log=use_log, clim=clim,
            cmap=cmap, plot_lim=plot_lim, ax=axes[case.name], colorbar=False, fontsize=fontsize,
        )
    cbar = _finish_group(fig, axes, bottom_row, mesh, _log_label(q.label, use_log), titles, fontsize)
    return fig, axes, cbar


def grouped_ratio_plot(
    cases: CaseSet,
    quantity: Union[str, Quantity] = "total_density",
    view: str = "XZ",
    log: Optional[bool] = None,
    clim: Optional[Sequence[float]] = None,
    cmap=None,
    plot_lim: Optional[float] = 300.0,
    slice_at: float = 0.0,
    titles: bool = True,
    figsize: Tuple[float, float] = (10.0, 6.0),
    fontsize: float = 12,
):
    """Ratio of each asymmetric case to the reference case (2 x 2 panels).

    Replaces the ``make_*_Ratio_Plot`` family.  Returns ``fig, axes, cbar``.
    """
    q = get_quantity(quantity)
    use_log = q.ratio_log if log is None else log
    if clim is None and use_log == q.ratio_log:
        clim = q.ratio_clim

    ref = cases.reference
    ref_values = compute(ref, q)
    fig = plt.figure(figsize=figsize)
    asym = cases.asymmetric
    axes, bottom_row = _grouped_axes(fig, asym, ratio=True)
    mesh = None
    for case in asym:
        values = compute(case, q)
        if values.shape != ref_values.shape:
            raise ValueError(f"Case {case.name} and the reference are on different grids")
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = values / ref_values
        _, _, mesh, _ = plot_slice(
            ratio, case.grids, view=view, slice_at=slice_at, log=use_log, clim=clim, cmap=cmap,
            plot_lim=plot_lim, ax=axes[case.name], colorbar=False, fontsize=fontsize,
        )
    label = q.ratio_label or f"{q.label} ratio"
    cbar = _finish_group(fig, axes, bottom_row, mesh, _log_label(label, use_log), titles, fontsize)
    return fig, axes, cbar
