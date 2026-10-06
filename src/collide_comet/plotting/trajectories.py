"""Test-particle trajectory figures (port of the 3D plot in ``plot_trajectories_v3.m``)."""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from .. import constants as c

__all__ = ["kinetic_energy_eV", "plot_trajectory_3d"]


def kinetic_energy_eV(V: np.ndarray, me: float = c.ME) -> np.ndarray:
    """Kinetic energy [eV] for velocities ``V`` [m/s], shape ``(n, 3)``."""
    return 0.5 * me * np.sum(np.asarray(V, dtype=float) ** 2, axis=1) / c.QE


def plot_trajectory_3d(
    X: np.ndarray,
    V: Optional[np.ndarray] = None,
    color_by: Optional[np.ndarray] = None,
    color_label: str = "Total Energy [eV]",
    origin: Sequence[float] = (0.0, 0.0, 0.0),
    scale: float = 1.0,
    comet_radius: Optional[float] = None,
    cmap: str = "jet",
    view: Tuple[float, float] = (10.0, -135.0),
    ax=None,
    fontsize: float = 12,
):
    """3D trajectory coloured by energy (or any per-step quantity).

    Parameters
    ----------
    X : (n_steps, 3) positions, as returned by
        :func:`collide_comet.io.load_trajectories`.
    V : (n_steps, 3) velocities [m/s]; used to colour by kinetic energy when
        ``color_by`` is not given.
    origin, scale : plotted position is ``(X - origin) * scale``.  The MATLAB
        script uses the comet's grid index as origin and the cell size in km
        as scale (e.g. ``origin=(134, 113, 113), scale=10``).
    comet_radius : draw a sphere of this radius [km] at the origin.

    Returns ``fig, ax, collection``.
    """
    X = (np.asarray(X, dtype=float) - np.asarray(origin, dtype=float)) * scale
    if color_by is None:
        if V is None:
            color_by, color_label = np.arange(X.shape[0], dtype=float), "Step"
        else:
            color_by = kinetic_energy_eV(V)
    color_by = np.asarray(color_by, dtype=float)

    if ax is None:
        fig = plt.figure(figsize=(7.0, 6.5))
        ax = fig.add_subplot(projection="3d")
    else:
        fig = ax.figure
    segments = np.stack([X[:-1], X[1:]], axis=1)
    lines = Line3DCollection(segments, cmap=cmap, linewidth=0.9)
    lines.set_array(0.5 * (color_by[:-1] + color_by[1:]))
    ax.add_collection3d(lines)
    cbar = fig.colorbar(lines, ax=ax, shrink=0.7, pad=0.1)
    cbar.set_label(color_label, fontsize=fontsize)

    if comet_radius:
        u, v = np.mgrid[0 : 2 * np.pi : 30j, 0 : np.pi : 15j]
        ax.plot_surface(comet_radius * np.cos(u) * np.sin(v), comet_radius * np.sin(u) * np.sin(v),
                        comet_radius * np.cos(v), color="0.5", alpha=0.6, linewidth=0)
    lo, hi = X.min(axis=0), X.max(axis=0)
    centre, half = 0.5 * (lo + hi), 0.5 * (hi - lo).max() or 1.0
    ax.set_xlim(centre[0] - half, centre[0] + half)
    ax.set_ylim(centre[1] - half, centre[1] + half)
    ax.set_zlim(centre[2] - half, centre[2] + half)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("X [km]", fontsize=fontsize)
    ax.set_ylabel("Y [km]", fontsize=fontsize)
    ax.set_zlabel("Z [km]", fontsize=fontsize)
    ax.view_init(*view)
    return fig, ax, lines
