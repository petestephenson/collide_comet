"""Collision cross-section figures.

Ports of the plots in ``inelasticXSplots.m`` that only need the collision
table (``Coll_XS.nc``): cross sections per process and energy-loss rates.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .. import constants as c
from ..io import CrossSections, load_cross_sections

__all__ = ["loss_rates", "plot_cross_sections", "plot_loss_rates"]

_IX_DOUBLE_IONI = 11  # O++ (MATLAB ix_ioni2 = 12)


def _names(xs: CrossSections, names: Optional[Sequence[str]]) -> Sequence[str]:
    if names is not None:
        return names
    if xs.n_process == len(c.PROCESS_NAMES):
        return c.PROCESS_NAMES
    return [f"Process {i + 1}" for i in range(xs.n_process)]


def loss_rates(xs: CrossSections) -> np.ndarray:
    """Energy-loss cross sections [eV cm^2], shape ``(n_energy, n_process)``.

    Port of the ``LR`` array in ``inelasticXSplots.m`` for the 12-process
    water table: excitations lose their threshold energy, single ionisation
    ``(E_th + E) / 2`` and double ionisation ``(E_th + 2E) / 3``.
    """
    if xs.n_process != 12:
        raise ValueError("loss_rates expects the 12-process water table (Coll_XS.nc)")
    E = xs.energy[:, None]
    lr = np.zeros_like(xs.xs)
    exc = slice(0, c.IX_IONI)
    ioni = slice(c.IX_IONI, _IX_DOUBLE_IONI)
    lr[:, exc] = xs.xs[:, exc] * xs.e_th[exc]
    lr[:, ioni] = 0.5 * xs.xs[:, ioni] * (xs.e_th[ioni] + E)
    lr[:, _IX_DOUBLE_IONI] = xs.xs[:, _IX_DOUBLE_IONI] * (xs.e_th[_IX_DOUBLE_IONI] + 2 * xs.energy) / 3
    return lr


def plot_cross_sections(
    xs: Optional[CrossSections] = None,
    names: Optional[Sequence[str]] = None,
    figsize: Tuple[float, float] = (11.0, 4.8),
    fontsize: float = 12,
):
    """Excitation (left) and ionisation (right) cross sections versus energy.

    Returns ``fig, axes``.
    """
    xs = load_cross_sections() if xs is None else xs
    names = _names(xs, names)
    split = min(c.IX_IONI, xs.n_process)
    groups = [range(0, split), range(split, xs.n_process)]
    limits = [((1e-2, 3e3), (1e-18, 1e-12)), ((10, 3e2), (1e-20, 3e-16))]
    titles = ("Elastic and excitation", "Ionisation")
    groups = [g for g in groups if len(g)]

    fig, axes = plt.subplots(1, len(groups), figsize=figsize, constrained_layout=True, squeeze=False)
    axes = axes[0]
    for ax, group, (xlim, ylim), title in zip(axes, groups, limits, titles):
        for i in group:
            positive = xs.xs[:, i] > 0
            ax.plot(xs.energy[positive], xs.xs[positive, i], linewidth=2.2, label=names[i])
        ax.set_xscale("log")
        ax.set_yscale("log")
        if xs.n_process == 12:
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
        ax.set_xlabel("Energy [eV]", fontsize=fontsize)
        ax.set_ylabel(r"Cross section [cm$^2$]", fontsize=fontsize)
        ax.set_title(title, fontsize=fontsize)
        ax.grid(True, which="major", alpha=0.3)
        ax.tick_params(labelsize=fontsize)
        ax.legend(fontsize=fontsize - 2)
    return fig, axes


def plot_loss_rates(
    xs: Optional[CrossSections] = None,
    figsize: Tuple[float, float] = (8.0, 4.8),
    fontsize: float = 12,
):
    """Energy-loss cross sections, with all ionisation channels combined.

    Returns ``fig, ax``.
    """
    xs = load_cross_sections() if xs is None else xs
    lr = loss_rates(xs)
    labels = ("Elastic", "Rotation", "Vib 010", "Vib Stretch", "Electronic", "Dissociation")
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    for i, label in enumerate(labels):
        positive = lr[:, i] > 0
        if not positive.any():  # e.g. elastic collisions: no threshold energy
            continue
        ax.plot(xs.energy[positive], lr[positive, i], linewidth=2.2, label=label)
    ioni = lr[:, c.IX_IONI:].sum(axis=1)
    ax.plot(xs.energy[ioni > 0], ioni[ioni > 0], linewidth=2.2, label="Ionization")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-2, 3e2)
    ax.set_ylim(5e-19, 3e-14)
    ax.set_xlabel("Energy [eV]", fontsize=fontsize)
    ax.set_ylabel(r"Loss rate [eV cm$^2$]", fontsize=fontsize)
    ax.grid(True, which="major", alpha=0.3)
    ax.tick_params(labelsize=fontsize)
    ax.legend(loc="lower left", fontsize=fontsize - 2)
    return fig, ax
