"""collide_comet: read and plot outputs of the collisional test-particle model.

Typical use::

    import collide_comet as cc

    cases = cc.CaseSet("E:/TestPlSimData", pe_runs=[...], sw_runs=[...])
    fig, axes, cbar = cc.plotting.grouped_slice_plot(cases, "cold_density", view="XZ")
    cc.save_figure(fig, "cold_density_XZ")

Sub-modules
-----------
``io``         netCDF readers (run files, cross sections, fields, trajectories)
``physics``    derived quantities (densities, ionisation rates, temperature ...)
``analysis``   profiles, Monte-Carlo sampling, electron-depth analysis
``cases``      :class:`CaseSet` / :class:`CaseData` -- groups of PE + SW runs
``quantities`` registry of named, mappable 3D quantities
``plotting``   figures
``figures``    saving to the ``Figures`` folder
``synthetic``  small made-up run files for demos and tests
"""

from . import analysis, constants, io, physics, plotting, quantities, synthetic
from .cases import CaseData, CaseSet
from .figures import figure_dir, save_figure, set_figure_dir
from .io import (
    CrossSections,
    Grids,
    Outgassing,
    describe,
    load_cross_sections,
    load_electron_depth_dist,
    load_outgassing,
    load_run_data,
    load_velocity_arrays,
)
from .quantities import QUANTITIES, list_quantities

__version__ = "0.1.0"

__all__ = [
    "analysis",
    "constants",
    "io",
    "physics",
    "plotting",
    "quantities",
    "synthetic",
    "CaseData",
    "CaseSet",
    "CrossSections",
    "Grids",
    "Outgassing",
    "QUANTITIES",
    "describe",
    "figure_dir",
    "list_quantities",
    "load_cross_sections",
    "load_electron_depth_dist",
    "load_outgassing",
    "load_run_data",
    "load_velocity_arrays",
    "save_figure",
    "set_figure_dir",
]
