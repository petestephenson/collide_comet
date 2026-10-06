"""Plotting functions.  Everything here is also available as ``collide_comet.plotting.<name>``."""

from .depth import (
    el_depth_angle_maps_plot,
    el_depth_distribution_plot,
    el_depth_summary_plot,
    grouped_el_depth_2d_plot,
    plot_el_depth_2d,
)
from .edf import grouped_edf_plot, grouped_edf_ratio_plot, plot_edf, plot_edf_ratio
from .profiles import (
    angle_shell_profiles_plot,
    axis_cone_profiles_plot,
    cold_production_histogram,
    cold_production_summary_plot,
    r_angle_maps_plot,
    radial_profiles_plot,
    temperature_vs_radius_plot,
    timestep_histogram,
)
from .slices import VIEWS, grouped_ratio_plot, grouped_slice_plot, plot_quantity, plot_slice
from .trajectories import plot_trajectory_3d
from .xsections import plot_cross_sections, plot_loss_rates

__all__ = [
    "VIEWS",
    "plot_slice",
    "plot_quantity",
    "grouped_slice_plot",
    "grouped_ratio_plot",
    "plot_edf",
    "plot_edf_ratio",
    "grouped_edf_plot",
    "grouped_edf_ratio_plot",
    "radial_profiles_plot",
    "angle_shell_profiles_plot",
    "axis_cone_profiles_plot",
    "r_angle_maps_plot",
    "cold_production_summary_plot",
    "cold_production_histogram",
    "temperature_vs_radius_plot",
    "timestep_histogram",
    "plot_el_depth_2d",
    "grouped_el_depth_2d_plot",
    "el_depth_angle_maps_plot",
    "el_depth_summary_plot",
    "el_depth_distribution_plot",
    "plot_cross_sections",
    "plot_loss_rates",
    "plot_trajectory_3d",
]
