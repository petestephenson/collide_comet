"""Small, made-up run files with the same layout as the real model output.

These exist so that the notebooks and tests can run without access to the
simulation data.  **The numbers are not physical results** -- they are smooth
analytic shapes (a 1/r-like coma that is denser towards the outgassing axis,
Maxwellian energy spectra, log-normal electron-depth distributions) chosen
only to exercise every reader and plot.

>>> import collide_comet as cc
>>> cases = cc.synthetic.demo_case_set()
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional, Sequence, Union

import netCDF4
import numpy as np

from . import constants as c
from .cases import DEFAULT_CASE_NAMES, CaseSet
from .physics import exobase_radius, local_outgassing_from_cos

__all__ = ["write_synthetic_run", "write_demo_cases", "demo_case_set", "DEMO_PE_RUNS", "DEMO_SW_RUNS"]

DEMO_PE_RUNS = ("9000001", "9000002", "9000003", "9000004", "9000005")
DEMO_SW_RUNS = ("9100001", "9100002", "9100003", "9100004", "9100005")
_DEMO_QVECS = {"-X": (-1, 0, 0), "+Y": (0, 1, 0), "+Z": (0, 0, 1), "-Z": (0, 0, -1), "Ref": (1, 0, 0)}


def _energy_edges() -> np.ndarray:
    return np.r_[np.arange(0.0, 10.0, 0.5), np.arange(10.0, 60.0, 5.0), np.arange(60.0, 401.0, 20.0)]


def _maxwellian(e_mid: np.ndarray, dE: np.ndarray, temperature: np.ndarray) -> np.ndarray:
    """Unit-density Maxwellian energy spectrum, shape ``T.shape + (nE,)``."""
    T = np.asarray(temperature, dtype=float)[..., None]
    f = np.sqrt(e_mid) * np.exp(-e_mid / T)
    return f / (f * dE).sum(axis=-1, keepdims=True)


def _put(ds: netCDF4.Dataset, name: str, array, dim_prefix: str) -> None:
    """Write ``array`` (in this package's axis order) the way Fortran would."""
    data = np.asarray(array)
    if data.ndim == 0:
        ds.createVariable(name, "f8")[...] = float(data)
        return
    data = data.T  # reverse axes: netCDF (C) order is the reverse of Fortran order
    dims = []
    for i, size in enumerate(data.shape):
        dim = f"{dim_prefix}_{i}_{size}"
        if dim not in ds.dimensions:
            ds.createDimension(dim, size)
        dims.append(dim)
    dtype = "f4" if data.ndim >= 4 else "f8"
    var = ds.createVariable(name, dtype, dims, zlib=data.ndim >= 4, complevel=1)
    var[...] = data


def write_synthetic_run(
    path: Union[str, Path],
    kind: str = "pe",
    q_vec: Sequence[float] = (1.0, 0.0, 0.0),
    asymmetric: bool = True,
    n_cells: int = 30,
    half_width: float = 300.0,
    q_total: float = 1.0e27,
    seed: int = 0,
) -> Path:
    """Write one synthetic run file.

    Parameters
    ----------
    kind : {"pe", "sw"}
        Photoelectron-seeded or solar-wind-seeded run.
    q_vec : direction of maximum outgassing.
    asymmetric : if False the coma is spherically symmetric.
    n_cells : cells per dimension (even, so that 0 is a cell boundary).
    half_width : the grid spans ``[-half_width, half_width]`` km.
    q_total : total outgassing rate [s^-1].
    """
    if kind not in ("pe", "sw"):
        raise ValueError("kind must be 'pe' or 'sw'")
    rng = np.random.default_rng(seed)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    q_vec = np.asarray(q_vec, dtype=float)
    q_hat = q_vec / np.linalg.norm(q_vec)
    q_uniform, q_dayside = (0.25 * q_total, 0.75 * q_total) if asymmetric else (q_total, 0.0)

    edges = np.linspace(-half_width, half_width, n_cells + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    X, Y, Z = np.meshgrid(centres, centres, centres, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2)
    cos_q = (X * q_hat[0] + Y * q_hat[1] + Z * q_hat[2]) / r
    q_local = local_outgassing_from_cos(cos_q, q_uniform, q_dayside)
    r_exo = exobase_radius(q_local)  # km
    activity = q_local / q_total
    collisional = np.exp(-r / (1.5 * r_exo))  # ~1 well inside the exobase

    e_edges = _energy_edges()
    e_mid, dE = 0.5 * (e_edges[1:] + e_edges[:-1]), np.diff(e_edges)

    # --- energy-resolved densities per generation --------------------------------
    if kind == "pe":
        n_primary = 6.0 * activity * 60.0 / (r + 30.0)
        T_primary = 8.0 + 4.0 * (1.0 - collisional)
        n_secondary = 25.0 * activity * collisional * 40.0 / (r + 20.0)
        T_secondary = 0.25 + 2.0 * (1.0 - collisional)
    else:
        n_primary = 4.0 * (1.0 + 2.5 * np.exp(-r / 120.0)) * (1.0 - 0.6 * collisional)
        T_primary = 12.0 + 30.0 * np.exp(-r / 200.0) * (1.0 - collisional)
        n_secondary = 8.0 * activity * collisional * 40.0 / (r + 20.0)
        T_secondary = 0.3 + 3.0 * (1.0 - collisional)
    noise = 1.0 + 0.03 * rng.standard_normal(r.shape)
    plas_dens_E = np.stack(
        [
            (n_primary * noise)[..., None] * _maxwellian(e_mid, dE, T_primary),
            (n_secondary * noise)[..., None] * _maxwellian(e_mid, dE, T_secondary),
        ],
        axis=-1,
    )  # (x, y, z, E, gen)

    # --- velocity moments (comp, x, y, z, gen), as MATLAB saw them ---------------
    unit_r = np.stack([X, Y, Z]) / r
    bulk = np.stack([2.0e4 * unit_r * (1 - collisional), 5.0e3 * unit_r * (1 - collisional)], axis=-1)
    v_th2 = c.QE / c.ME * np.stack([T_primary, T_secondary], axis=-1)  # (x,y,z,gen), per component
    square = bulk**2 + v_th2[None, ...]

    # --- electron depth distribution (n_depth, n_xs, n_r, n_theta, 3) ------------
    n_depth, n_xs, n_r, n_theta = 91, 3, 30, 18
    depth_grid = np.logspace(-5, 4, n_depth)
    rad_edges = np.linspace(0.0, 300.0, n_r + 1)
    r_mid = 0.5 * (rad_edges[1:] + rad_edges[:-1])
    cos_bounds = np.linspace(-1.0, 1.0, n_theta + 1)
    cos_mid = 0.5 * (cos_bounds[1:] + cos_bounds[:-1])
    log_tau = np.log10(depth_grid)
    depth_dist = np.zeros((n_depth, n_xs, n_r, n_theta, 3))
    for i_dim in range(3):
        # azimuthal average of the angle to q_hat for a polar axis e_dim
        q_ring = local_outgassing_from_cos(cos_mid * q_hat[i_dim], q_uniform, q_dayside)
        exo_ring = exobase_radius(q_ring)
        median = (exo_ring[None, :] / r_mid[:, None]) ** 1.5  # (n_r, n_theta)
        weight = (q_ring[None, :] / q_total) * np.exp(-r_mid[:, None] / 150.0) * 1e20
        for i_xs in range(n_xs):
            centre = np.log10(median * (0.3, 0.6, 1.0)[i_xs])
            pdf = np.exp(-0.5 * ((log_tau[:, None, None] - centre[None]) / 0.8) ** 2)
            pdf /= pdf.sum(axis=0, keepdims=True)
            depth_dist[:, i_xs, :, :, i_dim] = pdf * weight[None]
    depth_dist[0] = 0.0  # underflow bin

    theta_dist = np.tile((1.0 + 0.5 * cos_mid)[:, None], (1, 3)) * 1e20
    dt_edges = np.logspace(-8, -3, 41)
    dt_count = 1e6 * np.exp(-0.5 * ((np.log10(np.sqrt(dt_edges[1:] * dt_edges[:-1])) + 5.5) / 0.6) ** 2)
    pe_eng = _maxwellian(e_mid, dE, np.array(10.0)) * dE

    with netCDF4.Dataset(str(path), "w", format="NETCDF4") as ds:
        ds.description = "SYNTHETIC demo data written by collide_comet.synthetic -- not model output"
        _put(ds, "X Boundaries", edges, "x")
        _put(ds, "Y Boundaries", edges, "y")
        _put(ds, "Z Boundaries", edges, "z")
        _put(ds, "Energy Boundaries", e_edges, "e")
        _put(ds, "Plasma Density", plas_dens_E, "pd")
        _put(ds, "Bulk Velocity", bulk, "bv")
        _put(ds, "Square Velocity", square, "sv")
        _put(ds, "Electron Mass", c.ME, "me")
        _put(ds, "OGR", q_total, "q")
        _put(ds, "OGR Uniform", q_uniform, "q")
        _put(ds, "OGR Dayside", q_dayside, "q")
        _put(ds, "Q max vector", q_hat, "qv")
        _put(ds, "Electron Depth Dist", depth_dist, "ed")
        _put(ds, "elDepthRadGrid", rad_edges, "edr")
        _put(ds, "elDepthGrid", depth_grid, "edg")
        _put(ds, "elDepthThetaGrid", np.arccos(cos_bounds), "edt")
        _put(ds, "Particle Theta Dist", theta_dist, "ptd")
        _put(ds, "dT Boundaries", dt_edges, "dtb")
        _put(ds, "dT Count", dt_count, "dtc")
        _put(ds, "Photoel Eng Dist", pe_eng, "ped")
    return path


def write_demo_cases(directory: Union[str, Path], n_cells: int = 30, overwrite: bool = False):
    """Write the ten demo run files (five cases x PE/SW) into ``directory``.

    Returns ``(pe_runs, sw_runs, names)``.
    """
    directory = Path(directory)
    for i, name in enumerate(DEFAULT_CASE_NAMES):
        for kind, runs in (("pe", DEMO_PE_RUNS), ("sw", DEMO_SW_RUNS)):
            path = directory / f"plas_dens{runs[i]}.nc"
            if overwrite or not path.exists():
                write_synthetic_run(path, kind, _DEMO_QVECS[name], asymmetric=(name != "Ref"),
                                    n_cells=n_cells, seed=10 * i + (kind == "sw"))
    return list(DEMO_PE_RUNS), list(DEMO_SW_RUNS), list(DEFAULT_CASE_NAMES)


def demo_case_set(directory: Optional[Union[str, Path]] = None, n_cells: int = 30, **kwargs) -> CaseSet:
    """A :class:`CaseSet` of synthetic demo runs (written on first use).

    Files go to ``directory`` (default: ``<temp dir>/collide_comet_demo_<n_cells>``)
    and are reused on later calls (about 75 MB for the default size; delete
    the folder to reclaim the space).  Remember that these are **not** model
    results.
    """
    if directory is None:
        directory = Path(tempfile.gettempdir()) / f"collide_comet_demo_{n_cells}"
    pe_runs, sw_runs, names = write_demo_cases(directory, n_cells)
    return CaseSet(directory, pe_runs, sw_runs, names, **kwargs)
