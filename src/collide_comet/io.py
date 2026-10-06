"""Readers for the netCDF files used by the collisional test-particle model.

Array layout
------------
The run files are written by Fortran, so the netCDF (C-order) dimensions are
the reverse of what MATLAB's ``netcdf.getVar`` returned.  Every reader here
reverses the axes on load, so arrays come back in the same, physically
natural, order the MATLAB code worked with:

=====================  ==========================================
variable               shape returned by this module
=====================  ==========================================
plasma density         ``(nx, ny, nz, nE, n_gen)``
bulk / square velocity ``(nx, ny, nz, 3, n_gen)``
electron depth dist    ``(n_depth, n_xs, n_r, n_theta, 3)``
cross sections         ``(n_energy, n_process)``
=====================  ==========================================

All spatial arrays are indexed ``[ix, iy, iz]``.  Unlike MATLAB there is no
``permute(.., [2 1 3])`` anywhere: meshes are built with ``indexing="ij"``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Optional, Sequence, Tuple, Union

import netCDF4
import numpy as np

PathLike = Union[str, Path]

__all__ = [
    "Grids",
    "Outgassing",
    "CrossSections",
    "describe",
    "load_run_data",
    "load_outgassing",
    "load_velocity_arrays",
    "load_electron_mass",
    "load_electron_depth_dist",
    "load_electron_theta_dist",
    "read_timestep_grid",
    "load_photoelectron_energy_dist",
    "load_cross_sections",
    "load_eb_fields",
    "load_trajectories",
]


# ---------------------------------------------------------------------------
# low-level helpers
# ---------------------------------------------------------------------------
def _open(path: PathLike) -> netCDF4.Dataset:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"netCDF file not found: {path}")
    return netCDF4.Dataset(str(path), "r")


def _get(ds: netCDF4.Dataset, name: str) -> np.ndarray:
    """Read variable ``name`` and reverse its axes (Fortran/MATLAB order)."""
    if name not in ds.variables:
        available = ", ".join(repr(k) for k in ds.variables)
        raise KeyError(
            f"Variable {name!r} not found in {ds.filepath()}. "
            f"Available variables: {available}"
        )
    var = ds.variables[name]
    var.set_auto_maskandscale(True)
    data = var[...]
    if np.ma.isMaskedArray(data):
        if np.issubdtype(data.dtype, np.floating):
            data = data.filled(np.nan)
        else:
            data = data.filled()
    return np.asarray(data).T


def _get_scalar(ds: netCDF4.Dataset, name: str) -> float:
    return float(np.asarray(_get(ds, name)).reshape(-1)[0])


def _get_vector(ds: netCDF4.Dataset, name: str) -> np.ndarray:
    return np.asarray(_get(ds, name), dtype=float).reshape(-1)


def describe(path: PathLike, print_out: bool = True) -> dict:
    """Summarise the variables in a netCDF file.

    Returns ``{name: shape}`` with shapes given in the axis order used by this
    package (i.e. reversed with respect to ``ncdump -h``).
    """
    with _open(path) as ds:
        info = {name: tuple(reversed(var.shape)) for name, var in ds.variables.items()}
    if print_out:
        print(Path(path).name)
        for name, shape in info.items():
            print(f"  {name!r:28s} {shape}")
    return info


# ---------------------------------------------------------------------------
# containers
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Grids:
    """Cell boundaries of a run: space [km] and energy [eV]."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    energy: np.ndarray

    # -- sizes -------------------------------------------------------------
    @property
    def shape(self) -> Tuple[int, int, int]:
        """Number of spatial cells ``(nx, ny, nz)``."""
        return (self.x.size - 1, self.y.size - 1, self.z.size - 1)

    @property
    def n_energy(self) -> int:
        return self.energy.size - 1

    # -- 1D helpers --------------------------------------------------------
    @property
    def edges(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return (self.x, self.y, self.z)

    @property
    def centres(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Cell-centre coordinates ``(xc, yc, zc)`` [km]."""
        return tuple(0.5 * (e[1:] + e[:-1]) for e in self.edges)

    @property
    def lower_edges(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Lower cell boundaries (what the MATLAB slice plots used as nodes)."""
        return tuple(e[:-1] for e in self.edges)

    @property
    def energy_mid(self) -> np.ndarray:
        return 0.5 * (self.energy[1:] + self.energy[:-1])

    @property
    def dE(self) -> np.ndarray:
        return np.diff(self.energy)

    # -- 3D helpers --------------------------------------------------------
    def mesh(self, where: str = "centres"):
        """3D coordinate arrays ``(X, Y, Z)`` indexed ``[ix, iy, iz]``.

        ``where`` is ``"centres"`` or ``"lower"`` (lower cell boundaries).
        """
        if where == "centres":
            axes = self.centres
        elif where == "lower":
            axes = self.lower_edges
        else:
            raise ValueError("where must be 'centres' or 'lower'")
        return np.meshgrid(*axes, indexing="ij")

    def radius(self, where: str = "centres") -> np.ndarray:
        """Cometocentric distance of every cell [km]."""
        X, Y, Z = self.mesh(where)
        return np.sqrt(X**2 + Y**2 + Z**2)

    @property
    def cell_volume(self) -> np.ndarray:
        """Cell volumes [km^3], shape ``(nx, ny, nz)``."""
        dx, dy, dz = (np.diff(e) for e in self.edges)
        return dx[:, None, None] * dy[None, :, None] * dz[None, None, :]

    def as_list(self):
        """``[x, y, z, energy]`` -- the MATLAB ``grids`` cell array."""
        return [self.x, self.y, self.z, self.energy]


@dataclass(frozen=True)
class Outgassing:
    """Outgassing description stored in a run file.

    The local outgassing rate towards direction ``n`` is
    ``q_uniform + 4 * q_dayside * max(n . q_vec, 0)``.
    """

    q_spherical: float
    q_dayside: float
    q_uniform: float
    q_vec: np.ndarray = field(default_factory=lambda: np.zeros(3))

    def symmetric(self) -> "Outgassing":
        """Spherically symmetric equivalent (used for the reference case).

        MATLAB: ``qDayside = 0; qUniform = qSpherical;``
        """
        return Outgassing(self.q_spherical, 0.0, self.q_spherical, self.q_vec)

    @property
    def q_hat(self) -> np.ndarray:
        """Unit vector towards maximum outgassing (NaN if ``q_vec`` is zero)."""
        norm = np.linalg.norm(self.q_vec)
        if norm == 0:
            return np.full(3, np.nan)
        return self.q_vec / norm


@dataclass(frozen=True)
class CrossSections:
    """Electron-impact collision cross sections.

    Attributes
    ----------
    energy : (n_energy,) incident electron energy [eV]
    xs : (n_energy, n_process) cross sections [cm^2]
    e_th : (n_process,) threshold energies [eV]
    scat_cdf, scat_theta, scat_energy : optional angular scattering tables
        (present in ``Adv_scat_Coll_XS.nc``).
    """

    energy: np.ndarray
    xs: np.ndarray
    e_th: np.ndarray
    scat_cdf: Optional[np.ndarray] = None
    scat_theta: Optional[np.ndarray] = None
    scat_energy: Optional[np.ndarray] = None

    @property
    def n_process(self) -> int:
        return self.e_th.size


# ---------------------------------------------------------------------------
# run files (plas_dens<run>.nc)
# ---------------------------------------------------------------------------
def load_grids(path: PathLike) -> Grids:
    """Read the spatial and energy cell boundaries of a run."""
    with _open(path) as ds:
        return Grids(
            x=_get_vector(ds, "X Boundaries"),
            y=_get_vector(ds, "Y Boundaries"),
            z=_get_vector(ds, "Z Boundaries"),
            energy=_get_vector(ds, "Energy Boundaries"),
        )


def load_run_data(path: PathLike) -> Tuple[Grids, np.ndarray]:
    """Read grids and the energy-resolved plasma density of a run.

    Port of ``asymPlots.loadRunData``.

    Returns
    -------
    grids : Grids
    plas_dens_E : ndarray, shape ``(nx, ny, nz, nE, n_gen)``
        Electron distribution [cm^-3 eV^-1] per test-particle generation.
    """
    with _open(path) as ds:
        grids = Grids(
            x=_get_vector(ds, "X Boundaries"),
            y=_get_vector(ds, "Y Boundaries"),
            z=_get_vector(ds, "Z Boundaries"),
            energy=_get_vector(ds, "Energy Boundaries"),
        )
        plas_dens_E = _get(ds, "Plasma Density")

    if plas_dens_E.ndim == 4:  # single generation written without a gen axis
        plas_dens_E = plas_dens_E[..., np.newaxis]
    expected = grids.shape + (grids.n_energy,)
    if plas_dens_E.shape[:4] != expected:
        raise ValueError(
            f"'Plasma Density' in {path} has shape {plas_dens_E.shape} "
            f"(after axis reversal) but the grids imply {expected} + (n_gen,)."
        )
    return grids, plas_dens_E


def load_outgassing(path: PathLike) -> Outgassing:
    """Read the outgassing properties.  Port of ``loadOutgassingProperties``."""
    with _open(path) as ds:
        return Outgassing(
            q_spherical=_get_scalar(ds, "OGR"),
            q_dayside=_get_scalar(ds, "OGR Dayside"),
            q_uniform=_get_scalar(ds, "OGR Uniform"),
            q_vec=_get_vector(ds, "Q max vector"),
        )


def _velocity_layout(arr: np.ndarray, name: str, path: PathLike) -> np.ndarray:
    """``(3, nx, ny, nz[, n_gen])`` -> ``(nx, ny, nz, 3, n_gen)``."""
    if arr.ndim == 4:
        arr = arr[..., np.newaxis]
    if arr.ndim != 5 or arr.shape[0] != 3:
        raise ValueError(
            f"{name!r} in {path} has unexpected shape {arr.shape}; expected "
            "(3, nx, ny, nz, n_gen) after axis reversal."
        )
    return np.moveaxis(arr, 0, 3)


def load_velocity_arrays(path: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """Read bulk and mean-square velocities.  Port of ``load_velocity_arrays``.

    Returns
    -------
    bulk_vel, square_vel : ndarray, shape ``(nx, ny, nz, 3, n_gen)``
        ``<v_i>`` [m/s] and ``<v_i^2>`` [m^2/s^2].  NaNs in ``square_vel`` are
        replaced by zero, as in the MATLAB loader.
    """
    with _open(path) as ds:
        bulk = _velocity_layout(_get(ds, "Bulk Velocity"), "Bulk Velocity", path)
        square = _velocity_layout(_get(ds, "Square Velocity"), "Square Velocity", path)
    square = np.where(np.isnan(square), 0.0, square)
    return bulk, square


def load_electron_mass(path: PathLike) -> float:
    """Electron mass used by the run [kg].  Port of ``get_electron_mass``."""
    with _open(path) as ds:
        return _get_scalar(ds, "Electron Mass")


def load_electron_depth_dist(
    path: PathLike, depth_grid: str = "logspace"
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read the electron (collisional) depth distribution.

    Port of ``asymPlots.load_electron_depth_dist``.

    Parameters
    ----------
    depth_grid : {"logspace", "file"}
        The MATLAB loader ignores the values of ``elDepthGrid`` in the file and
        replaces them with ``logspace(-5, 4, n)``; ``"logspace"`` (default)
        does the same.  Use ``"file"`` to get the stored values instead.

    Returns
    -------
    rad_edges : (n_r + 1,) radial bin boundaries [km]
    depth_grid : (n,) electron depth grid
    depth_dist : (n_depth, n_xs, n_r, n_theta, 3)
        Weight of electrons per depth bin, collision-process group, radius,
        polar angle and polar axis (x, y, z).
    """
    with _open(path) as ds:
        dist = _get(ds, "Electron Depth Dist")
        rad_edges = _get_vector(ds, "elDepthRadGrid")
        grid = _get_vector(ds, "elDepthGrid")
    if depth_grid == "logspace":
        grid = np.logspace(-5, 4, grid.size)
    elif depth_grid != "file":
        raise ValueError("depth_grid must be 'logspace' or 'file'")
    return rad_edges, grid, dist


def load_electron_theta_dist(path: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """Read the particle polar-angle distribution.

    Port of ``asymPlots.load_electron_theta_dist``.
    Returns ``(theta_grid, theta_dist)``.
    """
    with _open(path) as ds:
        dist = _get(ds, "Particle Theta Dist")
        grid = _get_vector(ds, "elDepthThetaGrid")
    return grid, dist


def read_timestep_grid(path: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """Read the histogram of particle time steps: ``(dt_edges, dt_count)``."""
    with _open(path) as ds:
        return _get_vector(ds, "dT Boundaries"), _get_vector(ds, "dT Count")


def load_photoelectron_energy_dist(path: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """Read the injected photoelectron energy distribution.

    Returns ``(energy_edges, counts)`` (variables ``Energy Boundaries`` and
    ``Photoel Eng Dist``).
    """
    with _open(path) as ds:
        return _get_vector(ds, "Energy Boundaries"), _get_vector(ds, "Photoel Eng Dist")


# ---------------------------------------------------------------------------
# cross sections
# ---------------------------------------------------------------------------
BUNDLED_XS_FILES = ("Coll_XS.nc", "Adv_scat_Coll_XS.nc", "Combined_XS.nc")


def bundled_xs_path(name: str = "Coll_XS.nc") -> Path:
    """Path of a cross-section file shipped inside the package."""
    if name not in BUNDLED_XS_FILES:
        raise ValueError(f"Unknown bundled file {name!r}; choose from {BUNDLED_XS_FILES}")
    return Path(str(resources.files("collide_comet").joinpath("data", name)))


def load_cross_sections(path: Optional[PathLike] = None) -> CrossSections:
    """Read a collision cross-section table.  Port of ``loadXSections``.

    Parameters
    ----------
    path : path-like, optional
        A cross-section netCDF file, or the name of a bundled one
        (``"Coll_XS.nc"``, ``"Adv_scat_Coll_XS.nc"``, ``"Combined_XS.nc"``).
        Defaults to the copy of ``Coll_XS.nc`` shipped with the package (the
        MATLAB code read this from a hard-coded OneDrive folder).
    """
    if path is None:
        path = bundled_xs_path("Coll_XS.nc")
    elif str(path) in BUNDLED_XS_FILES and not Path(path).exists():
        path = bundled_xs_path(str(path))

    with _open(path) as ds:
        xs = np.asarray(_get(ds, "XS"), dtype=float)
        energy = _get_vector(ds, "XS_Eng")
        e_th = _get_vector(ds, "E_th")
        extra = {}
        for key, name in (
            ("scat_cdf", "Scat_CDF"),
            ("scat_theta", "Scat_Theta"),
            ("scat_energy", "Scat_Eng"),
        ):
            if name in ds.variables:
                extra[key] = np.asarray(_get(ds, name), dtype=float)
    if xs.shape != (energy.size, e_th.size):
        raise ValueError(
            f"'XS' in {path} has shape {xs.shape}; expected {(energy.size, e_th.size)}"
        )
    return CrossSections(energy=energy, xs=xs, e_th=e_th, **extra)


# ---------------------------------------------------------------------------
# PIC fields and particle trajectories
# ---------------------------------------------------------------------------
def load_eb_fields(path: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """Read time-averaged PIC fields from an ``EB_fields_*.nc`` file.

    The MATLAB scripts take the first two variables of the file as E and B
    (``netcdf.getVar(ncid, 0)`` / ``(ncid, 1)``); this does the same.

    Returns
    -------
    E, B : ndarray, shape ``(nx, ny, nz, 3)``
    """
    with _open(path) as ds:
        names = list(ds.variables)
        if len(names) < 2:
            raise ValueError(f"{path} must contain at least two variables (E and B)")
        fields = []
        for name in names[:2]:
            arr = _get(ds, name)
            if arr.ndim != 4:
                raise ValueError(f"{name!r} in {path} has shape {arr.shape}; expected 4D")
            if arr.shape[0] == 3 and arr.shape[-1] != 3:
                arr = np.moveaxis(arr, 0, -1)
            fields.append(arr)
    return fields[0], fields[1]


def load_trajectories(
    path: PathLike, particles: Optional[Sequence[int]] = None, trim: bool = True
) -> dict:
    """Read test-particle trajectories from a ``trajectories<run>.nc`` file.

    Particles are stored as variable pairs ``X<i>`` / ``V<i>`` (1-based ``i``).

    Parameters
    ----------
    particles : sequence of int, optional
        Particle numbers to read (default: all in the file).
    trim : bool
        Drop the unused tail of each trajectory (positions equal to zero), as
        the MATLAB scripts do.

    Returns
    -------
    dict ``{i: (X, V)}`` with arrays of shape ``(n_steps, 3)``.
    """
    out = {}
    with _open(path) as ds:
        if particles is None:
            particles = sorted(
                int(name[1:]) for name in ds.variables if name[:1] == "X" and name[1:].isdigit()
            )
        for i in particles:
            X = np.asarray(_get(ds, f"X{i}"), dtype=float)
            V = np.asarray(_get(ds, f"V{i}"), dtype=float)
            # MATLAB transposes to (n_steps, 3)
            if X.shape[0] == 3 and X.shape[-1] != 3:
                X, V = X.T, V.T
            if trim:
                keep = np.zeros(X.shape[0], dtype=bool)
                keep[:-1] = X[1:, 0] != 0
                X, V = X[keep], V[keep]
            out[int(i)] = (X, V)
    return out
