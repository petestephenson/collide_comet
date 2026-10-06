"""Groups of simulation runs ("cases") and their derived quantities.

Every case of the asymmetric-outgassing study consists of two test-particle
runs -- one seeded with photoelectrons (PE) and one with solar-wind electrons
(SW) -- whose outputs are added together.  The MATLAB plotting functions each
took ``peCasesNo, swCasesNo, plasDensDir`` and re-read the files; here a
:class:`CaseSet` holds that information once and hands out :class:`CaseData`
objects that load lazily and cache what they compute.

Example
-------
>>> cases = CaseSet(
...     "E:/TestPlSimData",
...     pe_runs=["1944473", "1944689", "1944691", "1944692", "2857936"],
...     sw_runs=["2767385", "2777459", "2777466", "2777467", "2857233"],
...     names=["-X", "+Y", "+Z", "-Z", "Ref"],
... )
>>> ref = cases["Ref"]
>>> n_e, n_cold = ref.densities()
"""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from typing import Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np

from . import constants as c
from . import io, physics

__all__ = ["CaseSet", "CaseData", "DEFAULT_CASE_NAMES"]

DEFAULT_CASE_NAMES = ("-X", "+Y", "+Z", "-Z", "Ref")
_SOURCES = ("total", "pe", "sw")


def _is_reference_name(name: str) -> bool:
    low = str(name).strip().lower()
    return "ref" in low or "sym" in low


class CaseData:
    """One case: a PE run plus a SW run, with cached derived quantities.

    Obtain instances from a :class:`CaseSet` (``cases["-X"]``) or build one
    directly from two file paths.
    """

    def __init__(
        self,
        name: str,
        pe_path: Union[str, Path],
        sw_path: Union[str, Path],
        is_reference: Optional[bool] = None,
        xs: Optional[io.CrossSections] = None,
        _owner: Optional["CaseSet"] = None,
    ):
        self.name = str(name)
        self.pe_path = Path(pe_path)
        self.sw_path = Path(sw_path)
        self.is_reference = _is_reference_name(name) if is_reference is None else bool(is_reference)
        self._xs = xs
        self._owner = _owner
        self._raw = None  # large energy-resolved arrays
        self._small = {}  # grids, outgassing, electron mass ...
        self._derived = {}  # cached 3D results

    def __repr__(self) -> str:
        kind = "reference" if self.is_reference else "asymmetric"
        return f"CaseData({self.name!r}, pe={self.pe_path.name}, sw={self.sw_path.name}, {kind})"

    # -- raw data ----------------------------------------------------------
    def _load_raw(self) -> dict:
        if self._raw is None:
            grids, pe = io.load_run_data(self.pe_path)
            grids_sw, sw = io.load_run_data(self.sw_path)
            if grids.shape != grids_sw.shape or grids.n_energy != grids_sw.n_energy:
                raise ValueError(
                    f"Case {self.name}: PE and SW runs are on different grids "
                    f"({grids.shape}/{grids.n_energy} vs {grids_sw.shape}/{grids_sw.n_energy})."
                )
            self._small["grids"] = grids
            self._raw = {"pe": pe, "sw": sw}
        if self._owner is not None:
            self._owner._touch(self)
        return self._raw

    def release(self) -> None:
        """Free the large energy-resolved arrays (derived results are kept)."""
        self._raw = None

    def clear_cache(self) -> None:
        """Forget everything that has been loaded or computed."""
        self._raw = None
        self._small.clear()
        self._derived.clear()

    @property
    def grids(self) -> io.Grids:
        if "grids" not in self._small:
            self._small["grids"] = io.load_grids(self.pe_path)
        return self._small["grids"]

    @property
    def xs(self) -> io.CrossSections:
        if self._xs is None:
            self._xs = self._owner.xs if self._owner is not None else io.load_cross_sections()
        return self._xs

    @property
    def electron_mass(self) -> float:
        if "me" not in self._small:
            self._small["me"] = io.load_electron_mass(self.pe_path)
        return self._small["me"]

    @property
    def outgassing_raw(self) -> io.Outgassing:
        """Outgassing properties exactly as stored in the SW run file."""
        if "outgassing" not in self._small:
            self._small["outgassing"] = io.load_outgassing(self.sw_path)
        return self._small["outgassing"]

    @property
    def outgassing(self) -> io.Outgassing:
        """Outgassing used for the neutral coma.

        For the reference case the dayside term is switched off
        (``qDayside = 0; qUniform = qSpherical``), as in the MATLAB scripts.
        """
        raw = self.outgassing_raw
        return raw.symmetric() if self.is_reference else raw

    def plasma_density_E(self, source: str = "total", cometary_only: bool = False) -> np.ndarray:
        """Energy-resolved electron density ``(nx, ny, nz, nE, n_gen)``.

        Parameters
        ----------
        source : {"total", "pe", "sw"}
            Both runs summed, or the photoelectron / solar-wind run alone.
        cometary_only : bool
            Zero the primary solar-wind generation before summing (MATLAB
            ``plas_dens_E_Asym_sw(:,:,:,:,1) = 0``), leaving only electrons
            born in the coma.
        """
        if source not in _SOURCES:
            raise ValueError(f"source must be one of {_SOURCES}")
        raw = self._load_raw()
        if source == "pe":
            return raw["pe"]
        sw = raw["sw"]
        if cometary_only:
            sw = sw.copy()
            sw[..., c.GEN_PRIMARY] = 0.0
        return sw if source == "sw" else sw + raw["pe"]

    # -- cached derived quantities ------------------------------------------
    def _cached(self, key, func):
        if key not in self._derived:
            self._derived[key] = func()
        return self._derived[key]

    def densities(
        self,
        source: str = "total",
        gens: Sequence[int] = c.ALL_GENS,
        min_energy: float = 0.0,
        max_energy: float = 1.0,
        cometary_only: bool = False,
        smooth: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Total and energy-band (default 0-1 eV, "cold") densities [cm^-3]."""
        gens = tuple(int(g) for g in np.atleast_1d(gens))
        key = ("dens", source, gens, float(min_energy), float(max_energy), cometary_only, smooth)
        return self._cached(
            key,
            lambda: physics.plasma_densities(
                self.grids, self.plasma_density_E(source, cometary_only), gens,
                min_energy, max_energy, smooth,
            ),
        )

    def neutral_density(self) -> np.ndarray:
        """Water density at the cell centres [cm^-3]."""
        return self._cached(
            ("n_neutral",), lambda: physics.neutral_density_on_grid(self.grids, self.outgassing)
        )

    def eii_frequency(self, source: str = "total", gens: Sequence[int] = c.ALL_GENS) -> np.ndarray:
        """Electron-impact ionisation frequency [s^-1]."""
        gens = tuple(int(g) for g in np.atleast_1d(gens))
        return self._cached(
            ("eii_freq", source, gens),
            lambda: physics.eii_frequency(
                self.plasma_density_E(source), self.grids.energy, self.xs, gens
            ),
        )

    def eii_production(self, source: str = "total", gens: Sequence[int] = c.ALL_GENS) -> np.ndarray:
        """Electron-impact ionisation production rate [cm^-3 s^-1]."""
        return self.neutral_density() * self.eii_frequency(source, gens)

    def cold_electron_production(
        self,
        e_cold: float = 1.0,
        source: str = "total",
        gens: Sequence[int] = c.ALL_GENS,
        cometary_only: bool = False,
        legacy_speed_grid: bool = True,
    ) -> np.ndarray:
        """Production rate of electrons below ``e_cold`` eV [cm^-3 s^-1]."""
        gens = tuple(int(g) for g in np.atleast_1d(gens))
        return self._cached(
            ("cold_prod", float(e_cold), source, gens, cometary_only, legacy_speed_grid),
            lambda: physics.cold_electron_production(
                self.plasma_density_E(source, cometary_only), self.grids, self.outgassing,
                self.xs, e_cold, gens, legacy_speed_grid,
            ),
        )

    def temperature(
        self, mode: str = "combined", gens: Sequence[int] = c.ALL_GENS, sw_gen: int = c.GEN_PRIMARY
    ) -> np.ndarray:
        """Electron temperature [eV].

        Parameters
        ----------
        mode : {"combined", "sw"}
            ``"combined"``: density-weighted PE + SW population over ``gens``
            (MATLAB ``calcCombinedElectronTemperature``).  ``"sw"``: the
            ``sw_gen`` generation of the solar-wind run only.
        """
        gens = tuple(int(g) for g in np.atleast_1d(gens))

        def compute():
            me = self.electron_mass
            sw_bulk, sw_square = io.load_velocity_arrays(self.sw_path)
            if mode == "sw":
                return physics.temperature(sw_square[..., sw_gen], sw_bulk[..., sw_gen], me)
            pe_bulk, pe_square = io.load_velocity_arrays(self.pe_path)
            raw = self._load_raw()
            _, bulk, square = physics.combine_run_velocities(
                pe_bulk, pe_square, raw["pe"], sw_bulk, sw_square, raw["sw"], self.grids.dE, gens
            )
            return physics.temperature(square, bulk, me)

        if mode not in ("combined", "sw"):
            raise ValueError("mode must be 'combined' or 'sw'")
        return self._cached(("temperature", mode, gens, int(sw_gen)), compute)

    def electron_depth(self, source: str = "total", depth_grid: str = "logspace"):
        """Electron depth distribution: ``(rad_edges, depth_grid, depth_dist)``.

        ``source`` is ``"total"`` (PE + SW, as in the grouped MATLAB plots),
        ``"pe"`` or ``"sw"``.
        """
        if source not in _SOURCES:
            raise ValueError(f"source must be one of {_SOURCES}")

        def compute():
            if source == "sw":
                return io.load_electron_depth_dist(self.sw_path, depth_grid)
            rad, grid, dist = io.load_electron_depth_dist(self.pe_path, depth_grid)
            if source == "total":
                _, _, dist_sw = io.load_electron_depth_dist(self.sw_path, depth_grid)
                dist = dist + dist_sw
            return rad, grid, dist

        return self._cached(("el_depth", source, depth_grid), compute)


class CaseSet:
    """A set of cases sharing one data directory.

    Parameters
    ----------
    data_dir : path-like
        Folder containing the run files (MATLAB ``plasDensDir``).
    pe_runs, sw_runs : sequences of str/int
        Run numbers of the photoelectron and solar-wind runs, one per case
        (MATLAB ``peCasesNo`` / ``swCasesNo``).  Full paths ending in ``.nc``
        are also accepted.
    names : sequence of str
        Case labels.  Defaults to ``("-X", "+Y", "+Z", "-Z", "Ref")``, the
        ordering used throughout the MATLAB scripts.
    reference : str, optional
        Name of the spherically symmetric reference case.  Defaults to the
        case whose name contains "ref" or "sym".
    file_template : str
        File name pattern; ``{run}`` is replaced by the run number.
    xs : CrossSections or path-like, optional
        Collision cross sections (default: bundled ``Coll_XS.nc``).
    max_loaded : int
        How many cases may keep their large energy-resolved arrays in memory
        at once.  Derived 3D quantities are always kept.
    """

    def __init__(
        self,
        data_dir: Union[str, Path],
        pe_runs: Sequence[Union[str, int]],
        sw_runs: Sequence[Union[str, int]],
        names: Optional[Sequence[str]] = None,
        reference: Optional[str] = None,
        file_template: str = "plas_dens{run}.nc",
        xs: Union[None, io.CrossSections, str, Path] = None,
        max_loaded: int = 1,
    ):
        if len(pe_runs) != len(sw_runs):
            raise ValueError("pe_runs and sw_runs must have the same length")
        if names is None:
            if len(pe_runs) != len(DEFAULT_CASE_NAMES):
                raise ValueError("names must be given unless there are exactly 5 cases")
            names = DEFAULT_CASE_NAMES
        if len(names) != len(pe_runs):
            raise ValueError("names must have one entry per run pair")
        if len(set(names)) != len(names):
            raise ValueError("case names must be unique")

        self.data_dir = Path(data_dir)
        self.file_template = file_template
        self.max_loaded = max(1, int(max_loaded))
        self._xs = xs if isinstance(xs, io.CrossSections) or xs is None else io.load_cross_sections(xs)
        self._loaded: "OrderedDict[str, CaseData]" = OrderedDict()

        if reference is None:
            refs = [n for n in names if _is_reference_name(n)]
            reference = refs[0] if refs else None
        elif reference not in names:
            raise ValueError(f"reference {reference!r} is not one of the case names {list(names)}")
        self.reference_name = reference

        self._cases = OrderedDict(
            (
                str(name),
                CaseData(
                    name, self.run_path(pe), self.run_path(sw),
                    is_reference=(name == reference), _owner=self,
                ),
            )
            for name, pe, sw in zip(names, pe_runs, sw_runs)
        )

    # -- construction helpers ---------------------------------------------
    def run_path(self, run: Union[str, int, Path]) -> Path:
        """File path for a run number (or pass-through for a ``.nc`` path)."""
        text = str(run)
        if text.lower().endswith(".nc"):
            p = Path(text)
            return p if p.is_absolute() else self.data_dir / p
        return self.data_dir / self.file_template.format(run=text)

    @property
    def xs(self) -> io.CrossSections:
        if self._xs is None:
            self._xs = io.load_cross_sections()
        return self._xs

    def _touch(self, case: CaseData) -> None:
        """Keep at most ``max_loaded`` cases' large arrays in memory."""
        self._loaded.pop(case.name, None)
        self._loaded[case.name] = case
        while len(self._loaded) > self.max_loaded:
            _, old = self._loaded.popitem(last=False)
            old.release()

    # -- container protocol -------------------------------------------------
    @property
    def names(self) -> List[str]:
        return list(self._cases)

    def __len__(self) -> int:
        return len(self._cases)

    def __iter__(self) -> Iterator[CaseData]:
        return iter(self._cases.values())

    def __contains__(self, name) -> bool:
        return name in self._cases

    def __getitem__(self, key: Union[str, int]) -> CaseData:
        if isinstance(key, (int, np.integer)):
            return list(self._cases.values())[key]
        try:
            return self._cases[key]
        except KeyError:
            raise KeyError(f"No case named {key!r}; available: {self.names}") from None

    def __repr__(self) -> str:
        return f"CaseSet({str(self.data_dir)!r}, cases={self.names}, reference={self.reference_name!r})"

    @property
    def reference(self) -> CaseData:
        if self.reference_name is None:
            raise ValueError("This CaseSet has no reference case")
        return self._cases[self.reference_name]

    @property
    def asymmetric(self) -> List[CaseData]:
        """All cases except the reference, in order."""
        return [case for case in self if not case.is_reference]

    def missing_files(self) -> List[Path]:
        """Run files that do not exist on disk."""
        paths = [p for case in self for p in (case.pe_path, case.sw_path)]
        return [p for p in paths if not p.exists()]

    def clear_cache(self) -> None:
        for case in self:
            case.clear_cache()
        self._loaded.clear()
