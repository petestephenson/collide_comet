"""Named 3D quantities that can be mapped for a case.

Each entry replaces one family of near-identical MATLAB functions
(``make_totalelectrons_Plot``, ``make_coldelectrons_Plot``,
``make_grouped_EIIFreqPlot`` ...): the registry records how the quantity is
computed, its colour-bar label and the default colour limits those functions
hard-coded.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Optional, Tuple

import numpy as np

from . import constants as c
from .cases import CaseData

__all__ = ["Quantity", "QUANTITIES", "compute", "get_quantity", "list_quantities"]


@dataclass(frozen=True)
class Quantity:
    """Description of a mappable quantity.

    Attributes
    ----------
    func : callable ``CaseData -> (nx, ny, nz) array``
    label : colour-bar label (without any ``log10``)
    log : plot ``log10`` of the quantity by default
    clim : default colour limits (of ``log10(value)`` when ``log``)
    ratio_label : label for the asymmetric / reference ratio
    ratio_clim : default colour limits of the (log10) ratio
    ratio_log : plot ``log10`` of the ratio by default
    matlab : the MATLAB function(s) this entry replaces
    """

    func: Callable[[CaseData], np.ndarray]
    label: str
    log: bool = True
    clim: Optional[Tuple[float, float]] = None
    ratio_label: Optional[str] = None
    ratio_clim: Tuple[float, float] = (-1.0, 1.0)
    ratio_log: bool = True
    matlab: str = ""
    description: str = ""


_P, _S = c.GEN_PRIMARY, c.GEN_SECONDARY

QUANTITIES: Dict[str, Quantity] = {
    "total_density": Quantity(
        lambda case: case.densities("total", (_P, _S))[0],
        r"$n_e$ [cm$^{-3}$]", clim=(-0.5, 2.5),
        ratio_label=r"$n_{e,\mathrm{asym}} / n_{e,\mathrm{sym}}$",
        matlab="make_totalelectrons_Plot / make_totalelectrons_Ratio_Plot",
        description="All electrons (PE + SW runs, both generations)",
    ),
    "cold_density": Quantity(
        lambda case: case.densities("total", (_P, _S), 0.0, 1.0)[1],
        r"$n_{e,\mathrm{cold}}$ [cm$^{-3}$]", clim=(-1.5, 1.5),
        ratio_label=r"$n_{e,c,\mathrm{asym}} / n_{e,c,\mathrm{sym}}$", ratio_clim=(-2.0, 2.0),
        matlab="make_coldelectrons_Plot / make_coldelectrons_Ratio_Plot",
        description="Electrons below 1 eV (PE + SW runs, both generations)",
    ),
    "pe_density": Quantity(
        lambda case: case.densities("pe", (_P,))[0],
        r"$n_e$ [cm$^{-3}$]", clim=(-2.5, 1.5),
        ratio_label=r"$n_{e,\mathrm{asym}} / n_{e,\mathrm{sym}}$",
        matlab="make_pe_electrons_Plot / make_pe_electrons_Ratio_Plot",
        description="Primary photoelectrons",
    ),
    "sw_density": Quantity(
        lambda case: case.densities("sw", (_P,))[0],
        r"$n_e$ [cm$^{-3}$]", clim=(-0.5, 2.5),
        ratio_label=r"$n_{e,\mathrm{asym}} / n_{e,\mathrm{sym}}$",
        matlab="make_swelectrons_Plot / make_swelectrons_Ratio_Plot",
        description="Primary solar-wind electrons",
    ),
    "sw_warm_density": Quantity(
        lambda case: case.densities("sw", (_P,), 3.0, 350.0)[1],
        r"$n_e$ [cm$^{-3}$]", clim=(-0.5, 2.5),
        ratio_label=r"$n_{e,\mathrm{asym}} / n_{e,\mathrm{sym}}$",
        matlab="make_sw_warm_electrons_Plot",
        description="Primary solar-wind electrons between 3 and 350 eV",
    ),
    "sw_secondary_cold_density": Quantity(
        lambda case: case.densities("sw", (_S,), 0.0, 1.0)[1],
        r"$n_{e,c}$ [cm$^{-3}$]", clim=(-2.5, 1.5),
        ratio_label=r"$n_{e,c,\mathrm{asym}} / n_{e,c,\mathrm{sym}}$",
        matlab="make_sw_seccoldelectrons_Plot / make_sw_seccoldelectrons_Ratio_Plot",
        description="Cold secondaries produced by solar-wind electrons",
    ),
    "secondary_density": Quantity(
        lambda case: case.densities("total", (_S,))[0],
        r"$n_e$ [cm$^{-3}$]", clim=(-2.0, 2.5),
        ratio_label=r"$n_{e,\mathrm{asym}} / n_{e,\mathrm{sym}}$",
        matlab="make_secondaryelectrons_Plot / make_secondaryelectrons_Ratio_Plot",
        description="Secondary electrons (PE + SW runs)",
    ),
    "neutral_density": Quantity(
        lambda case: case.neutral_density(),
        r"$n_{\mathrm{H_2O}}$ [cm$^{-3}$]", clim=(0.0, 12.0),
        ratio_label=r"$n_{\mathrm{H_2O,asym}} / n_{\mathrm{H_2O,sym}}$",
        matlab="make_grouped_nDensityPlot",
        description="Neutral water density of the coma model",
    ),
    "eii_frequency": Quantity(
        lambda case: case.eii_frequency(),
        r"$\nu_e$ [s$^{-1}$]", clim=(-8.0, -5.0),
        ratio_label=r"$\nu_{e,\mathrm{asym}} / \nu_{e,\mathrm{sym}}$",
        matlab="make_grouped_EIIFreqPlot / make_grouped_EIIFreqRatio_Plot",
        description="Electron-impact ionisation frequency",
    ),
    "eii_production": Quantity(
        lambda case: case.eii_production(),
        r"$P_e$ [cm$^{-3}$ s$^{-1}$]", clim=(-3.0, 3.0),
        ratio_label=r"$P_{e,\mathrm{asym}} / P_{e,\mathrm{sym}}$",
        matlab="make_grouped_EIIProductionPlot / make_grouped_EIIProductionRatio_Plot",
        description="Electron-impact ionisation production rate",
    ),
    "cold_electron_production": Quantity(
        lambda case: case.cold_electron_production(),
        r"$P_{e,\mathrm{cold}}$ [cm$^{-3}$ s$^{-1}$]", clim=(-3.0, 3.0),
        ratio_label=r"$P_{e,\mathrm{cold,asym}} / P_{e,c,\mathrm{sym}}$",
        matlab="make_grouped_cold_e_prod_Plot",
        description="Production rate of electrons below 1 eV",
    ),
    "temperature": Quantity(
        lambda case: case.temperature("combined"),
        r"$T_e$ [eV]", log=False, clim=(0.0, 40.0),
        ratio_label=r"$T_{e,\mathrm{asym}} / T_{e,\mathrm{sym}}$", ratio_clim=(0.0, 2.0),
        ratio_log=False,
        matlab="make_grouped_Temperature_Plot (PlotMode='combined')",
        description="Temperature of the combined electron population",
    ),
    "sw_temperature": Quantity(
        lambda case: case.temperature("sw"),
        r"$T_e$ [eV]", log=False, clim=(0.0, 60.0),
        ratio_label=r"$T_{e,\mathrm{asym}} / T_{e,\mathrm{sym}}$", ratio_clim=(0.0, 2.0),
        ratio_log=False,
        matlab="make_SW_temp_plot / make_SW_tempRatioplot",
        description="Temperature of the primary solar-wind electrons",
    ),
}


def get_quantity(name) -> Quantity:
    """Look up a quantity by name (a :class:`Quantity` is passed through)."""
    if isinstance(name, Quantity):
        return name
    try:
        return QUANTITIES[name]
    except KeyError:
        raise KeyError(
            f"Unknown quantity {name!r}. Available: {', '.join(QUANTITIES)}"
        ) from None


def compute(case: CaseData, name) -> np.ndarray:
    """Evaluate quantity ``name`` for ``case``; returns ``(nx, ny, nz)``."""
    return get_quantity(name).func(case)


def list_quantities(print_out: bool = True) -> Dict[str, str]:
    """Names and descriptions of the registered quantities."""
    info = {name: q.description for name, q in QUANTITIES.items()}
    if print_out:
        for name, desc in info.items():
            print(f"{name:28s} {desc}")
    return info
