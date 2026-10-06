"""Shared settings for the notebooks -- edit this file to point at your data.

Every notebook starts with ``from run_config import get_cases`` so the data
folder and run numbers only have to be set here.
"""

from pathlib import Path

import collide_comet as cc

# --- where the run files (plas_dens<run>.nc) live -----------------------------
# Same folder as ``plasDensDir`` in asymmetricPlots.m.
DATA_DIR = Path("E:/TestPlSimData")

# --- which runs make up each case ---------------------------------------------
# Same run numbers and ordering as ``peCasesNo`` / ``swCasesNo`` in
# asymmetricPlots.m: one photoelectron (PE) and one solar-wind (SW) run per case.
CASE_NAMES = ["-X", "+Y", "+Z", "-Z", "Ref"]
PE_RUNS = ["1944473", "1944689", "1944691", "1944692", "2857936"]
SW_RUNS = ["2767385", "2777459", "2777466", "2777467", "2857233"]
REFERENCE = "Ref"

# --- figure output --------------------------------------------------------------
# Figures are written to <repository>/Figures/<subfolder>/<name>.<format>.
FIG_FORMATS = ["png"]  # e.g. ["png", "pdf"]
FIG_DPI = 300


def get_cases(allow_demo=True):
    """Return ``(cases, demo)``: the CaseSet to plot and whether it is demo data.

    If the run files in ``DATA_DIR`` cannot be found, a set of small synthetic
    runs is used instead so the notebooks still execute.  Demo figures are
    NOT model results; they are saved under ``Figures/demo`` to keep them apart.
    """
    cases = cc.CaseSet(DATA_DIR, PE_RUNS, SW_RUNS, CASE_NAMES, reference=REFERENCE)
    missing = cases.missing_files()
    if not missing:
        print(f"Using model output in {DATA_DIR}")
        return cases, False
    if not allow_demo:
        raise FileNotFoundError("Run files not found:\n  " + "\n  ".join(map(str, missing)))
    print(
        f"WARNING: {len(missing)} of {2 * len(cases)} run files were not found in {DATA_DIR}\n"
        f"         (first missing: {missing[0].name}).\n"
        "         Falling back to SYNTHETIC DEMO DATA -- the figures below are not model results.\n"
        "         Edit notebooks/run_config.py to point at your data."
    )
    return cc.synthetic.demo_case_set(), True


def make_saver(topic, demo):
    """Return ``save(fig, name)`` writing into ``Figures/<topic>`` (or ``Figures/demo/<topic>``)."""
    subfolder = Path("demo") / topic if demo else Path(topic)

    def save(fig, name):
        paths = cc.save_figure(fig, name, subfolder=subfolder, formats=FIG_FORMATS, dpi=FIG_DPI)
        for path in paths:
            print("saved", path)
        return paths

    return save
