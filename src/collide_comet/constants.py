"""Physical constants and model conventions.

The numerical values deliberately match the ones hard-coded in the original
MATLAB scripts (e.g. ``qe = 1.6e-19``) so that the Python port reproduces the
MATLAB numbers rather than differing in the third significant figure.
"""

#: Elementary charge [C] (value used throughout the MATLAB code).
QE = 1.6e-19
#: Electron mass [kg] (value used throughout the MATLAB code).
ME = 9.11e-31
#: Boltzmann constant [J/K].
KB = 1.38e-23
#: Elementary charge used for the eV -> J conversion of pressures [C].
QE_PRESSURE = 1.602e-19

#: Neutral outflow speed [km/s].
U_GAS = 1.0
#: Electron-neutral cross section used for exobase estimates [cm^2].
SIGMA_EXOBASE = 6.7e-17

#: Index (0-based) of the first ionisation process in the collision table.
#: Processes ``IX_IONI:`` are ionisation (MATLAB ``ix_ioni = 7``).
IX_IONI = 6

#: Names of the 12 electron-water collision processes in ``Coll_XS.nc``.
PROCESS_NAMES = (
    "Elastic",
    "Rotation",
    "Vibration 010",
    "Vibration Stretch",
    "Electronic",
    "Dissociation",
    r"H$_2$O$^+$",
    r"OH$^+$",
    r"O$^+$",
    r"H$_2^+$",
    r"H$^+$",
    r"O$^{2+}$",
)

#: Test-particle generations stored along the last axis of the run arrays.
#: MATLAB ``inc_gens = [1, 2]`` is ``gens=(0, 1)`` here.
GEN_PRIMARY = 0
GEN_SECONDARY = 1
ALL_GENS = (GEN_PRIMARY, GEN_SECONDARY)
