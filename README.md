# collide_comet

Python tools for reading and plotting the output of the collisional test-particle
model of electrons at a comet. This is the Python port of the MATLAB analysis code in
this repository (the `+asymPlots` package and `asymmetricPlots.m`).

* **`src/collide_comet/`** – the package: netCDF readers, derived quantities, analysis and plotting.
* **`notebooks/`** – Jupyter notebooks that make the figures and save them to **`Figures/`**.
* **`tests/`** – checks of the port against line-by-line transliterations of the MATLAB maths.
* The original MATLAB files are untouched and still sit at the top level.

## Install

```bash
cd collide_comet
python -m venv .venv            # optional but recommended
.venv\Scripts\activate          # Windows   (source .venv/bin/activate on macOS/Linux)
pip install -e ".[notebooks]"   # the package (editable) + JupyterLab
```

Requires Python ≥ 3.9 with numpy, scipy, matplotlib and netCDF4 (installed automatically).

## Using the notebooks

1. Edit **`notebooks/run_config.py`**: set `DATA_DIR` to the folder holding the
   `plas_dens<run>.nc` files and check the run numbers (they default to the ones in
   `asymmetricPlots.m`).
2. Start `jupyter lab` and open a notebook from `notebooks/`.

| Notebook | Contents |
| --- | --- |
| `01_single_run_quickstart` | inspect a run file, load it, map quantities for one case |
| `02_grouped_maps` | five-panel maps and asymmetric/reference ratios for every quantity and view; energy spectrograms |
| `03_profiles_and_statistics` | radial and angular profiles, radius–angle maps, cold-electron production, temperature vs r |
| `04_electron_depth` | electron-depth maps and summary, exobase estimates, depth distributions |
| `05_cross_sections` | the collision cross-section tables |

Figures are written to `Figures/<topic>/<name>.png` (formats and dpi are set in
`run_config.py`). If the run files cannot be found the notebooks fall back to small
**synthetic demo data** so they still run; those figures are not model results and go to
`Figures/demo/` (git-ignored).

## Using the package directly

```python
import collide_comet as cc
from collide_comet import plotting as ccp

cases = cc.CaseSet(
    "E:/TestPlSimData",
    pe_runs=["1944473", "1944689", "1944691", "1944692", "2857936"],   # photoelectron runs
    sw_runs=["2767385", "2777459", "2777466", "2777467", "2857233"],   # solar-wind runs
    names=["-X", "+Y", "+Z", "-Z", "Ref"],
)

fig, axes, cbar = ccp.grouped_slice_plot(cases, "cold_density", view="XZ")
cc.save_figure(fig, "cold_density_XZ")            # -> Figures/cold_density_XZ.png

n_e, n_cold = cases["-X"].densities()             # 3D arrays [cm^-3], cached
grids, plas_dens_E = cc.load_run_data("E:/TestPlSimData/plas_dens2857936.nc")
cc.describe("E:/TestPlSimData/plas_dens2857936.nc")   # list the variables in a file
cc.list_quantities()                               # everything that can be mapped
```

## Conventions (differences from MATLAB worth knowing)

* **Array order.** Arrays are indexed `[ix, iy, iz, iE, i_gen]`, the order MATLAB's
  `netcdf.getVar` gave. The readers reverse the netCDF axes for you, and meshes use
  `indexing="ij"`, so none of the `permute(..., [2 1 3])` calls are needed.
* **0-based indices.** `inc_gens = [1, 2]` is `gens=(0, 1)`; `iTau = 3` is `i_tau=2`;
  `iDim`/`iAngle` 1–3 become `0–2` or the names `"X"/"Y"/"Z"` and `"YZ"/"XZ"/"XY"`.
* **Cases are named.** `cases["-X"]`, `cases.reference`; the reference case automatically uses
  `qDayside = 0, qUniform = qSpherical` as the MATLAB scripts do for case 5.
* **Slices are drawn in 2D** with `pcolormesh` on the true cell boundaries instead of a 3D
  `slice` viewed along an axis. The default colour map is the same `redblue`.
* **Cross sections** default to the copy of `XS/Coll_XS.nc` bundled with the package rather than
  a hard-coded OneDrive path; pass `xs=` / a path to use another table.
* **Random numbers.** Monte-Carlo figures use NumPy's legacy Mersenne Twister, which produces the
  same stream as MATLAB's `rng(1); rand(...)`.
* **Memory.** A `CaseSet` keeps the large energy-resolved arrays of only one case in memory at a
  time (`max_loaded=`); the derived 3D quantities are cached.

Three deliberate deviations from the MATLAB numerics:

* `physics.pressure` returns pascals. MATLAB `calculatePressure` multiplies by an extra `k_B`
  that `calcAmbipolarfield` divides out again; the ambipolar field is unchanged.
* `physics.ambipolar_field` differentiates along the true x, y, z axes on cell centres (MATLAB's
  `gradient` call treated the second array dimension as x).
* Energy-spectrogram colour bars are labelled cm⁻³ eV⁻¹ (the MATLAB label also had s⁻¹).

One MATLAB quirk is kept so numbers match: `calculateColdeProduction` evaluates the electron
speed on energy nodes left over from the last collision process in its loop. Pass
`legacy_speed_grid=False` to use plain bin mid-points (affects two energy bins).

## MATLAB → Python

| MATLAB (`asymPlots.`) | Python |
| --- | --- |
| `loadRunData`, `loadOutgassingProperties`, `load_velocity_arrays`, `get_electron_mass`, `load_electron_depth_dist`, `load_electron_theta_dist`, `read_timestep_grid`, `loadXSections` | `cc.io.load_run_data`, `load_outgassing`, `load_velocity_arrays`, `load_electron_mass`, `load_electron_depth_dist`, `load_electron_theta_dist`, `read_timestep_grid`, `load_cross_sections` |
| `calculatePlasmaDensities`, `smooth_plas_dens_E` | `cc.physics.plasma_densities`, `smooth_plas_dens_E` |
| `calcNeutralDens`, `calc_local_q_from_pos`, `calc_local_q_from_angle` | `cc.physics.neutral_density`, `local_outgassing`, `local_outgassing_from_cos` |
| `calcEIoniFreq`, `calculateEIIProduction`, `calculateColdeProduction` | `cc.physics.eii_frequency`, `eii_production`, `cold_electron_production` |
| `calculateTemperature`, `combine_run_velocities`, `calculatePressure`, `calcAmbipolarfield` | `cc.physics.temperature`, `combine_run_velocities`, `pressure`, `ambipolar_field` |
| `makePlot` and `makePlasDensPlot`, `makeTemperaturePlot`, `makeEIIFreqPlot`, … | `ccp.plot_slice`, `ccp.plot_quantity(case, "<quantity>")` |
| `make_totalelectrons_Plot`, `make_coldelectrons_Plot`, `make_pe_electrons_Plot`, `make_swelectrons_Plot`, `make_sw_warm_electrons_Plot`, `make_secondaryelectrons_Plot`, `make_sw_seccoldelectrons_Plot`, `make_grouped_nDensityPlot`, `make_grouped_EIIFreqPlot`, `make_grouped_EIIProductionPlot`, `make_grouped_cold_e_prod_Plot`, `make_grouped_Temperature_Plot`, `make_SW_temp_plot` | `ccp.grouped_slice_plot(cases, "<quantity>", view)` |
| the matching `*_Ratio_Plot` / `make_SW_tempRatioplot` functions | `ccp.grouped_ratio_plot(cases, "<quantity>", view)` |
| `make_edf_plot`, `make_edfRatioplot`, `make_SW_edf_plot`, `make_SW_edfRatioplot` | `ccp.plot_edf`, `plot_edf_ratio`, `grouped_edf_plot`, `grouped_edf_ratio_plot` |
| `make_electrondens_radial_AllCases_Plot` | `ccp.radial_profiles_plot` |
| `make_electrondens_angleToQVec_shell_MC_AllCases_Plot_v2` (and the error-bar version, `error_style="bars"`) | `ccp.angle_shell_profiles_plot` |
| `make_electrondens_axisCone_MC_scaledExobase_Plot` | `ccp.axis_cone_profiles_plot` |
| `make_electrondens_rAngle_QVec_2d_MC_AllCases_Plot_v2`, and the exobase overlay / Ref-first layout of `make_electrondens_rAngle_QVec_2d_MC_AllCases_Plot` and `make_electrondens_rAngle_QVec_2d_AllCases_Plot` | `ccp.r_angle_maps_plot` (`show_exobase=`, `reference_first=`, `reference_axes=`, `method="mc"/"cells"`) |
| `make_grouped_cold_e_prod_Plot_v3` | `ccp.cold_production_summary_plot` |
| `make_cold_e_prod_volumeWeighted_Hist_AllCases` | `ccp.cold_production_histogram` |
| `makeTemperature_vs_r_Plot_v2` | `ccp.temperature_vs_radius_plot` |
| `make_2dElDepthPlot`, `make_grouped_ElDepth2d_Plot` | `ccp.plot_el_depth_2d`, `grouped_el_depth_2d_plot` |
| `compute_2dElDepthMedianAndExobase`, `compute_grouped_ElDepthMedianExobase` | `cc.analysis.median_depth_and_exobase` |
| `make_eDepth_QVec_2d_AllCases_Plot` | `ccp.el_depth_angle_maps_plot` |
| `make_eDepth_QVec_summaryPlot` | `ccp.el_depth_summary_plot` |
| `plot_eDepth_dist` | `ccp.el_depth_distribution_plot` |
| `inelasticXSplots.m` (cross sections, loss rates) | `ccp.plot_cross_sections`, `plot_loss_rates` |
| `plot_trajectories_v3.m` (file reading, 3D trajectory) | `cc.io.load_trajectories`, `load_eb_fields`, `ccp.plot_trajectory_3d` |

Where MATLAB had several versions of a function (`_v2`, `_v3`, non-MC variants) the newest one
called from `asymmetricPlots.m` was ported.

### Not ported

These scripts depend on `.mat` workspaces, VTK field dumps or data outside this repository and
were left as MATLAB: `Chapter5_LowOGR.m`, `FieldComp*.m`, `collvsNonColl_lowOGR.m`,
`ExBplots*.m`, `exobaseEstimatePlot*.m`, `AsymmetricFields_Jul23.m`, `SWupstreamboundary.m`,
`Traj_video_v2.m` and the drift analysis in `plot_trajectories_v2/v3.m`, `scatteringPlots.m`,
`photelectrondistPlots.m`, and the `+coldEPlots`, `+photoelectronProd` and `Cross_Sections`
folders.

## Tests

```bash
pip install -e ".[dev]"
pytest                                             # ~10 s
COLLIDE_COMET_TEST_NOTEBOOKS=1 pytest tests/test_notebooks.py   # also execute the notebooks
```

The tests run on synthetic files written by `collide_comet.synthetic`; the port has not yet been
compared against MATLAB output on real run files.
