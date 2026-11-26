# hsB–SM for Python

A Python implementation of the coupled hillslope-storage Boussinesq and soil-moisture model by Marcus Nóbrega Gomes Júnior. The model combines canopy interception, Penman–Monteith evapotranspiration, Milly infiltration, root-zone/transmission-zone Darcy exchange, lateral groundwater flow, a deep reservoir, and channel routing.

This repository extends the [MATLAB implementation](https://github.com/marcusnobrega-eng/1D_hsb-SM). Python runs independently of MATLAB and TopoToolbox using prepared geometry tables. The port includes explicit water-accounting corrections described in [Migration and validation](docs/MIGRATION.md). It is not a claim of identical full-model hydrographs or a new calibration for Wet Beaver Creek.

## Quick Start

This guide sets up a model for your catchment. Prepare the geometry, meteorological time series, parameters, and initial conditions before running it. Observed discharge is optional for a forward simulation and is needed to evaluate or calibrate the model.

### 1. Install the model

Install Python 3.11 or newer and Git, then download the repository.

```sh
git clone https://github.com/mariaandradeal/1D_hsb-SM-python.git
cd 1D_hsb-SM-python
```

On Windows PowerShell, create and activate a Python environment:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

On macOS or Linux:

```sh
python3 -m venv .venv
source .venv/bin/activate
```

Install the package in the active environment:

```sh
python -m pip install -e .
```

If PowerShell blocks activation, use `.\.venv\Scripts\python.exe` in place of `python` in the remaining commands. Install with:

```powershell
.\.venv\Scripts\python.exe -m pip install -e .
```

### 2. Prepare the input files

From the repository root, copy the [input template](templates/catchment) to a new directory:

```sh
python -c "from shutil import copytree; copytree('templates/catchment', 'my_model')"
```

The CSV files contain column headers only. Fill them with data for your catchment. The command requires a new destination directory and does not overwrite an existing model.

| File | Required contents |
|---|---|
| `geometry.csv` | `x_m`: cell-center distance from the outlet (m); `w_m`: planform width (m); `dx_m`: cell length (m) |
| `channel.csv` | `distance_m`: channel travel distance to the outlet (m); `weight`: nonnegative probability mass for each distance bin |
| `forcing.csv` | One row per time interval, with the fields and units below |
| `config.json` | File paths, timestep, channel slope, parameters, and initial conditions |

Order hillslope cells from the outlet toward the divide. Distances, widths, and lengths must be positive, and each center must lie within the cell edges defined by cumulative `dx_m`. Channel distances must increase; weights must have a positive sum and are normalized by the model. Use integrated probability masses for unequal channel-bin widths. Supply prepared geometry tables; the Python package does not delineate a catchment from a digital elevation model.

| Forcing column | Meaning and unit |
|---|---|
| `date` | Start of the interval, as `YYYY-MM-DD` or `YYYY-MM-DD HH:MM:SS` |
| `P` | Precipitation rate (m/s) |
| `Temp_C` | Air temperature (°C) |
| `Pres_kPa`, `VapP_kPa` | Air pressure and actual vapor pressure (kPa) |
| `Wind_ms` | Wind speed (m/s) |
| `SRad_Wm2`, `LRad_Wm2` | Incoming shortwave and longwave radiation (W/m²) |
| `Albd` | Surface albedo, from 0 to 1 |
| `LAI` | Leaf area index (m²/m²) |
| `fv` | Vegetated area fraction, from 0 to 1 |
| `SFC` | Optional snow-cover fraction, from 0 to 1; zero if the column is absent |

Use one consistent time basis. Dates must be unique, increasing, and evenly spaced, with no missing required values. If precipitation is expressed in mm/day, name the column `Prec_mm_day` instead of `P`; the reader converts it to m/s. The forcing is spatially uniform across the hillslope. For Excel input, install `python -m pip install -e ".[excel]"` and update the forcing filename in the configuration.

### 3. Set the configuration and parameters

Open `my_model/config.json`. Paths in this file are relative to its directory. Replace the three required `null` placeholders with numbers:

| Setting | Value to supply |
|---|---|
| `dt` | Forcing interval in seconds, equal to the spacing between timestamps |
| `channel.slope` | Positive channel slope (m/m) |
| `parameters.Hill.iota_rad` | Hillslope angle in radians |

Set `forcing.lai_basis` to `"ground"` if LAI is defined per unit ground area, or `"vegetated"` if it is defined per unit vegetated area. The reader divides ground-area LAI by `fv` once. Ground-area LAI must be zero wherever `fv` is zero.

The template lists the implementation defaults for the other parameters. Review them for your catchment before running. In particular, set soil-layer depths and hydraulic properties in `Soil.RZ` and `Soil.TZ`, aquifer conductivity, drainable porosity, and thickness in `Aq`, and channel velocity and routing parameters in `Route`. Conductivities are in m/s, depths are in m, and soil water contents are in m³/m³. Van Genuchten `alpha` is in m⁻¹ and `n` is dimensionless. The defaults are not calibrated catchment estimates.

Choose `Flags.interception` explicitly. `Deep.phi_d` controls the fraction of downward recharge entering the deep reservoir. Review the initial soil water contents (`theta0` in each soil layer). The optional `initial` object overrides the starting stores; accepted fields and defaults are listed under [Inputs and parameters](#inputs-and-parameters).

Six optional parameter values may remain `null`: `Veg.d`, `Veg.z0m`, and `Veg.z0h` use the canopy-height relations, while `Infil.K_mps` and `Infil.S_mps05` use the soil-based infiltration expressions. `Troch.alpha_inv_m` defaults to the inverse of the root-zone van Genuchten alpha; a positive override is a length in metres and affects Troch mode only.

### 4. Check the inputs and run

After filling the CSV files and configuration, check that the model can load them:

```sh
python -c "from hsbsm import load_config; c = load_config('my_model/config.json'); print(len(c.forcing), 'intervals;', len(c.grid.x), 'cells')"
```

Resolve any reported input errors, then run the model:

```sh
python -m hsbsm run my_model/config.json --output outputs/my_model
```

### 5. Inspect the first run

Check `outputs/my_model/diagnostics.json` for the water-balance residual and solver diagnostics. Inspect `timeseries.csv` and `hydrograph.png` for discharge, evapotranspiration, and soil water content. Confirm that their units, dates, and magnitudes are consistent with the inputs. A closed water balance checks numerical accounting; observed data are needed to assess predictive performance. Use a different output directory for each experiment. `--no-plots` skips figure generation.

## Outputs

A completed command-line run creates:

| Output | Contents |
|---|---|
| `timeseries.csv` | Outlet flow, ET, storage, and per-step water balance |
| `states_and_fluxes.npz` | Node states and fluxes, with geometry; load with `numpy.load` |
| `diagnostics.json` | Whole-run balance, routing tail, and solver diagnostics |
| `parameters.json`, `run_metadata.json` | Parameters and run provenance |
| `signatures.json` | Flow quantiles, runoff ratio, and drought duration |
| `monthly_climatology.csv`, `annual_water_totals.csv` | Monthly mean rates and simulated calendar-year volumes |
| `hydrograph.*`, `states_and_duration.*` | PNG previews and editable SVG/PDF figures |

Reusing an output directory replaces files with the same names.

## Wet Beaver Creek

The prepared example contains 6,847 daily forcing records (2005-01-01 through 2023-09-30), 44 hillslope cells, and channel-width weights. Start with one year:

```sh
python -m hsbsm run examples/wet_beaver_creek/config.json --steps 365 --output outputs/wbc-year --observations examples/wet_beaver_creek/observations.csv --warmup 30
```

Omit `--steps 365` to use the entire forcing record. These are the MATLAB parameter defaults, not parameters calibrated for this Python implementation. The active MATLAB driver disables interception and uses zero deep-recharge fraction. Those choices are preserved in the Wet Beaver Creek configuration.

The supplied observations end on 2018-01-01 (4,749 daily values). Later simulated dates have no observed-flow metrics. For calibration over the full observed period, use `--steps 4749`; a split extending beyond the available observations is rejected rather than filled artificially.

See [input provenance](examples/wet_beaver_creek/README.md) before interpreting field results. A snow-cover fraction reduces potential bare-soil evaporation, as in the source driver. There is no snow accumulation/melt model.

## Python API

```python
from hsbsm import load_config, run

cfg = load_config("my_model/config.json")
result = run(cfg)
result.save("outputs/my_model")
print(result.timeseries[["date", "Q_total_m3s"]].head())
print(result.diagnostics["relative_water_balance_error"])
```

The [setup notebook](examples/quickstart.ipynb) loads your completed configuration and runs the same workflow. It requires your input files and contains no demonstration simulation. Install `python -m pip install -e ".[notebook]"` and run `jupyter lab` if you prefer notebooks.

## Inputs and parameters

Configuration is JSON. All file paths inside it are relative to the configuration file. Start from `templates/catchment`, fill the input tables, and review the parameters.

- Geometry columns are `x_m` (cell center), `w_m` (planform width), and `dx_m` (cell length). Cells start at the outlet and progress toward the divide. Width and length are positive.
- Channel columns are `distance_m` and `weight`. Weights are normalized. Use bin probability masses for nonuniform distance bins; the bundled uniform-bin density weights are normalized equivalently.
- Forcing requires `date`, `P`, `Temp_C`, `Pres_kPa`, `Wind_ms`, `VapP_kPa`, `SRad_Wm2`, `LRad_Wm2`, `Albd`, `LAI`, and `fv`. `P` is m/s; `Prec`/`Prec_mm_day` are accepted aliases in mm/day. Temperature is °C, pressure and vapor pressure are kPa, radiation is W/m², wind is m/s. Optional `SFC` is a fraction in [0,1].
- Meteorology is spatially uniform; stores and resulting fluxes vary by hillslope cell. LAI can be specified over vegetated area (`lai_basis="vegetated"`) or ground area (`"ground"`, converted using `fv`).
- CSV is the default input format. To read XLSX forcing, install `python -m pip install -e ".[excel]"`.
- Dates must be increasing, unique, and equally spaced by `dt` seconds. Missing forcing values fail validation. Observations may have gaps and are matched by timestamp without extrapolation or index substitution.

Parameter groups retain the MATLAB names `Hill`, `Met`, `Veg`, `Soil.RZ`, `Soil.TZ`, `Infil`, `Aq`, `Deep`, `Route`, `Num`, and `Flags`. Unknown parameters are rejected. Optional initial states are `theta_rz`, `theta_tz`, `canopy_m`, `aquifer_head_m`, `deep_m`, and `infiltration_m`; each accepts a scalar or one value per cell. Soil water contents are in m³/m³; the other initial states are in m. Without overrides, soil water contents use each layer's `theta0`, canopy and aquifer storage start at zero, deep storage starts at 1 mm, and cumulative infiltration starts at zero. A nonzero initial canopy store requires interception to be enabled.

## Calibration and sensitivity

Calibration minimizes `1 - KGE` (2009) with SciPy differential evolution, a fixed random seed, and a chronological split. The optimization algorithm differs from MATLAB's genetic algorithm. The default 80% calibration window and 20% validation window are separate; `--warmup` excludes initial steps from metrics. Conductivities are searched in logarithmic space.

```sh
python -m hsbsm calibrate examples/wet_beaver_creek/config.json --steps 365 --observations examples/wet_beaver_creek/observations.csv --active Aq_k Aq_f --warmup 30 --maxiter 10 --popsize 5 --seed 42 --output outputs/calibration
python -m hsbsm sensitivity my_model/config.json --active Aq_k Aq_f RZ_Ks --output outputs/sensitivity
```

Calibration runs many simulations and takes longer than a forward run. `calibration.json` records bounds, trial failures, optimization status, and calibration/validation metrics. Reaching the iteration budget is reported explicitly and is not treated as proof of convergence. Use `hsbsm.calibration.calibrate(..., bounds={"Aq_k": [1e-5, 1e-3]})` to supply site-specific bounds. The registry is in `hsbsm.calibration.REGISTRY`.

## Documentation and checks

The [technical manual](docs/manual/main.pdf) extends the original LaTeX manual with Python installation, inputs, outputs, APIs, examples, and the corrected coupling. Its editable source is [main.tex](docs/manual/main.tex).

With a TeX distribution installed, rebuild it using `python tools/build_manual.py`. This runs pdfLaTeX and BibTeX directly and does not require Perl.

```sh
python -m pip install -e ".[dev]"
python -m pytest -q
```

Tests include reference outputs produced by executing unmodified MATLAB modules in GNU Octave, independent water-balance checks, routing volume, dry/saturated limits, input errors, and calibration isolation. See [validation details](docs/MIGRATION.md). A MATLAB/Octave runtime is unnecessary for ordinary tests.

`examples/synthetic` contains prescribed inputs used by the automated tests. It is not an observed catchment dataset or the starting point for the Quick Start.

DEM delineation and historical exploratory MATLAB scripts are outside this first Python release. Supply prepared tables, or use the original MATLAB/TopoToolbox preprocessing and export them. No MATLAB runtime is used by the Python model.

## Citation and license

Gomes Júnior, M. N. *hsB–SM: coupled hillslope-storage Boussinesq and soil-moisture model, Python implementation*, version 0.1.0 (2026). Cite the repository and the exact commit or release used in an analysis. The original MATLAB source is pinned in `docs/MIGRATION.md`; no DOI has been assigned here.

Code and documentation are distributed under the [MIT License](LICENSE). Example data retain their source attribution; see the data provenance notes.

The [MATLAB/Python/manual comparison](docs/MATLAB_COMPATIBILITY.md) records the active equations, corrected documentation, supported options, and intentional numerical differences.
