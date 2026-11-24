# hsB–SM for Python

A Python implementation of the coupled hillslope-storage Boussinesq and soil-moisture model by Marcus Nóbrega Gomes Júnior. The model combines canopy interception, Penman–Monteith evapotranspiration, Milly infiltration, root-zone/transmission-zone Darcy exchange, lateral groundwater flow, a deep reservoir, and channel routing.

This repository extends the [MATLAB implementation](https://github.com/marcusnobrega-eng/1D_hsb-SM). Python runs independently of MATLAB and TopoToolbox using prepared geometry tables. The port includes explicit water-accounting corrections described in [Migration and validation](docs/MIGRATION.md). It is not a claim of identical full-model hydrographs or a new calibration for Wet Beaver Creek.

## Quick Start

Install Python 3.11 or newer and Git. The commands below run a five-day synthetic storm example with all five stores enabled.

```sh
git clone https://github.com/marcusnobrega-eng/1D_hsb-SM-python.git
cd 1D_hsb-SM-python
python -m venv .venv
```

Activate the environment on Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

On macOS or Linux:

```sh
source .venv/bin/activate
```

Then install and run:

```sh
python -m pip install -e .
python -m hsbsm run examples/synthetic/config.json --output outputs/quickstart
```

If PowerShell blocks activation, use the environment directly without changing execution policy:

```powershell
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m hsbsm run examples/synthetic/config.json --output outputs/quickstart
```

The terminal reports progress and the relative water-balance error. Open `outputs/quickstart/hydrograph.png` to see rainfall, ET, discharge components, and soil moisture. A successful run creates:

| Output | Contents |
|---|---|
| `timeseries.csv` | Outlet flow, ET, storage, and per-step water balance |
| `states_and_fluxes.npz` | Node states and fluxes, with geometry; load with `numpy.load` |
| `diagnostics.json` | Whole-run balance, routing tail, and solver diagnostics |
| `parameters.json`, `run_metadata.json` | Parameters and run provenance |
| `signatures.json` | Flow quantiles, runoff ratio, and drought duration |
| `monthly_climatology.csv`, `annual_water_totals.csv` | Monthly mean rates and simulated calendar-year volumes |
| `hydrograph.*`, `states_and_duration.*` | PNG previews and editable SVG/PDF figures |

Use a different `--output` directory for each experiment. Reusing a directory replaces files with the same names. `--no-plots` skips figure generation.

## Wet Beaver Creek

The prepared example contains 6,847 daily forcing records (2005-01-01 through 2023-09-30), 44 hillslope cells, and channel-width weights. Start with one year:

```sh
python -m hsbsm run examples/wet_beaver_creek/config.json --steps 365 --output outputs/wbc-year --observations examples/wet_beaver_creek/observations.csv --warmup 30
```

Omit `--steps 365` to use the entire forcing record. These are the MATLAB parameter defaults, not parameters calibrated for this Python implementation. The active MATLAB driver disables interception and uses zero deep-recharge fraction. Those choices are preserved in the Wet Beaver Creek configuration; the synthetic configuration explicitly enables interception and sets `Deep.phi_d` to 0.2.

The supplied observations end on 2018-01-01 (4,749 daily values). Later simulated dates have no observed-flow metrics. For calibration over the full observed period, use `--steps 4749`; a split extending beyond the available observations is rejected rather than filled artificially.

See [input provenance](examples/wet_beaver_creek/README.md) before interpreting field results. A snow-cover fraction reduces potential bare-soil evaporation, as in the source driver. There is no snow accumulation/melt model.

## Python API

```python
from hsbsm import load_config, run

cfg = load_config("examples/synthetic/config.json")
cfg.params["Aq"]["k"] = 1e-3
result = run(cfg)
result.save("outputs/my_run")
print(result.timeseries[["date", "Q_total_m3s"]].head())
print(result.diagnostics["relative_water_balance_error"])
```

The [walkthrough notebook](examples/quickstart.ipynb) shows the same workflow. Install `python -m pip install -e ".[notebook]"` and run `jupyter lab` if you prefer notebooks.

## Inputs and parameters

Configuration is JSON. All file paths inside it are relative to the configuration file. Copy an example folder to create a new experiment, then edit the geometry, forcing, channel table, and parameter overrides.

- Geometry columns are `x_m` (cell center), `w_m` (planform width), and `dx_m` (cell length). Cells start at the outlet and progress toward the divide. Width and length are positive.
- Channel columns are `distance_m` and `weight`. Weights are normalized. Use bin probability masses for nonuniform distance bins; the bundled uniform-bin density weights are normalized equivalently.
- Forcing requires `date`, `P`, `Temp_C`, `Pres_kPa`, `Wind_ms`, `VapP_kPa`, `SRad_Wm2`, `LRad_Wm2`, `Albd`, `LAI`, and `fv`. `P` is m/s; `Prec`/`Prec_mm_day` are accepted aliases in mm/day. Temperature is °C, pressure and vapor pressure are kPa, radiation is W/m², wind is m/s. Optional `SFC` is a fraction in [0,1].
- Meteorology is spatially uniform; stores and resulting fluxes vary by hillslope cell. LAI can be specified over vegetated area (`lai_basis="vegetated"`) or ground area (`"ground"`, converted using `fv`).
- CSV is the default input format. To read XLSX forcing, install `python -m pip install -e ".[excel]"`.
- Dates must be increasing, unique, and equally spaced by `dt` seconds. Missing forcing values fail validation. Observations may have gaps and are matched by timestamp without extrapolation or index substitution.

Parameter groups retain the MATLAB names `Hill`, `Met`, `Veg`, `Soil.RZ`, `Soil.TZ`, `Infil`, `Aq`, `Deep`, `Route`, and `Num`. Unknown parameters are rejected. Optional initial states are `theta_rz`, `theta_tz`, `canopy_m`, `aquifer_head_m`, `deep_m`, and `infiltration_m`; each accepts a scalar or one value per cell. Initial deep storage defaults to 1 mm, matching the source.

## Calibration and sensitivity

Calibration minimizes `1 - KGE` (2009) with SciPy differential evolution, a fixed random seed, and a chronological split. The optimization algorithm differs from MATLAB's genetic algorithm. The default 80% calibration window and 20% validation window are separate; `--warmup` excludes initial steps from metrics. Conductivities are searched in logarithmic space.

```sh
python -m hsbsm calibrate examples/wet_beaver_creek/config.json --steps 365 --observations examples/wet_beaver_creek/observations.csv --active Aq_k Aq_f --warmup 30 --maxiter 10 --popsize 5 --seed 42 --output outputs/calibration
python -m hsbsm sensitivity examples/synthetic/config.json --active Aq_k Aq_f RZ_Ks --output outputs/sensitivity
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

DEM delineation and historical exploratory MATLAB scripts are outside this first Python release. Supply prepared tables, or use the original MATLAB/TopoToolbox preprocessing and export them. No MATLAB runtime is used by the Python model.

## Citation and license

Gomes Júnior, M. N. *hsB–SM: coupled hillslope-storage Boussinesq and soil-moisture model, Python implementation*, version 0.1.0 (2026). Cite the repository and the exact commit or release used in an analysis. The original MATLAB source is pinned in `docs/MIGRATION.md`; no DOI has been assigned here.

Code and documentation are distributed under the [MIT License](LICENSE). Example data retain their source attribution; see the data provenance notes.
