# MATLAB to Python migration and validation

## Source and implementation scope

The source is [1D_hsb-SM](https://github.com/marcusnobrega-eng/1D_hsb-SM), commit
`8e915a9d1669b37ff9c4a65fd740bcc1d67f551b`. The active chain is
`run_hsbsm_pipeline` → `hsbsm_config` → `hsbsm_run`. The current runner calls
`hsbsm_rz_tz_coupled` directly; older independent root-zone/TZ wrappers are not
the active soil coupling. The Python version follows that chain's process
equations and makes the switches explicit.

| MATLAB component | Python component |
|---|---|
| Prepared `Grid`, `Ch`, `Forc` | `config.Grid`, CSV inputs, `config.load_config` |
| `hsbsm_params`, `hsbsm_initial_states` | Parameter dictionary and explicit `initial` settings |
| `hsbsm_pm`, `hsbsm_canopy_pm` | `physics.penman_monteith`, `physics.canopy_step` |
| `hsbsm_infil_milly` | `physics.infiltration_step` (Milly and its Troch switch) |
| `hsbsm_rz_tz_coupled` | `physics.soil_step` |
| `hsB_solver`, including local `hsb_core` | `aquifer.fixed_step`, `aquifer.advance` |
| CIUH construction and convolution | `routing.unit_hydrograph`, `routing.route_flow` |
| `hsbsm_run`, pipeline | `model.run`, `Result.save`, command-line interface |
| Calibration, metrics, OAT sensitivity | `calibration`, `evaluation` |
| Hydrographs, spatial states, flow duration | `plotting` |

Historical alternative formulations, GUI/live plotting, every original
climatology dashboard, and DEM/TopoToolbox preprocessing have not been ported.
The core workflow covers simulation, bounded calibration, evaluation, sensitivity,
scientific plots, and portable outputs. Monthly mean rates and annual water
totals are exported automatically, together with annual performance metrics
when observations are supplied. See [recorded validation](VALIDATION.md) for
the executed examples and their measured errors.

## Deliberate differences

1. **Water shared by ET and exchange.** The original coupled module can debit
   the same available RZ water separately for transpiration and soil evaporation,
   then clip the updated store. Python debits transpiration first, limits soil
   evaporation to the remaining donor water, and limits Darcy exchange using
   that updated amount. Saturation excess from RZ/TZ is routed as explicit
   runoff instead of being discarded by clipping. The Darcy head, conductivity,
   and stress formulas are unchanged in their unconstrained regime.
2. **Aquifer units.** The hsB state is `S=f*w*h` in m². Its equivalent water
   depth is `S/w`, and its cell volume is `S*dx`. The source soil-coupling limit
   divides `S` by cell area, introducing an extra length factor. Python limits
   capillary withdrawal with `S/w`. Withdrawal occurs before lateral integration;
   nonnegative recharge is supplied during the aquifer substeps.
3. **Checked solver convergence.** The original relaxed Picard face-flux formula
   is retained, including its numerical film and drift stabilization. Python
   retries unconverged steps at a smaller timestep and checks mass residuals.
   Failure raises an error rather than silently producing zero baseflow.
   Supplied grid edges are used consistently. The one-cell grid has a full
   control volume from zero to twice the center coordinate.
4. **Overflow counted once.** The source runner adds aquifer overflow to surface
   routing inflow and adds it again to the final outlet hydrograph. Python adds
   it once, as surface inflow. Optional fill-spill returns node overflow to TZ,
   then RZ, and routes the remainder. The diagnostic overflow column is the
   generated overflow before fill-spill and must not be added to total discharge.
5. **Routing bins and retained water.** Python integrates the same ADE travel-time
   density over each timestep (inverse-Gaussian CDF), rather than resampling
   density values at widely separated timestamps. The kernel has unit mass.
   Undelivered water at the end of the run is reported as channel storage and
   routing-tail volume, so truncation cannot appear as a water loss.
6. **Explicit options.** The active MATLAB runner hard-codes interception off.
   It defaults `TZ.deep_frac` to zero even though `Deep.phi_d` is nonzero in the
   parameter file. The Python defaults reflect the active runner (both processes
   off); one canonical `Deep.phi_d` controls the split. The synthetic example
   explicitly activates both. The default initial deep store remains 1 mm.
7. **Inputs and observations.** Input units and array shapes are checked.
   Unknown overrides and non-finite forcing fail early. There is no fallback
   from calendar alignment to index matching. The bundled example records
   vegetation interpolation separately and reads snow cover from the source
   spreadsheet; the committed MATLAB preprocessed file omitted that field.
8. **Calibration.** SciPy differential evolution replaces MATLAB toolbox
   optimizers. The objective is the active MATLAB objective, `1-KGE` (2009).
   Calibration bounds are explicit and ordered, including the reversed RZ_Ks
   source bounds. Trial configurations and initial states are independent.
   No optimization reproducibility claim is made across MATLAB and Python.

These corrections change some coupled hydrographs. They are not hidden
roundoff differences. A fresh calibration is necessary before scientific use.
No new process equation, DEM algorithm, or snowpack model has been inferred from
the original diagrams.

## Reproducing the numerical checks

`tests/fixtures/matlab_reference.mat` was generated by executing the unmodified
source modules in GNU Octave 11.3.0. It includes PM potentials, canopy storage,
Milly infiltration, an unconstrained RZ/TZ exchange case, and ten hsB steps on a
nonuniform grid. Python tests compare those values at relative tolerances of
1e-11 to 1e-12 with absolute tolerances near numerical roundoff.

To regenerate with MATLAB or Octave, clone and check out the source commit, then
run from this repository:

```matlab
addpath('tools');
export_matlab_reference('/path/to/1D_hsb-SM', ...
                        'tests/fixtures/matlab_reference.mat');
```

The ordinary test suite requires only Python. It also tests canopy and soil
water balances, saturated/dry stores, aquifer overflow and capillary extraction,
routing impulse volume and tail, selected-component routing, fill-spill, invalid
forcing, timestamp alignment, metric degeneracies, and calibration-window
isolation. These tests address numerical implementation. They do not establish
field predictive skill, parameter uniqueness, or convergence of the physical
solution with respect to time and spatial resolution.

The complete legacy pipeline was not executed under MATLAB. Full-system
bitwise parity is neither demonstrated nor intended because of the corrections
above. The two-language comparison is limited to the named component cases.

## Figure conventions

Plots use the requested muted scientific palette and 1.5-point axes. The font
selection order is Helvetica Light, Helvetica Neue Light, Helvetica, Arial,
then DejaVu Sans, depending on installed fonts. SVG keeps text editable; PDF
and PNG are generated from the same numerical output. No external LaTeX
installation is needed to draw model figures.
