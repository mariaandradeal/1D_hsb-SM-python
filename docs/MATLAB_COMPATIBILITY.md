# MATLAB, Python, and manual comparison

This comparison uses MATLAB commit
[`8e915a9`](https://github.com/marcusnobrega-eng/1D_hsb-SM/tree/8e915a9d1669b37ff9c4a65fd740bcc1d67f551b).
It follows the active `run_hsbsm_pipeline` → `hsbsm_config` → `hsbsm_run`
path. A function present elsewhere in the repository is not necessarily called
by that path. The revised [manual](manual/main.pdf) describes the Python code
and identifies intentional differences from MATLAB.

## One-by-one comparison

| Item | Active MATLAB code | Python code | Earlier manual discrepancy | Resolution |
|---|---|---|---|---|
| 1. Infiltration-excess runoff | `hsbsm_infil_milly` returns `Pnet-I`; its saturation-excess output is zero. | `infiltration_step` does the same. | Surface-routing section defined both excess terms as `Pnet-I`, counting the same water twice. | Manual now separates infiltration excess from soil excess and aquifer overflow. No extra vegetation weighting is applied to ground-area fluxes. |
| 2. Aquifer overflow | `hsbsm_run` adds overflow to surface inflow and later adds it again to outlet discharge. | Overflow enters surface inflow once, after optional TZ/RZ refill. | A corrected note coexisted with inconsistent runoff equations. | Python's conservation correction is retained; equations and notes now agree. Diagnostic overflow must not be added to total discharge. |
| 3. Canopy drip | `hsbsm_canopy_pm` supports `drip_coef`; values below one leave excess that is discarded when storage is clipped. Interception is hard-coded off in the active runner. | `Flags.interception` defaults to false. When enabled, all excess drains, equivalent to MATLAB `drip_coef=1`. | Adjustable fractional drip was presented as conservative for all coefficients. | Manual specifies immediate complete drip and explains why fractional drip is unsupported. Turning interception off does not disable transpiration or soil evaporation. |
| 4. Infiltration and Troch option | Capacity is `Keff + S²/(2 max(F,Fmin))`. Troch mode can use `Troch.alpha_inv_m`. | Same capacity. This revision adds the capillary-length override. | Manual described an override that Python rejected; dry reset wording implied a reset to `F0` alone. | Override is now supported and checked against the unchanged MATLAB function. Its unit is m; effective alpha is its reciprocal. Dry counter resets when it exceeds the threshold: default fourth dry step, to `max(F0,Fmin)=0.01 m`. |
| 5. Routing coefficients | Normal parameter construction derives `c,d²` from `V,F,a0` and channel slope, overwriting direct coefficients when channel data are present. A standalone routing wrapper has a different interface. | Derives the same coefficients. Integrates the travel-time density over bins and retains undelivered channel water. | Manual offered direct `c,d²` configuration. | Manual lists the accepted inputs. Bin integration and channel-tail accounting remain explicit numerical differences. |
| 6. Coordinates and recharge | `x` increases from outlet to divide; negative outlet face flux is outward. The active split uses `TZ.deep_frac`, default zero. | Same signs; one canonical `Deep.phi_d`, default zero. | Reference table said downslope coordinate, negative inward flux, and `NhsB=gtz` without the deep split or upward return. | Table now gives `NhsB=(1-phi_d)q+ - q-` and `Ndeep=phi_d q+`, with the correct coordinate and signs. |
| 7. Aquifer water available to soil | The soil coupling divides `S` [m²] by node area when limiting upward transfer. | Uses `S/w` [m], consistent with `S=fwh`. | General descriptions did not distinguish depth from volume. | Python's dimensional correction is retained. Manual states depth `S/w` and cell volume `S*dx`. |
| 8. Soil ET and capillary return | Active runner calls `hsbsm_rz_tz_coupled`: signed Darcy exchange and pressure-head stress. Competing ET demands may independently use the same initial water, followed by clipping. | Same unconstrained equations, with shared remaining donor water and explicit soil excess runoff. | Older stage-2 evaporation and capillary-demand descriptions remained in tables and the diagram caption. | Active equations replace those descriptions. Python keeps conservative donor limits; legacy clipping behavior is not reproduced. |
| 9. Penman–Monteith calculation | Active runner uses common net radiation and aerodynamic resistance, emissivity 0.98, zero ground heat fraction, and prescribed surface resistances. | Same active formulas, verified with a component fixture. | Detailed equations suggested separate surface radiation, different constants, stability corrections, and optional Jarvis/stage-2 schemes. | Manual now gives the actual Tetens, density, radiation, resistance, and PM equations. No independent surface-temperature or sensible-heat closure is claimed. |
| 10. Variable names and outputs | MATLAB structures and module-local diagnostics differ from runner exports. | `Config`, `Result.timeseries`, and `Result.nodes`; intermediate canopy and infiltration quantities are local. | Some MATLAB-style names looked like available Python output fields. | Tables now identify Python parameters, exported fields, and local quantities separately. |
| 11. Parameter overrides | `run_hsbsm_pipeline` constructs `P` from `POverrides`, then constructs `Cfg` without assigning `Cfg.Params=P` before running. Other entry points may set configuration parameters directly. | JSON overrides are merged, validated, used, and exported. | No warning that the named pipeline can lose its constructed overrides. | Difference is documented. It is not evidence that every MATLAB calibration ignored its settings. |
| 12. Numerical controls and future processes | Source defaults include 30 Picard iterations, 10 s minimum step, and separate timestep controls. Constant drainable porosity and two soil buckets are active. | 60 Picard iterations, 0.001 s minimum step, checked convergence and retries; same constant-porosity/bucket structure. | General numerical and process options could be read as available settings. | Manual lists supported Python controls. Dynamic porosity and a vertically resolved Richards solver are not implemented. |

## Source locations

The MATLAB links below are pinned to the reviewed commit:

- [Canopy update](https://github.com/marcusnobrega-eng/1D_hsb-SM/blob/8e915a9d1669b37ff9c4a65fd740bcc1d67f551b/hsbsm_canopy_pm.m#L114),
  [infiltration overrides and capacity](https://github.com/marcusnobrega-eng/1D_hsb-SM/blob/8e915a9d1669b37ff9c4a65fd740bcc1d67f551b/hsbsm_infil_milly.m#L70).
- [Active soil call and deep split](https://github.com/marcusnobrega-eng/1D_hsb-SM/blob/8e915a9d1669b37ff9c4a65fd740bcc1d67f551b/hsbsm_run.m#L374),
  [overflow in surface inflow](https://github.com/marcusnobrega-eng/1D_hsb-SM/blob/8e915a9d1669b37ff9c4a65fd740bcc1d67f551b/hsbsm_run.m#L490),
  [outlet assembly](https://github.com/marcusnobrega-eng/1D_hsb-SM/blob/8e915a9d1669b37ff9c4a65fd740bcc1d67f551b/hsbsm_run.m#L594).
- [Parameter construction](https://github.com/marcusnobrega-eng/1D_hsb-SM/blob/8e915a9d1669b37ff9c4a65fd740bcc1d67f551b/hsbsm_params.m#L211),
  [pipeline handoff](https://github.com/marcusnobrega-eng/1D_hsb-SM/blob/8e915a9d1669b37ff9c4a65fd740bcc1d67f551b/run_hsbsm_pipeline.m#L98).

Python implementation: [surface and soil processes](../src/hsbsm/physics.py),
[runner and outputs](../src/hsbsm/model.py), [aquifer](../src/hsbsm/aquifer.py),
[routing](../src/hsbsm/routing.py), and [configuration](../src/hsbsm/config.py).

## Supported options and limits of agreement

`Troch.alpha_inv_m` is optional and positive, in metres. `null` uses the
root-zone van Genuchten alpha. It affects the internally estimated sorptivity
only in Troch mode. A supplied `Infil.S_mps05` takes precedence. The source's
`Troch.n` and `Troch.no` are read but do not affect this infiltration calculation;
they are not Python options. Alternative MATLAB floor controls (`fmax_mps`,
`Fmin_default_m`) and `Soil.Se_floor` overrides are also not exposed. Python
uses explicit positive `Fmin_m` and a fixed saturation floor of 1e-6.

Python accepts zero conductivity and zero sorptivity as explicit values. Some
MATLAB configuration paths replace these values with defaults or an internal
estimate. Use `null` in Python for automatic values. Unknown parameter names
are rejected, including legacy names that would otherwise have no effect.

The ordinary Milly defaults are unchanged by this revision. Component agreement
has been checked by executing unmodified MATLAB files in GNU Octave 11.3.0.
The test fixture now includes the Troch capillary-length case, in addition to
PM, canopy, Milly infiltration, unconstrained soil exchange, and aquifer cases.
The full pipeline has not been run in MATLAB. Conservation corrections and
different routing/time integration can change complete hydrographs, so this
comparison does not establish full-run numerical identity or field calibration.

For the reason and validation behind each conservation correction, see
[MIGRATION.md](MIGRATION.md) and [VALIDATION.md](VALIDATION.md).
