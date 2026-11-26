# Recorded validation

Numerical checks were executed locally on 24 November 2026, using Windows and
Python 3.12. The setup workflow and MATLAB/manual alignment were checked on 25 November 2026.
The numerical dependency versions are recorded in `requirements-tested.txt`.
The MATLAB baseline is commit `8e915a9d1669b37ff9c4a65fd740bcc1d67f551b`.

## Automated checks

`python -m pytest -q`: **35 passed**. The suite includes component comparisons
against the unmodified MATLAB functions executed in GNU Octave 11.3.0,
independent conservation checks, ADE-kernel integration against numerical
quadrature, input validation, and calibration/validation separation.
The complete original MATLAB pipeline was not executed. The 25 November
revision adds a MATLAB/Octave fixture for the Troch capillary-length override
and checks that it affects Troch mode only. Four invalid-override cases are
also rejected. The manual now gives the active process equations and options.

The revised Quick Start notebook was validated as nbformat 4. All four code
cells were executed in fresh Python kernels with both an unconfigured template
and completed catchment inputs. The unconfigured template produced no
simulation or output files. For the completed-input check, seven days of the
existing Wet Beaver Creek forcing and its 44-cell geometry were copied to an
isolated temporary directory. The documented command-line workflow and the
notebook produced matching discharge, with a relative tolerance of 1e-12. The
relative water-balance error was 3.7906e-14. These data were used only to check
the setup workflow; the distributed input tables contain headers only.
The template copy command was also checked to reject an existing destination
without changing its configuration.

## Forward simulations

| Case | Intervals | Relative water-balance error |
|---|---:|---:|
| Synthetic, all five stores | 120 hourly | 1.9113e-10 |
| Wet Beaver Creek, first year | 365 daily | 4.1593e-14 |
| Wet Beaver Creek, full forcing | 6,847 daily | 1.8754e-14 |

The full forcing period is 2005-01-01 through 2023-09-30. Its cumulative
residual is approximately −5.12e-5 m³ against 2.643e9 m³ of precipitation.
The synthetic and field configurations differ in interception and deep-recharge
settings, as documented with each example. Monthly and annual exports were
checked against the integrated whole-run volumes.

These checks establish numerical water accounting, not field predictive skill.
With the uncalibrated source defaults and a 365-day warmup, the full field run
has 4,384 overlapping observations, NSE −1.2425, KGE −0.0416, and PBIAS −55.20%.
The available observed record ends on 2018-01-01; there is no observation-based
evaluation after that date. These results require a fresh site calibration and
independent validation before predictive use.

Reproduce the full run with:

```sh
python -m hsbsm run examples/wet_beaver_creek/config.json --observations examples/wet_beaver_creek/observations.csv --warmup 365 --output outputs/wbc-full
```

## Calibration and sensitivity workflow

A deliberately short calibration check used the first 60 field days, `Aq_k`,
warmup 10, split 0.8, seed 42, `maxiter=1`, and `popsize=1`. All 11 evaluations
completed without trial failures. The objective improved from 1.16535 to
0.95750. The optimizer reported its iteration limit and `converged=false`, as
expected for this small budget. This is a workflow check, not a calibrated
parameter recommendation. The test suite separately checks that changing
validation observations does not change the fitted parameters.

The sensitivity command completed a baseline and ±10% perturbations of
`Aq_k`, `Aq_f`, and `RZ_Ks` on the synthetic case.

## Documentation and distribution

The revised 30-page manual was compiled with pdfLaTeX and BibTeX. The changed
pages were inspected after rendering; the remaining pages had been inspected
on 24 November. No unresolved references, duplicate labels, or overfull boxes
remain. Model plots were inspected in PNG and also exported as PDF/SVG.

The source distribution and wheel were built locally. GitHub Actions checks
Windows/Linux and Python 3.11–3.13; see the [workflow runs](https://github.com/mariaandradeal/1D_hsb-SM-python/actions/workflows/tests.yml)
for the status and logs of each published commit. The numerical results above
were recorded locally on Windows/Python 3.12. Minimum dependency versions have
not been tested as a separate matrix.

See [migration decisions](MIGRATION.md) for the intentional scientific and
numerical differences from the active MATLAB implementation.
