# Wet Beaver Creek example data

These files were derived from the public MATLAB repository
[marcusnobrega-eng/1D_hsb-SM](https://github.com/marcusnobrega-eng/1D_hsb-SM)
at commit `8e915a9d1669b37ff9c4a65fd740bcc1d67f551b`. SHA-256 checksums of the
source files are recorded in `source_hashes.json`.

- `geometry.csv` and `channel.csv` come from `Grid` and `Ch` in `pre_hsbsm.mat`.
  The source file and the separately committed width CSV files describe different
  prepared geometries. This example consistently uses the MAT-file geometry,
  with hillslope area 283,941,233.215 m². Channel area is not added a second time.
- `forcing_raw.csv` is a direct CSV export of
  `Input/Met_Forcing_Patrick_VegModis_SCF.xlsx`. It retains all original fields.
- `forcing.csv` differs in vegetation fields and one invalid albedo. Eight missing LAI entries,
  one missing FVC entry, and one zero FVC entry with positive LAI are replaced
  by linear interpolation between adjacent valid daily values. The albedo value
  -32.767 on 2016-11-23 is treated as invalid and interpolated in the same way.
  Every modified value is listed in `forcing_corrections.csv`. Precipitation and observed
  discharge are not interpolated. This preparation is an example-data choice,
  not an implicit behavior of the model reader.
- `LAI` in this source is ground-area LAI. The reader divides by FVC once to
  obtain vegetated-patch LAI. `SFC` is used only to reduce bare-soil potential
  evaporation; precipitation is assumed liquid water entering the model.
- `observations.csv` is the timestamp-preserving subset of
  `calibration/data/Streamflow_WetBeaverCreek_cms.csv` within the forcing period.
  Its discharge unit is m³/s. Missing observations are excluded from metrics.
  The supplied observation record ends on 2018-01-01, leaving 4,749 daily
  values in this example. The meteorological record continues through 2023-09-30.

The original repository does not provide a complete upstream data-license or
station/forcing-production record. These are research examples attributed to
the original repository, not a newly independently verified field dataset.
Check the original data provenance before reusing them in a publication.

Regenerate the files after installing the Excel extra:

```sh
python -m pip install -e ".[excel]"
python tools/prepare_examples.py /path/to/1D_hsb-SM
```

Use the pinned source commit above. The script does not download remote data or
run DEM processing. Default parameters are demonstration values. Calibration,
spinup selection, forcing interpretation, and independent validation are still
required for a site-specific scientific analysis.
