# Catchment input template

Copy this directory to a new model directory and fill its CSV tables with your
catchment data. The tables contain headers only.

In `config.json`, replace the required `null` values in `dt`, `channel.slope`,
and `parameters.Hill.iota_rad` with the forcing interval in seconds, channel
slope in m/m, and hillslope angle in radians. Choose `forcing.lai_basis` from
the definition of your leaf area index data. Review the remaining physical
parameters and initial conditions. Their supplied values are implementation
defaults, not estimates for your catchment.

Other `null` values in the default parameter dictionary have a different role.
`Veg.d`, `Veg.z0m`, and `Veg.z0h` use the canopy-height relations when null.
`Infil.K_mps` and `Infil.S_mps05` use the soil-based infiltration expressions
when null. These five optional values may remain null.

File paths are relative to `config.json`. Record the source of the input data
and parameter estimates in `metadata`. Keep an unchanged copy of the original
data and record any unit conversion or gap treatment in your project notes.

The repository [Quick Start](https://github.com/mariaandradeal/1D_hsb-SM-python#quick-start) defines the input
columns, units, configuration choices, and execution commands.
