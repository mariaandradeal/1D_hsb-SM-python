"""Rebuild bundled CSV examples from a pinned MATLAB checkout (no DEM needed)."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.io import loadmat

parser = argparse.ArgumentParser()
parser.add_argument("matlab_checkout", type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
source = args.matlab_checkout.resolve()
out = root/"examples"/"wet_beaver_creek"
out.mkdir(parents=True, exist_ok=True)
pre = loadmat(source/"pre_hsbsm.mat", simplify_cells=True)
g, ch = pre["Grid"], pre["Ch"]
pd.DataFrame({"x_m": g["x"], "w_m": g["w"], "dx_m": g["dx"]}).to_csv(out/"geometry.csv", index=False)
pd.DataFrame({"distance_m": ch["c_dist"], "weight": ch["w_norm"]}).to_csv(out/"channel.csv", index=False)
forcing = pd.read_excel(source/"Input"/"Met_Forcing_Patrick_VegModis_SCF.xlsx")
forcing.to_csv(out/"forcing_raw.csv", index=False, date_format="%Y-%m-%d")
corrections = []
for column in ("FVC", "LAI", "Albd_unitless"):
    original = forcing[column].copy()
    invalid = original.isna()
    if column == "FVC":
        invalid |= (original == 0) & (forcing.LAI > 0)
    if column == "Albd_unitless":
        invalid |= ~original.between(0, 1)
    prepared = original.mask(invalid).interpolate(method="linear", limit_area="inside")
    if prepared.isna().any():
        raise ValueError(f"Cannot interpolate unbounded gaps in {column}")
    for index in np.flatnonzero(invalid):
        corrections.append({"date": str(forcing.date.iloc[index].date()), "column": column,
                            "original": original.iloc[index], "prepared": prepared.iloc[index],
                            "rule": "Linear interpolation between adjacent valid daily values"})
    forcing[column] = prepared
pd.DataFrame(corrections).to_csv(out/"forcing_corrections.csv", index=False)
forcing.to_csv(out/"forcing.csv", index=False, date_format="%Y-%m-%d")
obs = pd.read_csv(source/"calibration"/"data"/"Streamflow_WetBeaverCreek_cms.csv")
obs = obs.rename(columns={"Time": "date", "Discharge_cms": "Q_m3s"})
obs.date = pd.to_datetime(obs.date)
obs = obs[(obs.date >= forcing.date.min()) & (obs.date <= forcing.date.max())]
obs.to_csv(out/"observations.csv", index=False, date_format="%Y-%m-%d")
spec = {"geometry": "geometry.csv", "channel": {"file": "channel.csv", "slope": float(ch["S0"])},
        "forcing": {"file": "forcing.csv", "lai_basis": "ground"}, "dt": 86400,
        "parameters": {"Hill": {"iota_rad": float(pre["meta"]["metrics"]["mean_slope_hs_rad"])}},
        "metadata": {"site": "Wet Beaver Creek", "source_commit": "8e915a9d1669b37ff9c4a65fd740bcc1d67f551b",
                     "calibration_status": "Default parameters; not a calibrated Python model",
                     "forcing_preparation": "See forcing_corrections.csv; raw source retained"}}
(out/"config.json").write_text(json.dumps(spec, indent=2)+"\n", encoding="utf-8")
hashes = {str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest() for p in (
    source/"pre_hsbsm.mat", source/"Input"/"Met_Forcing_Patrick_VegModis_SCF.xlsx",
    source/"calibration"/"data"/"Streamflow_WetBeaverCreek_cms.csv")}
(out/"source_hashes.json").write_text(json.dumps(hashes, indent=2)+"\n", encoding="utf-8")

out = root/"examples"/"synthetic"
out.mkdir(parents=True, exist_ok=True)
nx, nt = 12, 120
pd.DataFrame({"x_m": np.arange(nx)*25+12.5, "w_m": np.linspace(30, 12, nx), "dx_m": 25}).to_csv(out/"geometry.csv", index=False)
pd.DataFrame({"distance_m": [250, 750, 1250, 1750], "weight": [0.1, 0.35, 0.4, 0.15]}).to_csv(out/"channel.csv", index=False)
rain = np.zeros(nt)
rain[12:20] = 20
rain[55:59] = 35
hour = np.arange(nt)%24
daylight = np.maximum(0, np.sin((hour-6)*np.pi/12))
df = pd.DataFrame({"date": pd.date_range("2000-01-01", periods=nt, freq="h"),
    "P": rain/3.6e6, "Temp_C": 15+10*daylight, "Pres_kPa": 95., "Wind_ms": 2.,
    "VapP_kPa": 1.2, "SRad_Wm2": 600*daylight, "LRad_Wm2": 330., "Albd": 0.2,
    "LAI": 2., "fv": 0.6, "SFC": 0.})
df.to_csv(out/"forcing.csv", index=False)
spec = {"geometry": "geometry.csv", "channel": {"file": "channel.csv", "slope": 0.02},
    "forcing": {"file": "forcing.csv", "lai_basis": "vegetated"}, "dt": 3600,
    "parameters": {"Flags": {"interception": True}, "Deep": {"phi_d": 0.2},
                   "Hill": {"iota_rad": 0.1}},
    "metadata": {"site": "Synthetic demonstration", "precipitation_total_mm": float(rain.sum())}}
(out/"config.json").write_text(json.dumps(spec, indent=2)+"\n", encoding="utf-8")
print("Prepared synthetic and Wet Beaver Creek inputs")
