"""Portable forcing readers and non-pickle output files."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ALIASES = {
    "Temp_C": ("Temp", "Temp_degC"), "Pres_kPa": ("Pres",),
    "Wind_ms": ("Wind", "Wind_m_s"), "VapP_kPa": ("VapP",),
    "SRad_Wm2": ("SRad", "SRad_W_m2"), "LRad_Wm2": ("LRad", "LRad_W_m2"),
    "Albd": ("Albd_unitless",), "LAI": ("Lai", "Lai_m2_m2"),
    "fv": ("VegFrac_LAI", "Fv", "FVC"), "SFC": ("SCF", "snow_cover_fraction"),
}


def read_forcing(path, lai_basis="vegetated"):
    path = Path(path)
    df = pd.read_excel(path) if path.suffix.lower() in (".xls", ".xlsx") else pd.read_csv(path)
    df.columns = df.columns.str.strip()
    for target, aliases in ALIASES.items():
        if target not in df:
            for alias in aliases:
                if alias in df:
                    df = df.rename(columns={alias: target})
                    break
    if "date" not in df:
        for name in ("Date", "DateTime", "datetime", "Time"):
            if name in df:
                df = df.rename(columns={name: "date"})
                break
    if "date" not in df:
        raise ValueError("Forcing must include an explicit date column")
    df["date"] = pd.to_datetime(df.date, errors="raise")
    if "P" not in df:
        column = next((name for name in ("Prec", "Prec_mm_day") if name in df), None)
        if column is None:
            raise ValueError("Supply P [m/s] or Prec/Prec_mm_day [mm/day]")
        df["P"] = df[column].to_numpy(float)/86400000
    if lai_basis == "ground":
        # Ground LAI is zero on bare ground, so division is defined explicitly.
        if ((df.fv == 0) & (df.LAI != 0)).any():
            raise ValueError("Ground-area LAI must be zero wherever fv is zero")
        df["LAI"] = np.divide(df.LAI, df.fv, out=np.zeros(len(df)), where=df.fv.to_numpy() > 0)
    elif lai_basis != "vegetated":
        raise ValueError("lai_basis must be 'ground' or 'vegetated'")
    if "SFC" not in df:
        df["SFC"] = 0.0
    return df


def read_observations(path):
    df = pd.read_csv(path)
    dc = next((c for c in ("date", "Time", "DateTime", "Date") if c in df), None)
    qc = next((c for c in ("Q_m3s", "Discharge_cms", "Q_cms", "discharge_cms") if c in df), None)
    if dc is None or qc is None:
        raise ValueError("Observations need date and Q_m3s columns (legacy Time/Discharge_cms accepted)")
    t = pd.DatetimeIndex(pd.to_datetime(df[dc], errors="raise"))
    if t.hasnans or t.has_duplicates or not t.is_monotonic_increasing:
        raise ValueError("Observation dates must be unique, finite, and increasing")
    values = pd.to_numeric(df[qc], errors="raise").to_numpy(float, copy=True)
    values[values < 0] = np.nan
    return pd.Series(values, index=t, name="Q_obs_m3s")


def json_safe(value):
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer, np.bool_)):
        return value.item()
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (Path, pd.Timestamp)):
        return str(value)
    return value


def write_json(path, data):
    Path(path).write_text(json.dumps(json_safe(data), indent=2, allow_nan=False)+"\n", encoding="utf-8")
