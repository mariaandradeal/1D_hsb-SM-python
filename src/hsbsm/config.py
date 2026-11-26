"""Explicit SI geometry and validated configuration; paths are config-relative."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
import json
import numpy as np
import pandas as pd


def default_parameters() -> dict:
    """Active MATLAB defaults, with explicit canopy/deep switches (both off)."""
    return {
        "Hill": {"iota_rad": 0.1},
        "Met": {"zref_m": 10.0},
        "Veg": {"cw_per_LAI": 0.0012, "rc_min": 100.0, "h_c": 20.0,
                "d": None, "z0m": None, "z0h": None, "free_throughfall_f0": 0.0},
        "Soil": {
            "RZ": {"z": 0.65, "theta_min": 0.045, "theta_sat": 0.38, "theta0": 0.30,
                   "vG": {"alpha": 6.0, "n": 2.2, "Ks": 2e-6}},
            "TZ": {"z": 0.35, "theta_min": 0.045, "theta_sat": 0.38, "theta0": 0.30,
                   "vG": {"alpha": 6.0, "n": 2.2, "Ks": 2e-6}, "cap_frac_of_K": 0.30},
            "rs_soil": 157.70,
        },
        "Infil": {"model": "Milly", "K_mps": None, "S_mps05": None, "F0_m": 0.001,
                  "Fmin_m": 0.01, "resetDrySteps": 3, "max_c": 0.7},
        "Troch": {"alpha_inv_m": None},
        "Aq": {"k": 1.468e-3, "f": 0.09, "D": 0.35},
        "Deep": {"phi_d": 0.0, "k_out_s": 1 / (83.333 * 86400)},
        "Route": {"enabled": True, "F": 0.8, "V": 1.0, "a0": 2/3,
                  "route_baseflow": False, "route_deepflow": False,
                  "overflow_fill_hierarchy": False},
        "Num": {"theta": 1.0, "omega": 0.8, "picard_max": 60, "picard_tol": 1e-9,
                "dt_max": 3600.0, "dt_min": 0.001, "cfl": 0.3},
        "Flags": {"interception": False},
    }


def merge_parameters(base: dict, overrides: dict, prefix: str = "") -> dict:
    out = deepcopy(base)
    for key, value in overrides.items():
        if key not in out:
            raise ValueError(f"Unknown parameter: {prefix}{key}")
        if isinstance(out[key], dict):
            if not isinstance(value, dict):
                raise ValueError(f"{prefix}{key} must be an object")
            out[key] = merge_parameters(out[key], value, prefix + key + ".")
        else:
            out[key] = value
    return out


@dataclass
class Grid:
    """Centers x [m], planform widths w [m], and control-volume edges [m]."""
    x: np.ndarray
    w: np.ndarray
    edges: np.ndarray | None = None

    def __post_init__(self):
        self.x = np.asarray(self.x, dtype=float).reshape(-1)
        self.w = np.asarray(self.w, dtype=float).reshape(-1)
        if len(self.x) == 0 or self.w.shape != self.x.shape:
            raise ValueError("Grid x and w must have equal nonzero length")
        if not np.all(np.isfinite(self.x)) or not np.all(np.isfinite(self.w)):
            raise ValueError("Grid values must be finite")
        if np.any(self.w <= 0) or self.x[0] <= 0 or np.any(np.diff(self.x) <= 0):
            raise ValueError("Widths must be positive; centers must increase from x > 0")
        if self.edges is None:
            last = 2*self.x[0] if len(self.x) == 1 else self.x[-1] + np.diff(self.x)[-1]/2
            self.edges = np.r_[0.0, (self.x[1:] + self.x[:-1])/2, last]
        self.edges = np.asarray(self.edges, dtype=float).reshape(-1)
        if (len(self.edges) != len(self.x)+1 or not np.all(np.isfinite(self.edges))
                or np.any(np.diff(self.edges) <= 0) or self.edges[0] != 0
                or np.any(self.x <= self.edges[:-1]) or np.any(self.x >= self.edges[1:])):
            raise ValueError("Edges must increase from zero and enclose every center")

    @property
    def dx(self):
        return np.diff(self.edges)

    @property
    def area(self):
        return self.w * self.dx


FORCING_COLUMNS = ("P", "Temp_C", "Pres_kPa", "Wind_ms", "VapP_kPa", "SRad_Wm2",
                   "LRad_Wm2", "Albd", "LAI", "fv")


@dataclass
class Config:
    grid: Grid
    forcing: pd.DataFrame
    dt: float
    channel: dict = field(default_factory=dict)
    params: dict = field(default_factory=default_parameters)
    initial: dict = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    def validate(self):
        if not np.isfinite(self.dt) or self.dt <= 0:
            raise ValueError("dt must be finite and positive [s]")
        if self.forcing.empty:
            raise ValueError("Forcing is empty")
        missing = set(FORCING_COLUMNS) - set(self.forcing.columns)
        if missing:
            raise ValueError(f"Missing forcing columns: {sorted(missing)}")
        if not np.isfinite(self.forcing[list(FORCING_COLUMNS)].to_numpy(float)).all():
            raise ValueError("Forcing contains missing or non-finite values")
        for name in ("P", "Wind_ms", "VapP_kPa", "SRad_Wm2", "LRad_Wm2", "LAI"):
            if (self.forcing[name] < 0).any():
                raise ValueError(f"{name} must be nonnegative")
        for name in ("fv", "Albd", "SFC"):
            if name in self.forcing and not self.forcing[name].between(0, 1).all():
                raise ValueError(f"{name} must lie in [0,1]")
        if (self.forcing.Pres_kPa <= 0).any() or (self.forcing.Temp_C <= -100).any():
            raise ValueError("Check forcing temperature [C] and pressure [kPa]")
        if "date" in self.forcing:
            time = pd.DatetimeIndex(pd.to_datetime(self.forcing.date))
            if time.hasnans or not time.is_monotonic_increasing or time.has_duplicates:
                raise ValueError("Dates must be unique, finite, and increasing")
            if len(time) > 1 and not np.allclose(np.diff(time.as_unit("ns").asi8)/1e9, self.dt, atol=1e-6, rtol=0):
                raise ValueError("Forcing dates must have a constant cadence equal to dt")
        p = self.params
        # Also validate dictionaries modified through the Python API.
        merge_parameters(default_parameters(), p)
        def finite_scalars(group, prefix=""):
            nullable = {"Veg.d", "Veg.z0m", "Veg.z0h", "Infil.K_mps", "Infil.S_mps05",
                        "Troch.alpha_inv_m"}
            for key, value in group.items():
                name = prefix+key
                if isinstance(value, dict):
                    finite_scalars(value, name+".")
                elif value is None:
                    if name not in nullable:
                        raise ValueError(f"{name} cannot be null")
                elif isinstance(value, str):
                    if name != "Infil.model":
                        raise ValueError(f"{name} must be numeric or boolean")
                elif not np.isscalar(value) or not np.isfinite(value):
                    raise ValueError(f"{name} must be a finite scalar")
        finite_scalars(p)
        for group, names in (("Flags", ("interception",)), ("Route", (
                "enabled", "route_baseflow", "route_deepflow", "overflow_fill_hierarchy"))):
            for name in names:
                if not isinstance(p[group][name], (bool, np.bool_)):
                    raise ValueError(f"{group}.{name} must be true or false")
        v = p["Veg"]
        if not (v["cw_per_LAI"] >= 0 and v["rc_min"] >= 0 and v["h_c"] > 0
                and 0 <= v["free_throughfall_f0"] <= 1 and p["Met"]["zref_m"] > 0
                and p["Soil"]["rs_soil"] >= 0):
            raise ValueError("Invalid vegetation, aerodynamic height, or surface resistance")
        if (v["d"] is not None and v["d"] < 0) or any(v[k] is not None and v[k] <= 0 for k in ("z0m", "z0h")):
            raise ValueError("Displacement must be nonnegative and roughness lengths positive")
        for zone in ("RZ", "TZ"):
            z = p["Soil"][zone]
            if not (0 <= z["theta_min"] < z["theta_sat"] <= 1
                    and z["theta_min"] <= z["theta0"] <= z["theta_sat"] and z["z"] > 0
                    and z["vG"]["n"] > 1 and z["vG"]["alpha"] > 0 and z["vG"]["Ks"] >= 0):
                raise ValueError(f"Invalid {zone} hydraulic parameters")
        if not (0 < p["Aq"]["f"] <= 1 and p["Aq"]["k"] >= 0 and p["Aq"]["D"] > 0):
            raise ValueError("Invalid aquifer k, f, or D")
        if not (0 <= p["Hill"]["iota_rad"] < np.pi/2):
            raise ValueError("Hill.iota_rad must be in [0,pi/2)")
        if not (0 <= p["Deep"]["phi_d"] <= 1 and p["Deep"]["k_out_s"] >= 0):
            raise ValueError("Invalid deep-reservoir parameters")
        num = p["Num"]
        if not (0.5 <= num["theta"] <= 1 and 0 < num["omega"] <= 1
                and 0 < num["dt_min"] <= num["dt_max"] and num["picard_tol"] > 0
                and num["picard_max"] >= 1 and 0 < num["cfl"] <= 0.5):
            raise ValueError("Invalid numerical options")
        inf = p["Infil"]
        alpha_inv = p.get("Troch", {}).get("alpha_inv_m")
        if alpha_inv is not None and alpha_inv <= 0:
            raise ValueError("Troch.alpha_inv_m must be positive [m] or null")
        if (inf["model"].lower() not in ("milly", "troch") or inf["Fmin_m"] <= 0
                or inf["F0_m"] < 0 or inf["resetDrySteps"] < 0
                or (inf["K_mps"] is not None and inf["K_mps"] < 0)
                or (inf["S_mps05"] is not None and inf["S_mps05"] < 0)
                or not 0 <= inf["max_c"] <= 1):
            raise ValueError("Invalid infiltration parameters (Milly or Troch required)")
        if not 0 <= p["Soil"]["TZ"]["cap_frac_of_K"] <= 1:
            raise ValueError("cap_frac_of_K must lie in [0,1]")
        if p["Route"]["enabled"] and not self.channel:
            raise ValueError("Routing requires a channel width distribution")
        return self


def load_config(path: str | Path) -> Config:
    from .io import read_forcing
    path = Path(path).resolve()
    spec = json.loads(path.read_text(encoding="utf-8"))
    allowed = {"geometry", "channel", "forcing", "dt", "parameters", "initial", "metadata"}
    if set(spec) - allowed:
        raise ValueError(f"Unknown configuration keys: {sorted(set(spec)-allowed)}")
    g = pd.read_csv(path.parent / spec["geometry"])
    edges = None
    if "dx_m" in g:
        edges = np.r_[0.0, g.dx_m.to_numpy(float).cumsum()]
    grid = Grid(g.x_m.to_numpy(), g.w_m.to_numpy(), edges)
    ch = {}
    if spec.get("channel"):
        chspec = spec["channel"]
        c = pd.read_csv(path.parent / chspec["file"])
        ch = {"distance": c.distance_m.to_numpy(float), "weight": c.weight.to_numpy(float),
              "slope": chspec["slope"]}
    fspec = spec["forcing"]
    forcing = read_forcing(path.parent / fspec["file"],
                           lai_basis=fspec.get("lai_basis", "vegetated"))
    if fspec.get("start"):
        forcing = forcing.loc[forcing.date >= pd.Timestamp(fspec["start"])]
    if fspec.get("end"):
        forcing = forcing.loc[forcing.date <= pd.Timestamp(fspec["end"])]
    if fspec.get("steps") is not None:
        forcing = forcing.iloc[:int(fspec["steps"])]
    cfg = Config(grid, forcing.reset_index(drop=True), float(spec["dt"]), ch,
                 merge_parameters(default_parameters(), spec.get("parameters", {})),
                 spec.get("initial", {}), spec.get("metadata", {}))
    return cfg.validate()
