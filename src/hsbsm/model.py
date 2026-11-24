"""Five-store coupled model with auditable water balances."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
from .config import Config
from .physics import penman_monteith, canopy_step, infiltration_step, soil_step, deep_step
from .aquifer import advance
from .routing import unit_hydrograph, route_flow
from .io import write_json


@dataclass
class Result:
    timeseries: pd.DataFrame
    nodes: dict
    diagnostics: dict
    config: Config

    def save(self, directory, plots=True):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        self.timeseries.to_csv(directory/"timeseries.csv", index=False)
        np.savez_compressed(directory/"states_and_fluxes.npz", **self.nodes,
                            x_m=self.config.grid.x, w_m=self.config.grid.w,
                            edges_m=self.config.grid.edges)
        write_json(directory/"diagnostics.json", self.diagnostics)
        write_json(directory/"parameters.json", self.config.params)
        from .evaluation import climatology_tables
        monthly, annual = climatology_tables(self)
        monthly.to_csv(directory/"monthly_climatology.csv")
        annual.to_csv(directory/"annual_water_totals.csv")
        write_json(directory/"run_metadata.json", {"dt_s": self.config.dt,
            "initial": self.config.initial, "source": self.config.metadata,
            "dates": [str(self.timeseries.date.iloc[0]), str(self.timeseries.date.iloc[-1])],
            "n_steps": len(self.timeseries), "hsbsm_version": "0.1.0"})
        if plots:
            from .plotting import plot_result
            plot_result(self, directory)
        return directory


def run(config: Config, progress=None) -> Result:
    """Run constant-cadence forcing. State arrays are [time, node]."""
    config.validate()
    g, p, dt = config.grid, config.params, config.dt
    forcing = config.forcing.reset_index(drop=True)
    area, nx, nt = g.area, len(g.x), len(forcing)
    total_area = float(area.sum())
    rz, tz = p["Soil"]["RZ"], p["Soil"]["TZ"]

    def initial(name, default, lower, upper):
        v = np.broadcast_to(np.asarray(config.initial.get(name, default), float), (nx,)).copy()
        if not np.isfinite(v).all() or np.any(v < lower) or np.any(v > upper):
            raise ValueError(f"Initial {name} lies outside its physical range")
        return v

    allowed = {"theta_rz", "theta_tz", "canopy_m", "aquifer_head_m", "deep_m", "infiltration_m"}
    if set(config.initial)-allowed:
        raise ValueError(f"Unknown initial states: {set(config.initial)-allowed}")
    tr = initial("theta_rz", rz["theta0"], rz["theta_min"], rz["theta_sat"])
    tt = initial("theta_tz", tz["theta0"], tz["theta_min"], tz["theta_sat"])
    canopy = initial("canopy_m", 0.0, 0, np.inf)
    if not p["Flags"]["interception"] and np.any(canopy):
        raise ValueError("A nonzero canopy store requires interception=True")
    aq = initial("aquifer_head_m", 0.0, 0, p["Aq"]["D"])*p["Aq"]["f"]*g.w
    deep = initial("deep_m", 0.001, 0, np.inf)
    cumulative = initial("infiltration_m", 0.0, 0, np.inf)
    dry = np.zeros(nx, int)

    def store_volume():
        return float(np.dot(canopy+tr*rz["z"]+tt*tz["z"]+deep, area)+np.dot(aq, g.dx))

    initial_volume = previous_volume = store_volume()
    names = ("theta_rz", "theta_tz", "canopy_m", "aquifer_head_m", "aquifer_storage_m2", "deep_m",
             "infiltration_mps", "surface_runoff_mps", "transpiration_mps", "soil_evaporation_mps",
             "canopy_evaporation_mps", "rz_to_tz_mps", "tz_to_aq_mps", "deep_recharge_mps")
    nodes = {name: np.empty((nt, nx)) for name in names}
    series_names = ("Q_surface_input_m3s", "Q_base_raw_m3s", "Q_deep_raw_m3s", "Q_overflow_m3s",
                    "ET_mps", "storage_m3", "land_residual_m3", "aquifer_substeps")
    series = {name: np.zeros(nt) for name in series_names}
    pm = penman_monteith(forcing, p)
    met = {name: forcing[name].to_numpy(float) for name in forcing.columns if name != "date"}
    total_iterations = 0
    for i in range(nt):
        rain, fv, lai = met["P"][i], met["fv"][i], met["LAI"][i]
        wet, dry_pet, bare_pet = (np.asarray(pm[name])[i] for name in ("wet", "dry", "bare"))
        if p["Flags"]["interception"]:
            pnet, ecan, canopy = canopy_step(rain, wet, fv, lai, canopy, p["Veg"], dt)
        else:
            pnet, ecan = np.full(nx, rain), np.zeros(nx)
        inf, rex, cumulative, dry, _ = infiltration_step(pnet, cumulative, dry, tr, p, dt)
        h = aq/(p["Aq"]["f"]*g.w)
        soil = soil_step(inf, fv*dry_pet, (1-fv)*bare_pet*(1-met.get("SFC", np.zeros(nt))[i]),
                         tr, tt, h, aq/g.w, p["Soil"], dt)
        tr, tt = soil["theta_rz"], soil["theta_tz"]
        exchange = soil["tz_to_aq"]
        to_deep = np.maximum(exchange, 0)*p["Deep"]["phi_d"]
        to_aq = exchange-to_deep
        aq, qb, overflow, aqdiag = advance(g, aq, to_aq, p["Aq"], p["Hill"]["iota_rad"], dt, p["Num"])
        total_iterations += aqdiag["iterations"]
        original_overflow = overflow.sum()
        if p["Route"]["overflow_fill_hierarchy"]:
            excess_depth = overflow*dt/area
            fill_t = np.minimum(excess_depth, np.maximum(0, (tz["theta_sat"]-tt)*tz["z"]))
            tt += fill_t/tz["z"]
            excess_depth -= fill_t
            fill_r = np.minimum(excess_depth, np.maximum(0, (rz["theta_sat"]-tr)*rz["z"]))
            tr += fill_r/rz["z"]
            overflow = (excess_depth-fill_r)*area/dt
        deep, qdeep = deep_step(deep, to_deep, p["Deep"]["k_out_s"], dt)
        runoff = rex+soil["runoff"]+overflow/area
        qs = float(np.dot(runoff, area))
        qd = float(np.dot(qdeep, area))
        evap = ecan+soil["transpiration"]+soil["soil_evaporation"]
        et = float(np.dot(evap, area)/total_area)
        volume = store_volume()
        residual = (rain-et)*total_area*dt-(qs+qb+qd)*dt-(volume-previous_volume)
        previous_volume = volume
        values = (tr, tt, canopy, aq/(p["Aq"]["f"]*g.w), aq, deep, inf, runoff,
                  soil["transpiration"], soil["soil_evaporation"], ecan,
                  soil["rz_to_tz"], exchange, to_deep)
        for name, value in zip(names, values):
            nodes[name][i] = value
        vals = (qs, qb, qd, original_overflow, et, volume, residual, aqdiag["substeps"])
        for name, value in zip(series_names, vals):
            series[name][i] = value
        if progress and (i == nt-1 or i % max(1, nt//20) == 0):
            progress(i+1, nt)
    route = p["Route"]
    kernel = unit_hydrograph(config.channel, route, dt) if route["enabled"] else np.ones(1)
    channel_storage = np.zeros(nt)
    tail_volume = 0.0
    for component, raw_name, enabled in (("surface", "Q_surface_input_m3s", True),
            ("base", "Q_base_raw_m3s", route["route_baseflow"]),
            ("deep", "Q_deep_raw_m3s", route["route_deepflow"])):
        if enabled:
            out, retained, tail = route_flow(series[raw_name], kernel, dt)
            channel_storage += retained
            tail_volume += float(np.sum(tail)*dt)
        else:
            out = series[raw_name].copy()
        series[f"Q_{component}_m3s"] = out
    series["Q_total_m3s"] = sum(series[f"Q_{c}_m3s"] for c in ("surface", "base", "deep"))
    series["channel_storage_m3"] = channel_storage
    series["rain_mm_day"] = met["P"]*86400000
    series["ET_mm_day"] = series["ET_mps"]*86400000
    series["Q_total_mm_day"] = series["Q_total_m3s"]*86400000/total_area
    if "date" in forcing:
        dates = pd.to_datetime(forcing.date)
    else:
        dates = pd.Timestamp("2000-01-01")+pd.to_timedelta(np.arange(nt)*dt, unit="s")
    table = pd.DataFrame({"date": dates, **series})
    rain_volume = float(met["P"].sum()*dt*total_area)
    et_volume = float(series["ET_mps"].sum()*dt*total_area)
    q_volume = float(series["Q_total_m3s"].sum()*dt)
    residual = rain_volume-et_volume-q_volume-(volume-initial_volume)-channel_storage[-1]
    diagnostics = {"area_m2": total_area, "initial_storage_m3": initial_volume,
        "final_land_storage_m3": volume, "final_channel_storage_m3": float(channel_storage[-1]),
        "rain_volume_m3": rain_volume, "et_volume_m3": et_volume, "outlet_volume_m3": q_volume,
        "water_balance_residual_m3": float(residual),
        "relative_water_balance_error": float(abs(residual)/max(rain_volume+initial_volume, 1e-30)),
        "max_step_land_residual_m3": float(np.max(np.abs(series["land_residual_m3"]))),
        "routing_tail_volume_m3": tail_volume, "routing_kernel_sum": float(kernel.sum()),
        "aquifer_substeps": int(series["aquifer_substeps"].sum()), "picard_iterations": total_iterations}
    if not np.isfinite(table.select_dtypes(include=[np.number]).to_numpy()).all():
        raise ArithmeticError("Simulation produced non-finite outputs")
    return Result(table, nodes, diagnostics, config)
