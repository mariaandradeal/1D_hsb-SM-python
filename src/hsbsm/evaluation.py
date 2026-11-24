"""Observed-flow alignment, KGE (2009), NSE, and hydrologic summaries."""
import numpy as np
import pandas as pd


def metrics(observed, simulated):
    obs, sim = np.asarray(observed, float), np.asarray(simulated, float)
    if obs.shape != sim.shape:
        raise ValueError("Observed and simulated arrays must have equal shapes")
    keep = np.isfinite(obs) & np.isfinite(sim)
    obs, sim = obs[keep], sim[keep]
    result = {"n": int(len(obs)), "NSE": np.nan, "KGE": np.nan, "RMSE_m3s": np.nan,
              "PBIAS_percent": np.nan, "correlation": np.nan}
    if not len(obs):
        return result
    result["RMSE_m3s"] = float(np.sqrt(np.mean((sim-obs)**2)))
    if abs(obs.sum()) > 1e-30:
        result["PBIAS_percent"] = float(100*(sim-obs).sum()/obs.sum())
    if len(obs) < 2 or np.std(obs) <= 1e-15:
        return result
    result["NSE"] = float(1-np.sum((sim-obs)**2)/np.sum((obs-obs.mean())**2))
    if np.std(sim) > 1e-15 and abs(obs.mean()) > 1e-30:
        r = float(np.corrcoef(obs, sim)[0, 1])
        result["correlation"] = r
        result["KGE"] = float(1-np.sqrt((r-1)**2+(np.std(sim)/np.std(obs)-1)**2
                                       +(sim.mean()/obs.mean()-1)**2))
    return result


def align_observations(result, observations):
    """Exact timestamp intersection; missing observations remain missing."""
    t = pd.DatetimeIndex(pd.to_datetime(result.timeseries.date))
    return observations.reindex(t).to_numpy(float)


def evaluate(result, observations, warmup_steps=0):
    obs = align_observations(result, observations)
    sim = result.timeseries.Q_total_m3s.to_numpy()
    if warmup_steps < 0 or warmup_steps >= len(sim):
        raise ValueError("warmup_steps must leave at least one modeled step")
    table = result.timeseries[["date", "Q_total_m3s"]].copy()
    table["Q_obs_m3s"] = obs
    table["used"] = np.arange(len(table)) >= warmup_steps
    usable = table.loc[table.used]
    summary = metrics(usable.Q_obs_m3s, usable.Q_total_m3s)
    if summary["n"] == 0:
        raise ValueError("No finite observations overlap the modeled dates after warmup")
    annual = []
    for year, part in usable.groupby(pd.to_datetime(usable.date).dt.year):
        annual.append({"year": int(year), **metrics(part.Q_obs_m3s, part.Q_total_m3s)})
    return summary, pd.DataFrame(annual), table


def signatures(result, drought_threshold_mm_day=1.5):
    q = result.timeseries.Q_total_mm_day.to_numpy()
    longest = current = 0
    for dry in q < drought_threshold_mm_day:
        current = current+1 if dry else 0
        longest = max(longest, current)
    d = result.diagnostics
    return {"mean_Q_mm_day": float(q.mean()), "mean_ET_mm_day": float(result.timeseries.ET_mm_day.mean()),
            "Q_exceedance_5_50_95_mm_day": np.quantile(q, [0.95, 0.5, 0.05]).tolist(),
            "runoff_ratio": d["outlet_volume_m3"]/d["rain_volume_m3"] if d["rain_volume_m3"] else None,
            "longest_drought_days": longest*result.config.dt/86400,
            "drought_threshold_mm_day": drought_threshold_mm_day}


def climatology_tables(result):
    """Monthly mean rates and calendar-year volumes for the simulated window."""
    table = result.timeseries.copy()
    dates = pd.to_datetime(table.date)
    monthly = table.groupby(dates.dt.month)[["rain_mm_day", "ET_mm_day", "Q_total_mm_day"]].mean()
    monthly.index.name = "month"
    monthly["steps"] = table.groupby(dates.dt.month).size()
    area, dt = result.diagnostics["area_m2"], result.config.dt
    annual = table.groupby(dates.dt.year)[["rain_mm_day", "ET_mm_day", "Q_total_mm_day"]].sum()
    annual *= dt/86400000*area
    annual.columns = ["rain_m3", "ET_m3", "Q_m3"]
    annual.index.name = "year"
    annual["simulated_days"] = table.groupby(dates.dt.year).size()*dt/86400
    annual["runoff_ratio"] = annual.Q_m3/annual.rain_m3.replace(0, np.nan)
    return monthly, annual
