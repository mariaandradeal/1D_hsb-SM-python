"""Reproducible bounded calibration and one-at-a-time sensitivity."""
from copy import deepcopy
import numpy as np
from scipy.optimize import differential_evolution
from .model import run
from .evaluation import metrics

# Deliberately ordered bounds; the source RZ_Ks bound was reversed.
REGISTRY = {
    "RZ_Ks": ("Soil.RZ.vG.Ks", 1e-9, 1e-5, True),
    "TZ_Ks": ("Soil.TZ.vG.Ks", 1e-9, 1e-5, True),
    "RZ_z": ("Soil.RZ.z", 0.1, 2.0, False),
    "TZ_z": ("Soil.TZ.z", 0.1, 3.0, False),
    "Aq_k": ("Aq.k", 1e-7, 1e-2, True),
    "Aq_f": ("Aq.f", 0.01, 0.4, False),
    "Route_V": ("Route.V", 0.05, 3.0, False),
    "deep_fraction": ("Deep.phi_d", 0.0, 1.0, False),
    "deep_k": ("Deep.k_out_s", 1/(1000*86400), 1/86400, True),
    "rc_min": ("Veg.rc_min", 20.0, 1000.0, False),
    "rs_soil": ("Soil.rs_soil", 20.0, 1000.0, False),
}


def get_parameter(params, path):
    for key in path.split("."):
        params = params[key]
    return params


def set_parameter(params, path, value):
    keys = path.split(".")
    for key in keys[:-1]:
        params = params[key]
    if keys[-1] not in params:
        raise ValueError(f"Unknown parameter {path}")
    params[keys[-1]] = float(value)


def calibrate(config, observations, active=("Aq_k", "Aq_f"), split=0.8, warmup_steps=0,
              maxiter=10, popsize=5, seed=42, bounds=None, progress=None):
    """Minimize 1-KGE on a chronological calibration window only.

    Validation observations are never passed to the objective. Every trial
    starts from the same specified initial conditions. Conductivities are
    optimized in log10 space. SciPy DE is not MATLAB's GA algorithm.
    """
    import pandas as pd
    if not 0 < split < 1 or len(active) == 0 or len(set(active)) != len(active):
        raise ValueError("Provide unique active parameters and 0 < split < 1")
    if maxiter < 0 or popsize < 1:
        raise ValueError("maxiter >= 0 and popsize >= 1 are required")
    n = len(config.forcing)
    cut = int(np.floor(split*n))
    if not 0 <= warmup_steps < cut < n:
        raise ValueError("Warmup and split must leave calibration and validation windows")
    obs = observations.reindex(pd.DatetimeIndex(pd.to_datetime(config.forcing.date))).to_numpy(float)
    if np.isfinite(obs[warmup_steps:cut]).sum() < 5 or np.isfinite(obs[cut:]).sum() < 5:
        raise ValueError("At least five overlapping observations are required in each window")
    if np.nanstd(obs[warmup_steps:cut]) <= 1e-15:
        raise ValueError("KGE calibration requires nonconstant observed discharge")
    specs, search_bounds = [], []
    for name in active:
        if name not in REGISTRY:
            raise ValueError(f"Unknown calibration parameter {name}; choose from {list(REGISTRY)}")
        path, lo, hi, log = REGISTRY[name]
        if bounds and name in bounds:
            lo, hi = bounds[name]
        if not (np.isfinite([lo, hi]).all() and lo <= hi and (not log or lo > 0)):
            raise ValueError(f"Invalid bounds for {name}")
        specs.append((name, path, log))
        search_bounds.append((np.log10(lo), np.log10(hi)) if log else (lo, hi))
    calibration_config = deepcopy(config)
    calibration_config.forcing = config.forcing.iloc[:cut].copy()
    history = []

    def trial_config(base, x):
        trial = deepcopy(base)
        for value, (_, path, log) in zip(x, specs):
            set_parameter(trial.params, path, 10**value if log else value)
        return trial

    def objective(x):
        try:
            out = run(trial_config(calibration_config, x))
            score = metrics(obs[warmup_steps:cut], out.timeseries.Q_total_m3s.to_numpy()[warmup_steps:])["KGE"]
            value = 1-score if np.isfinite(score) else 1e6
            error = None
        except (ValueError, ArithmeticError) as exc:
            value, error = 1e6, str(exc)
        history.append({"evaluation": len(history)+1, "objective": float(value), "error": error})
        if progress:
            progress(len(history), value)
        return value

    x0 = [get_parameter(config.params, path) for _, path, _ in specs]
    x0 = [np.log10(x) if log else x for x, (_, _, log) in zip(x0, specs)]
    x0 = [np.clip(x, lo, hi) for x, (lo, hi) in zip(x0, search_bounds)]
    baseline_cost = objective(x0)
    optimum = differential_evolution(objective, search_bounds, seed=seed, x0=x0,
                                    maxiter=maxiter, popsize=popsize, polish=False, workers=1)
    if optimum.fun >= 1e6:
        raise ArithmeticError("Every calibration trial failed or had undefined KGE")
    best = trial_config(config, optimum.x)
    result = run(best)
    sim = result.timeseries.Q_total_m3s.to_numpy()
    report = {"algorithm": "scipy differential_evolution", "seed": seed,
        "active": list(active), "bounds_in_optimizer_space": search_bounds,
        "best_parameters": {name: get_parameter(best.params, path) for name, path, _ in specs},
        "objective_initial": float(baseline_cost), "objective_best": float(optimum.fun),
        "converged": bool(optimum.success), "optimizer_message": str(optimum.message),
        "evaluations": len(history), "failed_evaluations": sum(h["error"] is not None for h in history),
        "warmup_steps": warmup_steps, "split_index": cut,
        "calibration": metrics(obs[warmup_steps:cut], sim[warmup_steps:cut]),
        "validation": metrics(obs[cut:], sim[cut:]), "history": history}
    return result, report


def sensitivity(config, active, fraction=0.1):
    import pandas as pd
    if not 0 < fraction < 1:
        raise ValueError("Sensitivity fraction must be between zero and one")
    baseline = run(config)
    rows = [{"parameter": "baseline", "factor": 1.0,
             "outlet_volume_m3": baseline.diagnostics["outlet_volume_m3"]}]
    for name in active:
        if name not in REGISTRY:
            raise ValueError(f"Unknown parameter {name}")
        path = REGISTRY[name][0]
        for factor in (1-fraction, 1+fraction):
            trial = deepcopy(config)
            set_parameter(trial.params, path, get_parameter(config.params, path)*factor)
            result = run(trial)
            rows.append({"parameter": name, "factor": factor,
                         "outlet_volume_m3": result.diagnostics["outlet_volume_m3"]})
    return pd.DataFrame(rows)
