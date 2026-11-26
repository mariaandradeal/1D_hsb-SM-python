from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
import pytest
from hsbsm import load_config, run, Grid
from hsbsm.evaluation import metrics, evaluate
from hsbsm.calibration import calibrate
from hsbsm.io import read_observations

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("value", [0, -1, float("nan"), "3"])
def test_invalid_troch_capillary_length_rejected(value):
    cfg = load_config(ROOT/"examples/synthetic/config.json")
    cfg.params["Troch"]["alpha_inv_m"] = value
    with pytest.raises(ValueError):
        cfg.validate()


@pytest.mark.parametrize("kind", ["nan", "date", "negative", "parameter"])
def test_bad_inputs_fail_early(kind):
    c = load_config(ROOT/"examples/synthetic/config.json")
    if kind == "nan": c.forcing.loc[0, "Temp_C"] = np.nan
    if kind == "date": c.forcing.loc[1, "date"] = c.forcing.loc[0, "date"]
    if kind == "negative": c.forcing.loc[0, "P"] = -1
    if kind == "parameter": c.params["Soil"]["RZ"]["vG"]["n"] = 1
    with pytest.raises(ValueError):
        run(c)


def test_metrics_and_missing_observations():
    x = np.arange(1, 11, dtype=float)
    m = metrics(x, x)
    assert m["NSE"] == pytest.approx(1)
    assert m["KGE"] == pytest.approx(1)
    assert np.isnan(metrics(np.ones(5), np.ones(5))["KGE"])
    assert metrics([1, np.nan, 3], [1, 2, 3])["n"] == 2


def test_cli_runs_from_another_working_directory(tmp_path):
    cmd = [sys.executable, "-m", "hsbsm", "run", str(ROOT/"examples/synthetic/config.json"),
           "--steps", "12", "--output", str(tmp_path/"out"), "--no-plots"]
    completed = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert (tmp_path/"out/timeseries.csv").exists()
    with np.load(tmp_path/"out/states_and_fluxes.npz", allow_pickle=False) as data:
        assert data["theta_rz"].shape == (12, 12)


def test_calibration_is_seeded_and_validation_does_not_affect_fit():
    cfg = load_config(ROOT/"examples/synthetic/config.json")
    cfg.forcing = cfg.forcing.iloc[:24].copy()
    truth = run(cfg)
    obs = pd.Series(truth.timeseries.Q_total_m3s.to_numpy(), index=pd.DatetimeIndex(cfg.forcing.date))
    _, report = calibrate(cfg, obs, active=["Aq_k"], split=.7, maxiter=0, popsize=1, seed=11)
    changed = obs.copy()
    changed.iloc[int(.7*len(obs)):] *= 100
    _, second = calibrate(cfg, changed, active=["Aq_k"], split=.7, maxiter=0, popsize=1, seed=11)
    assert report["best_parameters"] == second["best_parameters"]
    assert report["objective_best"] == second["objective_best"]
    assert report["calibration"]["KGE"] > .999


def test_no_observation_overlap_raises():
    cfg = load_config(ROOT/"examples/synthetic/config.json")
    cfg.forcing = cfg.forcing.iloc[:2]
    result = run(cfg)
    obs = pd.Series([1., 2.], index=pd.date_range("2020-01-01", periods=2))
    with pytest.raises(ValueError, match="overlap"):
        evaluate(result, obs)


def test_single_cell_grid_has_physical_volume():
    grid = Grid([5], [2])
    assert grid.area[0] == 20


def test_field_inputs_and_observations_are_readable():
    cfg = load_config(ROOT/"examples/wet_beaver_creek/config.json")
    obs = read_observations(ROOT/"examples/wet_beaver_creek/observations.csv")
    assert len(cfg.forcing) == 6847
    assert len(cfg.grid.x) == 44
    assert cfg.grid.area.sum() == pytest.approx(283941233.21503556)
    assert obs.notna().sum() == 4749
    assert obs.index.max() == pd.Timestamp("2018-01-01")


def test_invalid_boolean_and_nonfinite_parameter_rejected():
    cfg = load_config(ROOT/"examples/synthetic/config.json")
    cfg.params["Flags"]["interception"] = "false"
    with pytest.raises(ValueError):
        cfg.validate()
    cfg.params["Flags"]["interception"] = True
    cfg.params["Veg"]["cw_per_LAI"] = np.nan
    with pytest.raises(ValueError):
        cfg.validate()
