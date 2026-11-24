from copy import deepcopy
from pathlib import Path
import numpy as np
import pytest
from numpy.testing import assert_allclose
from hsbsm import load_config, run, default_parameters, Grid
from hsbsm.physics import canopy_step, soil_step, deep_step
from hsbsm.aquifer import advance
from hsbsm.routing import unit_hydrograph, route_flow

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("fv,lai,store,rain", [(0, 0, 0, 1e-5), (1, 3, 0.02, 0), (0.6, 2, 0, 1e-4)])
def test_canopy_conserves_water(fv, lai, store, rain):
    dt = 3600
    pn, e, end = canopy_step(rain, 1e-7, fv, lai, store, default_parameters()["Veg"], dt)
    assert_allclose(store+rain*dt, end+(pn+e)*dt, atol=1e-14)
    assert end >= 0


@pytest.mark.parametrize("wet", [False, True])
def test_soil_limits_do_not_create_or_delete_water(wet):
    soil = default_parameters()["Soil"]
    rng = np.random.default_rng(42)
    r, t = soil["RZ"], soil["TZ"]
    tr, tt = rng.uniform(.045, .38, (2, 100))
    if wet:
        tr[:] = tt[:] = .38
    aq = rng.uniform(0, .02, 100)
    inf = rng.uniform(0, 1e-4, 100)
    dt = 86400
    s = soil_step(inf, 3e-6, 2e-6, tr, tt, aq/.09, aq, soil, dt)
    change = (s["theta_rz"]-tr)*r["z"]+(s["theta_tz"]-tt)*t["z"]
    flux = (inf-s["transpiration"]-s["soil_evaporation"]-s["tz_to_aq"]-s["runoff"])*dt
    assert_allclose(change, flux, atol=2e-15)
    assert np.all(s["theta_rz"] >= .045-1e-14)
    assert np.all(s["theta_tz"] <= .38+1e-14)
    assert np.all(-s["tz_to_aq"]*dt <= aq+1e-14)


def test_deep_reservoir_large_timestep():
    end, q = deep_step(np.array([0., .1]), np.array([1e-6, 0]), 1, 86400)
    assert_allclose(end+q*86400, [.0864, .1])
    assert np.min(end) >= 0


@pytest.mark.parametrize("recharge,initial", [(0., 0.), (1e-5, .34), (-1e-7, .2)])
def test_aquifer_mass_including_overflow_and_extraction(recharge, initial):
    p = default_parameters()
    g = Grid([5., 15., 25.], [8., 7., 5.])
    s0 = initial*p["Aq"]["f"]*g.w
    s, qb, of, d = advance(g, s0, recharge, p["Aq"], .1, 3600., p["Num"])
    residual = (np.dot(s0-s, g.dx)+recharge*g.area.sum()*3600)-(qb+of.sum())*3600
    assert abs(residual) < 1e-5
    assert np.min(s) >= 0
    assert np.all(s <= p["Aq"]["f"]*g.w*p["Aq"]["D"]+1e-9)


def test_routing_conserves_full_impulse_and_reports_tail():
    p = default_parameters()["Route"]
    kernel = unit_hydrograph(dict(distance=[500, 1500], weight=[.3, .7], slope=.01), p, 60)
    out, storage, tail = route_flow(np.array([1., 0., 0.]), kernel, 60)
    assert_allclose((out.sum()+tail.sum())*60, 60, atol=1e-12)
    assert_allclose(storage[-1], tail.sum()*60, atol=1e-12)
    assert np.min(kernel) >= 0


def test_routing_bins_match_independent_density_quadrature():
    from scipy.integrate import quad
    route = default_parameters()["Route"]
    ch = dict(distance=[100.], weight=[1.], slope=.01)
    dt = 60
    kernel = unit_hydrograph(ch, route, dt)
    c = (1+route["a0"])*route["V"]
    d2 = route["V"]**3/(9.80665*ch["slope"]*route["F"]**2)*(1-route["a0"]**2*route["F"]**2)
    def density(t):
        return 100/np.sqrt(2*np.pi*d2*t**3)*np.exp(-(100-c*t)**2/(2*d2*t)) if t else 0
    expected = [quad(density, j*dt, (j+1)*dt, epsabs=1e-12)[0] for j in range(len(kernel))]
    assert_allclose(kernel, expected, rtol=1e-8, atol=1e-10)


@pytest.mark.parametrize("base,deep,fill", [(False, False, False), (True, True, False), (True, False, True)])
def test_coupled_storm_balance_and_component_sum(base, deep, fill):
    cfg = load_config(ROOT/"examples/synthetic/config.json")
    cfg.params["Route"].update(route_baseflow=base, route_deepflow=deep, overflow_fill_hierarchy=fill)
    result = run(cfg)
    t = result.timeseries
    assert_allclose(t.Q_total_m3s, t.Q_surface_m3s+t.Q_base_m3s+t.Q_deep_m3s, atol=1e-14)
    assert result.diagnostics["relative_water_balance_error"] < 1e-7
    assert result.nodes["theta_rz"].shape == (120, 12)
    assert np.isfinite(result.nodes["aquifer_head_m"]).all()


def test_dry_forcing_cannot_create_water():
    cfg = load_config(ROOT/"examples/synthetic/config.json")
    cfg.forcing = cfg.forcing.iloc[:24].copy()
    cfg.forcing["P"] = 0
    cfg.initial = dict(theta_rz=.045, theta_tz=.045, deep_m=0.)
    r = run(cfg)
    assert_allclose(r.timeseries.Q_total_m3s, 0, atol=1e-14)
    assert r.diagnostics["relative_water_balance_error"] < 1e-10
