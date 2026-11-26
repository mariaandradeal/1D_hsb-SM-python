"""Surface and vertical processes translated from the active MATLAB modules.

Internal fluxes are m/s over ground area; storage depths are m. See
docs/MIGRATION.md for explicit corrections to the MATLAB water accounting.
"""
from __future__ import annotations
import numpy as np


def penman_monteith(forcing, params):
    """Three patch-scale potential rates, matching hsbsm_pm.m."""
    f, v = forcing, params["Veg"]
    tc = np.asarray(f["Temp_C"])
    tk = tc + 273.15
    pressure = np.asarray(f["Pres_kPa"])
    wind = np.maximum(f["Wind_ms"], 0.1)
    lai = np.maximum(f["LAI"], 0.05)
    rn = (1-np.clip(f["Albd"], 0, 0.95))*np.maximum(f["SRad_Wm2"], 0)
    rn += 0.98*(np.asarray(f["LRad_Wm2"])-5.670374419e-8*tk**4)
    es = 0.6108*np.exp(17.27*tc/(tc+237.3))
    delta = 4098*es/(tc+237.3)**2
    gamma = 1005*pressure/(0.622*2.45e6)
    vpd = np.maximum(es-np.maximum(f["VapP_kPa"], 0), 0)
    d = 0.65*v["h_c"] if v["d"] is None else v["d"]
    z0m = max(0.001, 0.1*v["h_c"] if v["z0m"] is None else v["z0m"])
    z0h = max(0.0001, 0.1*z0m if v["z0h"] is None else v["z0h"])
    zr = params["Met"]["zref_m"]
    if zr <= d:
        zr = d+z0m+1
    ra = np.maximum(np.log((zr-d)/z0m)*np.log((zr-d)/z0h)/(0.41**2*wind), 5)
    rc = v["rc_min"]/np.maximum(lai, 0.1)
    rho = pressure*1000/(287*tk)
    numerator = delta*1000*rn + rho*1005*vpd*1000/ra

    def rate(rs):
        return np.maximum(numerator/(2.45e9*(delta*1000+gamma*1000*(1+rs/ra))), 0)

    return {"wet": rate(0), "dry": rate(rc), "bare": rate(params["Soil"]["rs_soil"]),
            "rn": rn, "ra": ra}


def canopy_step(rain, wet_pet, fv, lai, storage, veg, dt):
    intercepted = fv*(1-veg["free_throughfall_f0"])*rain*dt
    evaporation = np.minimum(fv*wet_pet*dt, storage+intercepted)
    provisional = storage+intercepted-evaporation
    capacity = fv*veg["cw_per_LAI"]*lai
    drip = np.maximum(provisional-capacity, 0)
    return rain-intercepted/dt+drip/dt, evaporation/dt, provisional-drip


def infiltration_step(pnet, cumulative, dry_count, theta, params, dt):
    soil, inf = params["Soil"]["RZ"], params["Infil"]
    k = soil["vG"]["Ks"] if inf["K_mps"] is None else inf["K_mps"]
    se = np.clip((theta-soil["theta_min"])/(soil["theta_sat"]-soil["theta_min"]), 1e-6, 1)
    alpha = soil["vG"]["alpha"]
    alpha_inv = params.get("Troch", {}).get("alpha_inv_m")
    if inf["model"].lower() == "troch" and alpha_inv is not None:
        alpha = 1/alpha_inv
    s2 = 2*np.maximum(soil["theta_sat"]-theta, 0)*k/alpha
    if inf["S_mps05"] is not None:
        s2 = np.full_like(theta, inf["S_mps05"]**2)
    keff = np.minimum(inf["max_c"], se)*k if inf["model"].lower() == "troch" else k
    capacity = keff+s2/(2*np.maximum(cumulative, inf["Fmin_m"]))
    infiltration = np.minimum(pnet, capacity)
    cumulative = cumulative+infiltration*dt
    dry_count = np.where(pnet > 0, 0, dry_count+1)
    reset = dry_count > inf["resetDrySteps"]
    cumulative = np.where(reset, max(inf["F0_m"], inf["Fmin_m"]), cumulative)
    dry_count = np.where(reset, 0, dry_count)
    return infiltration, pnet-infiltration, cumulative, dry_count, capacity


def hydraulic_properties(theta, zone, head_floor):
    se = (theta-zone["theta_min"])/(zone["theta_sat"]-zone["theta_min"])
    n, alpha, ks = (zone["vG"][k] for k in ("n", "alpha", "Ks"))
    m = 1-1/n
    sh = np.clip(se, 1e-9, 1-1e-9)
    psi = np.maximum(-((sh**(-1/m)-1)**(1/n))/alpha, head_floor)
    sk = np.clip(se, 1e-6, 1)
    conductivity = ks*np.sqrt(sk)*(1-(1-sk**(1/m))**m)**2
    return psi, conductivity


def soil_step(infiltration, pet_trans, pet_soil, theta_rz, theta_tz,
              aq_head, aq_depth, soil, dt):
    """Signed Darcy exchanges with shared donor water and explicit excess runoff.

    The pressure/conductivity closures and stress functions match
    hsbsm_rz_tz_coupled.m. Flux limits operate on actual remaining water,
    so saturation bounds never silently create or delete water.
    """
    rz, tz = soil["RZ"], soil["TZ"]
    psi_r, kr = hydraulic_properties(theta_rz, rz, -10)
    psi_t, kt = hydraulic_properties(theta_tz, tz, -20)
    avail_r = np.maximum(0, (theta_rz-rz["theta_min"])*rz["z"])
    avail_t = np.maximum(0, (theta_tz-tz["theta_min"])*tz["z"])
    transpiration = np.minimum(pet_trans*np.clip((psi_r+10)/9, 0, 1), avail_r/dt)
    evaporation = np.minimum(pet_soil*np.clip((psi_r+10)/9.5, 0, 1),
                             np.maximum(0, avail_r/dt-transpiration))
    wr = theta_rz*rz["z"]+(infiltration-transpiration-evaporation)*dt
    wt = theta_tz*tz["z"]
    rmin, rmax = rz["theta_min"]*rz["z"], rz["theta_sat"]*rz["z"]
    tmin, tmax = tz["theta_min"]*tz["z"], tz["theta_sat"]*tz["z"]
    distance = (rz["z"]+tz["z"])/2
    qrt = (kr+kt)/2*((psi_r-psi_t)/distance+1)
    qrt = np.clip(qrt, -avail_t/dt, np.maximum(wr-rmin, 0)/dt)
    wr -= qrt*dt
    wt += qrt*dt
    dw = np.minimum(np.maximum(aq_head, 0), tz["z"])
    qta = kt*(psi_t+tz["z"]/2-dw)/np.maximum(tz["z"]-dw, 1e-3)
    qta = np.where(qta < 0, qta*tz["cap_frac_of_K"], qta)
    qta = np.minimum(qta, np.maximum(wt-tmin, 0)/dt)
    qta = np.maximum(qta, -np.minimum(aq_depth, np.maximum(tmax-wt, 0))/dt)
    wt -= qta*dt
    excess_r, excess_t = np.maximum(wr-rmax, 0), np.maximum(wt-tmax, 0)
    wr -= excess_r
    wt -= excess_t
    # Tiny roundoff is harmless, but substantive deficits indicate a coding error.
    if np.any(wr < rmin-1e-12) or np.any(wt < tmin-1e-12):
        raise ArithmeticError("Vertical exchange depleted a store below residual water")
    return {"theta_rz": wr/rz["z"], "theta_tz": wt/tz["z"],
            "transpiration": transpiration, "soil_evaporation": evaporation,
            "rz_to_tz": qrt, "tz_to_aq": qta, "runoff": (excess_r+excess_t)/dt}


def deep_step(storage, recharge, coefficient, dt):
    """MATLAB's input-first, availability-limited linear reservoir."""
    available = storage+recharge*dt
    release = np.minimum(coefficient*available, available/dt)
    return available-release*dt, release
