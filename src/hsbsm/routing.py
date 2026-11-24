"""Bin-integrated geomorphic advection-diffusion unit hydrograph."""
import numpy as np
from scipy.signal import fftconvolve
from scipy.stats import invgauss


def unit_hydrograph(channel, route, dt):
    """Return dimensionless bin masses; sum(weights) = 1.

    Integrating the MATLAB ADE kernel over bins avoids undersampling travel
    times shorter than a daily forcing interval. Bin 0 is delivery within dt.
    """
    x = np.asarray(channel["distance"], float)
    w = np.asarray(channel["weight"], float)
    if (x.ndim != 1 or w.shape != x.shape or len(x) == 0 or not np.isfinite(x).all()
            or not np.isfinite(w).all() or np.any(x <= 0) or np.any(np.diff(x) <= 0)
            or np.any(w < 0) or w.sum() <= 0):
        raise ValueError("Channel distances must increase and weights must be nonnegative")
    # 'weight' is a density for uniform bins, or an already integrated bin mass.
    # The example uses uniform bins. For nonuniform bins supply integrated masses.
    w = w/w.sum()
    slope, velocity, froude, a0 = channel["slope"], route["V"], route["F"], route["a0"]
    if not (slope > 0 and velocity > 0 and froude > 0 and a0 >= 0 and (a0*froude)**2 < 1):
        raise ValueError("Routing requires S0,V,F>0 and 0 <= a0*F < 1")
    c = (1+a0)*velocity
    d2 = velocity**3/(9.80665*slope*froude**2)*(1-a0*a0*froude*froude)
    mu, scale = d2/(c*x), x*x/d2
    horizon = float(np.max(invgauss.ppf(1-1e-10, mu, scale=scale)))
    count = max(1, int(np.ceil(horizon/dt)))
    if count > 200000:
        raise ValueError("Routing kernel is too long; check channel units and parameters")
    # Chunking avoids a large time x distance temporary for long channels.
    cdf = np.zeros(count+1)
    for start in range(0, count+1, 1024):
        times = np.arange(start, min(start+1024, count+1))*dt
        cdf[start:start+len(times)] = invgauss.cdf(times[:, None], mu, scale=scale) @ w
    kernel = np.maximum(np.diff(cdf), 0)
    if kernel.sum() <= 0:
        raise ArithmeticError("Routing kernel has no mass")
    return kernel/kernel.sum()


def route_flow(inflow, weights, dt):
    full = np.maximum(fftconvolve(inflow, weights), 0)
    output = full[:len(inflow)]
    retained = np.cumsum((inflow-output)*dt)
    if np.min(retained) < -1e-8*max(1, np.sum(inflow)*dt):
        raise ArithmeticError("Routing delivered water before it entered")
    return output, np.maximum(retained, 0), full[len(inflow):]
