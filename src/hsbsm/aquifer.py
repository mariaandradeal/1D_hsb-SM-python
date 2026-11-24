"""Conservative FV hsB fluxes and a checked, substepped Picard integrator."""
from __future__ import annotations
import numpy as np


def face_fluxes(grid, storage, aquifer, slope):
    """Direct translation of compute_face_fluxes in hsB_solver.m.

    x increases from outlet toward divide. Q[0] < 0 denotes outflow.
    Storage has units m² per unit hillslope length, not cell volume.
    """
    k, f = aquifer["k"], aquifer["f"]
    h = storage/(f*grid.w)
    dc, adv = -k/f*np.cos(slope), -k/f*np.sin(slope)
    speed = 2*abs(adv)
    q = np.zeros(len(storage)+1)
    film = max(f*grid.w[0]*0.002, 1e-5)
    q[0] = dc*max(storage[0], film)*h[0]/(grid.x[0]-grid.edges[0])
    q[0] += (adv-0.5*speed)*storage[0]
    if len(storage) > 1:
        floor = np.maximum(f*(grid.w[1:]+grid.w[:-1])/2*0.002, 1e-5)
        left, right = np.maximum(storage[:-1], floor), np.maximum(storage[1:], floor)
        harmonic = 2*left*right/(left+right)
        q[1:-1] = dc*harmonic*np.diff(h)/np.diff(grid.x)
        q[1:-1] += adv*storage[:-1]+0.5*speed*(storage[:-1]-storage[1:])
    return q


def fixed_step(grid, old, recharge, aquifer, slope, dt, numerics):
    """One MATLAB-compatible fixed step; reports convergence and mass residual."""
    theta, omega = numerics["theta"], numerics["omega"]
    cap = aquifer["f"]*grid.w*aquifer["D"]
    current = old.copy()
    previous = np.inf
    converged = False
    for iteration in range(int(numerics["picard_max"])):
        stage = (1-theta)*old+theta*current
        q = face_fluxes(grid, stage, aquifer, slope)
        raw = old+dt*(-np.diff(q)/grid.dx+recharge*grid.w)
        target = np.clip(raw, 0, cap)
        relaxed = (1-omega)*current+omega*target
        increment = np.max(np.abs(relaxed-current))
        if increment > 1.05*previous:
            omega = max(0.2, 0.5*omega)
        current = relaxed
        if increment < numerics["picard_tol"]:
            converged = True
            break
        previous = increment
    stage = (1-theta)*old+theta*current
    q = face_fluxes(grid, stage, aquifer, slope)
    raw = old+dt*(-np.diff(q)/grid.dx+recharge*grid.w)
    excess = np.maximum(raw-cap, 0)
    residual = np.sum((current-old)*grid.dx) - (
        np.sum(recharge*grid.area)*dt+q[0]*dt-np.sum(excess*grid.dx))
    return current, -q[0], excess*grid.dx/dt, float(residual), converged, iteration+1, raw


def advance(grid, storage, recharge, aquifer, slope, duration, numerics):
    """Accept only converged substeps, splitting signed extraction explicitly."""
    storage = storage.copy()
    recharge = np.broadcast_to(recharge, storage.shape).astype(float)
    withdrawal = np.maximum(-recharge, 0)*grid.w*duration
    if np.any(withdrawal > storage+1e-10):
        raise ValueError("Aquifer capillary extraction exceeds available storage")
    storage = np.maximum(storage-withdrawal, 0)
    recharge = np.maximum(recharge, 0)
    remaining, out_volume = float(duration), 0.0
    overflow_volume = np.zeros_like(storage)
    mass_error, steps, iterations = 0.0, 0, 0
    k, f = aquifer["k"], aquifer["f"]
    while remaining > max(1e-10, duration*1e-13):
        h = storage/(f*grid.w)
        speed = max(2*k/f*np.sin(slope), 1e-30)
        diffusion = max(k/f*max(float(h.max()), 0.002)*np.cos(slope), 1e-30)
        # Include the outlet half-cell in the stability estimate.
        spacing = min(grid.dx.min(), grid.x[0]-grid.edges[0])
        bound = numerics["cfl"]*min(spacing/speed, spacing**2/diffusion)
        dt = min(remaining, numerics["dt_max"], max(numerics["dt_min"], bound))
        while True:
            new, out, overflow, residual, ok, nit, raw = fixed_step(
                grid, storage, recharge, aquifer, slope, dt, numerics)
            scale = max(1.0, np.dot(storage, grid.dx), np.dot(recharge, grid.area)*dt)
            mass_tolerance = max(5*numerics["picard_tol"]*grid.dx.sum(), 1e-11*scale)
            if ok and np.min(raw) >= -numerics["picard_tol"] and abs(residual) <= mass_tolerance:
                break
            dt /= 2
            if dt < min(numerics["dt_min"], remaining):
                raise ArithmeticError("Aquifer iteration did not converge; reduce dt_min or check parameters")
        storage = new
        remaining -= dt
        out_volume += out*dt
        overflow_volume += overflow*dt
        mass_error += residual
        steps += 1
        iterations += nit
    return storage, out_volume/duration, overflow_volume/duration, {
        "residual_m3": mass_error, "substeps": steps, "iterations": iterations}
