"""Headless, publication-oriented figures with consistent SI labels."""
from pathlib import Path
import numpy as np


def plot_result(result, directory, comparison=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib import dates as mdates
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    available = {font.name for font in font_manager.fontManager.ttflist}
    family = next((f for f in ("Helvetica Light", "Helvetica Neue Light", "Helvetica", "Arial", "DejaVu Sans")
                   if f in available), "DejaVu Sans")
    style = {"font.family": family,
             "font.size": 10, "axes.linewidth": 1.5, "xtick.major.width": 1.5,
             "ytick.major.width": 1.5, "lines.linewidth": 2, "svg.fonttype": "none"}
    t, nodes = result.timeseries, result.nodes
    with plt.rc_context(style):
        fig, ax = plt.subplots(3, 1, figsize=(9, 8), sharex=True, layout="constrained")
        ax[0].bar(t.date, t.rain_mm_day, width=result.config.dt/86400, color="#97A6C4", label="Rainfall")
        ax[0].plot(t.date, t.ET_mm_day, color="#E25759", label="Actual ET")
        ax[0].set_ylabel("Water flux (mm/day)")
        ax[0].legend(frameon=False, ncol=2)
        colors = {"total": "#0B81A2", "surface": "#E25759", "base": "#59A89C", "deep": "#7E4794"}
        for name, color in colors.items():
            ax[1].plot(t.date, t[f"Q_{name}_m3s"], color=color, label=name.capitalize(),
                       lw=2 if name == "total" else 1.5)
        if comparison is not None:
            ax[1].plot(comparison.date, comparison.Q_obs_m3s, color="#555555", lw=1, label="Observed")
        ax[1].set_ylabel("Discharge (m³/s)")
        ax[1].legend(frameon=False, ncol=3)
        weights = result.config.grid.area
        ax[2].plot(t.date, np.average(nodes["theta_rz"], axis=1, weights=weights), label="Root zone", color="#0B81A2")
        ax[2].plot(t.date, np.average(nodes["theta_tz"], axis=1, weights=weights), label="Transmission zone", color="#E25759")
        ax[2].set_ylabel("Water content (m³/m³)")
        ax[2].legend(frameon=False, ncol=2, loc="upper center")
        ax[2].margins(y=0.18)
        ax[2].xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax[2].xaxis.get_major_locator()))
        for i, axis in enumerate(ax):
            axis.text(0.5, 1.02, f"({chr(97+i)})", transform=axis.transAxes, ha="center")
            axis.spines[["top", "right"]].set_visible(False)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(directory/f"hydrograph.{ext}", dpi=180)
        plt.close(fig)
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
        image = ax[0].pcolormesh(result.config.grid.x, np.arange(len(t))*result.config.dt/86400,
                                nodes["aquifer_head_m"], shading="auto", cmap="viridis", rasterized=True)
        ax[0].set(xlabel="Distance from outlet (m)", ylabel="Elapsed time (days)")
        fig.colorbar(image, ax=ax[0], label="Aquifer head (m)")
        flow = np.sort(t.Q_total_mm_day.to_numpy())[::-1]
        exceedance = 100*np.arange(1, len(flow)+1)/(len(flow)+1)
        ax[1].plot(exceedance, flow, color="#0B81A2")
        ax[1].set(xlabel="Exceedance probability (%)", ylabel="Discharge (mm/day)")
        for i, axis in enumerate(ax):
            axis.text(0.5, 1.02, f"({chr(97+i)})", transform=axis.transAxes, ha="center")
        for ext in ("png", "pdf", "svg"):
            fig.savefig(directory/f"states_and_duration.{ext}", dpi=180)
        plt.close(fig)
