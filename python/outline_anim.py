"""
Outline-slide animation: random signal traces transporting across a chip plane.

For the AMG-TREM ITL talk (preso-SST), title/outline slide.  A field of many
wavy horizontal tracks (inverter lines) fills a die-shaped plane; along each, a
colored pulse ("comet") transports left-to-right at a random speed.  The intent
is to motivate that a chip is essentially thousands of transported signals that
together encode its functionality.

Seamless loop: each pulse makes an integer number of traversals per loop, so
frame T lands back on frame 0.  Frames are embedded via \animategraphics.

Palette: cmcrameri batlow (falls back to matplotlib).  Run with the tda-sst venv:
    ~/venvs/tda-sst/Scripts/python outline_anim.py

Author: Zach Grey / Claude Code.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
ANIMDIR = os.path.join(FIGDIR, "outline_anim")


def get_cmap(name="batlow"):
    try:
        from cmcrameri import cm as cmc
        if hasattr(cmc, name):
            return getattr(cmc, name)
    except Exception:
        pass
    return plt.get_cmap("plasma")


def main():
    os.makedirs(ANIMDIR, exist_ok=True)
    rng = np.random.default_rng(7)
    T = 24                       # frames per loop
    N = 34                       # number of signal tracks
    cmap = get_cmap("batlow")

    y0 = np.linspace(0.05, 0.95, N) + rng.uniform(-0.01, 0.01, N)
    speed = rng.integers(1, 4, N)          # integer traversals/loop -> seamless
    phase = rng.uniform(0, 1, N)
    amp = rng.uniform(0.0, 0.018, N)       # routing waviness
    freq = rng.integers(1, 4, N)
    colors = cmap(rng.uniform(0.1, 0.9, N))
    xs = np.linspace(0, 1, 200)
    ntrail = 9

    def y_at(k, x):
        return y0[k] + amp[k] * np.sin(2 * np.pi * freq[k] * x)

    for t in range(T):
        fig, ax = plt.subplots(figsize=(4.2, 3.4))
        # chip die
        ax.add_patch(plt.Rectangle((0, 0), 1, 1, facecolor="#0B1A24",
                                   edgecolor="#2A4A5A", lw=1.5, zorder=0))
        for k in range(N):
            yk = y_at(k, xs)
            ax.plot(xs, yk, color="#1E3A48", lw=0.8, zorder=1)   # faint wire
            pos = (phase[k] + speed[k] * t / T) % 1.0
            # comet trail (no wrap: fades in as it enters from the left)
            for j in range(ntrail):
                px = pos - j * 0.028
                if px < 0:
                    break
                a = 1.0 - j / ntrail
                ax.plot(px, y_at(k, px), marker="o", ms=4.5 * a + 1.0,
                        color=colors[k], alpha=0.9 * a, zorder=3,
                        markeredgecolor="none")
            ax.plot(pos, y_at(k, pos), marker="o", ms=6.0, color=colors[k],
                    markeredgecolor="w", markeredgewidth=0.4, zorder=4)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.axis("off")
        fig.savefig(os.path.join(ANIMDIR, f"outline_anim-{t}.png"),
                    dpi=120, bbox_inches="tight", pad_inches=0.02)
        plt.close(fig)
    print(f"wrote {T} frames to {ANIMDIR}")


if __name__ == "__main__":
    main()
