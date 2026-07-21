"""
Outline-slide animation: signals traversing a maze of grid-routed tracks.

For the AMG-TREM ITL talk (preso-SST), title/outline slide.  Many random
Manhattan (axis-aligned) paths route left-to-right across a chip-die plane over
a grid; along each, a colored pulse ("comet") transports and turns corners --
signals traversing a maze.  Motivates that a chip is essentially thousands of
transported signals that together encode its functionality.

Seamless loop: each pulse makes an integer number of traversals of its path per
loop, so frame T lands back on frame 0.  Frames embedded via \animategraphics.

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
from matplotlib.collections import LineCollection

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
ANIMDIR = os.path.join(FIGDIR, "outline_anim")

G = 16                 # grid resolution
YLO, YHI = 0.04, 0.96  # vertical extent of the routable plane


def get_cmap(name="batlow"):
    try:
        from cmcrameri import cm as cmc
        if hasattr(cmc, name):
            return getattr(cmc, name)
    except Exception:
        pass
    return plt.get_cmap("plasma")


def make_path(rng):
    """Random axis-aligned (Manhattan) route from the left edge to the right
    edge over a G x G grid, with vertical detours -- a maze-like wire."""
    xi, yi = 0, int(rng.integers(0, G + 1))
    vdir = int(rng.choice([-1, 1]))
    nodes = [(xi, yi)]
    guard = 0
    while xi < G and guard < 8 * G:
        guard += 1
        r = rng.random()
        go_right = r < 0.52 or yi in (0, G)      # bias rightward; leave the rails
        if go_right:
            xi += 1
        else:
            if not (0 <= yi + vdir <= G):
                vdir = -vdir
            yi += vdir
            if rng.random() < 0.3:               # occasional turn
                vdir = -vdir
        nodes.append((xi, yi))
    while xi < G:                                # guarantee we reach the right edge
        xi += 1
        nodes.append((xi, yi))
    P = np.array(nodes, float)
    P[:, 0] = P[:, 0] / G                          # x in [0,1]
    P[:, 1] = YLO + (P[:, 1] / G) * (YHI - YLO)    # y in [YLO,YHI]
    return P


def arclen(P):
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
    return np.concatenate([[0.0], np.cumsum(seg)])


def point_at(P, cum, s):
    """Point at arc length s along polyline P (cum = cumulative lengths)."""
    s = np.clip(s, 0.0, cum[-1])
    i = int(np.searchsorted(cum, s) - 1)
    i = max(0, min(i, len(P) - 2))
    seglen = cum[i + 1] - cum[i]
    f = 0.0 if seglen <= 0 else (s - cum[i]) / seglen
    return P[i] * (1 - f) + P[i + 1] * f


def main():
    os.makedirs(ANIMDIR, exist_ok=True)
    rng = np.random.default_rng(11)
    T = 24
    N = 26
    cmap = get_cmap("batlow")

    paths = [make_path(rng) for _ in range(N)]
    cums = [arclen(P) for P in paths]
    speed = rng.integers(1, 3, N)           # integer traversals/loop -> seamless
    phase = rng.uniform(0, 1, N)
    colors = cmap(rng.uniform(0.1, 0.9, N))
    ntrail, ds = 9, 0.022

    gridlines = [[(i / G, YLO), (i / G, YHI)] for i in range(G + 1)] + \
                [[(0, YLO + j / G * (YHI - YLO)), (1, YLO + j / G * (YHI - YLO))]
                 for j in range(G + 1)]

    for t in range(T):
        fig, ax = plt.subplots(figsize=(4.2, 3.4))
        ax.add_patch(plt.Rectangle((0, 0), 1, 1, facecolor="#0B1A24",
                                   edgecolor="#2A4A5A", lw=1.5, zorder=0))
        ax.add_collection(LineCollection(gridlines, colors="#14303C",
                                         linewidths=0.5, zorder=1))
        for k in range(N):
            P, cum = paths[k], cums[k]
            ax.plot(P[:, 0], P[:, 1], color="#33627A", lw=1.3, zorder=2)  # wire
            L = cum[-1]
            f = (phase[k] + speed[k] * t / T) % 1.0
            s0 = f * L
            for j in range(ntrail):                # comet trail (no wrap)
                s = s0 - j * ds
                if s < 0:
                    break
                a = 1.0 - j / ntrail
                p = point_at(P, cum, s)
                ax.plot(p[0], p[1], marker="o", ms=4.5 * a + 1.0,
                        color=colors[k], alpha=0.9 * a, zorder=4,
                        markeredgecolor="none")
            ph = point_at(P, cum, s0)
            ax.plot(ph[0], ph[1], marker="o", ms=6.0, color=colors[k],
                    markeredgecolor="w", markeredgewidth=0.4, zorder=5)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.axis("off")
        fig.savefig(os.path.join(ANIMDIR, f"outline_anim-{t}.png"),
                    dpi=120, bbox_inches="tight", pad_inches=0.02)
        plt.close(fig)
    print(f"wrote {T} frames to {ANIMDIR}")


if __name__ == "__main__":
    main()
