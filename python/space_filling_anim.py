"""
Boustrophedon space-filling animation frames + standalone length blow-up.

For the AMG-TREM ITL talk (preso-SST), space-filling-trap slide.  Produces:

  figs/sf_anim/sf_anim-0.png ... sf_anim-{K}.png
      a looped sequence in which the boustrophedon delta-sweep of the unit
      square is progressively refined (more passes, thinner tube) until its
      delta-tube fills the plane -- embedded in the deck via \animategraphics.

  figs/space_filling_blowup.pdf
      the companion static panel: required length L >~ delta^{1-n} vs delta,
      diverging as delta -> 0 and exploding with the intrinsic dimension n.

Palette matches sphere_amg_demo.GEO.  Run with the tda-sst venv:
    ~/venvs/tda-sst/Scripts/python space_filling_anim.py

Author: Zach Grey / Claude Code.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
ANIMDIR = os.path.join(FIGDIR, "sf_anim")

GEO = {
    "active":   "#D81159",
    "inactive": "#FFBC42",   # marigold tube
    "embed":    "#6A4C93",
    "mean":     "#0B3954",   # petrol centerline
}
plt.rcParams.update({"font.size": 13})


def boustrophedon(n_pass: int):
    """Serpentine centerline sweeping [0,1]^2 with n_pass horizontal passes;
    returns (x, y, delta) with delta = half the pass spacing."""
    ys = (np.arange(n_pass) + 0.5) / n_pass
    xl, xr = 0.04, 0.96
    x, y = [], []
    for k, yy in enumerate(ys):
        seq = (xl, xr) if k % 2 == 0 else (xr, xl)
        x += [seq[0], seq[1]]; y += [yy, yy]
        if k < n_pass - 1:
            x += [seq[1]]; y += [(ys[k] + ys[k + 1]) / 2]
    return np.array(x), np.array(y), 0.5 / n_pass


def render_frames(passes, figsize=4.2, dpi=120):
    os.makedirs(ANIMDIR, exist_ok=True)
    for i, n_pass in enumerate(passes):
        x, y, delta = boustrophedon(int(n_pass))
        fig, ax = plt.subplots(figsize=(figsize, figsize))
        tube_pts = (2 * delta) * (figsize / 1.0) * 72 * 0.92
        ax.plot(x, y, color=GEO["inactive"], lw=tube_pts, alpha=0.55,
                solid_capstyle="round", solid_joinstyle="round")
        ax.plot(x, y, color=GEO["mean"], lw=1.4)
        ax.add_patch(plt.Rectangle((0, 0), 1, 1, fill=False, ec="0.4", lw=1.2))
        ax.set_title(r"$\delta = %.3f$" % delta, fontsize=13)
        ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
        ax.set_aspect("equal"); ax.axis("off")
        fig.savefig(os.path.join(ANIMDIR, f"sf_anim-{i}.png"),
                    dpi=dpi, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
    print(f"wrote {len(passes)} frames to {ANIMDIR}")


def render_blowup():
    d = np.logspace(-3, 0, 200)
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    for n, c in [(2, GEO["mean"]), (3, GEO["embed"]), (5, GEO["active"])]:
        ax.loglog(d, d ** (1 - n), color=c, lw=2.6, label=f"$n={n}$")
    ax.set_xlabel(r"tube radius $\delta$")
    ax.set_ylabel(r"required length $L \gtrsim \delta^{\,1-n}$")
    ax.set_title(r"length blows up as $\delta\to0$")
    ax.invert_xaxis()
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, title="intrinsic dim.")
    ax.spines[["top", "right"]].set_visible(False)
    out = os.path.join(FIGDIR, "space_filling_blowup.pdf")
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def main():
    # refinement schedule: coarse -> fine (grid interpolation getting denser)
    passes = np.unique(np.round(np.linspace(4, 34, 20)).astype(int))
    render_frames(passes)
    render_blowup()


if __name__ == "__main__":
    main()
