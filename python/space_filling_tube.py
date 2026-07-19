"""
Space-filling (boustrophedon) tube over the plane, and its length blow-up.

Supplemental figure for the AMG-TREM ITL talk (preso-SST), illustrating the
ill-posedness of unconstrained single-curve nonlinear dimension reduction
(Remark rem:dichotomy / Appendix app:reach of "A Riemannian View on Active
Subspaces"):

  (left)  a boustrophedon delta-sweep of the unit square -- a 1-D curve whose
          delta-tube covers the 2-D domain, so it reproduces ANY continuous f
          within its modulus of continuity (reduction error -> 0 for every f);
  (right) the price: the tube length obeys  L >= C * delta^{1-n} - delta,
          which blows up as delta -> 0 (fixed n) and exponentially in the
          intrinsic dimension n (fixed delta).  Vanishing error therefore
          certifies no genuine 1-D structure.

Geodesic regularization (AMG) is the well-posed extreme: exp-images of straight
lines in one tangent space -- identifiable, low-complexity.

Palette matches sphere_amg_demo.GEO.  Run with the tda-sst venv:
    ~/venvs/tda-sst/Scripts/python space_filling_tube.py

Author: Zach Grey / Claude Code.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")

GEO = {
    "active":   "#D81159",
    "inactive": "#FFBC42",  # marigold (the delta-tube)
    "embed":    "#6A4C93",
    "mean":     "#0B3954",  # deep petrol (centerline)
}
plt.rcParams.update({"font.size": 14})


def boustrophedon(n_pass: int = 9):
    """Serpentine centerline sweeping [0,1]^2 with n_pass horizontal passes."""
    ys = (np.arange(n_pass) + 0.5) / n_pass
    xs_l, xs_r = 0.04, 0.96
    x, y = [], []
    for k, yy in enumerate(ys):
        seq = (xs_l, xs_r) if k % 2 == 0 else (xs_r, xs_l)
        x += [seq[0], seq[1]]
        y += [yy, yy]
        if k < n_pass - 1:                      # vertical connector at the end
            x += [seq[1]]
            y += [(ys[k] + ys[k + 1]) / 2]
    return np.array(x), np.array(y), 0.5 / n_pass  # delta = half the pass spacing


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6))

    # ---- left: boustrophedon delta-tube filling the plane ------------------
    n_pass = 9
    x, y, delta = boustrophedon(n_pass)
    # tube: fat semi-transparent marigold band (linewidth ~ 2*delta in axes pts)
    tube_pts = (2 * delta) * (fig.get_size_inches()[1] / 2) * 72 * 0.92
    axA.plot(x, y, color=GEO["inactive"], lw=tube_pts, alpha=0.55,
             solid_capstyle="round", solid_joinstyle="round")
    axA.plot(x, y, color=GEO["mean"], lw=1.6)        # centerline gamma
    axA.add_patch(plt.Rectangle((0, 0), 1, 1, fill=False, ec="0.4", lw=1.2))
    # a delta scale bar
    axA.annotate("", xy=(0.985, (0.5) / n_pass), xytext=(0.985, (1.5) / n_pass),
                 arrowprops=dict(arrowstyle="<->", color=GEO["mean"], lw=1.2))
    axA.text(0.955, 1.0 / n_pass, r"$\delta$", color=GEO["mean"],
             ha="right", va="center", fontsize=13)
    axA.text(0.5, -0.09, r"a $\delta$-tube covers $\mathcal{X}\subset\mathbb{R}^2$"
             + "\n" + r"$\Rightarrow$ reproduces any $f$",
             ha="center", va="top", fontsize=11, color="0.3")
    axA.set_title(r"boustrophedon $\delta$-sweep $\gamma$")
    axA.set_xlim(-0.02, 1.02); axA.set_ylim(-0.02, 1.02)
    axA.set_aspect("equal"); axA.axis("off")

    # ---- right: length blow-up L >= C delta^{1-n} --------------------------
    d = np.logspace(-3, 0, 200)
    C = 1.0
    for n, c in [(2, GEO["mean"]), (3, GEO["embed"]), (5, GEO["active"])]:
        L = C * d ** (1 - n)
        axB.loglog(d, L, color=c, lw=2.6, label=f"$n={n}$")
    axB.set_xlabel(r"tube radius $\delta$")
    axB.set_ylabel(r"required length $L \gtrsim \delta^{\,1-n}$")
    axB.set_title(r"length blows up as $\delta\to0$")
    axB.invert_xaxis()                               # delta -> 0 to the right
    axB.grid(True, which="both", alpha=0.25)
    axB.legend(frameon=False, title="intrinsic dim.")
    axB.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    out = os.path.join(FIGDIR, "space_filling_tube.pdf")
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out, " (delta = %.3f)" % delta)


if __name__ == "__main__":
    main()
