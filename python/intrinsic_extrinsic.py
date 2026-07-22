"""
Intrinsic vs. extrinsic: two upper hemispheres, identical tangent gradients.

For the AMG-TREM ITL talk (preso-SST), the "Intrinsic G_0 vs. extrinsic C_iota"
slide.  Two copies of the same upper-hemisphere response carrying the SAME
tangential gradient field, given two different treatments:

  left  (intrinsic): a thick BLACK WIREFRAME mesh of the manifold, with the
                     function shown as a filled colored scatter (lapaz) at the
                     gradient roots -- no ambient frame.
  right (extrinsic): no mesh, but a large orthogonal axes system at the origin
                     (the ambient R^3 embedding frame).

The function is a filled colored scatter at the roots of the gradients (point
evaluations), NOT a colored surface -- this keeps the mesh visible.  Same
gradients, different treatments.

Reuses the ambient function registry, tangential-gradient grid, and colormap
resolver from the sphere_amg demo.  Run with the tda-sst venv:
    ~/venvs/tda-sst/Scripts/python intrinsic_extrinsic.py

Author: Zach Grey / Claude Code.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

import sphere_amg as amg
from sphere_amg_demo import get_cmap, _sphere_grad_grid, _view_direction

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
ELEV, AZIM = 20, -55
LIM = 0.74


def _radial_lines(ax, phis):
    """Smooth radial meridian arcs (pole -> equator) at azimuths ``phis`` -- the
    lines pass through the gradient roots and show the sphere's constant
    curvature.  No latitude rings."""
    th = np.linspace(0, np.pi / 2, 80)                       # smooth arc
    for phi in phis:
        ax.plot(np.sin(th) * np.cos(phi), np.sin(th) * np.sin(phi), np.cos(th),
                color="k", lw=1.3, zorder=2)


def _common(ax, Pg, Gg, near, glen, colors_pg):
    ax.computed_zorder = False
    ax.scatter(Pg[near, 0], Pg[near, 1], Pg[near, 2], c=colors_pg[near],
               s=90, edgecolors="k", linewidths=0.6, depthshade=False, zorder=6)
    ax.quiver(Pg[near, 0], Pg[near, 1], Pg[near, 2],
              Gg[near, 0], Gg[near, 1], Gg[near, 2],
              length=glen, color="0.12", linewidth=1.5, normalize=False,
              alpha=0.95, zorder=8)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.view_init(elev=ELEV, azim=AZIM)
    ax.set_xlim(-LIM, LIM); ax.set_ylim(-LIM, LIM); ax.set_zlim(-LIM, LIM)


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    cmap = get_cmap("lapaz")
    func = amg.make_function("linear_random")

    Pg, Gg = _sphere_grad_grid(func, n=6)                    # SAME gradients both panels
    glen = 0.28 / max(np.linalg.norm(Gg, axis=1).max(), 1e-12)
    near = (Pg @ _view_direction(ELEV, AZIM)) > 0
    f_pg = func.f(Pg)
    colors_pg = cmap((f_pg - f_pg.min()) / (np.ptp(f_pg) + 1e-30))

    fig = plt.figure(figsize=(9.2, 4.7))

    # radial meridian azimuths = the gradient-root azimuths (pass through roots)
    ring = Pg[1:]                                            # drop the pole
    phis = np.unique(np.round(np.arctan2(ring[:, 1], ring[:, 0]), 5))

    # ---- left: intrinsic (radial meridian lines only, no latitude rings) ----
    axL = fig.add_subplot(1, 2, 1, projection="3d")
    _radial_lines(axL, phis)
    _common(axL, Pg, Gg, near, glen, colors_pg)
    axL.text2D(0.5, 0.04, "intrinsic", transform=axL.transAxes, ha="center",
               fontsize=16)

    # ---- right: extrinsic (no mesh, ambient axes at the origin) -------------
    axR = fig.add_subplot(1, 2, 2, projection="3d")
    _common(axR, Pg, Gg, near, glen, colors_pg)
    L = 0.7
    axcol = "#0B3954"
    for vec, lab in [((L, 0, 0), r"$e_1$"), ((0, L, 0), r"$e_2$"), ((0, 0, L), r"$e_3$")]:
        v = np.array(vec)
        axR.quiver(0, 0, 0, *v, color=axcol, lw=2.6, arrow_length_ratio=0.10, zorder=11)
        axR.text(*(v * 1.09), lab, color=axcol, fontsize=13, zorder=12)
    axR.scatter([0], [0], [0], color=axcol, s=24, zorder=11)
    axR.text2D(0.5, 0.04, "extrinsic", transform=axR.transAxes, ha="center",
               fontsize=16)

    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0, wspace=0.02)
    out = os.path.join(FIGDIR, "intrinsic_extrinsic.pdf")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.0)
    fig.savefig(out.replace(".pdf", ".png"), dpi=170, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
