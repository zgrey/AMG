"""
Response sphere: a scalar function f colored on S^2 (lapaz), no AMG geodesics.

Motivation figures for the AMG-TREM ITL talk (preso-SST).  Two modes:

  --hemisphere full  (default): the whole sphere colored by f, no overlays --
      the plain "response restricted to the domain" picture.
  --hemisphere upper --gradients: an opaque upper-hemisphere cap colored by f
      with the tangential gradient field (arrows), the equator ring, and NO
      colored active/inactive geodesics -- the domain-reconciliation slide.

Colormap defaults to lapaz for consistency with the other sphere figures.
Reuses the ambient function registry, the tangential-gradient grid, and the
colormap resolver from the sphere_amg demo.  Run with the tda-sst venv, e.g.:
    ~/venvs/tda-sst/Scripts/python sphere_response.py \
        --func linear_random --hemisphere upper --gradients

Author: Zach Grey / Claude Code.
"""
from __future__ import annotations

import argparse
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

import sphere_amg as amg
from sphere_amg_demo import get_cmap, _sphere_grad_grid, _view_direction


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--func", default="nonlinear_ridge")
    ap.add_argument("--cmap", default="lapaz")
    ap.add_argument("--hemisphere", choices=["full", "upper"], default="full")
    ap.add_argument("--gradients", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    figdir = os.path.join(os.path.dirname(__file__), "figs")
    os.makedirs(figdir, exist_ok=True)
    cmap = get_cmap(args.cmap)
    func = amg.make_function(args.func)

    nu, nv = 160, 120
    vmax = np.pi / 2 if args.hemisphere == "upper" else np.pi
    u = np.linspace(0, 2 * np.pi, nu)
    v = np.linspace(0, vmax, nv)
    X = np.outer(np.cos(u), np.sin(v))
    Y = np.outer(np.sin(u), np.sin(v))
    Z = np.outer(np.ones_like(u), np.cos(v))
    F = func.f(np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])).reshape(nu, nv)
    facec = cmap((F - F.min()) / (np.ptp(F) + 1e-30))

    elev, azim = 20, -55
    fig = plt.figure(figsize=(5.2, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.computed_zorder = False
    ax.plot_surface(X, Y, Z, facecolors=facec, rstride=1, cstride=1,
                    linewidth=0, antialiased=True, shade=False, zorder=2)

    if args.hemisphere == "upper":
        eu = np.linspace(0, 2 * np.pi, 300)          # equator ring
        ax.plot(np.cos(eu), np.sin(eu), np.zeros_like(eu), "k-", lw=2.5, zorder=1)

    if args.gradients:
        Pg, Gg = _sphere_grad_grid(func, n=8)         # tangential gradients (upper cap)
        gmax = np.linalg.norm(Gg, axis=1).max()
        glen = 0.28 / gmax if gmax > 1e-12 else 0.28
        view = _view_direction(elev, azim)
        near = (Pg @ view) > 0                         # near side only
        ax.quiver(Pg[near, 0], Pg[near, 1], Pg[near, 2],
                  Gg[near, 0], Gg[near, 1], Gg[near, 2],
                  length=glen, color="0.12", linewidth=1.4, normalize=False,
                  alpha=0.95, zorder=7)

    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    lim = 0.72
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)

    out = args.out or os.path.join(figdir, "sphere_response.pdf")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.0)
    fig.savefig(out.replace(".pdf", ".png"), dpi=170, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)
    print("wrote", out, "(func=%s, cmap=%s, %s%s)"
          % (args.func, args.cmap, args.hemisphere,
             ", gradients" if args.gradients else ""))


if __name__ == "__main__":
    main()
