"""
Plain response sphere: a scalar function f colored on S^2 (lapaz), no overlays.

Motivation figure for the AMG-TREM ITL talk (preso-SST), latent-space slide
(f: Theta -> X -> R): a clean "response restricted to the domain manifold"
picture with no AMG geodesics, gradients, samples, or tangent plane -- just the
colored field, using the lapaz colormap for consistency with the other sphere
figures.

Reuses the ambient function registry and colormap resolver from the sphere_amg
demo.  Run with the tda-sst venv:
    ~/venvs/tda-sst/Scripts/python sphere_response.py [--func nonlinear_ridge]

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
from sphere_amg_demo import get_cmap

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--func", default="nonlinear_ridge")
    ap.add_argument("--cmap", default="lapaz")
    args = ap.parse_args()

    os.makedirs(FIGDIR, exist_ok=True)
    cmap = get_cmap(args.cmap)
    func = amg.make_function(args.func)

    nu, nv = 160, 120
    u = np.linspace(0, 2 * np.pi, nu)
    v = np.linspace(0, np.pi, nv)
    X = np.outer(np.cos(u), np.sin(v))
    Y = np.outer(np.sin(u), np.sin(v))
    Z = np.outer(np.ones_like(u), np.cos(v))
    F = func.f(np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])).reshape(nu, nv)
    facec = cmap((F - F.min()) / (np.ptp(F) + 1e-30))

    fig = plt.figure(figsize=(5.2, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z, facecolors=facec, rstride=1, cstride=1,
                    linewidth=0, antialiased=True, shade=False)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.view_init(elev=20, azim=-55)
    # tight crop around the sphere
    ax.set_xlim(-0.72, 0.72); ax.set_ylim(-0.72, 0.72); ax.set_zlim(-0.72, 0.72)

    out = os.path.join(FIGDIR, "sphere_response.pdf")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.0)
    fig.savefig(out.replace(".pdf", ".png"), dpi=170, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)
    print("wrote", out, "(func=%s, cmap=%s)" % (args.func, args.cmap))


if __name__ == "__main__":
    main()
