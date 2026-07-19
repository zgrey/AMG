"""
Domain construction on a sphere: a tangent-space ellipse wrapped onto S^2.

Supplemental figure for the AMG-TREM ITL talk (preso-SST), motivating how a
manifold-valued domain is *constructed* -- the wrapped construction of
Def. riemannian_view in "A Riemannian View on Active Subspaces": choose a base
point p0, sample a distribution in the tangent space T_{p0}M (its covariance is
an ellipse), and push it onto the manifold by exp_{p0}.  The same picture serves
both provably-spherical domains of the talk (warping-function preshapes; shape
preshapes).

Self-contained closed-form S^2 geometry (no external deps beyond numpy/mpl).
Palette matches sphere_amg_demo.GEO.  Run with the tda-sst venv:
    ~/venvs/tda-sst/Scripts/python tangent_ellipse_sphere.py

Author: Zach Grey / Claude Code.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")

GEO = {
    "active":   "#D81159",   # ellipse / wrapped curve
    "inactive": "#FFBC42",   # samples
    "embed":    "#6A4C93",
    "mean":     "#0B3954",   # base point / tangent frame
}
plt.rcParams.update({"font.size": 13})


def tangent_basis(p0):
    """Orthonormal 3x2 basis of the tangent plane at unit vector p0."""
    a = np.array([1.0, 0.0, 0.0]) if abs(p0[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = a - (a @ p0) * p0
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(p0, e1)
    return np.column_stack([e1, e2])


def sphere_exp(p0, V):
    """exp_{p0}(v) for each row v of V (tangent vectors at p0)."""
    n = np.linalg.norm(V, axis=1, keepdims=True)
    out = np.cos(n) * p0 + np.sin(n) * np.divide(V, n, out=np.zeros_like(V), where=n > 0)
    return out


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    rng = np.random.default_rng(5)

    p0 = np.array([0.0, 0.0, 1.0])             # north pole: clean horizontal tangent plane
    E = tangent_basis(p0)                       # 3x2

    # ellipse (the tangent-space sampling covariance) in tangent coords
    a, b, ang = 1.05, 0.48, np.deg2rad(30)
    t = np.linspace(0, 2 * np.pi, 240)
    ex, ey = a * np.cos(t), b * np.sin(t)
    R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    ell2 = (R @ np.vstack([ex, ey])).T         # (240,2)
    ell_vecs = ell2 @ E.T                       # tangent vectors in R^3
    ell_plane = p0 + ell_vecs                   # drawn on the affine tangent plane
    ell_wrap = sphere_exp(p0, ell_vecs)         # wrapped onto the sphere

    # samples uniformly filling the ellipse
    m = 60
    r = np.sqrt(rng.uniform(0, 1, m))
    th = rng.uniform(0, 2 * np.pi, m)
    samp2 = 0.92 * (R @ np.vstack([a * r * np.cos(th), b * r * np.sin(th)])).T
    samp_vecs = samp2 @ E.T
    samp_plane = p0 + samp_vecs
    samp_wrap = sphere_exp(p0, samp_vecs)

    # tangent-plane patch corners
    s = 1.15
    sq = np.array([[-s, -s], [s, -s], [s, s], [-s, s]]) @ E.T + p0

    fig = plt.figure(figsize=(6.2, 5.6))
    ax = fig.add_subplot(111, projection="3d")

    # sphere
    uu = np.linspace(0, 2 * np.pi, 60)
    vv = np.linspace(0, np.pi, 40)
    xs = np.outer(np.cos(uu), np.sin(vv))
    ys = np.outer(np.sin(uu), np.sin(vv))
    zs = np.outer(np.ones_like(uu), np.cos(vv))
    ax.plot_surface(xs, ys, zs, color="0.9", alpha=0.30, rstride=2, cstride=2,
                    linewidth=0, shade=False)

    # tangent plane patch
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    ax.add_collection3d(Poly3DCollection([sq], facecolor=GEO["mean"], alpha=0.10,
                                         edgecolor=GEO["mean"], linewidths=0.8))

    # ellipse in tangent plane + wrapped image
    ax.plot(ell_plane[:, 0], ell_plane[:, 1], ell_plane[:, 2],
            color=GEO["active"], lw=2.4, label="ellipse in $T_{p_0}$")
    ax.plot(ell_wrap[:, 0], ell_wrap[:, 1], ell_wrap[:, 2],
            color=GEO["active"], lw=2.4, ls="--", label=r"wrapped: $\exp_{p_0}$")

    # samples
    ax.scatter(samp_plane[:, 0], samp_plane[:, 1], samp_plane[:, 2],
               color=GEO["inactive"], s=14, depthshade=False, alpha=0.9)
    ax.scatter(samp_wrap[:, 0], samp_wrap[:, 1], samp_wrap[:, 2],
               color=GEO["inactive"], s=16, depthshade=False, edgecolor="k", lw=0.2)

    # base point + exp arrow
    ax.scatter(*(p0 * 1.01), color="w", s=90, depthshade=False,
               edgecolor=GEO["mean"], linewidth=2.0)
    ax.text(*(p0 * 1.28), r"$p_0$", color=GEO["mean"], fontsize=13)
    k = 90
    ax.quiver(ell_plane[k, 0], ell_plane[k, 1], ell_plane[k, 2],
              *(ell_wrap[k] - ell_plane[k]), color="0.35", lw=1.4,
              arrow_length_ratio=0.25)

    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.view_init(elev=18, azim=-60)
    ax.legend(loc="upper left", frameon=False, fontsize=10, bbox_to_anchor=(0.02, 0.98))

    out = os.path.join(FIGDIR, "tangent_ellipse_sphere.pdf")
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
