"""
Intrinsic frame field w_i(p) vs. the projected (extrinsic) frame, on S^2.

For the AMG-TREM ITL talk (preso-SST).  Visualizes the contrast developed in
sec. "The AMG frame field" of Applied_AMG-arxiv:

  intrinsic  w_1(p) = P^{-1}_p[w_1]   (eq:AMG_field): inverse parallel transport
             of the central eigenvector along radial geodesics.  Smooth,
             orthonormal, and UNIT-NORM everywhere.

  projected  x_hat -> pi_{x_hat}[b_hat]: pointwise orthogonal projection of a
             fixed ambient direction b_hat onto the tangent space.  Smooth but
             DEGENERATES -- its magnitude is the cosine of the angle to the
             tangent space, vanishing where b_hat is normal to the manifold
             (on the sphere, at x_hat = +/- b_hat).

Both fields agree at the base point p_0 (b_hat is chosen so pi_{p_0} b_hat = w_1);
away from p_0 the projected field shrinks to zero at +/- b_hat while the intrinsic
field stays unit length.  Arrow length AND color encode |field|, shared scale.

Parallel transport reused from sphere_amg (EVIE).  Run with the tda-sst venv:
    ~/venvs/tda-sst/Scripts/python frame_field_compare.py

Author: Zach Grey / Claude Code.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

import sphere_amg as amg
from sphere_amg_demo import get_cmap, _view_direction

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
ELEV, AZIM = 22, -60
LIM = 0.74
GLEN = 0.34


def hemi_grid(n_rings=6, n_az=16):
    """Points on the upper hemisphere (pole + rings), for the field arrows."""
    pts = [np.array([0.0, 0.0, 1.0])]
    for th in np.linspace(np.pi / 2 / n_rings, np.pi / 2 * 0.99, n_rings):
        for ph in np.linspace(0, 2 * np.pi, n_az, endpoint=False):
            pts.append(np.array([np.sin(th) * np.cos(ph),
                                 np.sin(th) * np.sin(ph), np.cos(th)]))
    return np.array(pts)


def _surface(ax):
    nu, nv = 60, 30
    u = np.linspace(0, 2 * np.pi, nu)
    v = np.linspace(0, np.pi / 2, nv)
    X = np.outer(np.cos(u), np.sin(v))
    Y = np.outer(np.sin(u), np.sin(v))
    Z = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(X, Y, Z, color="0.85", alpha=0.16, rstride=2, cstride=2,
                    linewidth=0, shade=False, zorder=0)
    eu = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(eu), np.sin(eu), np.zeros_like(eu), "k-", lw=1.8, zorder=1)


def _finish(ax, title):
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.view_init(elev=ELEV, azim=AZIM)
    ax.set_xlim(-LIM, LIM); ax.set_ylim(-LIM, LIM); ax.set_zlim(-LIM, LIM)
    # panel titles intentionally omitted -- the slide's equations carry them


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    cmap = get_cmap("batlow").reversed()          # 1 -> bold/dark, 0 -> faint
    norm = Normalize(0.0, 1.0)

    p0 = np.array([0.0, 0.0, 1.0])                 # base point (pole)
    w1 = np.array([1.0, 0.0, 0.0])                 # central eigenvector (tangent)
    bhat = w1.copy()                               # ambient dir; pi_{p0} bhat = w1

    P = hemi_grid()
    view = _view_direction(ELEV, AZIM)
    near = (P @ view) > 0.02

    # intrinsic frame field: parallel transport of w1 from p0 to each p (unit)
    Wi = np.array([amg.parallel_transport(p0, p, w1) if np.linalg.norm(p - p0) > 1e-9
                   else w1 for p in P])
    # projected field: pi_p[bhat] = bhat - (p.bhat) p
    Wp = bhat[None, :] - (P @ bhat)[:, None] * P
    mag_i = np.linalg.norm(Wi, axis=1)
    mag_p = np.linalg.norm(Wp, axis=1)

    fig = plt.figure(figsize=(9.6, 4.9))

    def draw(ax, W, mag):
        _surface(ax)
        for k in np.where(near)[0]:
            c = cmap(norm(mag[k]))
            ax.quiver(P[k, 0], P[k, 1], P[k, 2], W[k, 0], W[k, 1], W[k, 2],
                      length=GLEN, normalize=False, color=c, linewidth=2.0,
                      arrow_length_ratio=0.32, zorder=6)
        ax.scatter(*(p0 * 1.02), color="w", s=70, edgecolor="#0B3954",
                   linewidth=1.8, depthshade=False, zorder=9)

    axL = fig.add_subplot(1, 2, 1, projection="3d")
    draw(axL, Wi, mag_i)
    _finish(axL, None)

    axR = fig.add_subplot(1, 2, 2, projection="3d")
    draw(axR, Wp, mag_p)
    # mark the degeneracy points +/- b_hat (projected field vanishes)
    for s in (+1, -1):
        d = s * bhat
        if d @ view > 0.02 and d[2] >= -1e-6:
            axR.scatter(*(d * 1.02), color="#D81159", marker="x", s=90,
                        linewidths=2.5, depthshade=False, zorder=10)
    _finish(axR, None)

    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cax = fig.add_axes([0.42, 0.16, 0.16, 0.028])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label(r"$\|\text{field}\|$", fontsize=11)
    cb.set_ticks([0, 1])

    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0, wspace=0.0)
    out = os.path.join(FIGDIR, "frame_field_compare.pdf")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.0)
    fig.savefig(out.replace(".pdf", ".png"), dpi=170, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)
    print("wrote", out, " (mag_i in [%.2f,%.2f], mag_p in [%.2f,%.2f])"
          % (mag_i.min(), mag_i.max(), mag_p.min(), mag_p.max()))


if __name__ == "__main__":
    main()
