"""
Dominant projection vs. projected dominance (ordering of projection), on S^2.

For the AMG-TREM ITL talk (preso-SST), supplemental slide.  On the ordered
quadratic (non-ridge) f = x^T diag(1,2,3) x over S^2, the extrinsic operator
C_iota has two comparable ambient eigenvalues, so the cos^2(theta_i) reweighting
of the tangential projection FLIPS the ordering.  The two candidate dominant
tangent directions then separate by a visible angle:

  dominant projection  pi_0[b_hat]         (decompose C_iota, THEN project)
  projected dominance  top eig( E0^T C_iota E0 )   (project, THEN decompose)

Both are drawn as frame fields over a geodesic cap around p_0:
  * dominant projection field  x_hat -> pi_{x_hat}[b_hat]        (rose)
  * projected dominance field  p     -> P^{-1}_p[w_star]         (petrol, unit)
with the angle phi between the two central directions annotated.

See REMARK-eig-projection-noncommute.md (Applied_AMG) for the theory.
Run with the tda-sst venv:
    ~/venvs/tda-sst/Scripts/python projection_order.py

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
from sphere_amg_demo import _view_direction

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
ROSE, PETROL = "#D81159", "#0B3954"


def cap_grid(p0, R=1.05, n_rings=5, n_az=16):
    E = amg.tangent_basis(p0)
    pts = [p0.copy()]
    for r in np.linspace(R / n_rings, R, n_rings):
        for th in np.linspace(0, 2 * np.pi, n_az, endpoint=False):
            tv = r * (np.cos(th) * E[:, 0] + np.sin(th) * E[:, 1])
            pts.append(amg.sphere_exp(p0, tv))
    return np.array(pts)


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    rng = np.random.default_rng(1)
    # Genuine non-ridge with a well-separated spectrum.  The ordered quadratic
    # only reorders at NEAR-DEGENERATE top eigenvalues (ill-defined dominant);
    # the nonlinear non-ridge reorders with BOTH the ambient and block dominants
    # well-separated (ratios ~1.56), so the angle is honest and reproducible.
    func = amg.make_function("nonlinear_nonridge")

    # base point + LARGE support: the discrepancy is an O(R^2) curvature effect,
    # so it needs a big cap; here the cos^2(theta) reweighting flips the ordering
    # -> phi ~ 60 deg between dominant projection and projected dominance.
    p0 = np.array([0.16, -0.76, 0.63]); p0 /= np.linalg.norm(p0)
    E = amg.tangent_basis(p0)

    # wrapped sample around p0 for the extrinsic operator C_iota
    m = 2500
    rr = 1.45 * np.sqrt(rng.random(m)); th = 2 * np.pi * rng.random(m)
    tv = (rr * np.cos(th))[:, None] * E[:, 0] + (rr * np.sin(th))[:, None] * E[:, 1]
    P = amg.sphere_exp_batch(p0, tv)

    res = amg.compute_amg(func, p0, P)
    C = res.Gt.T @ res.Gt / res.Gt.shape[0]        # ambient GOP C_iota (3x3)
    # dominant projection: project the ambient-dominant eigenvector of C
    val, vec = np.linalg.eigh(C); o = np.argsort(val)[::-1]
    b_hat = vec[:, o[0]]
    dp0 = res.W.copy()                              # = normalized pi_0[b_hat]
    dp0 /= np.linalg.norm(dp0)
    # projected dominance: top eigenvector of the projected block E^T C E
    B = E.T @ C @ E
    bv, bvec = np.linalg.eigh(B); bo = np.argsort(bv)[::-1]
    w_star = E @ bvec[:, bo[0]]; w_star /= np.linalg.norm(w_star)
    if dp0 @ w_star < 0:                            # sign-align for display
        w_star = -w_star
    phi = np.degrees(np.arccos(np.clip(dp0 @ w_star, -1, 1)))
    print("angle phi (deg) =", round(phi, 1),
          "| ambient eigvals:", np.round(val[o], 3))

    # ---- figure -------------------------------------------------------------
    elev, azim = 30, np.degrees(np.arctan2(p0[1], p0[0]))
    view = _view_direction(elev, azim)
    G = cap_grid(p0, R=1.0, n_rings=4, n_az=12)
    near = (G @ view) > 0.05

    fig = plt.figure(figsize=(6.6, 5.8))
    ax = fig.add_subplot(111, projection="3d")
    ax.computed_zorder = False

    nu, nv = 60, 30
    u = np.linspace(0, 2 * np.pi, nu); v = np.linspace(0, np.pi / 2, nv)
    ax.plot_surface(np.outer(np.cos(u), np.sin(v)), np.outer(np.sin(u), np.sin(v)),
                    np.outer(np.ones_like(u), np.cos(v)), color="0.85", alpha=0.14,
                    rstride=2, cstride=2, linewidth=0, shade=False, zorder=0)

    glen = 0.26
    for k in np.where(near)[0]:
        x = G[k]
        if np.linalg.norm(x - p0) < 0.12:                # keep a clear zone at p0
            continue
        dp = (np.eye(3) - np.outer(x, x)) @ b_hat        # dominant projection field
        if np.linalg.norm(dp) > 0.08:
            d = dp / np.linalg.norm(dp)
            ax.quiver(*x, *d, length=glen, normalize=False, color=ROSE,
                      lw=1.4, arrow_length_ratio=0.3, zorder=6)
        pd = amg.parallel_transport(p0, x, w_star)       # projected dominance field
        ax.quiver(*x, *pd, length=glen, normalize=False, color=PETROL,
                  lw=1.4, arrow_length_ratio=0.3, zorder=5, alpha=0.9)

    # base point marker only (no central dominant vectors, no angle overlay)
    ax.scatter(*(p0 * 1.02), color="w", s=70, edgecolor="k", linewidth=1.4,
               depthshade=False, zorder=11)

    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0], [0], color=ROSE, lw=3,
               label=r"dominant projection  $\pi_0\widehat{b}$  (decompose $\to$ project)"),
        Line2D([0], [0], color=PETROL, lw=3,
               label=r"projected dominance  (project $\to$ decompose)")],
        loc="lower center", frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, -0.02))

    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    lim = 0.7
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)

    out = os.path.join(FIGDIR, "projection_order.pdf")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.0)
    fig.savefig(out.replace(".pdf", ".png"), dpi=170, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
