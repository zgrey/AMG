"""
Normal-coordinate comparison gallery: intrinsic transport vs the two extrinsic
routes (Remark "Projected dominance versus dominant projection").

For each ambient test function (rows: preferential quadratic, nonlinear
non-ridge), three panels over the same normal-coordinate geodesic ball:

  1. **intrinsic (transport)** -- the G0 integrand: tangential gradients
     parallel-transported to T_{p0} (black quiver), with the transported
     dominant field w1(p) = P^{-1}_p[w1] (rose arrows).  In the display frame
     (below) the transported field is *constant*: one eigendecomposition,
     coherent everywhere.
  2. **dominant projection (project -> decompose)** -- the block integrand:
     centrally-projected gradients E0^T grad (black quiver), with the
     per-point dominant-projection field: at each grid point the top
     eigenvector of E_x^T C_iota E_x (petrol LINE segments -- an eigenvector
     field is only a line field, sign-ambiguous per point, and each point
     costs a fresh 2x2 eigensolve).
  3. **projected dominance (decompose -> project)** -- the same projected
     gradients, with the field pi_x b (purple arrows): the *fixed* dominant
     ambient eigenvector b of C_iota projected pointwise; a genuine vector
     field, but norm-decaying and carrying the ambient ordering.

Display-frame convention: every tangent vector at a point p is
parallel-transported to the central tangent space and expressed in the E
basis -- the same convention as the existing normal-coordinate gradient
quiver (the G0 integrand).  Axis traces: intrinsic active/inactive (rose /
marigold dashes), the central-block dominant axis (petrol dashes, panel 2),
and the projected ambient dominant axis pi_0 b (purple dots, panel 3).

Run (tda-sst venv):
    ~/venvs/tda-sst/Scripts/python normal_coords_compare.py
Options: --funcs, --N, --seed, --R, --cmap (defaults match the paper's
examples: N=3000, seed=47, R=0.7).
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

import sphere_amg as amg
from sphere_amg_demo import (
    CMAP, GEO, TRACE, FIGDIR, get_cmap, trace2d, _normal_to_sphere,
)

ROW_LABELS = {
    "quadratic_pref": "preferential quadratic",
    "nonlinear_nonridge": "nonlinear non-ridge",
}

# Per-row sampling configuration.  The nonlinear non-ridge row reuses the
# projection_order.py configuration -- a *scanned* center and large wrapped
# support where the cos^2(theta) reweighting genuinely flips the ambient
# ordering (at the paper-default Karcher-mean center the two extrinsic
# candidates differ by only ~3 deg; the reordering is center-dependent).
ROW_CONFIG = {
    "nonlinear_nonridge": dict(
        p0=np.array([0.16, -0.76, 0.63]) / np.linalg.norm([0.16, -0.76, 0.63]),
        wrap_radius=1.45, m=2500, rng_seed=1),
}


def sample_row(name, N, seed):
    """Samples + center for one row: paper defaults, or the ROW_CONFIG
    override (wrapped uniform disk in the tangent space at a fixed center)."""
    cfg = ROW_CONFIG.get(name)
    if cfg is None:
        rng = np.random.default_rng(seed)
        P = amg.sample_disk(N, rng=rng)
        return amg.karcher_mean(P), P
    rng = np.random.default_rng(cfg["rng_seed"])
    p0 = cfg["p0"]
    E = amg.tangent_basis(p0)
    rr = cfg["wrap_radius"] * np.sqrt(rng.random(cfg["m"]))
    th = 2 * np.pi * rng.random(cfg["m"])
    tv = (rr * np.cos(th))[:, None] * E[:, 0] \
        + (rr * np.sin(th))[:, None] * E[:, 1]
    return p0, amg.sphere_exp_batch(p0, tv)
COL_TITLES = [
    "intrinsic (transport)",
    "dominant projection\n(project $\\rightarrow$ decompose)",
    "projected dominance\n(decompose $\\rightarrow$ project)",
]


def _to_center_frame(P_pts, p0, E, vecs):
    """Parallel-transport tangent vectors (rows of ``vecs`` at points
    ``P_pts``) to T_{p0} and return their E-basis components, (M, 2)."""
    out = np.array([amg.parallel_transport(P_pts[i], p0, vecs[i])
                    for i in range(len(P_pts))])
    return out @ E


def _grid(R, n, fill=0.93):
    g = np.linspace(-R * fill, R * fill, n)
    T1, T2 = np.meshgrid(g, g)
    t = np.column_stack([T1.ravel(), T2.ravel()])
    return t[np.hypot(t[:, 0], t[:, 1]) <= R * fill]


def _background(ax, func, res, R, ngrid=170):
    """f o exp contour over the geodesic ball (shared by all three panels)."""
    p0, E = res.p0, res.E
    g = np.linspace(-R, R, ngrid)
    T1, T2 = np.meshgrid(g, g)
    RR = np.hypot(T1, T2)
    tnml = np.column_stack([T1.ravel(), T2.ravel()])
    F = func.f(_normal_to_sphere(p0, E, tnml)).reshape(ngrid, ngrid)
    F = np.where(RR <= R, F, np.nan)
    cs = ax.contourf(T1, T2, F, 60, cmap=CMAP)
    ax.contour(T1, T2, F, 14, colors="k", linewidths=0.3, alpha=0.35)
    th = np.linspace(0, 2 * np.pi, 240)
    ax.plot(R * np.cos(th), R * np.sin(th), "k-", lw=1.2)
    ax.set_aspect("equal")
    ax.set_xlim(-1.06 * R, 1.06 * R)
    ax.set_ylim(-1.06 * R, 1.06 * R)
    ax.set_xticks([]), ax.set_yticks([])
    return cs


def _grad_quiver(ax, t, pg, scale=None):
    """Black gradient quiver at base points ``t`` with components ``pg``."""
    ax.quiver(t[:, 0], t[:, 1], pg[:, 0], pg[:, 1], color="k",
              width=0.0045, alpha=0.6, zorder=3, scale=scale,
              scale_units="xy", angles="xy")


def _field_arrows(ax, t, v, color, L=0.16, lw=0.006):
    """Vector-field arrows (pivot mid, display length L per unit norm)."""
    ax.quiver(t[:, 0], t[:, 1], v[:, 0], v[:, 1], color=color, pivot="mid",
              width=lw, scale=1.0 / L, scale_units="xy", angles="xy",
              zorder=5)


def _field_lines(ax, t, v, color, L=0.16, lw=2.0):
    """Headless line segments (eigenvector *line* field, sign-ambiguous)."""
    for (t1, t2), (v1, v2) in zip(t, v):
        n = np.hypot(v1, v2)
        if n < 1e-12:
            continue
        d1, d2 = v1 / n * L / 2, v2 / n * L / 2
        ax.plot([t1 - d1, t1 + d1], [t2 - d2, t2 + d2], color=color,
                lw=lw, solid_capstyle="round", zorder=5)


def _axis_trace(ax, u2, R, key):
    """Central axis through the origin along E-components ``u2`` (styled)."""
    n = np.linalg.norm(u2)
    if n < 1e-12:
        return
    u = u2 / n * R
    trace2d(ax, [-u[0], u[0]], [-u[1], u[1]], key, z=7)


def compare_row(axes, func, res, C_iota, R=0.7, ngrid=170):
    """Fill one gallery row (three axes).  Returns the contourf handle."""
    p0, E = res.p0, res.E
    b_hat = np.linalg.eigh(C_iota)[1][:, -1]        # dominant ambient e-vec
    # central-block (dominant-projection) axes at p0
    B0 = E.T @ C_iota @ E
    val0, vec0 = np.linalg.eigh(B0)
    block_top = vec0[:, -1]                          # E-components, T_{p0}
    # sign-align display axes to the intrinsic active direction
    u1 = E.T @ res.U_active
    u2i = E.T @ res.U_inactive
    if block_top @ u1 < 0:
        block_top = -block_top
    pb0 = E.T @ ((np.eye(3) - np.outer(p0, p0)) @ b_hat)  # pi_0 b, E-comps

    # gradient quivers: transported (panel 1) vs centrally-projected (2, 3)
    tg = _grid(R, 15)
    xg = _normal_to_sphere(p0, E, tg)
    Gt = amg.project_tangent(func.grad(xg), xg)
    pg_transport = _to_center_frame(xg, p0, E, Gt)
    pg_project = Gt @ E

    # field grids
    tf = _grid(R, 9)
    xf = _normal_to_sphere(p0, E, tf)
    # panel 1: transported dominant field, in the center frame == constant u1
    v_intr = np.tile(u1, (len(tf), 1))
    # panel 2: per-point top eigenvector of E_x^T C_iota E_x (fresh eigensolve
    # per point; raw eigh output -- sign is genuinely ambiguous)
    v_block = np.zeros((len(tf), 2))
    for i, x in enumerate(xf):
        Ex = amg.tangent_basis(x)
        _, vec = np.linalg.eigh(Ex.T @ C_iota @ Ex)
        v_block[i] = _to_center_frame(x[None, :], p0, E,
                                      (Ex @ vec[:, -1])[None, :])[0]
    # panel 3: pi_x b, a genuine (norm-decaying) vector field
    v_proj = _to_center_frame(
        xf, p0, E, amg.project_tangent(np.tile(b_hat, (len(xf), 1)), xf))

    gscale = float(np.percentile(np.linalg.norm(pg_transport, axis=1), 95)) \
        / (0.12 * R) if len(pg_transport) else None

    cs = _background(axes[0], func, res, R, ngrid)
    _grad_quiver(axes[0], tg, pg_transport, scale=gscale)
    _field_arrows(axes[0], tf, v_intr, GEO["active"], L=0.14 * R / 0.7)
    _axis_trace(axes[0], u1, R, "active")
    _axis_trace(axes[0], u2i, R, "inactive")

    _background(axes[1], func, res, R, ngrid)
    _grad_quiver(axes[1], tg, pg_project, scale=gscale)
    _field_lines(axes[1], tf, v_block, GEO["mean"], L=0.14 * R / 0.7)
    _axis_trace(axes[1], block_top, R, "block")

    _background(axes[2], func, res, R, ngrid)
    _grad_quiver(axes[2], tg, pg_project, scale=gscale)
    _field_arrows(axes[2], tf, v_proj, GEO["embedding"], L=0.14 * R / 0.7)
    _axis_trace(axes[2], pb0, R, "embedding")

    # honest diagnostics for the caption / sanity checks
    lam_blk = np.sort(val0)[::-1]
    ang = lambda a, b: np.degrees(np.arccos(np.clip(
        abs(a @ b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-30), 0, 1)))
    print(f"  [{func.name}]")
    print(f"    intrinsic lam        = {np.round(res.lam, 4)}")
    print(f"    block lam (E0'CE0)   = {np.round(lam_blk, 4)}")
    print(f"    ambient lam (C_iota) = {np.round(res.lam_emb, 4)}")
    print(f"    angle(w1, block-top) = {ang(u1, block_top):6.2f} deg")
    print(f"    angle(w1, pi0 b)     = {ang(u1, pb0):6.2f} deg   "
          f"||pi0 b|| = {np.linalg.norm(pb0):.3f}")
    return cs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--funcs", nargs="+",
                    default=["quadratic_pref", "nonlinear_nonridge"])
    ap.add_argument("--N", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=47)
    ap.add_argument("--R", type=float, default=0.7)
    ap.add_argument("--cmap", default="lapaz")
    args = ap.parse_args()

    global CMAP
    CMAP = get_cmap(args.cmap)
    TRACE.setdefault("block", dict(color=GEO["mean"], ls="--", lw=3, halo=2.25))

    nrow = len(args.funcs)
    fig, axes = plt.subplots(nrow, 3, figsize=(12.6, 4.15 * nrow))
    axes = np.atleast_2d(axes)

    for r, name in enumerate(args.funcs):
        func = amg.make_function(name, seed=args.seed)
        p0, P = sample_row(name, args.N, args.seed)
        res = amg.compute_amg(func, p0, P)
        C_iota = (res.Gt.T @ res.Gt) / len(P)
        cs = compare_row(axes[r], func, res, C_iota, R=args.R)
        axes[r, 0].set_ylabel(ROW_LABELS.get(name, name))
        fig.colorbar(cs, ax=list(axes[r]), shrink=0.9, pad=0.015,
                     fraction=0.035, label=r"$f\circ\exp_{p_0}$")
    for c, title in enumerate(COL_TITLES):
        axes[0, c].set_title(title)

    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ("png", "pdf"):
        path = os.path.join(FIGDIR, f"normal_coords_compare.{ext}")
        fig.savefig(path, dpi=170, bbox_inches="tight")
        print(f"written {path}")


if __name__ == "__main__":
    main()
