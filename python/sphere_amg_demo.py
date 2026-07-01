"""
Active Manifold-Geodesics on the 2-sphere -- visualisations.

Python/matplotlib port of the figures in ``emb_mfld_AMG.m``.  Produces:

  1. normal-coordinate contour of the function with the projected-gradient
     quiver (the local tangent-plane view at the central point p0);
  2. the function on the sphere (3-D) with samples, tangential gradients,
     the Karcher mean, tangent plane, and the active / inactive / embedding
     manifold-geodesics;
  3. the "shadow plot": f composed with the exponential map along a sweep of
     geodesics, highlighting the active (AMG), inactive (IAMG), and embedding
     sweeps, over the scatter of samples projected onto the active coordinate;
  4. the subspace-distance convergence study (loglog vs sample size N).

Run:  python sphere_amg_demo.py --func linear_random
Manifold geometry is reused from EVIE (see sphere_amg.py).

Author: Zach Grey (MATLAB); Python port: Claude Code.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless; we save PNGs
import matplotlib.pyplot as plt
from matplotlib import cm

import sphere_amg as amg

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")


def get_cmap(name="lapaz"):
    """Resolve a colormap by name: Crameri scientific colour maps (cmcrameri)
    first -- e.g. 'lapaz', 'batlow', 'vik', 'davos', 'roma' -- then matplotlib
    built-ins.  Default 'lapaz'."""
    try:
        from cmcrameri import cm as cmc
        if hasattr(cmc, name):
            return getattr(cmc, name)
    except Exception:
        pass
    return plt.get_cmap(name)


CMAP = get_cmap("lapaz")

# Geodesic accent palette -- warm/vivid hues that read on both the cool lapaz
# sphere and a white background, and are mutually distinct.  Tweak freely.
GEO = {
    "active":    "#D81159",  # rose / magenta
    "inactive":  "#FFBC42",  # marigold
    "embedding": "#6A4C93",  # royal purple
    "mean":      "#0B3954",  # deep petrol (Karcher-mean marker edge)
}

# Line style for the active / inactive / embedding traces, shared across the
# sphere, normal-coordinate, and shadow figures.  ***Adjust here to restyle
# every figure consistently*** (color, dash pattern, width, white-halo width).
TRACE = {
    "active":    dict(color=GEO["active"],    ls="--",  lw=3, halo=2.25),
    "inactive":  dict(color=GEO["inactive"],  ls="--", lw=2.5, halo=1.75),
    "embedding": dict(color=GEO["embedding"], ls=":",  lw=3, halo=2.25),
}

# The embedding (extrinsic) trace is drawn only when it is meaningfully distinct
# from the intrinsic active trace, i.e. subspace distance d(U1, W) > this tol.
# Ridge / non-ridge cases are coincident (d ~ 1e-3..1e-2) and hide it; the
# degenerate aligned case (d ~ 0.66) shows it.  Adjust to taste.
EMBED_SHOW_TOL = 0.05


def trace2d(ax, x, y, key, z=6, halo_color="w"):
    """Plot a styled trace (TRACE[key]) with a white halo so it reads over a
    coloured field or a dense scatter.  Used by the 2-D shadow/normal figures."""
    st = TRACE[key]
    ax.plot(x, y, color=halo_color, lw=2.25, solid_capstyle="round",
            zorder=z, alpha=0)
    ax.plot(x, y, color=st["color"], ls=st["ls"], lw=2,
            solid_capstyle="round", zorder=z + 1)


# --------------------------------------------------------------------------- #
def _normal_to_sphere(p0, E, tnml):
    """Map normal coordinates ``tnml`` (M, 2) to sphere points via the true
    normal-coordinate map Exp_{p0}(E t)."""
    return amg.sphere_exp_batch(p0, tnml @ E.T)


def figure_normal_coords(func, res, R=0.7, ngrid=170, path=None):
    """Contour of the function in normal coordinates over a *geodesic ball* of
    radius R at p0, with the projected Riemannian-gradient quiver and the
    intrinsic active/inactive axes.  Uses the true normal-coordinate map
    Exp_{p0}(t) (not the MATLAB's Exp(1, .)), so the neighbourhood is a
    faithful geodesic ball -- the key picture for interpreting the metric."""
    p0, E = res.p0, res.E

    g = np.linspace(-R, R, ngrid)
    T1, T2 = np.meshgrid(g, g)
    RR = np.hypot(T1, T2)
    tnml = np.column_stack([T1.ravel(), T2.ravel()])
    F = func.f(_normal_to_sphere(p0, E, tnml)).reshape(ngrid, ngrid)
    F = np.where(RR <= R, F, np.nan)                 # restrict to the geodesic ball

    # projected tangential gradient on a coarse grid inside the ball
    gc = np.linspace(-R, R, 17)
    Q1, Q2 = np.meshgrid(gc, gc)
    tq = np.column_stack([Q1.ravel(), Q2.ravel()])
    tq = tq[np.hypot(tq[:, 0], tq[:, 1]) <= 0.97 * R]
    xq = _normal_to_sphere(p0, E, tq)
    pg = amg.project_tangent(func.grad(xq), xq) @ E

    fig, ax = plt.subplots(figsize=(5.7, 5.1))
    cs = ax.contourf(T1, T2, F, 60, cmap=CMAP)
    ax.contour(T1, T2, F, 14, colors="k", linewidths=0.35, alpha=0.4)
    ax.quiver(tq[:, 0], tq[:, 1], pg[:, 0], pg[:, 1], color="k",
              width=0.004, alpha=0.85)

    th = np.linspace(0, 2 * np.pi, 240)
    ax.plot(R * np.cos(th), R * np.sin(th), "k-", lw=1.4)  # geodesic-ball boundary

    u1 = E.T @ res.U_active; u1 = u1 / np.linalg.norm(u1) * R
    u2 = E.T @ res.U_inactive; u2 = u2 / np.linalg.norm(u2) * R
    trace2d(ax, [-u1[0], u1[0]], [-u1[1], u1[1]], "active", z=6)
    trace2d(ax, [-u2[0], u2[0]], [-u2[1], u2[1]], "inactive", z=6)

    ax.set_xlabel(r"$t_1$"); ax.set_ylabel(r"$t_2$")
    ax.set_aspect("equal")
    ax.set_xlim(-1.05 * R, 1.05 * R); ax.set_ylim(-1.05 * R, 1.05 * R)
    fig.colorbar(cs, ax=ax, shrink=0.85, label=r"$f\circ\exp_{p_0}$")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return fig


# --------------------------------------------------------------------------- #
def _geodesic(p0, direction, t):
    """Points along the geodesic Exp_{p0}(t * direction) for an array t."""
    return np.array([amg.sphere_exp(p0, direction, ti) for ti in t])


def figure_sphere_3d(func, res, n_show=45, arc=1.25, show_embedding=True, path=None):
    """The function on the sphere with samples, gradients, and geodesics."""
    p0, P, Gt = res.p0, res.P, res.Gt
    r = 1.01  # float overlays just above the surface

    def geo(curve, key, z):
        # white halo underneath so the curve reads over the coloured sphere
        st = TRACE[key]
        ax.plot(*curve.T, "-", color="w", lw=st["lw"] + st["halo"], zorder=z, alpha=0)
        ax.plot(*curve.T, color=st["color"], ls=st["ls"], lw=st["lw"], zorder=z + 1)

    fig = plt.figure(figsize=(6.6, 6.4))
    ax = fig.add_subplot(111, projection="3d")
    ax.computed_zorder = False  # respect our explicit zorder (curves over surface)

    # --- upper hemisphere: coloured by the function (the admissible domain) ---
    nu, nv = 72, 34
    u = np.linspace(0, 2 * np.pi, nu)
    v_up = np.linspace(0, np.pi / 2, nv)
    Xu = np.outer(np.cos(u), np.sin(v_up))
    Yu = np.outer(np.sin(u), np.sin(v_up))
    Zu = np.outer(np.ones_like(u), np.cos(v_up))
    Fs = func.f(np.column_stack([Xu.ravel(), Yu.ravel(), Zu.ravel()])).reshape(nu, nv)
    facec = CMAP((Fs - Fs.min()) / (np.ptp(Fs) + 1e-30))
    ax.plot_surface(Xu, Yu, Zu, facecolors=facec, rstride=1, cstride=1,
                    linewidth=0, antialiased=True, alpha=0.9, shade=False)

    # --- lower hemisphere: grid only (f is NOT evaluated here); kept apparent
    #     to emphasise the inadmissible region ---
    v_dn = np.linspace(np.pi / 2, np.pi, nv)
    Xd = np.outer(np.cos(u), np.sin(v_dn))
    Yd = np.outer(np.sin(u), np.sin(v_dn))
    Zd = np.outer(np.ones_like(u), np.cos(v_dn))
    ax.plot_wireframe(Xd, Yd, Zd, color="0.45", linewidth=0.5, alpha=0.6,
                      rstride=2, cstride=2)

    # --- equator: solid black boundary of the admissible domain ---
    eu = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(eu), np.sin(eu), np.zeros_like(eu), "k-", lw=2.2, zorder=4)

    # samples + tangential gradients (thinned; no marker edges, semi-transparent)
    idx = np.arange(min(n_show, P.shape[0]))
    ax.scatter(P[idx, 0] * r, P[idx, 1] * r, P[idx, 2] * r, c="k", s=8, alpha=0.4,
               edgecolors="none", depthshade=False)
    ax.quiver(P[idx, 0] * r, P[idx, 1] * r , P[idx, 2] * r,
              Gt[idx, 0], Gt[idx, 1], Gt[idx, 2],
              length=0.17, color="0.15", linewidth=2, normalize=False, alpha=0.55)

    # geodesics floated above the surface (arc < pi/2 keeps them on the visible cap)
    t = np.linspace(-arc, arc, 200)
    AMGc = _geodesic(p0, res.U_active, t) * r
    IAMGc = _geodesic(p0, res.U_inactive, t) * r
    EMBc = _geodesic(p0, res.W, t) * r
    geo(IAMGc, "inactive", 2)
    if show_embedding:
        geo(EMBc, "embedding", 7)
    geo(AMGc, "active", 1)

    # Karcher mean, on top
    ax.scatter(*(p0 * r), c="w", s=120, depthshade=False, edgecolor=GEO["mean"],
               linewidth=2.2, zorder=10)

    # view from the perpendicular bisector of the active/inactive azimuths so
    # both geodesics read as oblique arcs crossing at p0 (not edge-on); a low
    # elevation exposes the southern (inadmissible) hemisphere
    az = np.degrees(np.arctan2(res.U_active[1], res.U_active[0])) + 45
    ax.view_init(elev=24, azim=az)
    ax.set_box_aspect((1, 1, 1), zoom=1.45)   # fill the frame
    _hide_3d_panes(ax)
    ax.set_position([0, 0, 1, 1])              # no surrounding white space
    if path:
        fig.savefig(path, dpi=150); plt.close(fig)
    return fig


def _hide_3d_panes(ax):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_visible(False)
        a.line.set_color((1, 1, 1, 0))


# --------------------------------------------------------------------------- #
def figure_shadow(func, res, nsweep=40, ngeo=200, show_embedding=True, path=None):
    """Shadow plot: f along a sweep of geodesics from active to inactive,
    with the active/inactive/embedding sweeps highlighted and the sample
    scatter over the active coordinate."""
    p0, P = res.p0, res.P
    u1, u2, w = res.U_active, res.U_inactive, res.W

    # active-coordinate of each sample: <Log_{p0}(P), u1>
    Gy = np.array([amg.sphere_log(p0, q) @ u1 for q in P])
    Fp = func.f(P)
    Tr = 1.2 * np.max(np.abs(Gy)) if Gy.size else 1.5
    tg = np.linspace(-Tr, Tr, ngeo)

    fig, ax = plt.subplots(figsize=(6.2, 4.6))

    # sweep directions rotating u1 -> u2 (normalised)
    for s in np.linspace(0, 1, nsweep):
        d = (1 - s) * u1 + s * u2
        d /= np.linalg.norm(d)
        shade = abs(d @ u1)               # 1 at active, ~0 at inactive
        col = str(float(np.clip(0.85 * (1 - shade), 0.0, 1.0)))  # grey; clamp fp noise
        vals = func.f(_geodesic(p0, d, tg))
        ax.plot(tg, vals, "-", color=col, lw=0.6, zorder=1)

    # scatter underneath, then the intrinsic traces on top (haloed) so they are
    # not obscured; embedding only when it differs from the active trace
    sc = ax.scatter(Gy, Fp, c=Fp, cmap=CMAP, s=26, edgecolors="none",
                    alpha=0.5, zorder=2)
    if show_embedding:
        trace2d(ax, tg, func.f(_geodesic(p0, w, tg)), "embedding", z=6)
    trace2d(ax, tg, func.f(_geodesic(p0, u2, tg)), "inactive", z=8)
    trace2d(ax, tg, func.f(_geodesic(p0, u1, tg)), "active", z=10)

    ax.set_xlabel(r"$t$"); ax.set_ylabel(r"$(f\circ\exp_{p_0})(t\,w)$")
    fig.colorbar(sc, ax=ax, shrink=0.85, label=r"$f$")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return fig


# --------------------------------------------------------------------------- #
def figure_convergence(func, p0, N_values, nboot=20, path=None):
    """Subspace-distance convergence of the intrinsic active direction vs N."""
    mean_err, std_err = amg.convergence_study(func, p0, N_values, nboot=nboot)

    # least-squares rate on log10
    A = np.column_stack([np.ones_like(N_values, float), np.log10(N_values)])
    coef, *_ = np.linalg.lstsq(A, np.log10(mean_err), rcond=None)
    rate = coef[1]

    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    ax.fill_between(N_values, np.clip(mean_err - 3 * std_err, 1e-16, None),
                    mean_err + 3 * std_err, color="0.5", alpha=0.25)
    ax.loglog(N_values, 10 ** (coef[0]) * N_values ** rate, "--",
              color=GEO["active"], lw=1.4)
    ax.loglog(N_values, mean_err, "o-", color="k", lw=2, ms=6)
    ax.set_xlabel(r"$N$"); ax.set_ylabel("subspace distance")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return fig, rate


# --------------------------------------------------------------------------- #
def figure_ridge_recovery(func, res, R_max=0.9, npts=9000, deg=6, seed=1, path=None):
    """Demonstrate that an ambient ridge is recovered as a ridge *in normal
    coordinates* over geodesic balls.

    Left:  f collapses onto the single active normal coordinate
           s1 = <t, w1> -- points at the same s1 but different inactive
           coordinate s2 share the same value (a curve, not a cloud) for a
           ridge; a genuine 2-D function spreads into a band.
    Right: the unexplained variance  1 - R^2_ridge  (fraction of f-variance a
           1-D fit in s1 misses) as a function of the geodesic-ball radius R.
           For a ridge this vanishes as R -> 0 (curvature-limited, ~R^2);
           for a non-ridge it plateaus at O(1)."""
    p0, E = res.p0, res.E
    w1 = E.T @ res.U_active; w1 /= np.linalg.norm(w1)
    w2 = E.T @ res.U_inactive; w2 /= np.linalg.norm(w2)

    rng = np.random.default_rng(seed)
    rr = R_max * np.sqrt(rng.random(npts))          # uniform over the disk
    th = 2 * np.pi * rng.random(npts)
    tnml = np.column_stack([rr * np.cos(th), rr * np.sin(th)])
    s1, s2 = tnml @ w1, tnml @ w2
    fval = func.f(_normal_to_sphere(p0, E, tnml))

    def residual_fraction(mask):
        x, y = s1[mask], fval[mask]
        if y.size < deg + 2 or np.var(y) < 1e-30:
            return np.nan
        c = np.polyfit(x, y, deg)
        return np.var(y - np.polyval(c, x)) / np.var(y)

    Rs = np.linspace(0.12, R_max, 16)
    resid = np.array([residual_fraction(rr <= R) for R in Rs])

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.4, 4.5))

    order = np.argsort(np.abs(s2))[::-1]            # large |s2| behind
    sc = axL.scatter(s1[order], fval[order], c=s2[order], cmap=CMAP, s=10,
                     edgecolors="none", alpha=0.55)
    axL.set_xlabel(r"active coordinate  $s_1=\langle t,\,w_1\rangle$")
    axL.set_ylabel(r"$f\circ\exp_{p_0}$")
    fig.colorbar(sc, ax=axL, shrink=0.85, label=r"inactive  $s_2=\langle t,\,w_2\rangle$")

    resid_plot = np.clip(resid, 1e-16, None)
    axR.loglog(Rs, resid_plot, "o-", color=GEO["active"], lw=2, ms=6)
    axR.set_xlabel(r"geodesic-ball radius  $R$")
    axR.set_ylabel(r"unexplained variance  $1-R^2_{\mathrm{ridge}}$")
    axR.grid(True, which="both", alpha=0.3)
    good = np.isfinite(resid) & (resid > 1e-13)
    slope = None
    if good.sum() >= 2:
        c = np.polyfit(np.log10(Rs[good]), np.log10(resid[good]), 1)
        axR.plot(Rs, 10 ** c[1] * Rs ** c[0], "k--", lw=1)   # power-law guide
        slope = c[0]
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return fig, resid, slope


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--func", default="linear_random",
                    help="ambient function name (see sphere_amg.make_function)")
    ap.add_argument("--N", type=int, default=3000, help="samples for the field figures")
    ap.add_argument("--seed", type=int, default=47)
    ap.add_argument("--nboot", type=int, default=20)
    ap.add_argument("--no-convergence", action="store_true")
    ap.add_argument("--cmap", default="lapaz",
                    help="colormap: Crameri (lapaz, batlow, vik, davos, ...) or matplotlib")
    args = ap.parse_args()

    global CMAP
    CMAP = get_cmap(args.cmap)

    os.makedirs(FIGDIR, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    func = amg.make_function(args.func, seed=args.seed)
    P = amg.sample_disk(args.N, rng=rng)
    p0 = amg.karcher_mean(P)
    res = amg.compute_amg(func, p0, P)

    print(f"function            : {func.name} -- {func.description}")
    print(f"Karcher mean p0     : {np.round(p0, 4)}")
    print(f"intrinsic eigenvalues: {np.round(res.lam, 5)}")
    print(f"extrinsic eigenvalues: {np.round(res.lam_emb, 5)}")
    print(f"active  U1          : {np.round(res.U_active, 4)}")
    print(f"embedding W         : {np.round(res.W, 4)}")

    tag = func.name
    d_uw = amg.subspace_distance(res.U_active, res.W)
    show_emb = d_uw > EMBED_SHOW_TOL   # hide the coincident embedding trace
    show_emb = False # TURN THIS OFF
    print(f"d(U1, W)            : {d_uw:.4f}  -> embedding trace "
          f"{'shown' if show_emb else 'hidden (coincident with active)'}")
    figure_normal_coords(func, res, path=os.path.join(FIGDIR, f"01_normal_coords_{tag}.png"))
    figure_sphere_3d(func, res, show_embedding=show_emb,
                     path=os.path.join(FIGDIR, f"02_sphere3d_{tag}.png"))
    figure_shadow(func, res, show_embedding=show_emb,
                  path=os.path.join(FIGDIR, f"03_shadow_{tag}.png"))
    _, _, rr_slope = figure_ridge_recovery(
        func, res, path=os.path.join(FIGDIR, f"05_ridge_recovery_{tag}.png"))
    if rr_slope is not None:
        print(f"ridge-recovery slope: R^{rr_slope:.2f}")
    # Convergence needs a non-degenerate reference: the ridge direction must
    # have a non-vanishing tangential component at p0 (the aligned case is
    # normal at the pole -> skip it, and say so).
    ref_ok = func.a is not None and np.linalg.norm(
        amg.project_tangent(func.a[None, :], p0[None, :])) > 1e-6
    if not args.no_convergence and ref_ok:
        Nvals = 2 ** np.arange(2, 12)   # 4 .. 2048
        _, rate = figure_convergence(func, p0, Nvals, nboot=args.nboot,
                                     path=os.path.join(FIGDIR, f"04_convergence_{tag}.png"))
        print(f"convergence rate    : N^{rate:.2f}")
    elif not args.no_convergence:
        print("convergence         : skipped (no non-degenerate reference direction)")

    print(f"figures written to  : {FIGDIR}")


if __name__ == "__main__":
    main()
