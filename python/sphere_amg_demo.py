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

# --------------------------------------------------------------------------- #
# Universal font sizing.  ***Bump FONT_SIZE to enlarge the labels, ticks,
# annotations, and colorbars across every figure at once.***  Tinker freely.
FONT_SIZE = 15
plt.rcParams.update({
    "font.size":        FONT_SIZE,
    "axes.titlesize":   FONT_SIZE,
    "axes.labelsize":   FONT_SIZE,
    "xtick.labelsize":  FONT_SIZE - 1,
    "ytick.labelsize":  FONT_SIZE - 1,
    "legend.fontsize":  FONT_SIZE - 1,
    "figure.titlesize": FONT_SIZE + 1,
})


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
    "black": dict(color='k', ls=":", lw = 1),
    "inactive_dots": dict(color=GEO["inactive"],  ls=":", lw=2),
    "active_dots": dict(color=GEO["active"],    ls=":",  lw=2),
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


def _draw_normal_coords(ax, func, res, R=0.7, ngrid=170, quiver=True):
    """Draw the normal-coordinate contour, gradient quiver (tangential gradient
    parallel-transported to p0 -- the G0 integrand), and active/inactive axes
    into ``ax``.  Returns the contourf handle (for an optional colorbar)."""
    p0, E = res.p0, res.E
    g = np.linspace(-R, R, ngrid)
    T1, T2 = np.meshgrid(g, g)
    RR = np.hypot(T1, T2)
    tnml = np.column_stack([T1.ravel(), T2.ravel()])
    F = func.f(_normal_to_sphere(p0, E, tnml)).reshape(ngrid, ngrid)
    F = np.where(RR <= R, F, np.nan)                 # restrict to the geodesic ball
    cs = ax.contourf(T1, T2, F, 60, cmap=CMAP)
    ax.contour(T1, T2, F, 14, colors="k", linewidths=0.35, alpha=0.4)
    if quiver:
        gc = np.linspace(-R, R, 17)
        Q1, Q2 = np.meshgrid(gc, gc)
        tq = np.column_stack([Q1.ravel(), Q2.ravel()])
        tq = tq[np.hypot(tq[:, 0], tq[:, 1]) <= 0.97 * R]
        xq = _normal_to_sphere(p0, E, tq)
        Gt_q = amg.project_tangent(func.grad(xq), xq)
        pg = np.array([amg.parallel_transport(xq[i], p0, Gt_q[i])
                       for i in range(len(xq))]) @ E
        ax.quiver(tq[:, 0], tq[:, 1], pg[:, 0], pg[:, 1], color="k",
                  width=0.004, alpha=0.85)
    th = np.linspace(0, 2 * np.pi, 240)
    ax.plot(R * np.cos(th), R * np.sin(th), "k-", lw=1.4)  # geodesic-ball boundary
    u1 = E.T @ res.U_active; u1 = u1 / np.linalg.norm(u1) * R
    u2 = E.T @ res.U_inactive; u2 = u2 / np.linalg.norm(u2) * R
    trace2d(ax, [-u1[0], u1[0]], [-u1[1], u1[1]], "active", z=6)
    trace2d(ax, [-u2[0], u2[0]], [-u2[1], u2[1]], "inactive", z=6)
    ax.set_aspect("equal")
    ax.set_xlim(-1.05 * R, 1.05 * R); ax.set_ylim(-1.05 * R, 1.05 * R)
    return cs


def figure_normal_coords(func, res, R=0.7, ngrid=170, path=None):
    """Standalone normal-coordinate figure (axis labels + colorbar)."""
    fig, ax = plt.subplots(figsize=(5.7, 5.1))
    cs = _draw_normal_coords(ax, func, res, R=R, ngrid=ngrid)
    ax.set_xlabel(r"$t_1$"); ax.set_ylabel(r"$t_2$")
    fig.colorbar(cs, ax=ax, shrink=0.85, label=r"$f\circ\exp_{p_0}$")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return fig


# --------------------------------------------------------------------------- #
def _geodesic(p0, direction, t):
    """Points along the geodesic Exp_{p0}(t * direction) for an array t."""
    return np.array([amg.sphere_exp(p0, direction, ti) for ti in t])

def _sphere_grad_grid(func, n=8, rad_max=1.0):
    """Sparse regular (theta, phi) grid of upper-hemisphere points and their
    tangential gradients.  The sampled domain is the spherical cap
    0 <= theta <= arcsin(rad_max); its outer ring sits at
    z = cos(arcsin(rad_max)) = sqrt(1 - rad_max^2), so ``rad_max`` near 1
    fills the cap out to the equator (z -> 0)."""
    if not 0.0 < rad_max <= 1.0:
        raise ValueError("rad_max must lie in (0, 1]; theta = arcsin(rad_max) "
                         "is undefined otherwise.")

    theta_max = np.arcsin(rad_max)
    theta = np.linspace(0.0, theta_max, n)                   # rings: pole -> boundary
    phi = np.linspace(0.0, 2.0 * np.pi, 2*n, endpoint=False)   # azimuths

    # theta[0] == 0 is the north pole: a single point for every phi.  Emit it
    # once, then take the full (theta > 0) x phi product for the rings so the
    # output carries no duplicate rows (mirrors the disk grid's uniqueness).
    pole = np.array([[0.0, 0.0, 1.0]])
    T, Ph = np.meshgrid(theta[1:], phi, indexing="ij")
    rings = np.column_stack([
        (np.sin(T) * np.cos(Ph)).ravel(),
        (np.sin(T) * np.sin(Ph)).ravel(),
        np.cos(T).ravel(),
    ])
    Pg = np.vstack([pole, rings])

    Pg /= np.linalg.norm(Pg, axis=1, keepdims=True)  # unit by construction; guards fp error
    Gg = amg.project_tangent(func.grad(Pg), Pg)
    return Pg, Gg

def _view_direction(elev, azim):
    """Unit vector from the origin toward the camera for a matplotlib 3-D view
    (used to split artists into near / far for hand-managed occlusion)."""
    e, a = np.radians(elev), np.radians(azim)
    return np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])


def _draw_sphere_3d(ax, func, res, arc=1.25, n_grid=7, n_sample=12,
                    show_embedding=True, seed=0, zoom=1.45):
    """Draw the sphere, sparse gradient grid, ringed Monte-Carlo subset, and
    geodesics into a 3-D ``ax`` (behaviour matches the standalone figure)."""
    p0 = res.p0
    r = 1.01  # float the curves/markers just above the surface
    elev = 24
    az = np.degrees(np.arctan2(res.U_active[1], res.U_active[0])) + 45
    view = _view_direction(elev, az)
    ax.computed_zorder = False  # respect explicit zorder (hand-managed occlusion)

    def geo(curve, key, z):
        st = TRACE[key]
        ax.plot(*curve.T, "-", color="w", lw=st["lw"] + st["halo"], zorder=z, alpha=0)
        ax.plot(*curve.T, color=st["color"], ls=st["ls"], lw=st["lw"], zorder=z + 1)

    nu, nv = 72, 34
    u = np.linspace(0, 2 * np.pi, nu)

    # lower hemisphere: grid only (inadmissible region), behind everything
    v_dn = np.linspace(np.pi / 2, np.pi, nv)
    ax.plot_wireframe(np.outer(np.cos(u), np.sin(v_dn)),
                      np.outer(np.sin(u), np.sin(v_dn)),
                      np.outer(np.ones_like(u), np.cos(v_dn)),
                      color="0.45", linewidth=0.5, alpha=0.6, rstride=2, cstride=2,
                      zorder=0)

    Pg, Gg = _sphere_grad_grid(func, n=n_grid)
    near = (Pg @ view) > 0
    rng = np.random.default_rng(seed)
    sub = np.zeros(len(Pg), bool)
    sub[rng.choice(len(Pg), size=min(n_sample, len(Pg)), replace=False)] = True

    def draw_grads(mask, z):
        if not np.any(mask):
            return
        ax.quiver(Pg[mask, 0], Pg[mask, 1], Pg[mask, 2],
                  Gg[mask, 0], Gg[mask, 1], Gg[mask, 2],
                  length=0.20, color="0.12", linewidth=1.4, normalize=False,
                  alpha=0.9, zorder=z)
        ring = mask & sub                      # Monte Carlo sample: ring the feet
        if np.any(ring):
            ax.scatter(Pg[ring, 0], Pg[ring, 1], Pg[ring, 2], s=150,
                       facecolors="none", edgecolors="k", linewidths=2.0,
                       depthshade=False, zorder=z + 0.5)

    v_up = np.linspace(0, np.pi / 2, nv)
    Xu = np.outer(np.cos(u), np.sin(v_up))
    Yu = np.outer(np.sin(u), np.sin(v_up))
    Zu = np.outer(np.ones_like(u), np.cos(v_up))
    Fs = func.f(np.column_stack([Xu.ravel(), Yu.ravel(), Zu.ravel()])).reshape(nu, nv)
    facec = CMAP((Fs - Fs.min()) / (np.ptp(Fs) + 1e-30))
    ax.plot_surface(Xu, Yu, Zu, facecolors=facec, rstride=1, cstride=1,
                    linewidth=0, antialiased=True, alpha=1.0, shade=False, zorder=2)

    draw_grads(near, z=7)                       # near side (over the surface)

    eu = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(eu), np.sin(eu), np.zeros_like(eu), "k-", lw=3, zorder=1)

    t = np.linspace(-arc * 0.20, arc, 200)
    geo(_geodesic(p0, res.U_active, t) * r, "active", 3)
    geo(_geodesic(p0, res.U_inactive, -t) * r, "inactive", 3)
    #if show_embedding:
    #    geo(_geodesic(p0, res.W, t) * r, "embedding", 3)

    ax.scatter(*(p0 * r), c="w", s=120, depthshade=False, edgecolor=GEO["mean"],
               linewidth=2.2, zorder=12)

    ax.view_init(elev=elev, azim=az)
    ax.set_box_aspect((1, 1, 1), zoom=zoom)
    _hide_3d_panes(ax)


def figure_sphere_3d(func, res, arc=1.25, n_grid=7, n_sample=12,
                     show_embedding=True, seed=0, path=None):
    """Standalone 3-D sphere figure (fills the frame)."""
    fig = plt.figure(figsize=(6.6, 6.4))
    ax = fig.add_subplot(111, projection="3d")
    _draw_sphere_3d(ax, func, res, arc=arc, n_grid=n_grid, n_sample=n_sample,
                    show_embedding=show_embedding, seed=seed)
    ax.set_position([0, 0, 1, 1])
    if path:
        fig.savefig(path, dpi=150); plt.close(fig)
    return fig


def _hide_3d_panes(ax):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_visible(False)
        a.line.set_color((1, 1, 1, 0))


# --------------------------------------------------------------------------- #
def _draw_shadow(ax, func, res, nsweep=40, ngeo=200, show_embedding=True):
    """Draw the shadow sweeps, sample scatter, and haloed intrinsic traces into
    ``ax``.  Returns the scatter handle (for an optional colorbar)."""
    p0, P = res.p0, res.P
    u1, u2, w = res.U_active, res.U_inactive, res.W
    Gy = np.array([amg.sphere_log(p0, q) @ u1 for q in P])   # <Log_{p0}(P), u1>
    Fp = func.f(P)
    Tr = 1.2 * np.max(np.abs(Gy)) if Gy.size else 1.5
    tg = np.linspace(-Tr, Tr, ngeo)
    for s in np.linspace(0, 1, nsweep):                      # sweep u1 -> u2
        d = (1 - s) * u1 + s * u2
        d /= np.linalg.norm(d)
        shade = abs(d @ u1)
        col = str(float(np.clip(0.85 * (1 - shade), 0.0, 1.0)))
        ax.plot(tg, func.f(_geodesic(p0, d, tg)), "-", color=col, lw=0.6, zorder=1)
    sc = ax.scatter(Gy, Fp, c=Fp, cmap=CMAP, s=26, edgecolors="none",
                    alpha=0.5, zorder=2)
    if show_embedding:
        trace2d(ax, tg, func.f(_geodesic(p0, w, tg)), "embedding", z=6)
    trace2d(ax, tg, func.f(_geodesic(p0, u2, tg)), "inactive", z=8)
    trace2d(ax, tg, func.f(_geodesic(p0, u1, tg)), "active", z=10)
    return sc


def figure_shadow(func, res, nsweep=40, ngeo=200, show_embedding=True, path=None):
    """Standalone shadow figure (axis labels + colorbar)."""
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    sc = _draw_shadow(ax, func, res, nsweep=nsweep, ngeo=ngeo,
                      show_embedding=show_embedding)
    ax.set_xlabel(r"$t$"); ax.set_ylabel(r"$(f\circ\exp_{p_0})(t\,w)$")
    fig.colorbar(sc, ax=ax, shrink=0.85, label=r"$f$")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return fig


# --------------------------------------------------------------------------- #
PANEL_FUNCS = ["linear_random", "linear_aligned", "quadratic_ridge",
               "quadratic_pref", "nonlinear_ridge", "nonlinear_nonridge"]
PANEL_LABELS = {"linear_random": "linear ridge", "linear_aligned": "linear (aligned)",
                "quadratic_ridge": "quadratic ridge", "quadratic_pref": "quadratic",
                "nonlinear_ridge": "nonlinear ridge",
                "nonlinear_nonridge": "nonlinear non-ridge"}


def figure_panel(funcs=None, N=2000, seed=47, R=0.7, path=None):
    """Gallery of the toy examples: rows = functions, columns = normal
    coordinates / on the sphere / shadow plot.  Shared/minimal annotations
    (column headers + row labels only), no colorbars, tight layout."""
    funcs = PANEL_FUNCS if funcs is None else funcs
    heads = ["normal coordinates", "on the sphere", "shadow plot"]
    nrow = len(funcs)
    fig = plt.figure(figsize=(9.6, 3.05 * nrow))
    top = []
    for i, name in enumerate(funcs):
        fn = amg.make_function(name, seed=seed)
        rng = np.random.default_rng(seed)
        P = amg.sample_disk(N, rng=rng)
        p0 = amg.karcher_mean(P)
        res = amg.compute_amg(fn, p0, P)
        show_emb = amg.subspace_distance(res.U_active, res.W) > EMBED_SHOW_TOL

        ax0 = fig.add_subplot(nrow, 3, 3 * i + 1)
        _draw_normal_coords(ax0, fn, res, R=R)
        ax0.set_xticks([]); ax0.set_yticks([])
        ax0.set_ylabel(PANEL_LABELS.get(name, name))

        ax1 = fig.add_subplot(nrow, 3, 3 * i + 2, projection="3d")
        _draw_sphere_3d(ax1, fn, res, show_embedding=show_emb, zoom=1.7)

        ax2 = fig.add_subplot(nrow, 3, 3 * i + 3)
        _draw_shadow(ax2, fn, res, show_embedding=show_emb)
        ax2.set_xticks([]); ax2.set_yticks([])

        if i == 0:
            top = [ax0, ax1, ax2]

    fig.tight_layout(h_pad=0.4, w_pad=0.4)
    # column headers above the top row (fig.text avoids the 3-D title overlap)
    for ax, head in zip(top, heads):
        pos = ax.get_position()
        fig.text(pos.x0 + pos.width / 2, pos.y1 + 0.004, head,
                 ha="center", va="bottom", fontsize=FONT_SIZE)
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
def figure_ridge_recovery(func, res, R_max=1.0, npts=9000, seed=1, path=None):
    """Ridge recovery over a shrinking geodesic ball (two panels).

    Left: f against the active geodesic coordinate s1 = <exp^{-1}_{p0}(x), w1>,
    with the sample scatter in nested shells for decreasing ball radius R
    (lightening with R; dotted lines mark each ball's extent +/-R).  As R -> 0
    the extent contracts and the scatter collapses onto the 1-D active-geodesic
    sweep.  Right: the RMS deviation of the samples from that active-geodesic
    sweep as a function of R -- the convergence, curvature-limited as R -> 0.
    A genuine 2-D (non-ridge) function keeps a vertical spread / a nonzero RMS
    plateau even at small R."""
    p0, E = res.p0, res.E
    u1 = res.U_active
    u2 = res.U_inactive
    w1 = E.T @ u1; w1 /= np.linalg.norm(w1)

    rng = np.random.default_rng(seed)
    rr = R_max * np.sqrt(rng.random(npts))           # geodesic ball radius |log_{p0}|
    th = 2 * np.pi * rng.random(npts)
    tnml = np.column_stack([rr * np.cos(th), rr * np.sin(th)])
    s1 = tnml @ w1                                   # active geodesic coordinate
    fval = func.f(_normal_to_sphere(p0, E, tnml))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.4, 4.7))

    # --- left: shadow over nested shrinking balls ---
    Rs = [1.0, 0.6, 0.3, 0.15]                       # nested, decreasing
    greys = [0.80, 0.55, 0.35, 0.0]                 # light (large R) -> dark (small R)
    ytop = np.nanmin(fval)+0.75
    for R, g in zip(Rs, greys):
        m = rr <= R
        axL.scatter(s1[m], fval[m], s=7, color=str(g), edgecolors="none",
                    alpha=0.75, zorder=2)
        slocal = np.linspace(-R, R, 2)
        axL.scatter(slocal, func.f(_geodesic(p0, u2, slocal)), marker="o", color=str(g), alpha = 0.75, s=20)
        axL.scatter(slocal, func.f(_geodesic(p0, u1, slocal)), marker="o", color=str(g), alpha = 0.75, s=20)
        for sgn in (-1, 1):                          # extent of the R-ball in s1
            axL.axvline(sgn * R, color=str(g), lw=0.9, ls=":", zorder=2)
        axL.text(R, ytop, fr"$R\!=\!{R:g}$", color="0.25", fontsize=FONT_SIZE - 4,
                 ha="right", va="top", rotation=90)
    sgrid = np.linspace(-R_max, R_max, 300)
    trace2d(axL, sgrid, func.f(_geodesic(p0, u2, sgrid)), "inactive", z=0)
    trace2d(axL, sgrid, func.f(_geodesic(p0, u1, sgrid)), "active", z=1)
    axL.set_xlabel(r"active geodesic coordinate  $\langle\exp^{-1}_{p_0}(\hat{x}),\,w_1\rangle$")
    axL.set_ylabel(r"$f\circ\exp_{p_0}$")
    axL.set_xlim(-1.15 * R_max, 1.15 * R_max)

    # --- right: RMS deviation from the active-geodesic sweep vs ball radius ---
    f_on_curve = func.f(amg.sphere_exp_batch(p0, s1[:, None] * u1[None, :]))
    dev = fval - f_on_curve
    Rgrid = np.linspace(0.1, R_max, 16)
    rms = np.array([np.sqrt(np.mean(dev[rr <= R] ** 2)) if np.any(rr <= R) else np.nan
                    for R in Rgrid])
    axR.loglog(Rgrid, np.clip(rms, 1e-16, None), "o-", color="k", lw=2, ms=6, zorder=1)
    good = np.isfinite(rms) & (rms > 1e-13)
    if good.sum() >= 2:
        c = np.polyfit(np.log10(Rgrid[good]), np.log10(rms[good]), 1)
        axR.loglog(Rgrid, 10 ** c[1] * Rgrid ** c[0], "--", color=GEO["active"], lw=3)
        axR.text(0.05, 0.92, fr"$\propto R^{{{c[0]:.1f}}}$", transform=axR.transAxes,
                 color=GEO["active"], fontsize=FONT_SIZE - 1, va="top")
    axR.set_xlabel(r"geodesic-ball radius  $R$")
    axR.set_ylabel(r"RMS deviation from active geodesic")
    axR.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return fig


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
    ap.add_argument("--panel", action="store_true",
                    help="render the full example gallery (all toy functions) and exit")
    args = ap.parse_args()

    global CMAP
    CMAP = get_cmap(args.cmap)

    os.makedirs(FIGDIR, exist_ok=True)

    if args.panel:
        figure_panel(seed=args.seed, path=os.path.join(FIGDIR, "00_panel.png"))
        print(f"gallery panel written to {os.path.join(FIGDIR, '00_panel.png')}")
        return

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
    print(f"d(U1, W)            : {d_uw:.4f}  -> embedding trace "
          f"{'shown' if show_emb else 'hidden (coincident with active)'}")
    figure_normal_coords(func, res, path=os.path.join(FIGDIR, f"01_normal_coords_{tag}.png"))
    figure_sphere_3d(func, res, show_embedding=show_emb,
                     path=os.path.join(FIGDIR, f"02_sphere3d_{tag}.png"))
    figure_shadow(func, res, show_embedding=show_emb,
                  path=os.path.join(FIGDIR, f"03_shadow_{tag}.png"))
    figure_ridge_recovery(
        func, res, path=os.path.join(FIGDIR, f"05_ridge_recovery_{tag}.png"))
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
