"""
Illustrative TREM line-signature signals as preshapes on a hypersphere.

Supplemental visualization for the AMG-TREM ITL seminar deck (preso-SST).
Several synthetic trend-filtered "characteristic curves" (piecewise-linear
photon-count-vs-time signatures with distinct slopes/delays, one per inverter
line) are length-normalized into unit-norm *preshapes* living on S^{n-1}.  We
render the preshape cloud on a 2-sphere via its three leading coordinates and
overlay a Frechet mean and geodesic arcs -- echoing the sphere figures in the
"A Riemannian View on Active Subspaces" paper.

The signal data are SYNTHETIC and purely illustrative of the geometry; they are
NOT measured TREM data.  Swap in real trend-filtered signatures when available.

Palette matches sphere_amg_demo.GEO (marigold = the TREM characteristic-curve
color).  Run with the tda-sst venv:
    ~/venvs/tda-sst/Scripts/python signals_hypersphere.py

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

# palette (kept in sync with sphere_amg_demo.GEO)
GEO = {
    "active":   "#D81159",
    "inactive": "#FFBC42",  # marigold
    "embed":    "#6A4C93",
    "mean":     "#0B3954",
}
plt.rcParams.update({"font.size": 14, "axes.labelsize": 14})


def make_signals(n_lines: int = 6, n: int = 200, seed: int = 7):
    """Synthetic piecewise-linear line signatures f_k(t) on t in [0,1].

    Each line has a distinct onset delay and rise slope (fast lines rise early
    and steeply); a mild plateau follows.  Returns (t, F) with F[k] the k-th
    signal sampled at n points.
    """
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, 1.0, n)
    F = np.zeros((n_lines, n))
    delays = np.linspace(0.10, 0.45, n_lines)      # onset (transport delay)
    slopes = np.linspace(6.0, 2.2, n_lines)        # rise rate (fast -> slow)
    for k in range(n_lines):
        d, s = delays[k], slopes[k]
        rise = np.clip(s * (t - d), 0.0, None)
        plateau = 1.0 - np.exp(-np.clip(t - d, 0.0, None) * 3.0)
        f = 0.5 * rise + 0.5 * plateau
        f = f + 0.01 * rng.standard_normal(n)      # faint measurement noise
        F[k] = f
    return t, F, delays


def to_preshape(F: np.ndarray) -> np.ndarray:
    """Center (remove translation) and scale (unit norm) each row -> S^{n-1}."""
    Fc = F - F.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(Fc, axis=1, keepdims=True)
    return Fc / norms


def leading_three(P: np.ndarray):
    """Project the preshape cloud onto its 3 leading principal axes and
    re-normalize to the unit 2-sphere (an honest low-dim view of S^{n-1})."""
    mu = P.mean(axis=0)
    U, S, Vt = np.linalg.svd(P - mu, full_matrices=False)
    coords = (P - mu) @ Vt[:3].T
    coords = coords / np.linalg.norm(coords, axis=1, keepdims=True)
    return coords


def slerp(a: np.ndarray, b: np.ndarray, m: int = 60) -> np.ndarray:
    """Geodesic (great-circle) arc between two unit vectors."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    om = np.arccos(np.clip(a @ b, -1.0, 1.0))
    if om < 1e-8:
        return np.tile(a, (m, 1))
    ts = np.linspace(0.0, 1.0, m)
    return (np.sin((1 - ts)[:, None] * om) * a + np.sin(ts[:, None] * om) * b) / np.sin(om)


def frechet_mean_sphere(X: np.ndarray, iters: int = 100) -> np.ndarray:
    """Intrinsic (Karcher) mean of unit vectors on S^2."""
    mu = X.mean(axis=0)
    mu = mu / np.linalg.norm(mu)
    for _ in range(iters):
        # log at mu
        d = X @ mu
        d = np.clip(d, -1.0, 1.0)
        ang = np.arccos(d)
        tang = X - d[:, None] * mu
        nrm = np.linalg.norm(tang, axis=1, keepdims=True)
        nrm[nrm < 1e-12] = 1.0
        logs = (ang[:, None]) * tang / nrm
        v = logs.mean(axis=0)
        nv = np.linalg.norm(v)
        if nv < 1e-10:
            break
        mu = np.cos(nv) * mu + np.sin(nv) * v / nv
        mu = mu / np.linalg.norm(mu)
    return mu


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    n_lines = 6
    t, F, delays = make_signals(n_lines=n_lines)
    P = to_preshape(F)
    X = leading_three(P)
    mu = frechet_mean_sphere(X)

    # color lines by delay (fast=rose, slow=marigold) via a 2-color blend
    from matplotlib.colors import to_rgb
    c0, c1 = np.array(to_rgb(GEO["active"])), np.array(to_rgb(GEO["inactive"]))
    fr = (delays - delays.min()) / (np.ptp(delays) + 1e-12)
    colors = [tuple((1 - a) * c0 + a * c1) for a in fr]

    fig = plt.figure(figsize=(11, 4.6))

    # ---- left: the trend-filtered signatures -------------------------------
    axL = fig.add_subplot(1, 2, 1)
    for k in range(n_lines):
        axL.plot(t, F[k], color=colors[k], lw=2.4,
                 label=f"line {k+1}" if k in (0, n_lines - 1) else None)
    axL.set_xlabel("time (normalized)")
    axL.set_ylabel("photon count (normalized)")
    axL.set_title("Trend-filtered line signatures")
    axL.spines[["top", "right"]].set_visible(False)
    axL.text(0.02, 0.95, "fast (low delay)", color=GEO["active"],
             transform=axL.transAxes, va="top", fontsize=11)
    axL.text(0.02, 0.86, "slow (high delay)", color=GEO["inactive"],
             transform=axL.transAxes, va="top", fontsize=11)
    axL.text(0.98, 0.03, "illustrative / synthetic", transform=axL.transAxes,
             ha="right", va="bottom", fontsize=9, style="italic", color="0.5")

    # ---- right: preshapes on the 2-sphere ----------------------------------
    axR = fig.add_subplot(1, 2, 2, projection="3d")
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(0, np.pi, 40)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    axR.plot_surface(xs, ys, zs, color="0.9", alpha=0.35,
                     rstride=2, cstride=2, linewidth=0, shade=False)
    # geodesic arcs from the mean to each preshape
    for k in range(n_lines):
        arc = slerp(mu, X[k])
        axR.plot(arc[:, 0], arc[:, 1], arc[:, 2], color=colors[k], lw=1.6, alpha=0.8)
        axR.scatter(*(X[k] * 1.02), color=colors[k], s=55, depthshade=False,
                    edgecolor="k", linewidth=0.4)
    axR.scatter(*(mu * 1.03), color="w", s=130, depthshade=False,
                edgecolor=GEO["mean"], linewidth=2.0, marker="o")
    axR.text(*(mu * 1.28), r"$\bar{p}$", color=GEO["mean"], fontsize=13)
    axR.set_title(r"Preshapes on $S^{n-1}$  (3 leading coords)")
    axR.set_box_aspect((1, 1, 1))
    axR.set_xticks([]); axR.set_yticks([]); axR.set_zticks([])
    axR.grid(False)
    for pane in (axR.xaxis, axR.yaxis, axR.zaxis):
        pane.pane.set_visible(False)
    axR.view_init(elev=22, azim=35)

    fig.tight_layout()
    out = os.path.join(FIGDIR, "signals_hypersphere.pdf")
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
