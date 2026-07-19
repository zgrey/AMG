"""
Two-manifold schematic for the TREM -> AMG motivation slide (preso-SST deck).

A single composite figure making the load-bearing point of the talk: the TREM
measurement produces TWO manifold-valued data streams.

  (left)  spatiotemporal emission PATTERN  ->  segmented boundary curve
          ->  Separable-Shape-Tensor representative on the Grassmannian Gr(2,n).
  (right) univariate photon-count SIGNAL    ->  trend-filtered functional data
          ->  length-normalized preshape on a hypersphere S^{n-1}.

Everything shown is a SCHEMATIC / illustrative cartoon (synthetic data), meant
to set up "the domain of analysis is a manifold, not R^n" before the geometry.

Palette matches sphere_amg_demo.GEO.  Run with the tda-sst venv:
    ~/venvs/tda-sst/Scripts/python trem_manifold_schematic.py

Author: Zach Grey / Claude Code.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")

GEO = {
    "active":   "#D81159",
    "inactive": "#FFBC42",  # marigold (segmented curves / signals)
    "embed":    "#6A4C93",
    "mean":     "#0B3954",
    "pattern":  "#123C5A",   # deep blue "zeros" background of a binarized panel
}
plt.rcParams.update({"font.size": 13})


def arrow(ax, x0, x1, y=0.5, text=None):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>",
                                 mutation_scale=18, lw=2.0, color="0.35"))
    if text:
        ax.text((x0 + x1) / 2, y + 0.06, text, ha="center", va="bottom",
                fontsize=10, color="0.35")


def draw_pattern(ax):
    """Cartoon binarized space-time TREM panel with a segmented boundary."""
    rng = np.random.default_rng(3)
    H, W = 60, 44
    img = np.zeros((H, W))
    # a slanted "characteristic" band (ones) inside a marigold boundary
    xx, yy = np.meshgrid(np.linspace(0, 1, W), np.linspace(0, 1, H))
    band = np.abs(yy - (0.25 + 0.5 * xx)) < (0.12 + 0.03 * np.sin(6 * xx))
    img[band] = 1.0
    ax.imshow(img, cmap="cividis", origin="lower", extent=[0, 1, 0, 1],
              aspect="auto", vmin=0, vmax=1)
    # marigold segmented boundary (the closed curve)
    theta = np.linspace(0, 2 * np.pi, 200)
    cx, cy = 0.5, 0.5
    r = 0.34 + 0.05 * np.sin(3 * theta) + 0.03 * np.cos(5 * theta)
    bx = cx + r * np.cos(theta) * 0.62
    by = cy + r * np.sin(theta) * 0.9
    ax.plot(bx, by, color=GEO["inactive"], lw=2.6)
    ax.set_xlabel("space (along line)")
    ax.set_ylabel("time")
    ax.set_title("emission pattern")
    ax.set_xticks([]); ax.set_yticks([])


def draw_curve(ax):
    """The extracted closed boundary curve, standardized."""
    theta = np.linspace(0, 2 * np.pi, 200)
    r = 1.0 + 0.16 * np.sin(3 * theta) + 0.09 * np.cos(5 * theta)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    ax.plot(x, y, color=GEO["inactive"], lw=2.8)
    ax.fill(x, y, color=GEO["inactive"], alpha=0.12)
    ax.scatter(x[::14], y[::14], color=GEO["mean"], s=14, zorder=3)
    ax.set_aspect("equal")
    ax.set_title("boundary curve")
    ax.axis("off")


def draw_grassmann(ax):
    """Cartoon Gr(2,n): a plane (2-frame) through the origin in R^3."""
    ax.set_title(r"$Gr(2,n)$")
    # draw a tilted plane patch
    gx, gy = np.meshgrid(np.linspace(-1, 1, 2), np.linspace(-1, 1, 2))
    gz = 0.35 * gx - 0.2 * gy
    ax.plot_surface(gx, gy, gz, color=GEO["active"], alpha=0.30,
                    linewidth=0, shade=False)
    # two frame vectors
    for vv, c in [((1, 0, 0.35), GEO["active"]), ((0, 1, -0.2), GEO["mean"])]:
        vv = np.array(vv, float); vv = vv / np.linalg.norm(vv)
        ax.quiver(0, 0, 0, *vv, color=c, lw=2.2, arrow_length_ratio=0.15)
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.grid(False)
    for p in (ax.xaxis, ax.yaxis, ax.zaxis):
        p.pane.set_visible(False)
    ax.view_init(elev=20, azim=-60)


def draw_signal(ax):
    """Noisy photon-count trace + trend-filtered signal overlay."""
    rng = np.random.default_rng(11)
    t = np.linspace(0, 1, 300)
    clean = 0.5 * np.clip(4.0 * (t - 0.28), 0, None) + \
        0.5 * (1 - np.exp(-np.clip(t - 0.28, 0, None) * 3))
    noisy = clean + 0.06 * rng.standard_normal(t.size)
    ax.plot(t, noisy, color="0.7", lw=0.9, label="raw counts")
    ax.plot(t, clean, color=GEO["inactive"], lw=2.6, label="trend-filtered")
    ax.set_title("time signal")
    ax.set_xlabel("time"); ax.set_yticks([])
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper left")


def draw_hypersphere(ax):
    """Cartoon preshape hypersphere with a couple of points + geodesic."""
    ax.set_title(r"$S^{n-1}$")
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 30)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(xs, ys, zs, color="0.9", alpha=0.35,
                    rstride=2, cstride=2, linewidth=0, shade=False)
    a = np.array([0.6, 0.3, 0.74]); a /= np.linalg.norm(a)
    b = np.array([-0.2, 0.7, 0.68]); b /= np.linalg.norm(b)
    om = np.arccos(np.clip(a @ b, -1, 1))
    ts = np.linspace(0, 1, 40)
    arc = (np.sin((1 - ts)[:, None] * om) * a + np.sin(ts[:, None] * om) * b) / np.sin(om)
    ax.plot(arc[:, 0], arc[:, 1], arc[:, 2], color=GEO["active"], lw=2.2)
    for p, c in [(a, GEO["inactive"]), (b, GEO["mean"])]:
        ax.scatter(*(p * 1.02), color=c, s=60, depthshade=False, edgecolor="k", lw=0.4)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.grid(False)
    for p in (ax.xaxis, ax.yaxis, ax.zaxis):
        p.pane.set_visible(False)
    ax.view_init(elev=22, azim=35)


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    fig = plt.figure(figsize=(13, 5.2))
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 0.08], hspace=0.05, wspace=0.45)

    # top row of 6 mini-panels: pattern | curve | Gr  ||  signal | (arrow) | S
    axP = fig.add_subplot(gs[0, 0]); draw_pattern(axP)
    axC = fig.add_subplot(gs[0, 1]); draw_curve(axC)
    axG = fig.add_subplot(gs[0, 2], projection="3d"); draw_grassmann(axG)
    axS = fig.add_subplot(gs[0, 3]); draw_signal(axS)
    axH = fig.add_subplot(gs[0, 5], projection="3d"); draw_hypersphere(axH)

    # arrows between stages (use a thin overlay axis spanning the figure)
    ov = fig.add_axes([0, 0, 1, 1]); ov.axis("off"); ov.set_xlim(0, 1); ov.set_ylim(0, 1)
    for x0, x1 in [(0.175, 0.205), (0.34, 0.37)]:
        ov.add_patch(FancyArrowPatch((x0, 0.55), (x1, 0.55), arrowstyle="-|>",
                     mutation_scale=16, lw=1.8, color="0.4"))
    ov.add_patch(FancyArrowPatch((0.665, 0.55), (0.70, 0.55), arrowstyle="-|>",
                 mutation_scale=16, lw=1.8, color="0.4"))

    # section labels
    ov.text(0.25, 0.96, "PATTERNS  (shape manifold)", ha="center",
            fontsize=13, weight="bold", color=GEO["mean"])
    ov.text(0.75, 0.96, "SIGNALS  (functional data)", ha="center",
            fontsize=13, weight="bold", color=GEO["mean"])
    ov.axvline(0.51, ymin=0.12, ymax=0.9, color="0.8", lw=1.2, ls="--")
    ov.text(0.5, 0.02, "TREM measurement  ->  two manifold-valued domains  "
            "(illustrative schematic)", ha="center", fontsize=10,
            style="italic", color="0.5")

    out = os.path.join(FIGDIR, "trem_manifold_schematic.pdf")
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
