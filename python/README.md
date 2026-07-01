# AMG on the 2-sphere (Python)

Python / matplotlib port of [`../emb_mfld_AMG.m`](../emb_mfld_AMG.m) — Active
Manifold-Geodesics on the sphere, treated as `S^2 = Gr(1, 3)`.

Manifold geometry is **reused from EVIE** (`evie.grassmannian`, in the TDA-SST
repo): `Gr_exp`, `Gr_log`, `Gr_intrinsicmean`, and `Gr_parallel_trans` (the
Edelman Thm. 2.4 parallel transport, added to EVIE for this work). The p=1
Grassmannian maps were verified to reproduce the closed-form sphere maps to
machine precision within the injectivity radius (angle < π/2).

## Files

- **`sphere_amg.py`** — core. Sphere/Grassmannian ops (via EVIE), an ambient
  test-function registry (linear/quadratic/nonlinear ridge and non-ridge) with
  analytic gradients (finite-difference checked), Karcher mean, tangential
  (Riemannian) gradient, parallel transport of gradients to the central tangent
  space, intrinsic active/inactive directions (`G₀` eigendecomposition), the
  extrinsic/embedding comparison (Thm. `normal_eigvals`), and a
  subspace-distance convergence study.
- **`sphere_amg_demo.py`** — matplotlib figures: normal-coordinate geodesic-ball
  view, function on the sphere (upper hemisphere colored, southern grid, active/
  inactive/embedding geodesics), shadow plot, convergence, and the
  ridge-recovery demonstration.

## Run

Uses the `tda-sst` venv (numpy / scipy / matplotlib + evie):

```bash
~/venvs/tda-sst/Scripts/python sphere_amg_demo.py --func linear_random
```

Options: `--func {linear_random, linear_aligned, quadratic_pref,
quadratic_ridge, nonlinear_ridge, nonlinear_nonridge}`, `--cmap
{lapaz, batlow, vik, ...}` (Crameri via `cmcrameri`, or any matplotlib map),
`--nboot N`, `--no-convergence`. Figures are written to `figs/`.

Self-tests: `python sphere_amg.py` (analytic-gradient checks + end-to-end AMG).

## Styling

Colours and line styles are centralised at the top of `sphere_amg_demo.py`:
`CMAP` (colormap), `GEO` (geodesic accent palette), `TRACE` (per-trace
color/dash/width/halo, shared across figures), and `EMBED_SHOW_TOL` (the
embedding trace is drawn only when it is meaningfully distinct from the
intrinsic active direction).
