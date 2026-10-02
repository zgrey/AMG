"""Can a non-ridge with NO critical point near the centre show R^2 decay *and* persistent scatter?

Supports the v1.5 Fig. 6 caption discussion (Applied_AMG, notes/RESEARCH_NOTES.md, Round 3,
Step 4 item 2).  Uses the figure pipeline itself (sphere_amg; same hemisphere samples and
seeds as sphere_amg_demo.figure_ridge_recovery).  Test function: the quadratic of Fig. 6
tilted so that its centre is no longer critical,

    f(x) = x^T diag(1,2,3) x + t * b.x,     b a tangent direction at e3.

For each tilt t it reports the fitted log-log slope of the RMS deviation from the
active-geodesic profile over R in [0.1, 1], and the relative deviation (RMS / std of f over
the ball) at R = 0.1 and R = 1, with the AMG direction w1 either fixed from the hemisphere
(as in the figure) or recomputed on each ball.

Finding (2026-10-01): R^2 decay together with persistent relative scatter requires a
(near-)critical point within the plotted ball.  In the local quadratic model the relative
deviation is ~ |H| R / (|g| + |H| R), which stays order one only while R >~ |g| / |H|; with the
figure's fixed hemisphere direction, a non-critical centre adds a first-order misalignment
term that pulls the slope toward 1.

Run:  ~/venvs/tda-sst/Scripts/python check_noncritical_nonridge.py
"""
import numpy as np

import sphere_amg as amg

rng = np.random.default_rng(47)
P = amg.sample_disk(3000, rng=rng)
p0 = amg.karcher_mean(P)

H = np.diag([1.0, 2.0, 3.0])
b = np.array([1.0, 0.6, 0.0]); b /= np.linalg.norm(b)          # tangent tilt at e3


def tilted(t):
    return amg.AmbientFunction(f"quad+{t}b", lambda X: np.sum((X @ H) * X, 1) + t * X @ b,
                               lambda X: 2.0 * X @ H + t * b, None, "tilted quadratic")


def study(func):
    res = amg.compute_amg(func, p0, P)                           # hemisphere AMG (figure)
    E = res.E
    g0 = amg.project_tangent(func.grad(p0[None, :]), p0[None, :]).ravel()
    r2 = np.random.default_rng(1)                                # figure's seed
    rr = np.sqrt(r2.random(9000)); th = 2 * np.pi * r2.random(9000)
    tn = np.column_stack([rr * np.cos(th), rr * np.sin(th)])
    X = amg.sphere_exp_batch(p0, tn @ E.T)
    fval = func.f(X)
    Rg = np.linspace(0.1, 1.0, 16)
    out = {}
    for mode in ("fixed", "local"):
        rms, rel = [], []
        for R in Rg:
            m = rr <= R
            u1 = res.U_active if mode == "fixed" else amg.compute_amg(func, p0, X[m]).U_active
            w1 = E.T @ u1; w1 /= np.linalg.norm(w1)
            s1 = tn[m] @ w1
            dev = fval[m] - func.f(amg.sphere_exp_batch(p0, s1[:, None] * u1[None, :]))
            rms.append(np.sqrt(np.mean(dev ** 2))); rel.append(rms[-1] / np.std(fval[m]))
        slope = np.polyfit(np.log(Rg), np.log(rms), 1)[0]
        out[mode] = (slope, rel[0], rel[-1])
    return np.linalg.norm(g0), out


def main():
    cases = [(f"quadratic + {t} tilt", tilted(t)) for t in (0.0, 0.3, 1.0, 3.0)]
    cases += [(nm, amg.make_function(nm, seed=47))
              for nm in ("linear_random", "nonlinear_ridge", "nonlinear_nonridge")]
    print(f"{'case':24s} {'|grad f(p0)|':>12s} | {'fixed w1: slope, rel@0.1, rel@1':>34s} | "
          f"{'local w1: slope, rel@0.1, rel@1':>34s}")
    for name, fn in cases:
        gn, o = study(fn)
        fs = "{:6.2f}  {:7.3f}  {:7.3f}".format(*o["fixed"])
        ls = "{:6.2f}  {:7.3f}  {:7.3f}".format(*o["local"])
        print(f"{name:24s} {gn:12.3f} | {fs:>34s} | {ls:>34s}")


if __name__ == "__main__":
    main()
