"""Numerical verification of the local-agreement results of the AMG paper (arXiv v2).

Deterministic tensor quadrature on a tangent ball B_R(0) of T_{p0}M, so the intrinsic G0 and the
extrinsic central representation E0^T C_iota E0 are integrated against the *same* wrapped measure
(no Monte Carlo noise, and no shared-randomness shortcut: each matrix is assembled from its own
integrand at the quadrature nodes).  Closed-form sphere and cylinder geometry; no EVIE dependency.

Checks (labels follow RESEARCH_NOTES.md, v2 cycle, Round 0):

  A  Prop. 3 (P1).  ||Delta||_2 <= kappa^2 int d^2 |grad f|^2 dmu, with kappa the largest normal
     curvature; on the unit sphere also the sharper 2 int (1 - cos d) |grad f|^2 dmu.  Also the trace
     identity  tr G0 - tr E0^T C E0 = int |(I - pi0) dio[grad f]|^2 dmu.
     Run on S^2, S^3 (kappa = 1) and on a cylinder of radius rho (kappa = 1/rho, intrinsically flat).
  B  Trailing spectrum and r >= 2 (P2, the remark).  Fitted log-log rates in R of
     |lambda_i(E0^T C E0) - lambda_i(G0)|, lambda_i(G0) and the top-r subspace distance, against the
     Davis--Kahan / Yu--Wang--Samworth bounds, for a symmetric (uniform) and a skewed density.
  C  Jacobi comparison (P3, Lemma).  |grad F - P grad f| <= (sinh d / d - 1) |grad f| on the unit
     sphere (|sec| <= Lambda = 1), F = f o exp_{p0}.
  D  Manifold ridge bound (P3, Prop. 4) on S^2 with the indicator density, r = 1:
     ||f - h(w1^T log_{p0})||_{L^2(mu)}  <=  sqrt(C_P) [ lambda_2^{1/2} + beta(R) (lambda_1+lambda_2)^{1/2} ],
     with h the exact conditional mean on the inactive fibres and C_P = (2R/pi)^2 (Payne--Weinberger).

Run:  ~/venvs/tda-sst/Scripts/python verify_local_agreement.py  [--seeds 5]
"""
import argparse

import numpy as np

RS = np.array([0.4, 0.2, 0.1, 0.05, 0.025])


# --------------------------------------------------------------------------- quadrature
def ball_quadrature(n, R, nr=48, nang=72):
    """Nodes x (N, n) and weights w (N,) for Lebesgue measure on B_R(0) subset R^n, n in {2, 3}."""
    r, wr = np.polynomial.legendre.leggauss(nr)
    r = 0.5 * R * (r + 1.0)
    wr = 0.5 * R * wr
    if n == 2:
        th = 2 * np.pi * np.arange(nang) / nang
        dirs = np.stack([np.cos(th), np.sin(th)], 1)
        wd = np.full(nang, 2 * np.pi / nang)
    elif n == 3:
        c, wc = np.polynomial.legendre.leggauss(nang // 2)
        ph = 2 * np.pi * np.arange(nang) / nang
        s = np.sqrt(1 - c ** 2)
        dirs = np.stack([np.outer(s, np.cos(ph)), np.outer(s, np.sin(ph)),
                         np.outer(c, np.ones_like(ph))], -1).reshape(-1, 3)
        wd = np.outer(wc, np.full(nang, 2 * np.pi / nang)).ravel()
    else:
        raise ValueError(n)
    x = (r[:, None, None] * dirs[None, :, :]).reshape(-1, n)
    w = (wr[:, None] * r[:, None] ** (n - 1) * wd[None, :]).ravel()
    return x, w


def density(x, R, skew):
    """Uniform (symmetric under x -> -x) or skewed, rho ~ 1 + 0.9 x_1 / R (positive on the ball)."""
    return (1 + 0.9 * x[:, 0] / R) if skew else np.ones(len(x))


# --------------------------------------------------------------------------- sphere geometry
class Sphere:
    """Unit S^n in R^{n+1}, p0 = e_{n+1}, E0 = first n columns of I; ambient test function
    f(x) = b.x + 0.5 x^T A x (random b, A), so grad f = b + A x."""
    kappa = 1.0

    def __init__(self, n, rng):
        self.n, m = n, n + 1
        self.b = rng.standard_normal(m); self.b /= np.linalg.norm(self.b)
        A = rng.standard_normal((m, m)); self.A = 0.5 * (A + A.T)
        self.p0 = np.zeros(m); self.p0[-1] = 1.0
        self.E0 = np.eye(m)[:, :n]

    def f(self, X):
        return X @ self.b + 0.5 * np.sum((X @ self.A) * X, 1)

    def exp(self, x):
        d = np.linalg.norm(x, axis=1)
        xb = np.divide(x, d[:, None], out=np.zeros_like(x), where=d[:, None] > 0)
        e = xb @ self.E0.T
        return np.cos(d)[:, None] * self.p0 + np.sin(d)[:, None] * e

    def fields(self, x):
        """Central-coordinate transported gradient gt, projected gradient q, |grad f|, d, and the
        normal-to-center part |(I - pi0) ghat| for nodes x in T_{p0}."""
        d = np.linalg.norm(x, axis=1)
        xb = np.divide(x, d[:, None], out=np.zeros_like(x), where=d[:, None] > 0)
        e = xb @ self.E0.T
        xh = np.cos(d)[:, None] * self.p0 + np.sin(d)[:, None] * e
        u = -np.sin(d)[:, None] * self.p0 + np.cos(d)[:, None] * e      # radial dir at xh
        amb = self.b + xh @ self.A
        gh = amb - np.sum(amb * xh, 1)[:, None] * xh                       # tangential gradient
        a = np.sum(gh * u, 1)
        Pg = gh + a[:, None] * (e - u)                                     # transport to p0
        nrm_out = np.abs(gh @ self.p0)                                     # |(I - pi0) ghat|
        return Pg @ self.E0, gh @ self.E0, np.linalg.norm(gh, axis=1), d, xb, nrm_out


class Cylinder:
    """Cylinder of radius rho in R^3, p0 = (rho, 0, 0); intrinsic orthonormal coordinates
    (u, v) = (rho*theta, z) are normal coordinates (the cylinder is flat), transport is the
    identity in the frame (e_theta, e_z), kappa = 1/rho, injectivity radius pi*rho."""
    n = 2

    def __init__(self, rho, rng):
        self.rho, self.kappa = rho, 1.0 / rho
        self.c = rng.standard_normal(4)

    def grad_uv(self, x):
        u, v = x[:, 0], x[:, 1]
        c = self.c
        fu = c[0] + c[2] * np.cos(c[2] * u + c[3] * v)
        fv = c[1] + c[3] * np.cos(c[2] * u + c[3] * v)
        return np.stack([fu, fv], 1)

    def fields(self, x):
        g = self.grad_uv(x)                       # components in (e_theta(theta), e_z)
        th = x[:, 0] / self.rho
        gt = g                                    # transported: same components at p0
        q = np.stack([g[:, 0] * np.cos(th), g[:, 1]], 1)   # pi0 e_theta(theta) = cos(theta) e_theta(0)
        d = np.linalg.norm(x, axis=1)
        nrm_out = np.abs(g[:, 0] * np.sin(th))
        return gt, q, np.linalg.norm(g, axis=1), d, None, nrm_out


# --------------------------------------------------------------------------- core assembly
def assemble(geom, R, skew=False):
    x, w = ball_quadrature(geom.n, R)
    w = w * density(x, R, skew); w = w / w.sum()                          # probability measure mu
    gt, q, gn, d, xb, nrm_out = geom.fields(x)
    G = (gt * w[:, None]).T @ gt
    C = (q * w[:, None]).T @ q
    return dict(x=x, w=w, gt=gt, q=q, gn=gn, d=d, xb=xb, nrm_out=nrm_out, G=G, C=C)


def eig_desc(M):
    lam, V = np.linalg.eigh(M)
    return lam[::-1], V[:, ::-1]


def subspace_dist(V1, V2, r):
    return np.linalg.norm(V1[:, :r] @ V1[:, :r].T - V2[:, :r] @ V2[:, :r].T, 2)


def slope(R, y):
    ok = np.isfinite(y) & (y > 1e-15)
    return np.polyfit(np.log(R[ok]), np.log(y[ok]), 1)[0] if ok.sum() >= 2 else np.nan


# --------------------------------------------------------------------------- checks
def check_A(geom, label, skew=False):
    print(f"\n[A] Prop. 3 bound + trace identity -- {label}{' (skewed mu)' if skew else ''}, "
          f"kappa = {geom.kappa:g}")
    print("     R     ||Delta||   2int(1-cos d)|g|^2   kappa^2 int d^2|g|^2   ratio   trace-id err")
    worst = 0.0
    for R in RS:
        a = assemble(geom, R, skew)
        nD = np.linalg.norm(a["C"] - a["G"], 2)
        bnd = geom.kappa ** 2 * np.sum(a["w"] * a["d"] ** 2 * a["gn"] ** 2)
        sph = (np.sum(a["w"] * 2 * (1 - np.cos(a["d"])) * a["gn"] ** 2)
               if isinstance(geom, Sphere) else np.nan)
        tr_gap = np.trace(a["G"]) - np.trace(a["C"])
        tr_id = np.sum(a["w"] * a["nrm_out"] ** 2)
        worst = max(worst, nD / bnd)
        print(f"  {R:6.3f}  {nD:10.3e}  {sph:15.3e}  {bnd:20.3e}  {nD / bnd:7.3f}"
              f"   {abs(tr_gap - tr_id) / max(tr_id, 1e-300):9.1e}")
    status = "OK" if worst <= 1 + 1e-12 else "VIOLATED"
    print(f"     max ||Delta|| / bound = {worst:.3f}  -> {status}")
    return worst <= 1 + 1e-12


def check_B(n, seeds, skew):
    tag = "skewed" if skew else "uniform (symmetric)"
    rs = list(range(1, n)) if n > 2 else [1]
    print(f"\n[B] Trailing spectrum / r >= 2 on S^{n}, {tag} mu: fitted rates in R over R in "
          f"{RS.tolist()}  (median over {seeds} random f)")
    rows = []
    for s in range(seeds):
        geom = Sphere(n, np.random.default_rng(100 + s))
        dl = np.zeros((len(RS), n)); lg_ = np.zeros((len(RS), n))
        sT = np.zeros((len(RS), len(rs))); yws = np.zeros((len(RS), len(rs)))
        dk_ok = np.zeros((len(RS), len(rs)), bool); gnorm0 = None
        for k, R in enumerate(RS):
            a = assemble(geom, R, skew)
            lG, VG = eig_desc(a["G"]); lC, VC = eig_desc(a["C"])
            nD = np.linalg.norm(a["C"] - a["G"], 2)
            dl[k] = np.abs(lC - lG); lg_[k] = lG
            for j, r in enumerate(rs):
                sT[k, j] = subspace_dist(VG, VC, r)
                eta = lG[r - 1] - lG[r]
                yws[k, j] = 2 * np.sqrt(r) * nD / eta            # Yu--Wang--Samworth Thm 2
                dk_ok[k, j] = eta > nD                           # classical DK applicable
        row = [slope(RS, dl[:, i]) for i in range(n)] + [slope(RS, lg_[:, i]) for i in range(1, n)]
        row += [slope(RS, sT[:, j]) for j in range(len(rs))] + [slope(RS, yws[:, j]) for j in range(len(rs))]
        row += [dk_ok[:, j].mean() for j in range(len(rs))]
        rows.append(row)
    rows = np.array(rows)
    med = np.nanmedian(rows, 0); lo = np.nanmin(rows, 0); hi = np.nanmax(rows, 0)
    names = ([f"|dlam{i + 1}|" for i in range(n)] + [f"lam{i + 1}" for i in range(1, n)]
             + [f"sinT(r={r})" for r in rs] + [f"YWS bnd(r={r})" for r in rs]
             + [f"frac(eta_r>||D||, r={r})" for r in rs])
    for nm, m_, a_, b_ in zip(names, med, lo, hi):
        print(f"     {nm:24s} median {m_:6.2f}   [min {a_:6.2f}, max {b_:6.2f}]")


def one_minus_sinc(d):
    """1 - sin(d)/d without cancellation (series below d = 0.1)."""
    d2 = d * d
    small = d2 / 6 - d2 ** 2 / 120 + d2 ** 3 / 5040 - d2 ** 4 / 362880
    return np.where(d < 0.1, small, 1 - np.sin(d) / np.where(d > 0, d, 1))


def sinhc_minus_one(d):
    """sinh(d)/d - 1 without cancellation (series below d = 0.1)."""
    d2 = d * d
    small = d2 / 6 + d2 ** 2 / 120 + d2 ** 3 / 5040 + d2 ** 4 / 362880
    return np.where(d < 0.1, small, np.sinh(d) / np.where(d > 0, d, 1) - 1)


def check_C(n, seeds):
    print(f"\n[C] Jacobi comparison on S^{n}: max over nodes of |grad F - P grad f| / "
          f"((sinh d/d - 1)|grad f|)  (must be <= 1)")
    worst = 0.0
    for s in range(seeds):
        geom = Sphere(n, np.random.default_rng(200 + s))
        for R in RS:
            a = assemble(geom, R)
            xb, gt, d = a["xb"], a["gt"], a["d"]
            rad = np.sum(gt * xb, 1)[:, None] * xb
            # On the unit sphere d exp_{p0} is exact radially and scales the transverse part by
            # sin(d)/d, so grad F - P grad f = -(1 - sin d / d) (transverse part of P grad f).
            lhs = one_minus_sinc(d) * np.linalg.norm(gt - rad, axis=1)
            rhs = sinhc_minus_one(d) * a["gn"]
            ok = rhs > 1e-14
            worst = max(worst, np.max(lhs[ok] / rhs[ok]))
    print(f"     max ratio = {worst:.4f}  -> {'OK' if worst <= 1 + 1e-9 else 'VIOLATED'}")
    return worst <= 1 + 1e-9


def check_D(seeds, ny=96, nz=64):
    print("\n[D] Prop. 4 manifold ridge bound on S^2, indicator density, r = 1 "
          "(error = exact conditional-mean error)")
    print("     f        R     error      bound     ratio   rel.err = error / std(f)")
    worst = 0.0
    cases = [(f"rand{s}", Sphere(2, np.random.default_rng(300 + s))) for s in range(seeds)]
    crit = Sphere(2, np.random.default_rng(0))                    # critical-point centre:
    crit.b = np.zeros(3); crit.A = np.diag([2.0, 4.0, 6.0])       # f = x^T diag(1,2,3) x at p0 = e3
    cases.append(("crit", crit))
    for name, geom in cases:
        for R in RS[:4]:
            a = assemble(geom, R)                                  # uniform density
            lG, VG = eig_desc(a["G"])
            w1, w2 = VG[:, 0], VG[:, 1]
            y, wy = np.polynomial.legendre.leggauss(ny)
            y = R * y; wy = R * wy
            zt, wz = np.polynomial.legendre.leggauss(nz)
            half = np.sqrt(np.maximum(R ** 2 - y ** 2, 0))
            Z = half[:, None] * zt[None, :]; WZ = half[:, None] * wz[None, :]
            X = y[:, None, None] * w1 + Z[:, :, None] * w2          # (ny, nz, 2) tangent nodes
            F = geom.f(geom.exp(X.reshape(-1, 2))).reshape(ny, nz)
            h = np.sum(F * WZ, 1) / np.sum(WZ, 1)                   # conditional mean on each fibre
            mass = np.pi * R ** 2
            err = np.sqrt(np.sum(wy[:, None] * WZ * (F - h[:, None]) ** 2) / mass)
            fbar = np.sum(wy[:, None] * WZ * F) / mass
            std = np.sqrt(np.sum(wy[:, None] * WZ * (F - fbar) ** 2) / mass)
            beta = float(sinhc_minus_one(np.array(R)))              # Lambda = 1
            bound = (2 * R / np.pi) * (np.sqrt(max(lG[1], 0)) + beta * np.sqrt(lG.sum()))
            worst = max(worst, err / bound)
            print(f"     {name:6s} {R:6.3f}  {err:9.3e}  {bound:9.3e}  {err / bound:6.3f}   {err / std:8.3e}")
    print(f"     max error / bound = {worst:.3f}  -> {'OK' if worst <= 1 + 1e-9 else 'VIOLATED'}")
    return worst <= 1 + 1e-9


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, default=5, help="random test functions per study")
    args = ap.parse_args()
    ok = True
    ok &= check_A(Sphere(2, np.random.default_rng(7)), "S^2")
    ok &= check_A(Sphere(3, np.random.default_rng(8)), "S^3", skew=True)
    ok &= check_A(Cylinder(0.5, np.random.default_rng(9)), "cylinder rho = 0.5 (flat)")
    for n in (2, 3):
        for skew in (False, True):
            check_B(n, args.seeds, skew)
    ok &= check_C(3, args.seeds)
    ok &= check_D(args.seeds)
    print(f"\nAll inequality checks {'PASSED' if ok else 'FAILED'}.")


if __name__ == "__main__":
    main()
