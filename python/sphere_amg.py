"""
Active Manifold-Geodesics (AMG) on the 2-sphere.

Python port of ``emb_mfld_AMG.m`` (Z. Grey).  The manifold geometry is reused
from EVIE (``evie.grassmannian``), treating the sphere as ``S^2 = Gr(1, 3)``:
a point is a unit vector (a line through the origin), tangent vectors are its
orthogonal complement, and the log/exp/parallel-transport/Karcher-mean routines
are the generic Grassmannian ones specialised to ``p = 1``.  This was verified
to reproduce the closed-form sphere maps to machine precision for base points
and targets within the injectivity radius (angle < pi/2 -- the Gr(1, 3) cut
locus, which is also the regime the AMG construction assumes).

The AMG construction (see "A Riemannian View on Active Subspaces"):
    * sample a neighbourhood around a central point p0 (the Karcher mean);
    * evaluate the ambient gradient and project it to each tangent space
      (the Riemannian gradient of the restricted function);
    * parallel-transport every gradient back to the central tangent space
      T_{p0} along the connecting geodesic;
    * the eigen-decomposition of the resulting average outer product G0 gives
      the intrinsic active (largest eigenvalue) and inactive directions, whose
      exponentials are the active / inactive manifold-geodesics.

The extrinsic (embedding) comparison follows Thm. "normal_eigvals": the
ambient second-moment C_iota projected onto T_{p0}.  Its top eigenvector is the
*extrinsic* active direction; where the ambient gradient leaks into the normal
space it is conflated by this projection -- the phenomenon the sphere example
is designed to expose.

Author: Zach Grey (MATLAB); Python port: Claude Code.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

from evie.grassmannian import (
    Gr_exp,
    Gr_log,
    Gr_intrinsicmean,
    Gr_parallel_trans,
)

# --------------------------------------------------------------------------- #
# Sphere (Gr(1, 3)) primitives -- thin wrappers over EVIE's Grassmannian ops.
# Points and vectors are plain length-3 arrays; we reshape to (3, 1) Stiefel
# representatives at the EVIE boundary and ravel back.
# --------------------------------------------------------------------------- #

def _col(v: np.ndarray) -> np.ndarray:
    return np.asarray(v, dtype=float).reshape(3, 1)


def sphere_exp(base: np.ndarray, direction: np.ndarray, t: float = 1.0) -> np.ndarray:
    """Exp_base(t * direction) on S^2.  ``direction`` need not be unit; the
    geodesic distance travelled is ``t * ||direction||``."""
    return Gr_exp(t, _col(base), _col(direction)).ravel()


def sphere_exp_batch(base: np.ndarray, dirs: np.ndarray) -> np.ndarray:
    """Vectorised Exp_base over many tangent directions (rows of ``dirs``).

    The ``p = 1`` specialisation of :func:`evie.grassmannian.Gr_exp`, used for
    dense grid / geodesic evaluation in the visualisations where per-point
    EVIE calls would be needlessly slow.  Equals ``Gr_exp(1, base, d)`` row by
    row (verified to machine precision); a zero direction maps to ``base``."""
    base = np.asarray(base, float)
    dirs = np.atleast_2d(np.asarray(dirs, float))
    nrm = np.linalg.norm(dirs, axis=1, keepdims=True)
    unit = np.divide(dirs, nrm, out=np.zeros_like(dirs), where=nrm > 1e-15)
    return np.cos(nrm) * base + np.sin(nrm) * unit


def sphere_log(base: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Log_base(target): tangent vector at ``base`` pointing to ``target``
    (norm = geodesic distance).  ``align=False`` keeps the genuine S^2 branch
    (no antipodal identification), valid while the angle is < pi/2."""
    H, _ = Gr_log(_col(base), _col(target), align=False)
    return H.ravel()


def parallel_transport(base_from: np.ndarray, base_to: np.ndarray,
                       vec: np.ndarray) -> np.ndarray:
    """Parallel-transport ``vec`` (tangent at ``base_from``) to ``base_to``
    along the connecting geodesic (Edelman Thm. 2.4, via EVIE)."""
    H = sphere_log(base_from, base_to)
    return Gr_parallel_trans(1.0, _col(base_from), _col(H), _col(vec)).ravel()


def karcher_mean(P: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Intrinsic (Karcher) mean of sphere points ``P`` (N, 3), oriented to the
    upper hemisphere (z > 0) to fix the Gr(1, 3) sign ambiguity."""
    X = np.ascontiguousarray(P.T).reshape(3, 1, -1)
    x0 = Gr_intrinsicmean(X, tol=tol, verbose=False).ravel()
    s = np.sign(x0[2]) if abs(x0[2]) > 1e-12 else 1.0
    return x0 * s


def tangent_basis(p0: np.ndarray) -> np.ndarray:
    """Orthonormal (3, 2) basis of the tangent space at ``p0`` (columns span
    {v : p0 . v = 0}).  This is the intended object behind the MATLAB
    ``E = eye(3) + p0 p0'`` eigen-trick, computed here directly and
    unambiguously as the orthogonal complement of ``p0``."""
    # SVD of the (1, 3) row p0^T: right-singular vectors 2 and 3 span the null.
    _, _, Vt = np.linalg.svd(p0.reshape(1, 3))
    return Vt[1:].T  # (3, 2)


def project_tangent(G: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Project ambient vectors ``G`` (N, 3) onto the tangent spaces at points
    ``P`` (N, 3): the Riemannian gradient of the restricted function."""
    return G - np.sum(G * P, axis=1, keepdims=True) * P


# --------------------------------------------------------------------------- #
# Ambient function registry.  Each function is defined on R^3 (evaluated on the
# sphere) and carries an *analytic* gradient; correctness is checked against
# finite differences in the __main__ self-test, so we never rely on a stale
# hand-derived expression.
# --------------------------------------------------------------------------- #

@dataclass
class AmbientFunction:
    name: str
    f: Callable[[np.ndarray], np.ndarray]        # (N, 3) -> (N,)
    grad: Callable[[np.ndarray], np.ndarray]     # (N, 3) -> (N, 3)
    a: Optional[np.ndarray] = None               # ridge direction, when defined
    description: str = ""


def make_function(name: str, seed: int = 47) -> AmbientFunction:
    """Build one of the ambient test functions from ``emb_mfld_AMG.m``.

    Names
    -----
    ``linear_random``   : f = a.x, a a random unit vector -- a clean 1-D active
                          direction (the MATLAB's first, commented option).
    ``linear_aligned``  : f = z (a = e_3).  Degenerate: the ambient gradient is
                          *normal* at the north-pole mean, so the tangent metric
                          is rotationally symmetric (no distinct active
                          direction) -- the extrinsic-conflation illustration.
    ``quadratic_pref``  : f = x^T diag(1..3) x -- preferential quadratic.
    ``quadratic_ridge`` : f = x^T H x with a rank-1 H block -- a quadratic ridge.
    ``nonlinear_ridge`` : f = sin(2*pi a.x) + cos(pi/2 a.x) -- 1-D but curved.
    ``nonlinear_nonridge`` : sinusoids across an orthonormal frame -- no
                          low-dimensional active subspace (the "space-filling"
                          cautionary case).
    """
    rng = np.random.default_rng(seed)
    m = 3

    if name == "linear_random":
        a = rng.uniform(-1, 1, m); a /= np.linalg.norm(a)
        return AmbientFunction(
            name, lambda X: X @ a, lambda X: np.tile(a, (X.shape[0], 1)), a,
            "linear f = a . x, generic a")

    if name == "linear_aligned":
        a = np.array([0.0, 0.0, 1.0])
        return AmbientFunction(
            name, lambda X: X @ a, lambda X: np.tile(a, (X.shape[0], 1)), a,
            "linear f = z; gradient normal at the pole (degenerate)")

    if name == "quadratic_pref":
        H = np.diag(np.linspace(1.0, m, m))
        return AmbientFunction(
            name, lambda X: np.sum((X @ H) * X, axis=1),
            lambda X: 2.0 * X @ H, None, "quadratic f = x^T diag(1..3) x")

    if name == "quadratic_ridge":
        H = np.zeros((m, m)); H[0, 0] = 1.0
        return AmbientFunction(
            name, lambda X: np.sum((X @ H) * X, axis=1),
            lambda X: 2.0 * X @ H, None, "rank-1 quadratic ridge f = (a.x)^2")

    if name == "nonlinear_ridge":
        a = rng.uniform(-1, 1, m); a /= np.linalg.norm(a)
        def f(X):
            s = X @ a
            return np.sin(2 * np.pi * s) + np.cos(np.pi / 2 * s)
        def g(X):
            s = X @ a
            gp = 2 * np.pi * np.cos(2 * np.pi * s) - (np.pi / 2) * np.sin(np.pi / 2 * s)
            return gp[:, None] * a
        return AmbientFunction(name, f, g, a, "ridge g(a.x), curved level sets")

    if name == "nonlinear_nonridge":
        A, _ = np.linalg.qr(rng.uniform(-1, 1, (m, m)))  # orthonormal frame
        freqs = np.array([1.0, 1.7, 2.3]) * np.pi
        def f(X):
            S = X @ A
            return np.sum(np.sin(freqs * S), axis=1)
        def g(X):
            S = X @ A
            dS = (freqs * np.cos(freqs * S))          # (N, 3) wrt frame coords
            return dS @ A.T                            # chain rule back to R^3
        return AmbientFunction(name, f, g, None,
                               "multi-directional sinusoid: no low-D active subspace")

    raise ValueError(f"unknown function name: {name!r}")


# --------------------------------------------------------------------------- #
# Core AMG computation.
# --------------------------------------------------------------------------- #

@dataclass
class AMGResult:
    p0: np.ndarray            # central point (Karcher mean), (3,)
    E: np.ndarray             # tangent basis at p0, (3, 2)
    P: np.ndarray             # samples on the sphere, (N, 3)
    Gt: np.ndarray            # tangential (Riemannian) gradients at samples, (N, 3)
    Vlog: np.ndarray          # gradients parallel-transported to T_{p0}, (N, 3)
    U_active: np.ndarray      # intrinsic active direction at p0, (3,)
    U_inactive: np.ndarray    # intrinsic inactive direction at p0, (3,)
    lam: np.ndarray           # intrinsic eigenvalues of G0 (mean-sq dir. deriv.), (2,)
    W: np.ndarray             # extrinsic (Mukherjee) direction, projected to T_{p0}, (3,)
    lam_emb: np.ndarray       # ambient eigenvalues of the extrinsic GOP C_iota
    W_proj_norm: float = 1.0  # survival of the projection of the dominant DR dir. (~0 => conflated)


def compute_amg(func: AmbientFunction, p0: np.ndarray, P: np.ndarray) -> AMGResult:
    """Intrinsic AMG at ``p0`` from samples ``P`` (N, 3), plus the extrinsic
    (embedding) comparison via the projected ambient second moment."""
    N = P.shape[0]
    E = tangent_basis(p0)

    # Riemannian gradients at the samples (ambient grad projected to tangent).
    G = func.grad(P)
    Gt = project_tangent(G, P)

    # Parallel-transport each gradient back to the central tangent space.
    Vlog = np.array([parallel_transport(P[i], p0, Gt[i]) for i in range(N)])

    # Intrinsic G0 = (1/N) sum vlog vlog^T; eigenvectors via SVD of (1/sqrt N) Vlog^T.
    U, S, _ = np.linalg.svd(Vlog.T / np.sqrt(N), full_matrices=False)
    U_active, U_inactive = U[:, 0], U[:, 1]
    lam = S[:2] ** 2  # eigenvalues of G0 = squared singular values

    # Extrinsic (Mukherjee) comparison.  C_iota is the ambient p x p gradient
    # outer product of the *tangential* gradients -- eq. emb_opg, i.e. the
    # Wu-Mukherjee manifold GOP E[dphi(grad_M f) (x) dphi(grad_M f)].  Its
    # dominant eigenvector is an *ambient* direction (the DR direction); the
    # extrinsic estimate at p0 is its projection to the central tangent space
    # (eigendecompose-then-project).  The projection is the partial isometry of
    # Thm normal_eigvals: where the DR direction is normal at p0 it collapses.
    C_iota = (Gt.T @ Gt) / N                        # ambient 3x3 GOP
    val, vec = np.linalg.eigh(C_iota)
    order = np.argsort(val)[::-1]
    lam_emb = val[order]                            # ambient eigenvalues
    a_hat = vec[:, order[0]]                        # dominant ambient DR direction
    proj = (np.eye(3) - np.outer(p0, p0)) @ a_hat   # project onto T_{p0}
    W_proj_norm = float(np.linalg.norm(proj))       # survival (~0 => conflation)
    W = proj / W_proj_norm if W_proj_norm > 1e-12 else proj

    # Orient the active direction for display, then sign-align W to it (order
    # matters: flipping U_active after aligning W would un-align them).
    if U_active[np.argmax(np.abs(U_active))] < 0:
        U_active = -U_active
    if U_active @ W < 0:
        W = -W

    return AMGResult(p0=p0, E=E, P=P, Gt=Gt, Vlog=Vlog,
                     U_active=U_active, U_inactive=U_inactive, lam=lam,
                     W=W, lam_emb=lam_emb, W_proj_norm=W_proj_norm)


# --------------------------------------------------------------------------- #
# Neighbourhood sampling on the upper hemisphere (matching the MATLAB disk).
# --------------------------------------------------------------------------- #

def sample_disk(N: int, rad_max: float = 0.99, theta_max: float = 2 * np.pi,
                rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """Sample ``N`` points on the upper hemisphere by drawing (r, theta) in the
    parametrisation disk and lifting through z = sqrt(1 - x^2 - y^2)."""
    rng = np.random.default_rng() if rng is None else rng
    r = rad_max * rng.random(N)
    th = theta_max * rng.random(N)
    x, y = r * np.cos(th), r * np.sin(th)
    z = np.sqrt(np.clip(1 - x**2 - y**2, 0.0, 1.0))
    P = np.column_stack([x, y, z])
    return P / np.linalg.norm(P, axis=1, keepdims=True)


# --------------------------------------------------------------------------- #
# Subspace-distance convergence study (bootstrap over sample size N).
# --------------------------------------------------------------------------- #

def subspace_distance(u: np.ndarray, v: np.ndarray) -> float:
    """Distance between the 1-D subspaces spanned by unit vectors u, v:
    ||u u^T - v v^T||_2."""
    return np.linalg.norm(np.outer(u, u) - np.outer(v, v), 2)


def convergence_study(func: AmbientFunction, p0: np.ndarray,
                      N_values: np.ndarray, nboot: int = 20,
                      rad_max: float = 0.99, seed: int = 0):
    """Mean subspace distance between the estimated intrinsic active direction
    and the analytic reference (tangential projection of the ridge direction
    ``a`` at ``p0``), as a function of sample size ``N``.  Returns
    ``(mean_err, std_err)`` each of shape ``(len(N_values),)``."""
    if func.a is None:
        raise ValueError("convergence reference needs a ridge direction `a`")
    ref = project_tangent(func.a[None, :], p0[None, :]).ravel()
    if np.linalg.norm(ref) < 1e-10:
        raise ValueError("reference direction is normal at p0 (degenerate case)")
    ref /= np.linalg.norm(ref)

    err = np.zeros((len(N_values), nboot))
    for j in range(nboot):
        rng = np.random.default_rng(seed + j)
        for i, N in enumerate(N_values):
            P = sample_disk(int(N), rad_max=rad_max, rng=rng)
            res = compute_amg(func, p0, P)
            err[i, j] = subspace_distance(res.U_active, ref)
    return err.mean(axis=1), err.std(axis=1) / np.sqrt(nboot)


# --------------------------------------------------------------------------- #
# Self-test: analytic gradients vs finite differences, and end-to-end sanity.
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    rng = np.random.default_rng(0)

    print("== analytic gradient vs finite difference (ambient R^3) ==")
    Xt = rng.standard_normal((200, 3))
    for nm in ["linear_random", "linear_aligned", "quadratic_pref",
               "quadratic_ridge", "nonlinear_ridge", "nonlinear_nonridge"]:
        fn = make_function(nm)
        g_an = fn.grad(Xt)
        eps = 1e-6
        g_fd = np.zeros_like(Xt)
        for k in range(3):
            e = np.zeros(3); e[k] = eps
            g_fd[:, k] = (fn.f(Xt + e) - fn.f(Xt - e)) / (2 * eps)
        rel = np.linalg.norm(g_an - g_fd) / (np.linalg.norm(g_fd) + 1e-30)
        print(f"  {nm:20s} rel.grad.err = {rel:.2e}")

    print("\n== end-to-end AMG on the sphere (linear_random) ==")
    fn = make_function("linear_random")
    P = sample_disk(2000, rng=rng)
    p0 = karcher_mean(P)
    res = compute_amg(fn, p0, P)
    ref = project_tangent(fn.a[None, :], p0[None, :]).ravel()
    ref /= np.linalg.norm(ref)
    print(f"  p0 (Karcher mean)     = {np.round(p0, 4)}")
    print(f"  intrinsic eigenvalues = {np.round(res.lam, 5)}")
    print(f"  active direction U1   = {np.round(res.U_active, 4)}")
    print(f"  analytic reference    = {np.round(ref, 4)}")
    print(f"  subspace dist (U1,ref)= {subspace_distance(res.U_active, ref):.3e}")
