# Research Plan — Active Flows

**Working title:** *Active Flows: Localized Activity Fields and Their Integral Curves on
Riemannian Manifolds* (name continues the lineage: active subspaces → active
manifold-geodesics → active flows).

**Status:** plan drafted 2026-07-25 from the `python/integral_curves.ipynb` experiments
(this repo). To be broken out into its own Overleaf doc/repo (Phase 0); this file then
seeds that project's `RESEARCH_NOTES.md`. Code home stays **this repo** (AMG), venv
`tda-sst`, per the workspace convention (Overleaf repos are documents, never code).

---

## 1. Motivation and positioning

The AMG paper (Applied_AMG, arXiv pending) establishes the two-task framing:
**integrate** a gradient second moment and eigendecompose to identify activity, then
**extend** that activity over the manifold. Its extensions are *global*: one central
decomposition, extended by parallel transport (intrinsic) or — illegitimately — by
projected dominance. The integral-curves experiments show a third posture: **localize**
the intrinsic construction (re-wrapped measure `G_x^δ`, or kernel weights `G_x^h`) and
**flow** along the resulting dominant/inactive line fields. Empirically (nonlinear
non-ridge on S², scanned-center cap):

- the localized intrinsic flow is genuinely response-adapted (bends with the response,
  organizes into nested dominant loops; inactive family tracks response bands like
  approximate leaves);
- the trail system is robust to the localization mechanism (top-hat δ=0.8 vs Gaussian
  h=0.6 nearly coincide);
- **the localized intrinsic field visibly starts to recover the dominant-projection
  field** `top-eig(E_x^T C_ι E_x)` — the author's flagship observation, with a
  conjectured quadratic agreement rate for ambient ridge functions.

This is the response-driven analog of **principal flows** (Panaretos–Pham–Yao 2014,
which integrate a localized tangent-PCA direction) and mechanically identical to
**tractography** (integrating dominant eigenvector fields of diffusion tensors). Nobody
appears to have done it with *gradient* second moments; the AMG paper's agreement
theory gives it quantitative teeth the PCA version lacks.

## 2. Objects

Over Riemannian n-manifold (M, g), sufficiently smooth f, and point x ∈ M:

- **Re-wrapped local tensor** `G_x^δ` — Def. `riemannian_view` re-centered at x:
  wrapped uniform measure on the δ-ball of T_xM, gradients transported to x.
- **Kernel-localized tensor** `G_x^h` — fixed global samples, weights K(d(q,x)/h)
  (Gaussian etc.); top-hat K recovers a sampled `G_x^δ`.
- **Dominant / inactive line fields** — eigenvector fields of either tensor
  (sign-ambiguous: line fields, not vector fields).
- **Active flows** — integral curves of the dominant field (sign-continued);
  **inactive flows** — integral curves of the trailing field(s): leaf surrogates.
- **Singular set** Σ — points where the local spectral gap η(x) closes
  (line-field umbilics, index ±1/2).
- **Extrinsic comparators** — dominant-projection field `top-eig(E_x^T C_ι E_x)`
  (fixed global C_ι, per-block eigensolve) and projected-dominance field `π_x b̂`
  (the pencil through ±b̂; the "never do this" baseline).

## 3. Theory targets

**T1 (flagship mechanism) — compression is implicit localization.**
Claim (verified sketch on S², 2026-07-25, TO BE FORMALIZED): by Lemma
`sphere_proj_trans`, for q at geodesic distance d from x,
`E_x^T dι[∇f(q)] = P_{q→x}[∇f(q)] − (1−cos d)·(radial component)`,
i.e. the block integrand is the *transported* integrand with radial component
contracted by cos d. Hence
`E_x^T C_ι E_x = ∫ M_q[P∇f] ⊗ M_q[P∇f] dμ`, `M_q = I − (1−cos d) e_r e_r^T`:
the dominant-projection tensor **is** the intrinsic construction under an
anisotropic radial cos-d shrinkage kernel (exact on hyperspheres; general M to
second order via the second fundamental form). This explains the observed recovery:
localizing the intrinsic field with bandwidth matched to the cos²-effective kernel
approximates the dominant-projection field. Deliverable: a lemma ("compression as
anisotropic localization") + the matched-bandwidth corollary.

**T2 (author's conjecture) — quadratic localized-agreement rate.**
For ambient ridge functions restricted to ι(M) (then general f):
(i) *matched-support version*: angle(top-eig G_x^δ, top-eig E_x^T C_ι^δ E_x) =
O(δ²/η_x) — likely a pointwise application of Prop. `local_agreement` at center x
with support radius δ; check whether it is literally a corollary or needs the
kernel-weighted extension. (ii) *global-block version*: quantify
angle(top-eig G_x^δ, top-eig E_x^T C_ι E_x) as a function of δ via T1 — conjecture:
minimized near a cos²-matched δ*, with quadratic behavior in the mismatch.
Numerics first (E2), then proof.

**T3 — well-posedness of active flows.** Existence/uniqueness of integral curves off
Σ (Picard–Lindelöf for the unit line field where η > 0; regularity of eigenvector
fields of smooth tensor fields). Behavior at Σ: classification of generic umbilics
(index ±1/2, trisector/wedge patterns — tensor-field-topology literature), and the
Poincaré–Hopf constraint forcing singularities on caps/spheres. Flows must treat Σ
as structural (response umbilics), not failure.

**T4 — variational characterization.** Principal flows arise as stationary curves of
a localized-variance functional; formulate the response-driven analog (e.g. curves
maximizing ∫ λ₁(G_γ(t)^δ) dt subject to unit speed, or an EL equation whose solutions
are the active flows). Payoffs: existence through Σ, identifiability, and a principled
answer to `rem:dichotomy` — does second-moment regularization bound curve complexity
(no space-filling)? A complexity/length bound is REQUIRED before claiming "principled
nonlinear dimension reduction."

**T5 — statistics.** Monte Carlo / quadrature consistency of `G_x^δ` and `G_x^h`
fields (uniform-in-x rates on compacts off Σ); error propagation along traced curves
(Grönwall along the flow, gap-dependent); bandwidth selection (bias–variance;
matched-mass calibration h ↔ δ; plug-in or Lepski-style rules). Data-driven version
with learned gradients (Mukherjee/Wu OPG rates in intrinsic dimension).

**T6 — foliation link.** Inactive flows as leaf surrogates: if trailing eigenvalues
vanish identically on a neighborhood, inactive flows ARE leaves (Frobenius, as in the
paper); quantify leaf-approximation error for near-ridges (O(δ²)? tie to T2). The
n>2 case: inactive *distribution* integrability vs bracket obstructions — flows
generalize to foliations only when the Frobenius condition (approximately) holds.

**T7 — ridge tracing connection.** Active flows over ambient ridges should trace the
ridge structure — connect to the `ridge_tracing` project (the latent
ridge-functions ↔ active-subspaces link flagged in the workspace CLAUDE.md becomes
concrete here). Compare quad-trace/SRVF machinery to gap-colored active flows.

## 4. Literature review targets (verify primary sources; pull real theorem numbers)

- **Principal flows/curves on manifolds:** Panaretos, Pham & Yao (JASA 2014); Yao &
  collaborators' follow-ons (principal sub-manifolds, principal boundaries); Hastie &
  Stuetzle 1989 (principal curves — the Euclidean ancestor + its own ill-posedness
  fixes); Kenobi/Dryden principal nested spheres (adjacent); Fletcher PGA; Pennec
  barycentric subspaces (already cited in AMG paper).
- **Tractography / tensor-field topology:** Basser et al. (1994; 2000 in vivo fiber
  tractography) — dominant-eigenvector integration + gap (FA) stopping criteria;
  Delmarcelle & Hesselink (1994) tensor-field topology (umbilics, index ±1/2,
  trisector/wedge); vector/direction-field design literature (N-symmetry fields,
  Ray et al.) for line-field singularity handling.
- **Gradient-outer-product learners:** Mukherjee, Wu, Zhou (2010); Wu et al. (2010)
  (already load-bearing in the AMG paper); Trivedi–Wang kernel GOP variants; recent
  neural active manifolds (Zanoni 2025 — the fully nonlinear rival to position against).
- **Kernel localization on manifolds:** kernel density/regression on Riemannian
  manifolds (Pelletier; Henry–Rodriguez), bandwidth theory to steal for T5.
- **Ridge tracing / ridge sets:** Eberly's ridge definitions (height ridges),
  Ozertem & Erdogmus (2011) subspace-constrained mean shift (SCMS — VERY close
  mechanically: integrating projected gradient flows to find ridges; must be cited
  and differentiated: SCMS uses density Hessians, we use response-gradient second
  moments + transport), Genovese et al. ridge estimation theory.
- Verify each against the actual papers before citing (research-integrity rule 3).

## 5. Numerical experiments (all in this repo; tda-sst venv)

- **E1 — bandwidth study:** sweep δ and h (matched effective mass), field angles vs
  bandwidth at fixed probe points; the δ→0 gradient-line limit and δ→∞ global limit;
  a principled matched-mass calibration curve h(δ).
- **E2 — rate study for T2 (the conjecture; PRIORITY):** ambient linear + nonlinear
  ridges on S²: measure angle(top G_x^δ, top E_x^T C_ι^δ E_x) vs δ at a grid of
  probe points; log-log slopes (expect 2, gap-normalized); repeat for the global-block
  version to locate δ*; repeat for the non-ridge (does the rate survive?). Quadrature
  vs Monte Carlo error bars separated (nboot).
- **E3 — sanity on the preferential quadratic:** flows should organize around the
  quadratic's tangential axis; near-geodesic active flows expected at small curvature
  of the response.
- **E4 — singularity census:** η(x) heat map over the cap; locate/classify umbilics
  (index via loop winding); gap-colored trails (stop criteria à la tractography FA
  thresholds); verify robustness of Σ to bandwidth.
- **E5 — leaf accuracy (T6):** for near-ridges, Hausdorff distance between inactive
  flows and true level sets vs δ.
- **E6 — higher dimensions:** S^n via the same Gr(1, n+1) EVIE ops (code already
  generic in p=1); then a Grassmannian Gr(2, n) / SPD example (separable shape
  tensors) — the shape-space payoff; expect line fields → distributions, tie to T6.
- **E7 — data-driven gradients:** replace analytic ∇f with OPG-learned gradients
  (ties to Paper-2 infrastructure; App-A airfoil data as an anchor when available).
- **E8 — application teaser:** TREM transport curves / hail or airfoil sensitivity
  flows — ONE compelling real-data figure for the intro (defer heavy application to
  Paper 2 pipeline).

## 6. Visualization plan

- Gap-colored trails (η(x) along curves; tractography-style stopping visualized).
- LIC-style dense line-field rendering for the four-field comparison (talk-quality).
- 3D sphere trail figures (exists; polish) + the normal-coords comparison gallery
  (exists in this repo, deliberately kept out of the AMG paper — it migrates HERE).
- Animations for talks: manim scenes live in the zgrey-github-io venv per convention.
- The "trail visualization showing integral lines from the various dominant fields"
  (author's original food-for-thought) becomes the signature figure.

## 7. Infrastructure & repo plan

- **Phase 0:** create the Overleaf doc/repo (suggested name: `Active_Flows`;
  author to create on Overleaf, then git-clone alongside the other Overleaf
  projects). Seed its `RESEARCH_NOTES.md` from this plan. Add a row to
  `Overleaf/CLAUDE.md`'s project table. Decide venue early (candidates: SIMODS —
  coordinates with the AMG submission; JMLR for the statistical/T4–T5 weight;
  SIAGA if the geometry/topology (T3) dominates).
- **Code:** refactor `integral_curves.ipynb` into `python/active_flows.py`
  (fields, tracer, quadratures) + pytest coverage (transport identity vs EVIE,
  quadrature exactness, tracer reversibility, known-answer tests on the aligned
  linear case); notebooks stay thin drivers. Keep the vectorized S² transport
  (verified 1.2e-15 vs EVIE).
- **Naming discipline:** projected dominance (decompose→project, never extend) vs
  dominant projection (project→decompose) — the paper's assignment is canonical;
  the ITL deck legend and `projection_order.py` still carry the REVERSED labels and
  must be reconciled before any of this circulates.

## 8. Milestones (with stopping points, per the Overleaf workflow)

1. **M0 (setup):** Overleaf repo created; plan migrated; notebook refactored to
   module + tests. ⏹
2. **M1 (evidence):** E1 + E2 complete — the conjectured quadratic rate measured
   (or honestly refuted) with clean log-log evidence. ⏹ *Go/no-go: T2 slope.*
3. **M2 (core theory):** T1 lemma + T2 statement/proof drafted; E3/E4 supporting
   figures. ⏹
4. **M3 (flows are objects):** T3 well-posedness + singular-set classification;
   gap-colored trail figures. ⏹
5. **M4 (principled DR claim):** T4 variational formulation + identifiability
   bound (or an honest counterexample). This gates the title/framing. ⏹
6. **M5 (statistics + higher-dim):** T5 rates, E6. ⏹
7. **M6 (paper):** full draft → author review loop (purple convention) → submit.

## 9. Risks & honesty checklist

- The "localization recovers dominant projection" observation is ONE function, ONE
  center, ONE bandwidth pair so far — E2 must test across the registry, centers,
  and bandwidths before it headlines anything.
- Bandwidth selection can look arbitrary — T5 needs a defensible rule, else report
  sensitivity honestly.
- Line-field singularities may fragment flows in higher dimensions; dominant
  *direction* requires a gap — state η-dependence everywhere (the AMG paper's
  discipline carries over).
- Identifiability: do NOT claim "principled nonlinear dimension reduction" until
  T4's complexity control exists — `rem:dichotomy` applies to us too.
- SCMS (Ozertem–Erdogmus) is mechanically close; the differentiation (response
  gradients + transport + wrapped measures vs density Hessians) must be argued, not
  assumed.
- Cost: re-wrapped fields need fresh f/∇f evaluations per quadrature node — fine for
  analytic toys, expensive for simulations; the kernel variant over fixed data is the
  practical path and must be developed as a first-class citizen, not a fallback.
