# K_G(3) ≤ 1.4494177506: proof objects

This bundle contains the exact certificate behind the main result of the paper *An improved upper bound on the
Grothendieck constant of order three from symmetric local models*:

> **K_G(3) ≤ 1.4494177506**, so the two-qubit Werner state has a local hidden-variable model for all projective
> measurements at visibility v ≤ 0.6899322156, and for all POVMs at v ≤ 0.45995.

The previous bound was K_G(3) ≤ 1.4546336567 (Designolle *et al.*, Phys. Rev. Research 5, 043059 (2023)). The bundle
does not ask you to trust the search that found the certificate: the claim can be re-checked from these files alone,
in exact arithmetic, by a short standard-library program.

## What is certified

The certificate `certificate/KG3_upper_1.4494177506_rad721.json` lists:

* 721 rational vectors u_1, …, u_721 ∈ ℚ³ (Bob uses the same set, `"bob": "alice"`). Their lengths are within 8·10⁻⁵ of 1.
* A rational scale λ = 3447864100660712344387584 / 4980863930023830797932435 (≈ 0.692222102249682).
* 19,318 strategies (n_k, s_k, t_k): positive integer weights over the denominator D = 70368744254637 (≈ 2⁴⁶), and sign
  vectors s_k, t_k ∈ {±1}⁷²¹ as hexadecimal bit masks (bit i set ⇔ sign +1).
* 60 signed index permutations (the icosahedral rotations acting on the 721 antipodal pairs), and the transposition flag.

The local model is L = 1/(G(1+tr)) Σ_h [h·N + (h·N)ᵀ] / D, with N = Σ_k n_k s_k t_kᵀ, G = 60, tr = 1, and
(h·X)_ij = σ^A_h(i) σ^B_h(j) X[π^A_h(i), π^B_h(j)]. Every term is a weighted deterministic strategy, so L is a local
correlation matrix with total weight ω = Σ_k n_k / D.

**Theorem 1 of the paper.** If ω + ‖λUUᵀ − L‖_F ≤ 1 then λUUᵀ is local, by Lemma 1 (the Frobenius unit ball lies in the
local polytope of full correlators). If in addition conv{±u_i} contains the ball of radius η, then the Werner state is
local for all projective measurements at visibility λη², and K_G(3) ≤ 1/(λη²).

For this certificate ω = 0.999999998769…, ‖E‖_F ≤ 2.31·10⁻¹⁰ and η² ≥ 0.996691977035244 (an exact rational with a
388-digit denominator). Hence K_G(3) ≤ 1.4494177506 after rounding up.

## Contents

| path | content |
|---|---|
| `certificate/KG3_upper_1.4494177506_rad721.json` | the certificate (8 MB, exact JSON); board post 241 of the project's run, 2026-10-15 |
| `check_release.py` | standalone checker (Python 3 standard library only; shares no code with the project verifier) |
| `negative_controls.py` | applies 15 corruptions (format, model, target, symmetry data and manifest tampering) and requires the checker to reject each |
| `expected_output.txt` | output of `check_release.py` and `negative_controls.py` on this bundle (timings will differ) |
| `verifier/` | the frozen project verifier v1.1 (`verify.py`, `hull3.py`, `accum.c`, tests, pinned hashes; `accum.c` needs a C compiler) |
| `data/records.json`, `seeds/` | the record table the verifier compares against, and the polyhedra its tests use |
| `MANIFEST.sha256` | SHA-256 of every other file |

## How to check

```bash
python3 -I check_release.py .
```

This takes about half a minute and needs no packages. The checker:

0. verifies the SHA-256 manifest;
1. parses every number as an exact rational (JSON floats are rejected), and checks the weights, masks, permutations and
   flags;
2. accumulates N exactly. Bob's index is packed into 64-bit lanes of Python big integers (the lane sums are bounded by
   Σ n_k ≤ D < 2⁶³). It then applies the 60 signed permutations and the transposition;
3. computes every entry of E = λUUᵀ − L as an exact rational and rounds |E_ij| up to a multiple of 2⁻⁸⁰. It then decides
   ω + ‖E‖_F ≤ 1 exactly on integers;
4. builds the convex hull of the 1442 points ±u_i by gift wrapping. Floating point only proposes the candidate points;
   every orientation test is exact. It then checks, exactly:
   * (a) the 2880 facets form a closed oriented surface: every directed edge occurs once and its reverse once;
   * (b) every facet (a, b, c) has det(a, b, c) > 0, i.e. the origin lies strictly inside its plane;
   * (c) five random generic rays from the origin each meet exactly one facet;
   * (d) every point lies on the closed inner side of every facet plane (all 4,152,960 pairs);
   * (e) the Euler characteristic is 2.

   Then η² is the minimum squared distance of a facet plane from the origin;
5. checks that 1/(λη²), rounded up to 10 digits, equals 1.4494177506, lies below the published 1.4546336567, and lies
   above the proven lower bound 1.4367003433 (Designolle, Vértesi, Pokutta, PRA 113, 022401 (2026));
6. runs the frozen project verifier as a separate program and requires the same bound (skip with `--no-verifier`, e.g.
   without a C compiler).

```bash
python3 -I negative_controls.py .
```

This takes about four minutes.

**Why (a)–(c) suffice for the inradius bound** (paper, Appendix A).
* By (b) the radial projection x ↦ x/|x| is well defined on the facets. By (a) it maps a closed oriented surface to the
  sphere, so it has a degree.
* For a generic direction, the degree is the signed number of facets whose cone contains it. By (b) every sign is +1,
  so by (c) the degree is 1.
* A map of nonzero degree is onto, so every unit direction u lies in the cone over some facet T. The ray along u meets T
  at distance at least dist(0, plane T) ≥ η, and that point lies in conv{±u_i}. Hence conv{±u_i} ⊇ η·Ball.
* Check (d) is not needed for this lower bound. It shows that η is exactly the inradius.

## Provenance

The certificate was found on 2026-10-15 by a multi-agent language-model system (Claude Sonnet 5.5 agents, coordinated by
Claude Opus 5.5) working against the frozen project verifier. Its ingredients are:
* symmetry reduction over the icosahedral rotations;
* Frank–Wolfe pools, interior-point linear programming and primal-dual column generation;
* the 721-direction polished geodesic sphere GP(12,0), with a small per-orbit radius and position deformation.

Every accepted certificate of the run (41 in all, results table in the paper's Appendix B) was re-checked by an
independent program. All 41 were re-verified with verifier v1.1, which performs every hull test exactly; v1.0 skipped
point/facet tests that a floating-point screen found far inside. All bounds were unchanged.
