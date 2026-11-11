# K_G(3) ≤ 1.4494177506: proof objects

This bundle (release tag `grothendieck-kg3-v3`) contains the exact certificate behind the main result of the paper
*An improved upper bound on the Grothendieck constant of order three from symmetric local models*, and the paper's
40 other certificates:

> **K_G(3) ≤ 1.4494177506**, so the two-qubit Werner state has a local hidden-variable model for all projective
> measurements at visibility v ≤ 0.6899322156.

The result is about projective measurements only. For general measurements (POVMs) the Werner state is known to be
local up to v = 1/2 (Renner, Phys. Rev. Lett. 132, 250202 (2024); Zhang and Chitambar, Phys. Rev. Lett. 132, 250201
(2024)), which is stronger than the 2v/3 = 0.45995 that follows from the bound above.

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
| `all_certificates/` | the other 40 certificates accepted during the runs (paper, Appendix B), LZMA-compressed (24 MB; 215 MB uncompressed), and `INDEX.json`, which lists all 41 with their file hashes and the certified λ, v and K_G(3) of the paper's Tables II and III |
| `negative_controls.py` | applies 17 corruptions (format, model, target, symmetry data, manifest tampering, index tampering) and requires the checker to reject each |
| `expected_output.txt` | output of `check_release.py` (default mode and `--all`) and of `negative_controls.py` on this bundle (timings will differ) |
| `verifier/` | the frozen project verifier v1.1 (`verify.py`, `hull3.py`, `accum.c`, tests, pinned hashes; `accum.c` needs a C compiler) |
| `data/records.json`, `seeds/` | the record table the verifier compares against, and the polyhedra its tests use |
| `MANIFEST.sha256` | SHA-256 of every other file, except the two maintainer scripts |
| `make_manifest.py`, `build_all_certificates.py` | maintainer scripts (the second needs the project's `candidates/` directory, which is not in the bundle) |

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

### The other certificates

```bash
python3 -I check_release.py . --table2      # the eight certificates of the paper's Table II (about 3 minutes)
python3 -I check_release.py . --all         # all 41 certificates of Table III (about 20 minutes)
```

For every certificate listed in `all_certificates/INDEX.json` the checker decompresses the file (standard-library
`lzma`), compares the SHA-256 of the JSON with the index, runs steps 1–4 above, and requires that the number of
directions, the number of atoms, λ (rounded down to 7 digits), v = λη² (rounded down to 10 digits) and the bound
K_G(3) ≤ 1/(λη²) (rounded up to 10 digits) equal the index. These are the values printed in the paper. Unless
`--no-verifier` is given, the frozen project verifier is also run on each certificate and must report the same bound
(its verdict is `BELOW_RECORD` for all but the last, because the record it compares against is the final bound).

```bash
python3 -I check_release.py . --cert FILE [--claim 1.45...]     # any certificate, .json or .json.xz
```

prints the certified bound of one certificate and compares it with `--claim` if given.

Six of the 40 certificates carried absolute directory names of the search machine in their free-text metadata field
`script`. That prefix was removed, by replacing this one JSON string; nothing else was changed. `INDEX.json` flags
these files (`script_sanitized`) and gives, for every file, the SHA-256 prefix of the file as submitted during the
run (`run_sha256_16`, also part of the file name) and the verifier's canonical hash, which ignores the metadata
fields and is the same before and after. The certificate of the main theorem, `certificate/KG3_upper_1.4494177506_rad721.json`,
is byte-identical to the file submitted during the run (SHA-256 `9f35fedc…da9d`, as printed in the paper); its
`script` field therefore still contains such directory names. Neither checker reads that field.

### Pinned hashes

`verifier/PINS.json` pins the verifier's sources, the record table and also the `accum` binary compiled on the
machine of the runs. `accum` is not shipped; the verifier compiles `accum.c` on first use, and with another compiler or
platform its report then says `pins: MODIFIED: accum`. That is expected and harmless here: step 6 does not use the
pin report, the sources themselves are covered by `MANIFEST.sha256`, and steps 1–5 do not use the verifier at all.

### Rays

Check (c) needs one direction in general position. The checker draws its directions from a pseudo-random generator
with a fixed seed (2026), so a run is reproducible; whether a direction is generic and which facets it meets is
decided exactly. Randomness plays no role in correctness.

**Why (a)–(c) suffice for the inradius bound** (paper, Lemma 2 in Appendix A).
* By (b) the three vertices of every facet are linearly independent, so the cone over a facet is a simplicial cone.
  Call a direction generic if it lies in no cone over an edge, and let c(u) be the number of facet cones containing u.
* c is the same for all generic directions. Move u along a path that avoids the vertex directions and the finitely
  many points where edge cones in different planes meet. When u crosses the cone over an edge {a, b}, it leaves one
  of the two facets at that edge and enters the other: by (a) these are a facet with the directed edge (a, b) and one
  with (b, a), and by (b) they lie on opposite sides of the plane through 0, a, b. So c does not change.
* By (c), c = 1 for one generic direction, hence for all. Generic directions are dense and the cones are closed, so
  every unit direction u lies in the cone over some facet T. The ray along u meets T at distance at least
  dist(0, plane T) ≥ η, and that point lies in conv{±u_i}. Hence conv{±u_i} ⊇ η·Ball.
* In topological terms c is the degree of the radial projection of the facet list, a closed oriented surface.
* Check (d) is not needed for this lower bound. It shows that η is exactly the inradius.

## Provenance

The certificate was found on 2026-10-15 by a multi-agent language-model system (Claude Sonnet 5.5 agents, coordinated by
Claude Opus 5.5) working against the frozen project verifier. Its ingredients are:
* symmetry reduction over the icosahedral rotations;
* Frank–Wolfe pools, interior-point linear programming and primal-dual column generation;
* the 721-direction polished icosahedral geodesic sphere {3,5+}_(12,0) (the dual of the Goldberg polyhedron
  GP(12,0)), with a small per-orbit radius and position deformation.

Every accepted certificate of the run (41 in all, results table in the paper's Appendix B) was re-checked by an
independent program. All 41 were re-verified with verifier v1.1, which performs every hull test exactly; v1.0 skipped
point/facet tests that a floating-point screen found far inside. All bounds were unchanged.

## Changes since tag `grothendieck-kg3-v2`

This bundle is the release `grothendieck-kg3-v3` (2026-11-10), the one named in the paper's Data availability
statement. The bundle published under the tag `grothendieck-kg3-v2` (2026-10-18) contained the theorem's certificate
only. The changes since then, made in response to the referee report of 2026-11-09:

* **Added `all_certificates/`:** the other 40 certificates (LZMA-compressed) and `INDEX.json`. Six of them have a
  sanitized `script` metadata field, as described above.
* **`check_release.py`:** new modes `--all` and `--table2`; `--cert FILE` now accepts `.json.xz`, and for a
  certificate other than the theorem's it prints the bound instead of comparing with the theorem's value
  (`--claim X` compares with X). Steps 0–6 of the default mode are unchanged, and so is the function `check` used by
  the negative controls, apart from reading `.xz` files.
* **`negative_controls.py`:** two more controls, for `INDEX.json` (a wrong file hash, a wrong bound): 17 in all.
* **`README.md`:** the POVM statement was removed (it is superseded by the known POVM threshold 1/2); the direction
  set is named as a geodesic polyhedron, the dual of GP(12,0), not as GP(12,0); the inradius argument is given in the
  elementary form of the paper's revised Lemma 2; notes on the pinned `accum` binary and on the rays were added.
* **`make_manifest.py`:** skips the new maintainer script `build_all_certificates.py` (new file).
* **`MANIFEST.sha256`, `expected_output.txt`:** regenerated.
* **Unchanged:** `certificate/KG3_upper_1.4494177506_rad721.json` (same SHA-256), `verifier/`, `data/`, `seeds/`.
