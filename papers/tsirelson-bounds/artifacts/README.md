# Certified Tsirelson bounds: certificates, catalog, verifier and exact replay

This folder is the artifact release (`tsirelson-bounds-v2`) for the paper *Certified Tsirelson bounds for
229 tight bipartite Bell inequalities*. It contains every certificate behind the paper's tables, the
inequality catalog, the frozen verifier with its tests, an independently written exact checker, and a
script that re-checks every advertised window from this folder alone.

## Quick start

```bash
cd release                                   # this folder, copied anywhere
shasum -a 256 -c SHA256SUMS                  # Linux: sha256sum -c SHA256SUMS
python3 replay/replay.py --tests --controls
```

The replay needs Python ≥ 3.9 and a C compiler (`cc`). It uses no third-party packages. With eight
processes on a 10-core laptop (Apple M-series, shared with other jobs) the full run took 12 minutes, and
it ends with `REPLAY PASS`. Most of the time goes into exact Bareiss elimination at N = 97. It writes `replay_report.json` to the current directory. Use `--only J97,A19` to replay
selected inequalities, and `--jobs N` to set the number of processes.

Single certificates can be checked directly:

```bash
python3 verifier/verify.py certificates/J97/upper.json            # frozen verifier: prints the verdict
python3 replay/exactcheck.py certificates/J97/lower.json --catalog data/catalog.json
python3 case_studies/A19/check_A19.py                              # the A19 case study
```

## What the replay checks

1. **Integrity.** Every file has the hash listed in `SHA256SUMS`. The verifier sources and frozen data
   have the hashes pinned in `verifier/PINS.json`.
2. **Floating-point environment.** `replay/fpenv.c` is built with the verifier's compiler flags. It checks
   binary64 evaluation, round-to-nearest, gradual underflow, the absence of fused multiply-add
   contraction, `nextafter`, a correctly rounded `sqrt` and overflow to infinity. As a negative control it
   is also built with `-ffp-contract=fast`. On hardware with FMA (arm64) that build reports the fusion.
3. **Every certificate, twice.** For each of the 238 strategies and 229 SOS certificates in `manifest.json`:
   * The frozen verifier `verifier/verify.py` checks it. It computes the exact violation two ways, checks
     the exact SOS identity, and obtains the rigorous shift τ from the residual kernel `verifier/rescheck.c`.
   * The independent checker `replay/exactcheck.py` checks it again. That code was written separately
     and shares nothing with `verifier/`. It has its own exact strategy evaluation, its own word reduction
     and identity check, and a proof that `Y + τ_stated I ⪰ 0` at the stated shift **in exact integer
     arithmetic**. For N ≤ 120 words it also runs exact Bareiss elimination at the stated shift.
   * Both results are compared with the manifest: the exact violation, λ, N, σ and τ. They are also
     compared with the printed window: printed lower ≤ exact violation, and exact bound ≤ printed upper.
4. **Tables.** The counts in the paper (239 inequalities, 229 two-sided, 223 windows below 1e-7, 146
   below 1e-8) are recomputed. No certified lower bound may exceed a certified upper bound. Appendix A is
   regenerated from the manifest and must equal `appendix_A.tex`, which is the LaTeX source of the paper's
   table. Table I (`table_I.tsv`) must agree with the manifest wherever it shows a certified value.
5. **Ladders and A19.** The certified rungs in `ladders/ladders.tsv` are re-evaluated exactly, and
   `case_studies/A19/check_A19.py` must reproduce the case study.
6. `--tests` runs the verifier's own test suite. `--controls` runs 15 negative controls
   (`replay/controls.py`), each of which must be rejected. The controls are:
   * a tampered SOS identity, and a tampered λ;
   * an SOS certificate whose identity still holds but whose Y is indefinite;
   * an understated shift for a genuine certificate, proven false by Bareiss elimination;
   * a tampered strategy, tampered inequality coefficients, tampered printed windows, and a tampered file;
   * indefinite and overflowing inputs to the residual kernel;
   * a 900-case randomized kernel soundness battery (singular, indefinite, subnormal and large matrices),
     in which every accepted shift is re-proved by exact elimination.

## Trust base

* **Lower bounds** use exact rational arithmetic only (Python integers and `fractions.Fraction`). Each
  projector is `V (V^dag V)^-1 V^dag`, built exactly, so it is an orthogonal projector by construction.
* **Upper bounds** rest on Theorem 1 of the paper: if `sum_ij Y_ij w_i^* w_j = λ 1 - W` holds in the
  universal algebra of projectors and `Y + τ I ⪰ 0`, then every commuting-operator strategy satisfies
  `<W> ≤ λ + N τ`. The identity is checked in exact integer arithmetic. Positive semidefiniteness is proven
  in two ways:
  * by the verifier's kernel `rescheck.c`, which needs the IEEE 754 assumptions that `fpenv.c` tests, and
  * by `exactcheck.py`, which needs no floating-point assumption at all. Floating point only produces a
    candidate factor L. The proof `Y + (σ + f) I = L L^T + (F + f I) ⪰ 0`, with `F` and
    `f = ceil(||F||_F)` computed exactly in integers, holds for any L whatsoever.

  The replay requires the exact shift to be at most the stated one, so the stated bound is proven
  without trusting the floating-point kernel.

## Layout

| path | contents |
|---|---|
| `manifest.json` | One entry per inequality, in the order of Appendix A. Each entry gives the coefficients (Collins–Gisin form), the published value, the printed window, and for each certificate: file, SHA-256, board submission, exact violation (lower) or exact λ, N, σ, residual, τ and bound (upper; binary64 values also in hex). |
| `appendix_A.tex`, `appendix_A.tsv` | Appendix A of the paper (LaTeX rows identical to the paper; plain-text form), regenerated by `generate/make_appendix.py` |
| `table_I.tsv` | Table I of the paper (open cases) |
| `certificates/<name>/lower.json`, `upper.json` | 238 strategies and 229 SOS certificates, byte-identical to the frozen submissions accepted by the verifier |
| `data/` | `catalog.json` (the 239 inequalities), `records.json` (reference values), `targets.json`, `make_records.py` (generates them from `data/sources/`: the transcribed tables of Pál–Vértesi 2009/2010 and Brunner–Gisin 2008, and the A-list files of Naceur et al., MIT license) |
| `verifier/` | frozen verifier `verify.py`, `qexact.py`, `rescheck.c`, `PINS.json`, `freeze.sh`, tests and fixtures |
| `replay/` | `replay.py`, `exactcheck.py`, `controls.py`, `fpenv.c` |
| `generate/` | `make_appendix.py`; `build_release.py`, the script that built this folder from the project repository (provenance; it needs the repository) |
| `case_studies/A19/` | the A19 case study: README, `check_A19.py`, inputs |
| `ladders/` | raw dimension-versus-value data of the ladder strategies, their certified rungs, the agents' notebooks |
| `search/` | the programs that found the certificates (provenance only; they need numpy/scipy and are not needed to verify anything) |
| `SHA256SUMS` | SHA-256 of every file |

## Certificate formats

Numbers are exact: integers, `"p/q"` strings or decimal strings, never JSON floats.

**Strategy.** `{"type": "strategy", "inequality": NAME, "field": "real"|"complex", "dA": dA, "dB": dB,
"A": [...], "B": [...], "psi": [...]}`.
* `A` has one entry per Alice setting x; each entry is a list of spanning vectors (length dA) of the
  projector A_x. An empty list means A_x = 0. `B` is the same for Bob.
* Complex entries are `[re, im]` pairs.
* `psi` has dA·dB entries, indexed `a*dB + b`, and need not be normalized.
* The certified quantity is the violation `<psi|W|psi>/<psi|psi> - b0`.

**SOS certificate.** `{"type": "sos", "inequality": NAME, "words": [...], "Y": {"denominator": D,
"upper": [...]}, "lambda": "p/q"}`.
* Words are strings such as `"A1A3B2"`, with 1-based settings; `""` is the identity. Alice's letters
  commute with Bob's, and repeated letters collapse.
* `Y` has integer numerators over the positive integer `D`. It is given by its upper triangle (row i
  lists the entries j = i..N-1) or sparsely (`"entries": [[i, j, num], ...]`).
* The identity `sum_ij Y_ij w_i^* w_j = λ 1 - W` must hold exactly. λ is the coefficient of the identity.

**Inequality** (in `data/catalog.json` and `manifest.json`):
`W = sum_x bA[x] A_x + sum_y bB[y] B_y + sum_xy bAB[x][y] A_x B_y`, where `A_x`, `B_y` are the projectors
onto outcome 1. The violation is `<W> - b0`, and the local bound is 0 in this normalization.

## Provenance

Every certificate was submitted during the project's multi-agent runs, frozen when submitted, and checked
by the pinned verifier. `board_file` and `board_id` in the manifest identify the submission. The two
exceptions are the exact J97 qubit strategy (a verifier test fixture) and the copies under `ladders/` and
`case_studies/`. Third-party data: `data/sources/naceur2025_alist/` is from
github.com/nininaceur/CertifiedQuantumBounds (MIT license, `LICENSE` included).
