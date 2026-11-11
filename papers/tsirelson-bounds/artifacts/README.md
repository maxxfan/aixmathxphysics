# Certified Tsirelson bounds: certificates, catalog, verifier and exact replay

This folder is the artifact release for the paper *Certified Tsirelson bounds for
238 tight bipartite Bell inequalities*. It contains every certificate behind the paper's tables, the
inequality catalog, the frozen verifier with its tests, a separately written exact checker, and a
script that re-checks every advertised window from this folder alone.

**Version.** The published release is the tag `tsirelson-bounds-v2`. This folder is the candidate for
`tsirelson-bounds-v3`: it differs from v2 as listed under "Changes since tsirelson-bounds-v2" at the end of
this file, and it has not been tagged or published.

## Quick start

```bash
cd release                                   # this folder, copied anywhere
shasum -a 256 -c SHA256SUMS                  # Linux: sha256sum -c SHA256SUMS
python3 replay/replay.py --tests --controls
```

The replay needs Python ≥ 3.9 and a C compiler (`cc`). It uses no third-party packages. With eight
processes on a 10-core laptop (Apple M-series, shared with other jobs) the full run took 12 to 16 minutes, and
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
3. **Every certificate, twice.** For each of the 238 strategies and 238 SOS certificates in `manifest.json`:
   * The frozen verifier `verifier/verify.py` checks it. It computes the exact violation two ways, checks
     the exact SOS identity, and obtains the rigorous shift τ from the residual kernel `verifier/rescheck.c`.
   * The second checker `replay/exactcheck.py` checks it again. That code was written separately (by the
     same system, with the same mathematical method) and shares no code with `verifier/`. It has its own exact strategy evaluation, its own word reduction
     and identity check, and a proof that `Y + τ_stated I ⪰ 0` at the stated shift **in exact integer
     arithmetic**. For N ≤ 120 words it also runs exact Bareiss elimination at the stated shift.
   * Both results are compared with the manifest: the exact violation, λ, N, σ and τ. They are also
     compared with the printed window: printed lower ≤ exact violation, and exact bound ≤ printed upper.
4. **Tables.** The counts in the paper (239 inequalities, 238 two-sided, 225 windows below 1e-7, 147
   below 1e-8) are recomputed. No certified lower bound may exceed a certified upper bound. Appendix A is
   regenerated from the manifest and must equal `appendix_A.tex`, which is the LaTeX source of the paper's
   table. The certified columns of Table I (`table_I.tsv`) must equal the manifest windows, and its width
   column must be certified upper minus certified lower.
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
| `table_I.tsv` | Table I of the paper (open cases): the numerical window of Pál–Vértesi 2010, the uncertified numerical order-4 entry of the Naceur et al. repository file (orientation only), and the certified window with its width |
| `certificates/<name>/lower.json`, `upper.json` | 238 strategies and 238 SOS certificates. 228 SOS certificates and all strategies but one are byte-identical to frozen board submissions accepted by the verifier; see "Provenance" for the others |
| `provenance/` | `board_submissions.tsv`: export of the board's submission table (550 submissions: time, agent, target, file, verdict, value); `postboard_index.json` and `postboard_*.report.json`: index and frozen-verifier reports of the certificates computed after the runs |
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
onto outcome 1. The violation is `<W> - b0`, where `b0` is the local bound: 0 for every inequality of the
catalog except `I4422_7`, which has `b0 = 1`.

## Provenance

Every strategy, and 228 of the 238 SOS certificates, were submitted during the project's multi-agent runs,
frozen when submitted, and checked by the pinned verifier. `board_file` and `board_id` in the manifest
identify the submission, and `provenance/board_submissions.tsv` lists all 550 submissions. The exceptions:

* the exact J97 qubit strategy is a verifier test fixture;
* ten SOS certificates were computed after the runs, on 2026-11-10, with the same solver
  (`search/work/certify/`): A14, A67 and A80 at full NPA level 3, and A47, A62, A64, A68, A82, A84, A89 at
  level 2 plus the words aab and abb (276 words). They were checked by the same pinned verifier in a local
  run and are **not** board submissions: their manifest entries have `board_id: null` and a `provenance`
  note, and the verifier's reports are in `provenance/`. The replay checks them exactly like the others;
* the copies under `ladders/` and `case_studies/`.

Third-party data: `data/sources/naceur2025_alist/` is from github.com/nininaceur/CertifiedQuantumBounds
(MIT license, `LICENSE` included; repository commit c5c3b60, files last changed in commit d32a3db). The file
`maximal_violations_exact.txt` in that folder has no header. Our reading of its columns (reference value,
numerical NPA value at order 3, at order 4) is an inference, and its values are **not certified**.
`data/records.json` used the order-4 entries as the bar for upper-bound submissions during the runs; the
paper's Table I shows them in a separate column for orientation only.

The header comment of `verifier/verify.py` says the records come from "Pal-Vertesi 2009". The file is
frozen; the records are in fact Pál–Vértesi 2010 Table I with the overlays described in
`data/make_records.py`.

## Changes since tsirelson-bounds-v2

State of 2026-11-10, after the referee exam of 2026-11-09 (`reviews/exam-20261109.md` in the project
repository). **Not tagged, not published.** `SHA256SUMS` and `manifest.json` differ from the v2 tag.

* **Ten new upper certificates** (see "Provenance"): A14, A47, A62, A64, A68, A80, A82, A84, A89 had no SOS
  certificate in v2, and the A67 certificate (276 words, bound 0.3993386504) is replaced by a level-3 one.
  Counts: upper certificates 229 → 238, two-sided windows 229 → 238, windows below 1e-7 223 → 225 (A67 at
  7.6e-9 and A80 at 6.9e-8), below 1e-8 146 → 147, certificates 467 → 476. No lower bound changed.
* **Table I (`table_I.tsv`) has new columns.** v2 printed the uncertified order-4 numbers of the Naceur et al.
  repository as "literature upper" for the A list and mixed them with certified values. Now: `pv2010_lower`,
  `pv2010_upper` (published, numerical), `order4_numerical_uncertified`, `certified_lower`, `certified_upper`,
  `certified_width`. The width is certified upper minus certified lower, nothing else.
* **Replay** (`replay/replay.py`): the expected counts are the new ones; the Table I check now requires the
  certified columns to equal the manifest windows and the width to be their difference. The checks of the
  certificates themselves are unchanged.
* **Manifest**: `release` is `tsirelson-bounds-v3`; the `convention` text states the local bound b0 correctly
  (0 except for I4422_7, where b0 = 1); entries of certificates computed after the runs carry `board_id: null`
  and a `provenance` note; the paper title counts 238 two-sided windows.
* **New folder `provenance/`**: export of the board's submission table and the verifier reports of the ten
  new certificates.
* **Unchanged**: `verifier/` (same pins), `data/`, all 238 strategies, the other 228 SOS certificates,
  `replay/exactcheck.py`, `replay/controls.py`, `replay/fpenv.c`, the ladder data, the A19 certificates.
  `case_studies/A19/README.md` gained a paragraph on the repository file's order-3 entry for A19.
* Wording: `replay/exactcheck.py` is described as "separately written" (same system, same method, no shared
  code), not "independent".

