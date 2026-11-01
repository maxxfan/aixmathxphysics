# Release for "Absolutely maximally entangled states with local dimensions 8 and 16, a coupling rule for composite dimensions, and a confirmation that no qubit-pair stabilizer AME(8,4) state exists" (revision 4)

This folder holds the data and the checkers behind every computational claim of the paper, in a form that runs from a
fresh copy with no network access. Published copy:
https://github.com/maxxfan/aixmathxphysics/tree/ame-states-v4/papers/ame-states/artifacts

```bash
shasum -a 256 -c SHA256SUMS      # or: sha256sum -c SHA256SUMS
./replay.sh                      # = ./replay.sh quick
./replay.sh full
```

`replay.sh` checks `SHA256SUMS`, runs the five parts below and prints one PASS/FAIL line per part and a verdict
(`RELEASE: PASS`, exit 0). Each part can be run on its own (`<part>/run.sh`); each part's README maps claims to commands
and expected outputs.

## Requirements and runtimes

* bash, a C compiler called as `cc` (or set `CC`), Python 3.9 or later.
* numpy for 24 of the 60 checks in `appendixB/` (they use `$PYTHON_NUMPY`, or `python3` if it imports numpy). Every
  other script uses the Python standard library only. Without numpy those 24 checks are skipped and the verdict is
  `RELEASE: INCOMPLETE` (exit 2).
* quick: about 40 minutes of CPU. `JOBS=4 ./replay.sh quick` took 14.6 minutes of wall time from a fresh copy
  (states 54 s, coupling 9 s, frames 410 s, appendixB 310 s, ancillary 94 s; 10-core laptop shared with other jobs);
  with `JOBS=1` expect about 45 minutes. `frames/` runs 4 checks in parallel unless `JOBS` is set.
* full: about 2.5 hours of CPU (appendixB full alone is about 100 minutes), so about an hour with `JOBS=4`; not timed
  as a whole (appendixB full took 30 minutes of wall time with `JOBS=4`).
* Memory: below 1.5 GB per process.

Measured times are in the part READMEs; they were taken on a 10-core Apple-silicon laptop shared with other jobs.

## Claims and where they are checked

| paper | claim | folder | replay |
|---|---|---|---|
| Thms 1, 2, Table III | AME(15,8), AME(17,8), AME(19,16), AME(21,16) (circulant graph states), AME(18,16), AME(20,16) (Hermitian self-dual MDS codes over GF(256)) | `states/` | printed coefficients and field conventions; field-level check with its own GF(2^m) arithmetic; conversion to the certified canonical rows; every balanced cut checked by the standalone checker |
| Sec. V, Lemma 3, Table IV | 26 new cells AME(11,d), AME(13,d), AME(15,d) by coupling (plus calibrations AME(5,10), AME(11,10)) | `coupling/`, `checker/` | 28 certificates, the 27 distinct ingredients in canonical form and at field level (generator matrices, fields, the bordered [[16,0,8]]_4 and the exceptional AME factors), ingredient-to-cell map, exact rank checks of every required set in C and in pure Python, four negative controls that must fail (one of them violates only a <= b) |
| Sec. VI B, Lemma 4, Table V (confirming Thm 3) | no F_2-linear (qubit-pair) stabilizer AME(8,4) state: frame search | `frames/` | the reduction and completeness argument (`frames/PROOF.md`), the ten gauge representatives from both implementations with their correspondence, both enumeration codes, positive controls, count hashes |
| App. B, Props 2 and 3, Table VI | every symmetry sector of a qubit-pair stabilizer AME(8,4) is empty (cross-check) | `appendixB/` | 60 checks: every computation of the appendix, counting steps made executable where possible, a per-row list of what remains a hand derivation |
| Sec. III, Table II | GF(8) circulant counts (336 / 48 hits, 28 / 4 orbits), first-moment columns (exact denominator q^n+1), zero-hit scans | `ancillary/circulants/` | exhaustive rescans with their own field tables and a symmetry-free Python re-check of the hits |
| Sec. IV | the shortened [17,8,10]_256 code of AME(18,16) is not GRS (19 vs 21 quadrics) | `ancillary/quadric/` | rank of the quadratic monomials, with a rational-normal-curve control |
| Sec. IV | the shadow inequality excludes no cell left uncovered by the snapshot | `ancillary/shadow/` | purity-form shadow test, calibrated on the table's exclusions |
| Secs IV, V, Table I | what "known", "nonexistent" and "not covered by the snapshot" mean | `snapshot/` | the three status files pinned by the verifier (SHA-256 prefixes checked), their sources and the program that built them; status of every source (refereed / preprint / repository / announcement) |

## Field conventions

GF(8) = F_2[x]/(x^3+x+1), GF(16) = F_2[x]/(x^4+x+1), GF(256) = F_2[x]/(x^8+x^4+x^3+x^2+1); a field element is written
as the integer whose bit i is the coefficient of x^i. Odd-characteristic and other fields used by the coupling
ingredients are listed in `coupling/fieldlevel/README.md`. Canonical form of a stabilizer state (n parties,
d = p^m): N = m n qudits over F_p, party j = qudits jm, ..., jm+m-1, N rows (a | b) for X(a)Z(b) (`checker/README.md`).

## What the replay does not cover, and known discrepancies

* The zero-hit circulant scans at (n,q) = (19,9) and (19,11) take about 10 minutes and 1 hour and are not in
  `replay.sh`; the commands and expected outputs are in `ancillary/README.md`.
* `frames/README.md` Sec. 9: one raw line count in a relaxed control differs by 30 (an edit of `vsearch8.c` after the
  original run added 30 near-miss lines with the same prefix); all search counts are reproduced exactly. Some expected
  values are documented numbers rather than logs, some random seeds were not recorded, and a few tools not needed for
  exhaustiveness are not packaged.
* `appendixB/README.md` "Discrepancies and caveats": the "1.4×10^6" solutions at threshold 4 in Table VI come from a
  capped search and are reproduced under a solution cap; the positive control of the (2^4) sector and the executable
  versions of several counting steps (`appendixB/extra/`) were written for this release; the (7,1) sector has one
  implementation, as the paper says.
* Scope (paper Prop. 1 and conclusion): the AME(8,4) results concern stabilizer states in the group-type
  nice-error-basis model (F_2-linear stabilizers with two qubits per party, which include the F_4-linear ones).
  Non-stabilizer AME(8,4) states are not covered.

## Provenance

All files were copied or generated from the project repository (quantum-information, `problems/p10-ame-states/`) at
revision 4 of the paper. Scripts copied from the original runs are marked as such in each part, with every change marked
`# [release]`; scripts written for the release say so in their header or in the part README.
