# Release for "Absolutely maximally entangled states with local dimensions 8 and 16, a coupling rule for composite dimensions, and a confirmation that no qubit-pair stabilizer AME(8,4) state exists" (revision 5)

This folder holds the data and the checkers behind every computational claim of the paper, in a form that runs from a
fresh copy with no network access. The last published copy is tag `ame-states-v4.1`,
https://github.com/maxxfan/aixmathxphysics/tree/ame-states-v4.1/papers/ame-states/artifacts
**This folder differs from that tag** (see "Changes since tag ame-states-v4.1" below) and has not been tagged or published;
the paper's Data availability section names the tag `ame-states-v5`, which does not exist until a human creates it.

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
* full: about 2.5 hours of CPU (appendixB full alone is about 100 minutes). `JOBS=4 ./replay.sh full` took 39 minutes
  of wall time from a fresh copy (states 40 s, coupling 57 s, frames 326 s, appendixB 1419 s, ancillary 475 s).
* Both modes of this release were then replayed from a fresh copy with `JOBS=4`, on the same laptop under heavier load
  from other jobs: `RELEASE: PASS` in 65 minutes (`full`) and 22 minutes (`quick`).
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
| Sec. III ("Status in the code tables") | the [[15,0,8]]_8 and [[17,0,9]]_8 codes that codetables.de shows since 2026-10-21 are symmetric circulant graph codes contained in our scan lists; the n = 17 code equals the state of Theorem 1 up to local Cliffords and a party permutation; the n = 15 code does not (it is another of the 4 classes of circulant AME(15,8) states, which are pairwise inequivalent under local Cliffords and permutations); the 4 orbits for n = 17 form one class | `ancillary/codetables/`, `snapshot/sources/codetables_2026-11-10/` | parses the cached pages, own GF(8) arithmetic, every cut re-checked, explicit equivalence map checked row by row; a local-Clifford invariant computed for all 28 orbit representatives |
| Sec. V (instances of Lemma 3) | the Lean-verified AME(9,10) state of 2026-08-28 is the case K = 4, a = 2, b = 5 of Lemma 3 | `ancillary/ame910/` | both ingredients re-checked (252 cuts mod 5; 127 common stabilizers, minimum weight 5) |
| Secs IV, V, Table I | what "known", "nonexistent" and "not covered by the snapshot" mean | `snapshot/` | the three status files pinned by the verifier (SHA-256 prefixes checked), their sources and the program that built them; status of every source (refereed / preprint / repository / announcement); what the snapshots miss (AME(9,10); the code-table change of 2026-10-21); the dated prior-art check of 2026-11-10 |

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

## Changes since tag ame-states-v4.1

Made on 2026-11-10 in answer to the exam of 2026-11-09 (`reviews/exam-20261109.md` in the project repository). **Not
tagged, not published.** No state, certificate, search, checker or expected count of v4.1 changed; every file of v4.1
other than those listed under "changed" is byte-identical.

* added `snapshot/sources/codetables_2026-11-10/` (41 pages of codetables.de fetched 2026-11-10, `INDEX.json`),
  `snapshot/prior_art_check_2026-11-10.md`;
* added `ancillary/codetables/` (comparison of the code-table codes [[15,0,8]]_8 and [[17,0,9]]_8 with Theorem 1 and the
  scan lists; an invariant that separates the four classes of n = 15 solutions) and `ancillary/ame910/` (ingredients of the Lean-verified AME(9,10)); both run in `ancillary/run.sh`, in
  both modes;
* changed `README.md`, `SHA256SUMS`, `ancillary/README.md`, `ancillary/run.sh`, `snapshot/README.md` (the snapshots are
  historical: AME(9,10) is known, the code tables changed on 2026-10-21, Ball-Moreno-Simoens is published in IEEE
  Trans. Inf. Theory);
* not changed although out of date: `frames/REVIEW.md` (a dated review; it says the journal version of
  Ball-Moreno-Simoens was "not confirmed by us") and `snapshot/literature_log_2026-10-19.md` (the log of the original
  run).

The record of the replays of this state of the folder (both modes, from a fresh copy) is kept outside the release, in
the project repository (`paper/README.md`, revision 5), so that the replayed files are exactly the released files. The
runtimes quoted above are those measured for v4.1; the three added checks take about 90 s more.

## Changes since tag ame-states-v4

Tag `ame-states-v4` lacked `ancillary/circulants/expected/hits_full.txt`, so `./replay.sh full` ended in
`RELEASE: FAIL` on every machine, although the run itself re-checked all 336 and 48 hits with 0 failures. This release
adds that file, records the measured time of the full mode and points to the new tag. No state, certificate, search,
checker or result changed.

## Provenance

All files were copied or generated from the project repository (quantum-information, `problems/p10-ame-states/`) at
revision 5 of the paper (the files of tag ame-states-v4.1 at revision 4.1). Scripts copied from the original runs are marked as such in each part, with every change marked
`# [release]`; scripts written for the release say so in their header or in the part README.
