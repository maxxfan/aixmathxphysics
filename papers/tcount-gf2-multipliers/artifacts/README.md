# GF(2^n) multiplier T-count records: proof objects

This bundle contains the exact objects behind the thirteen T-count records of the paper *Certified
low-T-count multipliers for binary fields from Chinese remaindering*. These are
n = 5–8, 10–16, 32 and 64. The bundle does not ask you to trust the authors' search: every
claim below can be re-checked from these files alone.

## What is certified

For each n, the circuit `circuits/gf2_<n>_T<t>.qasm` implements the in-place multiplier
|a,b,c⟩ ↦ |a,b,c ⊕ a·b mod P⟩ on 3n qubits, with no ancillas. P is the modulus in the matching
record file, and the circuit uses exactly t T gates. It is an explicit OpenQASM 2.0 gate list made
of four parts, in order:

1. Hadamards on the c register.
2. One phase gadget per factor: a CNOT ladder, one `t`, and the reversed ladder.
3. A diagonal Clifford correction (`s`, `z`, `sdg`, `cz`).
4. Hadamards on the c register again.

The correction is derived deterministically from the factor list (paper, Lemma 1). It is also
stored in each record file under `clifford_correction`.

**Model.** The T-count here is the number of T gates in circuits of this form: Hadamards only on
the output register, and a phase polynomial built from CNOT, T and diagonal Cliffords. The
"exact" small-algebra values in the paper (17 for GF(4), 15 for F2[y]/(y²), 13 for F2×F2) are
exhaustive optima *within this model*. The value 17 for GF(4) was proved optimal by Ruiz et al.
(AlphaTensor-Quantum, with the Z3 solver) and obtained again by Mazumder (arXiv:2610.01024); our
enumeration confirms it independently. The value 13 for F2×F2 (two disjoint CCZ gates) was reported
by Campbell and Howard and confirmed by Mazumder. The value 15 for F2[y]/(y²) is from this work.
The enumeration is in `exact/` and is repeated by the checker. These values are not lower bounds
for circuits with intermediate Hadamards, measurements or classical feed-forward. A lower T-count does not by itself imply a
lower runtime or physical-qubit cost. The CNOT counts reported here come from naive ladders and
are not optimized.

## Contents

| path | content |
|---|---|
| `records/gf2_<n>_T<t>.json` | modulus, the t factor vectors (character i = variable i: a_0..a_{n-1}, b_0..b_{n-1}, c_0..c_{n-1}), Clifford correction, gate counts, construction summary, source file in the project repository |
| `circuits/gf2_<n>_T<t>.qasm` | the full circuit |
| `blocks/` | small-algebra decompositions used as building blocks: GF(4) 17, GF(8) 29, GF(16) 39 (AlphaTensor-Quantum, CC-BY 4.0, mapped to the polynomial basis), F2[y]/(y²) 15, F2[y]/(y³) 25, GF(4)[ε]/(ε²) 35 |
| `construction/curve64.json` | genus-2 construction behind n = 64: the curve y² + y = x⁵ over GF(16), the shifted point P_a, the degree-16 place Q, the 33 evaluation functionals λ′_k (as 33×16 matrices over GF(16)), the coefficients ω_k, and the identity a·b = Σ_k λ′_k(a) λ′_k(b) ω_k |
| `exact/` | exhaustive enumeration of all 2²² decompositions for GF(4), F2[y]/(y²) and F2×F2: a standalone program (standard library), the auditor's original numpy program, and a README |
| `check_release.py` | standalone checker (Python 3 standard library only; shares no code with the project verifier) |
| `negative_controls.py` | applies 25 corruptions to copies of the bundle (circuit, record, block, format and manifest tampering) and requires the checker to reject each |
| `construction/regen_n64.py` | regenerates the raw 1287-gate n = 64 network from `curve64.json` and the GF(16) block alone, and checks it against the digest in `curve64.json` |
| `expected_output.txt` | output of the four commands below on this bundle |
| `verifier/`, `data/records.json` | the frozen project verifier (`verify.py`, `phasecheck.c`, pinned hashes) and the record table it compares against |
| `published_audit/` | the check of published values against their released data: method and results (`AUDIT.md`), scripts, the mapping code `qlib/basis_map.py` used for the positive certification, hashes of the audited files, raw output |
| `MANIFEST.sha256` | SHA-256 of every file |

## How to check

```bash
python3 -I check_release.py .
```

For every record, the checker:

1. tests the modulus for irreducibility;
2. checks the signature-tensor criterion on the factor list;
3. proves that the QASM gate list equals the multiplier for **all** inputs. It tracks the linear
   form carried by each wire through the CNOTs, adds the exact Z₈ phase polynomial of every
   T/S/Z/CZ gate, and requires the sum to equal 4·c·(ab mod P) coefficient by coefficient, with
   every wire returned to its input. This is a proof, not a sample;
4. cross-checks by simulating the circuit on basis states: all 2^(3n) inputs for n ≤ 8, and 2048
   random inputs for larger n.

It also checks the blocks and the n = 64 curve data, repeats the exhaustive enumeration for the
three two-dimensional algebras (minima 17, 15 and 13), verifies `MANIFEST.sha256`, and requires
exactly the thirteen records and six blocks the paper claims. The run takes under a minute, and
its output should match the first block of `expected_output.txt`.

There are four replay commands. Run all of them from the bundle root:

```bash
python3 -I check_release.py .                 # the checker (all records, blocks, curve data, exact minima)
python3 -I negative_controls.py .             # 25 corruptions, each must be rejected (several minutes)
python3 -I construction/regen_n64.py .        # rebuilds the raw 1287-gate n = 64 network from the curve data
python3 -I exact/enumerate_small.py .         # second, plain enumeration of the 2^22 decompositions (about 15 s)
```

`negative_controls.py` applies 25 corruptions and shows that the checker rejects each one. Most
come from red-teaming the checker:
a statement hidden on a header line or behind a non-newline line separator, another gate library,
non-ASCII digits or leading zeros, malformed moduli and factor strings, a missing or duplicated
record, an empty bundle, wrong block counts, a ring relabelled as a field, wrong gate-count,
Clifford-correction or comparison fields, and a stale or incomplete manifest.

The checker certifies the circuits ("This work" in the paper). It does not certify the published
values used for comparison. Each record file lists them under `comparison`:

- `reported` is the lowest published value, with `reported_status`:
  `ok`; `invalid` when its released data is not a GF(2^n) multiplication tensor (see
  `published_audit/`; printed with `*`); `unverified` when the value was never released and so
  could not be checked (n = 64: only the released start of the same pipeline is shown defective;
  printed with `?`). `reported_invalid` is true exactly for status `invalid`.
- `checked` is the lowest published value whose released decomposition we verified.
- `vandaele2025_table3` is the value in Table 3 of Vandaele, Quantum 9, 1860 (2025), column
  "Ref. [61] & FastTODD" (n = 5–8, 10, 16, 32, 64; for n = 16, 32, 64 these are 425, 1255, 3817).
  The optimized circuits are not released (the generator and the optimizer are public; we did not
  rerun them at these sizes), so we checked neither their tensors nor their circuit model.

The project verifier can be run on any record too:
`python3 verifier/verify.py records/gf2_8_T113.json` (it reads `data/records.json`). It compiles
`verifier/phasecheck` next to its source; `check_release.py` ignores that one build product, so the
manifest check still passes afterwards. On another machine the verifier prints
`pins: MODIFIED: phasecheck`, because the pinned hash of the compiled binary is machine-specific;
the pins of `verify.py`, `phasecheck.c` and `records.json` must match. `data/records.json` is the
frozen record table: its "by" text for n = 2 still credits the optimality proof to this project,
which is out of date (see **Model** above); the file is pinned by hash and was left unchanged.

## How the records were obtained

Each record file has a `construction` field: the composed start (subfield interpolation, CRT over
places with local-ring blocks, or the genus-2 curve) and the reduction or annealing step. The
starts are deterministic. The reductions and annealing runs are seeded but time-limited, so the
exact final vectors are the record of the result; re-running a search is not expected to reproduce
them bit for bit.

## Changes since tag `tcount-gf2-v3`

This bundle is the one of the git tag `tcount-gf2-v4` (repository `maxxfan/aixmathxphysics`, path
`papers/tcount-gf2-multipliers/artifacts`); it accompanies the revision of the paper dated
10 November 2026. It replaces the bundle published under the tag `tcount-gf2-v3`, which went with
the third version of the paper and remains available under that tag.

Unchanged, byte for byte: all thirteen `circuits/*.qasm`, the factor lists, moduli, Clifford
corrections and gate counts in `records/*.json`, the factor lists in `blocks/*.json`,
`construction/`, `verifier/`, `data/records.json`, `negative_controls.py`. No T-count changed.

Changed:

- `records/*.json`: the `comparison` field. It now separates a published value whose released data
  is invalid from one that was never released (n = 64: 3729, now `unverified`), and it adds
  Vandaele's Table 3 values (425, 1255, 3817 for n = 16, 32, 64), which the earlier bundle omitted.
- `blocks/*.json`: the `source` text. The optimality of 17 for GF(4) is credited to Ruiz et al.
- `published_audit/`: `AUDIT.md` shows the current records (57 at n = 5 and 179 at n = 11, where
  it showed 58 and 181) and lists the n = 64 and n = 128 starts; `qlib/basis_map.py` and
  `qlib/gf2mult.py` are shipped, so `certify_npy.py` runs from the bundle; `output.txt` has the
  n = 128 line and the output of `certify_npy.py`.
- `exact/`: new. The enumeration programs behind the exact values.
- `check_release.py`: checks the new comparison fields, repeats the exact enumeration, and no
  longer counts the verifier's compiled binary `verifier/phasecheck` as an unlisted file.
- `README.md`, `expected_output.txt`, `MANIFEST.sha256`.
