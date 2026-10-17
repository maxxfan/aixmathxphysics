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
"optimal" small-algebra values in the paper (17 for GF(4), 15 for F2[y]/(y²), 13 for F2×F2) are
exhaustive optima *within this model*. They are not lower bounds for circuits with intermediate
Hadamards, measurements or classical feed-forward. A lower T-count does not by itself imply a
lower runtime or physical-qubit cost. The CNOT counts reported here come from naive ladders and
are not optimized.

## Contents

| path | content |
|---|---|
| `records/gf2_<n>_T<t>.json` | modulus, the t factor vectors (character i = variable i: a_0..a_{n-1}, b_0..b_{n-1}, c_0..c_{n-1}), Clifford correction, gate counts, construction summary, source file in the project repository |
| `circuits/gf2_<n>_T<t>.qasm` | the full circuit |
| `blocks/` | small-algebra decompositions used as building blocks: GF(4) 17, GF(8) 29, GF(16) 39 (AlphaTensor-Quantum, CC-BY 4.0, mapped to the polynomial basis), F2[y]/(y²) 15, F2[y]/(y³) 25, GF(4)[ε]/(ε²) 35 |
| `construction/curve64.json` | genus-2 construction behind n = 64: the curve y² + y = x⁵ over GF(16), the shifted point P_a, the degree-16 place Q, the 33 evaluation functionals λ′_k (as 33×16 matrices over GF(16)), the coefficients ω_k, and the identity a·b = Σ_k λ′_k(a) λ′_k(b) ω_k |
| `check_release.py` | standalone checker (Python 3 standard library only; shares no code with the project verifier) |
| `negative_controls.py` | applies 25 corruptions to copies of the bundle (circuit, record, block, format and manifest tampering) and requires the checker to reject each |
| `construction/regen_n64.py` | regenerates the raw 1287-gate n = 64 network from `curve64.json` and the GF(16) block alone, and checks it against the digest in `curve64.json` |
| `expected_output.txt` | output of the checker on this bundle |
| `verifier/`, `data/records.json` | the frozen project verifier (`verify.py`, `phasecheck.c`, pinned hashes) and the record table it compares against |
| `published_audit/` | the check of published values against their released data (method, scripts, hashes of the audited files, output) |
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

It also checks the blocks and the n = 64 curve data, verifies `MANIFEST.sha256`, and requires
exactly the thirteen records and six blocks the paper claims. The run takes about ten seconds, and
its output should match `expected_output.txt`. `python3 -I negative_controls.py .` applies 25
corruptions and shows that the checker rejects each one. Most come from red-teaming the checker:
a statement hidden on a header line or behind a non-newline line separator, another gate library,
non-ASCII digits or leading zeros, malformed moduli and factor strings, a missing or duplicated
record, an empty bundle, wrong block counts, a ring relabelled as a field, wrong gate-count,
Clifford-correction or comparison fields, and a stale or incomplete manifest.

The checker certifies the circuits ("This work" in the paper). It does not certify the published
values used for comparison: each record file lists them under `comparison`. There, `reported` is
the best published value, flagged `reported_invalid` when its released data is not a GF(2^n)
multiplication tensor (see `published_audit/`), and `checked` is the best published value whose
released decomposition we verified. The project verifier can be run on any record too:
`python3 verifier/verify.py records/gf2_8_T113.json` (it reads `data/records.json`).

## How the records were obtained

Each record file has a `construction` field: the composed start (subfield interpolation, CRT over
places with local-ring blocks, or the genus-2 curve) and the reduction or annealing step. The
starts are deterministic. The reductions and annealing runs are seeded but time-limited, so the
exact final vectors are the record of the result; re-running a search is not expected to reproduce
them bit for bit.
