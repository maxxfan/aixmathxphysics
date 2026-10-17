# Audit of published GF(2^n)-multiplier T-counts against released data (2026-10-17)

Prompted by an external review of our paper, which asked us to separate the best *implied* bound
(for example 3 × T(5) = 147 at n = 10, from the published T(5) = 49) from certified circuits. To
use the published 49-gate block we downloaded the released data, with the user's permission. It
turned out not to be GF(32) multiplication. We then checked every released matrix. Data
provenance is in `data/fisher26/README.md`. The scripts are in `tools/published_audit/` and the
raw output is in `results/published_audit_output.txt`.

## Method

A released parity matrix (rows = T gates, columns = the qubits a, b, c) determines its signature
tensor exactly. It decomposes the GF(2^n) multiplier, possibly composed with CNOT circuits on the
registers, if and only if three things hold:

1. every single and every pair of variables appears in an even number of rows;
2. every odd triple has the form (a_i, b_j, c_k);
3. the bilinear map M(a,b) defined by those triples is isotopic to multiplication in GF(2^n).

Isotopy means equality after invertible linear maps on a, b and c. These are exactly the
freedoms of basis, modulus and output CNOT network. A field has no zero divisors, and isotopy
preserves that property. So if some a ≠ 0 makes the map b ↦ M(a,b) singular, M is not GF(2^n)
multiplication in any basis.

Condition 3 is checked in two ways:

- **Zero-divisor test:** test every a for n ≤ 16, and 4000 random a for n = 32.
- **Positive certification:** map the matrix to our polynomial basis (`qlib/basis_map.py`) and run
  the frozen verifier.

## Results

| source | files | conditions 1–2 | zero divisors | verdict |
|---|---|---|---|---|
| Fisher et al. (F26), Khoruzhii-start finals | n = 8 (123), 9 (139), 10 (157) | yes | none | **certified** GF(2^n) multipliers (mapped, verified) |
| F26 / VarTodd, Khoruzhii starts | n = 3..10 | yes | none | field-like |
| F26 / VarTodd, Vandaele starts (`gf_mult_Vandaele_wo_ancilla`) | n = 3..16, 32 (40, 42, 56, 98, 132, 164, 189, 320/324, 334, 376, 396, 460, 492, 536, 567, 1701 columns) | yes | **yes**, for every n: for example 2 of the 7 nonzero a at n = 3, 12/31 at n = 5, 852/2047 at n = 11, 32762/65535 at n = 16 | **not** GF(2^n) multiplication |
| F26, Vandaele-start finals | n = 5 (49), 6 (79), 9 (165), 10 (187), 11 (215), 12 (243), 13 (293), 14 (319), 15 (363), 16 (385; Table I says 383), 32 (1177) | yes | **yes** | **not** GF(2^n) multiplication |

The defect is already present in the starting matrices, and the evolution preserves it: at n = 5
the start (98 columns) and the final (49) encode the identical map. In each defective start, the
products x^0 … x^{n-1} map correctly, but some high-degree reductions are wrong. Examples:

- n = 3, P = x³+x+1: x⁴ ↦ 1+x+x², where the correct value is x+x².
- n = 4, P = x⁴+x+1: x⁶ ↦ 1+x²+x³, where the correct value is x²+x³.
- n = 5, P = x⁵+x²+1: x⁶ ↦ 1+x+x³, x⁷ ↦ x+x²+x⁴ and x⁸ ↦ 1+x³, where the correct values are x+x³, x²+x⁴ and 1+x²+x³.

We could not identify the exact cause; no simple CNOT-convention change reproduces it. The other
start folder (`gf_mult_Vandaele_standart_wo_ancilla`) holds the unoptimized benchmark circuits
merged into one phase polynomial. Its tensors contain triples outside the (a, b, c) form, so they
are not the multiplier tensor either.

**Consequence.** The F26 Table I values from the Vandaele start are not T-counts of GF(2^n)
multipliers. This covers n = 3: 27, 4: 39 (equal to ATQ's valid 39), 5: 49, 11: 215, 12: 243,
13: 293, 14: 319, 15: 363, 16: 383 and 32: 1177. The value n = 64: 3729 was not released and comes
from the same pipeline. The Khoruzhii-start values (n = 7: 101, 8: 123, 9: 139, 10: 157) are
correct where we could check them (n = 8, 9, 10).

## Vandaele's n = 3 value (23)

His circuits and code are public (github.com/VivienVandaele, Rust, no license; used locally only).
The GF(2^n) benchmark circuits (`circuits/inputs/gf2^n_mult.qc`, n ≥ 4) declare only a and b as
inputs, so the output register starts in |0⟩. Simulation of gf2^5 shows it is correct
out-of-place (1024/1024). In-place it maps c ↦ Ñc ⊕ ab for an invertible CNOT map Ñ, which is
Clifford-equivalent to the in-place multiplier.

The synthesis code from Ref. [61] (`cargo run -r 3 2 0`) writes a 12-qubit GF(8) circuit with 8
CCZ gates and 3 clean ancillas. The arguments `3 1 0` produce the identical file, which is wrong
for x³+x+1. His optimizer (default pipeline: FastTMerge, InternalHOpt, FastTODD) brings this
circuit to **31** T gates, still with 3 ancillas. We simulated the result: it is correct in-place
for all c, with the ancillas at 0. The published 23 started from a 6-CCZ circuit that is not in
either repository, so it could not be reproduced. On the benchmark circuits for n = 4–10 the same
pipeline gives 49, 81, 113, 155, 205, 257 and 315 T gates.

## Best published values that we could verify (for the in-place, ancilla-free multiplier)

| n | best verifiable published | ours |
|---|---|---|
| 3 | 29 (ATQ; also Khoruzhii et al.) | 29 |
| 4 | 39 (ATQ) | 39 |
| 5 | 59 (ATQ) | **58** |
| 6 | 77 (ATQ) | **69** |
| 7 | 101 (Khoruzhii et al.; F26 K-start; not checked by us) | **97** |
| 8 | 123 (ATQ, F26 K-start certified) | **113** |
| 9 | 139 (F26 K-start, certified) | 139 |
| 10 | 157 (F26 K-start, certified) | **151** |
| 11–15 | none known | 181, 191, 219, 229, 253 |
| 16, 32, 64 | none with released data. V25 Table 3 lists 425, 1255 and 3817 for its ancilla-using circuits | 267, 571, 1255 |

The verifier's record table (`data/records.json`) now uses 59 for n = 5 and keeps every other
published value as a conservative bar.
