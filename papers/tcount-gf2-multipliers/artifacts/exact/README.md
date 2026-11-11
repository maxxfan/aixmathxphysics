# Exhaustive enumeration for the two-dimensional algebras

The paper calls three block values exact:

| algebra | minimum T-count | minimizers | status |
|---|---|---|---|
| GF(4) | 17 | 1008 | proved by Ruiz et al. (AlphaTensor-Quantum, Z3 solver); also Mazumder (exact coset decoding); confirmed here |
| F2[y]/(y²) | 15 | 112 | this work |
| F2 × F2 (two disjoint CCZ gates) | 13 | 14 | value reported by Campbell and Howard, confirmed by Mazumder; reproduced here |

Each value is a minimum over **all** decompositions of the algebra's cubic form on N = 6 variables.
The decompositions are the solutions of 41 linear equations in 63 unknowns over F2 (paper, Lemma 1),
so they form a coset of a code of dimension 22, and the minimum is found by visiting all 2²²
elements. These are optima within the phase-polynomial model of the paper (Hadamards only on the
output register). They are not lower bounds for circuits with intermediate Hadamards, measurements
or feed-forward.

Three programs do this enumeration. They share no code.

| program | method | needs |
|---|---|---|
| `enumerate_small.py` | plain Gray-code walk over the coset; prints the complete weight distribution, a minimizer, and checks that the released blocks `blocks/block_alpha_n2.json` and `blocks/block_ring_e2.json` attain the minimum | Python 3 standard library |
| `../check_release.py` (function `exact_small`) | bit-sliced: all 2²² coset elements at once in big integers | Python 3 standard library |
| `exact_alg2.py P0 P1` | the project auditor's original program (numpy arrays); `1 1` = GF(4), `0 0` = F2[y]/(y²), `0 1` = F2 × F2 | numpy |

```bash
python3 -I exact/enumerate_small.py .       # about 15 s; output is the last block of expected_output.txt
python3 exact/exact_alg2.py 1 1             # coset size 4194304, min weight 17, minimisers 1008
python3 exact/exact_alg2.py 0 0             # min weight 15, minimisers 112
python3 exact/exact_alg2.py 0 1             # min weight 13, minimisers 14
```

`check_release.py` runs its own enumeration on every run and fails unless the minima are 17, 15
and 13 with 1008, 112 and 14 minimizers.
