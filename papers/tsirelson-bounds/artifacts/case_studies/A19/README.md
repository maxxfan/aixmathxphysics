# Case study: the published maximum of A19 is too low

Reproduce everything below with

```bash
python3 case_studies/A19/check_A19.py
```

(standard library only, about 2 seconds; run from any directory).

## The inequality

A19 is inequality 19 of the Ito–Imai–Avis list of tight 4522 Bell inequalities (four binary settings for
Alice, five for Bob). We use the Collins–Gisin form of Pál and Vértesi: with projectors `A_x`, `B_y`,

```
W = sum_x bA[x] A_x + sum_y bB[y] B_y + sum_{x,y} bAB[x][y] A_x B_y,        violation = <psi|W|psi> - b0,

bA  = [0, -1, -2, 0]
bB  = [-3, -3, 0, 0, 0]
bAB = [[ 2, 1, -1, -1,  0],
       [ 2, 1,  1,  0, -1],
       [ 2, 1,  0,  1,  1],
       [-2, 2,  0,  0,  0]]
b0  = 0
```

The coefficients come from `data/naceur2025_A19.txt`, a file in the repository of Naceur, Wang, Magron and
Acín (github.com/nininaceur/CertifiedQuantumBounds, MIT license). That file states the inequality as a
polynomial to be minimized. Our `W` is minus that polynomial. Under this sign the local bound is 0, and
`check_A19.py` confirms it by enumerating all 2^4 · 2^5 deterministic strategies. The same coefficients are
in `manifest.json` and `data/catalog.json`.

## The earlier source value

* Pál and Vértesi, Phys. Rev. A 79, 022120 (2009), Table II, row `A_19`. It gives violation 0.6226300 from
  a complex qubit strategy (C^2, projector ranks 1111 / 11111) and an NPA bound at level 2a equal to that
  value, so the case is listed as tight. The row is in `data/pv2009_table_row.tex`.
* The benchmark data of Naceur et al. (`maximal_violations_exact.txt`; the A19 line is in
  `data/naceur2025_maximal_violations_A19.txt`) list the same value, -0.6226300 in their sign convention.

## The certificates (all exact)

| object | file | value |
|---|---|---|
| complex qubit strategy (the published optimum) | `data/complex_qubit_lower.json` (board #403) | 0.62262997438… |
| real strategy, local dimension 4, all nine projectors of rank 2 | `../../certificates/A19/lower.json` (board #544) | 0.62263044335… |
| SOS certificate, NPA level 2 plus the words aab, abb, N = 202 | `data/upper_N202.json` (board #466) | ≤ 0.62263045521 |
| SOS certificate, NPA level 3, N = 318 | `../../certificates/A19/upper.json` (board #543) | ≤ 0.62263044795 |

The strategies list the spanning vectors of each projector and the state vector `psi`, all as exact
rationals or Gaussian rationals. The projector is `V (V^dag V)^-1 V^dag`, so it is an exact orthogonal
projector by construction. The SOS certificates list the words `w_i` and the exact rational matrix `Y`
(integer numerators over `Y.denominator`). The identity `sum Y_ij w_i^* w_j = lambda 1 - W` holds in the
universal projector algebra. The certified bound is `lambda + N tau` (Theorem 1 of the paper). For the
level-3 certificate the verifier needed no Cholesky shift (sigma = 0). Its stated shift is
tau = 2.77e-13, and `check_A19.py` proves `Y + tau I ⪰ 0` at that shift in exact integer arithmetic
(the exact residual gives 7.1e-15). `board_reports.json` holds the verifier's original reports.

## Conclusion

The certified window is [0.6226304433, 0.6226304480]. The exact real dimension-4 strategy exceeds the
published value 0.6226300 by 4.43e-7, far more than rounding to seven decimals allows (5e-8). Every NPA
level is an upper bound on the quantum maximum, so the published level-2a value of 0.6226300 cannot be a
correct upper bound. It must reflect solver inaccuracy. The published strategy value is the complex-qubit
local optimum 0.6226299744. The maximum is attained, to within 4.6e-9, by the real four-dimensional
strategy, and it is 4.7e-7 higher.

What the case study does not show: the window does not identify the exact maximum, and we do not know which
part of the original computation was inaccurate.
