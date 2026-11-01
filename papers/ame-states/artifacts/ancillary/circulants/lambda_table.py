#!/usr/bin/env python3
"""Recompute the first-moment columns of Table II exactly (stdlib).

E = sum_{w=1}^{floor(n/2)} C(n,w) (q^2-1)^w / (q^n + 1) is the expected number of nonzero stabilizers of party weight at
most floor(n/2) in a uniformly random F_q-linear Lagrangian subspace of F_q^(2n) (each fixed nonzero vector lies in it
with probability 1/(q^n+1)). lambda = E / (2n(q-1)) divides by scalar multiples and dihedral images; that division is a
heuristic (it treats orbits as independent). The column "circulants" is (q^floor(n/2) - 1)/(q - 1).
Printed as in the paper: E/(q-1) to 3 significant digits (comma-grouped above 100), lambda to one decimal.
"""
from fractions import Fraction
from math import comb

ROWS = [(15, 8), (17, 8), (16, 8), (19, 8), (16, 7), (19, 7), (16, 5), (19, 5), (19, 9), (19, 11), (13, 4), (15, 4)]
for n, q in ROWS:
    E = Fraction(sum(comb(n, w) * (q * q - 1) ** w for w in range(1, n // 2 + 1)), q ** n + 1)
    Eq = E / (q - 1)
    lam = E / (2 * n * (q - 1))
    eqf = float(Eq)
    eqs = f"{eqf:,.0f}" if eqf >= 100 else f"{eqf:.3g}"
    print(f"({n},{q})  E/(q-1) = {eqs:>8s}  lambda = {float(lam):5.1f}  circulants = {(q ** (n // 2) - 1) // (q - 1):,}")
