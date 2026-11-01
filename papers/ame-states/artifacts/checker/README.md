# Standalone checker

`ame_check.py` (Python 3.9+, standard library) and `cutrank.c` (compiled at run time with `cc -O2`) were written for this
release and share no code with the project verifier.

It checks a stabilizer state in canonical F_p form (n parties, d = p^m, N = m*n qudits, party j = qudits j*m..j*m+m-1,
N rows (a | b) for X(a)Z(b)):
1. the rows are isotropic for a.b' - b.a' (mod p) and have rank N (a pure stabilizer state);
2. k-uniformity: for every set B of n-k parties the N x 2m(n-k) block of the X and Z columns of B has rank N (paper
   Lemma 1: a nonzero stabilizer supported on the complement of B exists exactly when this rank drops). AME(n,d) is
   k = floor(n/2).
For a coupling certificate it checks that a is k-uniform on 2k+2 parties, b is AME(2k+2, b) and a <= b (Lemma 3).

```bash
python3 ame_check.py state FILE...          # AME check
python3 ame_check.py uniform K FILE         # K-uniformity
python3 ame_check.py coupling FILE...       # coupling certificates
```

By default every set is tested by the C enumerator and a seeded random sample of 200 sets is re-tested in Python;
`--python-only` tests every set in Python. Each line reports the number of sets tested, the number of full rank and the
number deficient; exit status 0 iff all checks pass.
