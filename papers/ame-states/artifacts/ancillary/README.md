# Ancillary claims

| claim (paper) | folder | command | expected |
|---|---|---|---|
| Sec. III: all symmetric circulants over GF(8) (first coefficient 1): 299,593 for n = 15 with 336 AME, 2,396,745 for n = 17 with 48 AME; under multipliers and Frobenius (scaling absorbed) 28 and 4 orbits, all of size 12 | `circulants/` | `python3 circulants/rescan.py hits quick` (about 2 min; `full` re-checks all 384 hits in Python, about 5 min) | `circulants/expected/hits_quick.txt`, hit lists and orbit representatives in `circulants/expected/hits.json` |
| Table II: first-moment columns E/(q-1), lambda (exact; denominator q^n+1, see the script's docstring) and family sizes | `circulants/` | `python3 circulants/lambda_table.py` | `circulants/expected/lambda_table.txt` |
| Table II: zero-hit exhaustive scans | `circulants/` | `python3 circulants/rescan.py zero 16,8 16,7 16,5 19,5 13,4 15,4 19,7 19,8` (in `./run.sh full`); `python3 circulants/rescan.py zero 19,9` (about 10 min) and `python3 circulants/rescan.py zero 19,11` (68 min of CPU when measured) by hand | `circulants/expected/zero_full.txt`, `circulants/expected/zero_19_9.txt`, `circulants/expected/zero_19_11.txt` |
| Sec. IV: the shortened [17,8,10]_256 code of the AME(18,16) code is not GRS: its 17 columns impose rank 17 on the 36 quadratic monomials (19 independent quadrics), a GRS code of length 2k+1 = 17 gives rank 15 (21 quadrics; Glynn, Ball-Pepe as quoted in Ball arXiv:1907.04391); control: 17 points of the rational normal curve give 21 | `quadric/` | `python3 quadric/grs_quadric_test.py 92,74,70,173,117,59,208,212,182 [4]` (shortened at coordinate 0 or 4) | `quadric/expected/quadric.txt` |
| Sec. IV: the shadow inequality in purity form excludes no cell left uncovered by the snapshot (in particular not AME(8,4), AME(22,16), AME(23,16)); it reproduces every shadow exclusion of the table and rejects AME(4,2), AME(8,3) | `shadow/` | `python3 shadow/shadow_check.py` and `python3 shadow/shadow_check.py 8 4 22 16 23 16 4 2 8 3` | `shadow/expected/` |

The circulant kernel `circulants/circ.c` is the search code used in the run (unchanged). It tests one set per dihedral
orbit (valid for symmetric circulants). The Python re-check in `rescan.py` tests every set without using the symmetry.
`rescan.py` writes the field tables itself (GF(4) = F_2[x]/(x^2+x+1), GF(8) = F_2[x]/(x^3+x+1), GF(9) = F_3[x]/(x^2+x+2)).

`./run.sh` replays everything except the two long scans (about 3 minutes; `./run.sh full` about 10 minutes).
