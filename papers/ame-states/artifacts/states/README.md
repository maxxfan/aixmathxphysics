# The six explicit states (paper Theorems 1 and 2, Table III)

AME(15,8) and AME(17,8) (symmetric circulant graph states over GF(8)), AME(19,16) and AME(21,16) (circulant graph states
over GF(16); n = 21 with the twisted symmetry c_(8j mod 21) = c_j^4), and AME(18,16), AME(20,16) (Hermitian self-dual MDS
codes over GF(256) with generator [I | A]).

| file | content |
|---|---|
| `printed.json` | the coefficients printed in the paper, with the field conventions (GF(8) = F_2[x]/(x^3+x+1), GF(16) = F_2[x]/(x^4+x+1), GF(256) = F_2[x]/(x^8+x^4+x^3+x^2+1); element = integer with bit i the coefficient of x^i) |
| `candidates/AME_n_d.json` | the certified states in the verifier's canonical F_2 form (m qubits per party), exactly as submitted to the board |
| `check_states.py` | field-level checks with its own GF(2^m) arithmetic: Lemma 2 for the graphs (every set of floor(n/2) parties, one representative per rotation orbit), A conj(A)^T = I and all C(n,n/2) column selections for the Hermitian codes, the twist of AME(21,16); `--convert` rebuilds the canonical rows from the printed data with `qlib_ame.py` (the converter that produced the candidates) and compares the F_2 row spaces |
| `qlib_ame.py` | the project's converter library (copied unchanged) |
| `run.sh` | replay (about 1 minute); `expected/` holds the expected outputs |

The canonical check (`../checker/ame_check.py state candidates/*.json`) tests every balanced cut of the certified states:
6,435, 24,310, 48,620, 92,378, 184,756 and 352,716 sets for n = 15, 17, 18, 19, 20, 21.
