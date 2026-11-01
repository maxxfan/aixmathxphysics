# Certificates and checkers: exact quantum ordered search

Release for the paper *Exact certificates for quantum ordered search: seven queries for 700,000 elements and six for
91,000* (`../ordered-search.pdf`). A certificate for a pair (k, n) proves that k queries search an ordered list of n
elements exactly. Applied recursively, such an algorithm searches N elements with k⌈log_n N⌉ queries, about
(k / log2 n) log2 N. A certificate is a chain of symmetric Laurent polynomials q_1, ..., q_{k-1} in the
Farhi–Goldstone–Gutmann–Sipser characterization, with integer coefficients, together with spectral factors h_t that
make nonnegativity a diagonal-dominance check.

## Certificates

| k | n | k / log2 n | file | where |
|---|---|---|---|---|
| 3 | 56 | 0.516587 | `cert_k3_n56.json` | `results/certificates/` |
| 4 | 605 | 0.432864 | `cert_k4_n605.json` | `results/certificates/` |
| 5 | 5000 | 0.406911 | `cert_k5_n5000.json.gz` | `results/certificates/` |
| 6 | 45000 | 0.388158 | `cert_k6_n45000.json.gz` | `results/certificates/` |
| 6 | 80000 | 0.368376 | `cert_k6_n80000.json.gz` | `results/certificates/` |
| 6 | 90000 | 0.364573 | `cert_k6_n90000.json.gz` | `results/certificates/` |
| 6 | 91000 | 0.364220 | `cert_k6_n91000.json.gz` | `results/certificates/` |
| 7 | 650000 | 0.362505 | `cert_k7_n650000.json.gz` (64.4 MB) | release asset |
| 7 | 700000 | 0.360509 | `cert_k7_n700000.json.gz` (68.8 MB) | release asset |

The two seven-query certificates are too large for the repository. Download them from the assets of the release
<https://github.com/maxxfan/aixmathxphysics/releases/tag/ordered-search-v1> into `results/certificates/`.
`results/certificates/SHA256SUMS` pins all nine files:

```bash
cd results/certificates && shasum -a 256 -c SHA256SUMS
```

## Re-checking

The frozen verifier (`verifier/verify.py`, version 1.1, hash pinned in `verifier/PINS.json`) uses only the Python
standard library and exact integer arithmetic:

```bash
python3 verifier/tests/test_verifier.py
python3 verifier/verify.py results/certificates/cert_k6_n91000.json.gz
python3 verifier/verify.py results/certificates/cert_k7_n700000.json.gz
```

The independent checker shares no code with the verifier and checks every condition by a different route in FLINT
arithmetic. It needs Python with python-flint:

```bash
python3 results/check_independent.py results/certificates/cert_k7_n700000.json.gz
```

`results/audits.json` records the audits of every certificate: both checkers' verdicts and positivity margins.
