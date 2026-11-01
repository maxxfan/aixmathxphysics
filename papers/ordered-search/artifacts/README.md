# Certificates and checkers: exact quantum ordered search (version 2)

Release for the paper *Exact certificates for quantum ordered search: seven queries for 700,000 elements and six for
91,000* (`../ordered-search.pdf`). This revision answers the external review of 31 October 2026. Each certificate is
an exact feasibility certificate for a pair (k, n). It is an explicit chain of symmetric Laurent polynomials with
integer coefficients that, through the Farhi–Goldstone–Gutmann–Sipser characterization, proves that an exact k-query
algorithm for ordered search on n elements exists. The unitaries are not constructed. Applied recursively, such an
algorithm searches N elements with k⌈log_n N⌉ ≤ (k / log2 n)·log2 N + k queries.

## Certificates

| k | n | bound on queries for N elements | file | where |
|---|---|---|---|---|
| 3 | 56 | 0.516587·log2 N + 3 | `cert_k3_n56.json` | `results/certificates/` |
| 4 | 605 | 0.432864·log2 N + 4 | `cert_k4_n605.json` | `results/certificates/` |
| 5 | 5000 | 0.406911·log2 N + 5 | `cert_k5_n5000.json.gz` | `results/certificates/` |
| 6 | 45000 | 0.388158·log2 N + 6 | `cert_k6_n45000.json.gz` | `results/certificates/` |
| 6 | 80000 | 0.368376·log2 N + 6 | `cert_k6_n80000.json.gz` | `results/certificates/` |
| 6 | 90000 | 0.364573·log2 N + 6 | `cert_k6_n90000.json.gz` | `results/certificates/` |
| 6 | 91000 | 0.364220·log2 N + 6 | `cert_k6_n91000.json.gz` | `results/certificates/` |
| 7 | 650000 | 0.362505·log2 N + 7 | `cert_k7_n650000.json.gz` (64,438,009 bytes) | release asset |
| 7 | 700000 | 0.360509·log2 N + 7 | `cert_k7_n700000.json.gz` (68,811,411 bytes) | release asset |

The coefficients are rounded up. The two seven-query certificates are too large for the repository. Download them
from the assets of <https://github.com/maxxfan/aixmathxphysics/releases/tag/ordered-search-v2> into
`results/certificates/`.

## Checking the files

`MANIFEST.sha256` pins all 23 released files: every checker, test, fixture and data file, and all nine certificates.
With the two assets in place, from this folder:

```bash
shasum -a 256 -c MANIFEST.sha256
```

## Re-checking the certificates

The frozen verifier (`verifier/verify.py`, version 1.1, pin `7d6fc293c48a7479` in `verifier/PINS.json`) uses only the
Python standard library and exact integer arithmetic. The two seven-query certificates take about 8 minutes each.

```bash
python3 verifier/tests/test_verifier.py
python3 verifier/verify.py results/certificates/cert_k7_n700000.json.gz
```

The independent checker shares no code with the verifier, and it checks every condition by a different route in
FLINT arithmetic. It needs Python with python-flint. It rejects bad input with explicit exceptions rather than
`assert`, so it behaves the same under `python -O`. It checks certificates and is not meant as a service for hostile
input.

```bash
python3 results/test_check_independent.py
python3 results/check_independent.py results/certificates/cert_k7_n700000.json.gz
```

`results/audits.json` records the audits of every certificate. Version 1 of this release is kept at tag
`ordered-search-v1`.
