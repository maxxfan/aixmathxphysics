# Certificates and checkers: exact quantum ordered search (version 5)

Release `ordered-search-v5` for the paper *Exact certificates for quantum ordered search: seven queries for 850,000
elements and six for 91,000* (`../ordered-search.pdf`). Each certificate is an exact feasibility certificate for a
pair (k, n). It is an explicit chain of symmetric Laurent polynomials with integer coefficients that, through the
Farhi–Goldstone–Gutmann–Sipser characterization, proves that an exact k-query algorithm for ordered search on n
elements exists. The unitaries are not constructed. Applied recursively, such an algorithm searches N elements with
k⌈log_n N⌉ ≤ (k / log2 n)·log2 N + k queries.

## Changes since ordered-search-v4

No certificate changed, and the verifier and the independent checker are the same code as in versions 2 to 4.

- The paper was revised after a referee exam of 9 November 2026 (text, tables and references; no new result).
- This README gave "about 8 minutes each" as the verifier time for 650,000 and 700,000 elements. The audit record
  has 7.1 and 10.0 minutes; the times below are now the ones in `results/audits.json`.
- New in the release: the search and certificate-building programs (`qlib/`, `work/lead/`) and the 45 run logs
  (`results/logs/`). They are not needed to check a certificate.
- `MANIFEST.sha256` now pins 80 files instead of 25: the same 25, the 10 program files and the 45 logs. Of the 25,
  three changed: `results/audit.py` (docstring only), `data/records.json` (notes added) and
  `results/make_manifest.py` (longer file list).

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
| 7 | 775000 | 0.357804·log2 N + 7 | `cert_k7_n775000.json.gz` (75,503,378 bytes) | release asset |
| 7 | 850000 | 0.355383·log2 N + 7 | `cert_k7_n850000.json.gz` (82,519,946 bytes) | release asset |

The coefficients are rounded up. The four seven-query certificates are too large for the repository. Download them
from the assets of <https://github.com/maxxfan/aixmathxphysics/releases/tag/ordered-search-v5> into
`results/certificates/`. They are byte-identical to the assets of the earlier releases.

## Checking the files

`MANIFEST.sha256` pins 80 released files: every checker, test, fixture and data file, all eleven certificates, the
programs and the logs. With the four assets in place, from this folder:

```bash
shasum -a 256 -c MANIFEST.sha256
```

Without the assets the same command reports the four missing files and checks the other 76.
`results/certificates/SHA256SUMS` lists the eleven certificates alone.

## Re-checking the certificates

The frozen verifier (`verifier/verify.py`, version 1.1, pin `7d6fc293c48a7479` in `verifier/PINS.json`) uses only the
Python standard library and exact integer arithmetic. On the laptop of the paper it took 7.1, 10.0, 12.2 and 12.0
minutes for the seven-query certificates for 650,000, 700,000, 775,000 and 850,000 elements (425, 601, 730 and
719 s in `results/audits.json`; the runs for 700,000 and 775,000 shared the machine with another search). The seven
smaller certificates take 21 s or less each.

```bash
python3 verifier/tests/test_verifier.py
python3 verifier/verify.py results/certificates/cert_k7_n850000.json.gz
```

The independent checker shares no code with the verifier, and it checks every condition by a different route in
FLINT arithmetic. It needs Python with python-flint. It took 12 to 23 s for each seven-query certificate. It rejects
bad input with explicit exceptions rather than `assert`, so it behaves the same under `python -O`. It checks
certificates and is not meant as a service for hostile input.

```bash
python3 results/test_check_independent.py
python3 results/check_independent.py results/certificates/cert_k7_n850000.json.gz
```

`results/audits.json` records the audits of every certificate.

## Programs and logs

`qlib/` (FGGS chain arithmetic, interpolation, the certificate builder and its tests) and `work/lead/` (the linear
programs for small n and the smoothing search for large n) are the programs that found the chains and built the
certificates. They need numpy and scipy. `python3 qlib/test_certify.py` tests the builder. The programs write their
working files to `work/lead/out/`, which is not released.

`results/logs/*.log.txt` are copies of the logs of the searches, certificate builds and audits. The absolute path of
the checkout was removed from them and nothing else was changed. The pass counts, times, spike counts and residuals
quoted in the paper come from these logs. The logs number passes from 0 and the paper counts them from 1.

Versions 1 to 4 of this release are kept at tags `ordered-search-v1` to `ordered-search-v4`.
