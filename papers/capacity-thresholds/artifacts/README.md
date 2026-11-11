# Certificates, evaluators and replay: certified thresholds for positive quantum capacity (artifact version 2)

Release for the paper *Certified thresholds for positive quantum capacity of the depolarizing, independent X-Z and
two-Pauli channels* (revision of 10 November 2026). Each certificate is an explicit rank-two permutation-symmetric
input on n qubits: two rows of amplitudes in the Dicke basis, taken as exact double-precision numbers, together with
an exact rational noise level p. A certificate proves Q > 0 by certifying that the coherent information I_c of the
input is positive.

**Status.** This is version 2, tag `capacity-thresholds-v2` of the site repository
(https://github.com/maxxfan/aixmathxphysics/tree/capacity-thresholds-v2/papers/capacity-thresholds/artifacts).
It supersedes the artifact published under tag `capacity-thresholds-v1`; see "Changes since capacity-thresholds-v1" at
the end.

## What is certified, and by what

| value | n | certified by | also evaluated by (not rigorous) |
|---|---|---|---|
| depolarizing, Q > 0 on [0, 0.06679861] | 480 | frozen verifier (board #30) and 128-bit evaluator | kernel evaluator (float64) |
| X-Z, Q > 0 on [0, 0.12464123] | 80 | frozen verifier (board #28) and 128-bit evaluator | kernel evaluator (float64 and 30 digits), master-equation float evaluator |
| two-Pauli, Q > 0 at 0.12367528 | 80 | frozen verifier (board #27) and 128-bit evaluator | kernel evaluator (float64 and 30 digits), master-equation float evaluator |
| **depolarizing, Q > 0 on [0, 0.06688201]** | 480 | **128-bit evaluator only** | kernel evaluator (float64) |
| **X-Z, Q > 0 on [0, 0.12770564]** | 320 | **128-bit evaluator only** | kernel evaluator (float64) |
| **two-Pauli, Q > 0 at 0.12677472** | 320 | **128-bit evaluator only** | kernel evaluator (float64) |

- The three values in bold are the headline values of the paper. They rest on one certifying program, `arb/arb_ic.py`
  with `arb/rigcore.py` and `arb/rigeig.py`. Its soundness argument is Appendix B of the paper; its validation suite
  is in `arb/validation/`.
- The frozen verifier cannot decide the headline inputs, for two different reasons. The X-Z and two-Pauli inputs
  (n = 320) are beyond its size limit (n <= 80 for these families). The depolarizing input (n = 480) is within its size
  limit (n <= 480) but below its resolution: the verifier returns **UNDECIDED, I_c in [-4.0218e-9, 2.8387e-9] bits**
  (`./replay.sh extra`, step `verifier_dep480_headline`), an interval that contains the 128-bit value 4.852e-10.
- The kernel evaluator (`second/`) and the master-equation float evaluator (`indep/`) return numbers, not enclosures.
  They corroborate the certified signs; they are not proofs.
- Depolarizing and X-Z are semigroups in p, so a certificate at p proves Q > 0 on [0, p]. A two-Pauli certificate
  holds at its own p. Parameter domains: depolarizing p in [0, 1/4], X-Z p in [0, 1/2], two-Pauli p in [0, 1/2].
- The depolarizing value is not the largest published lower bound: arXiv:2609.39747 claims 0.0676387 per Pauli
  (0.2029161 total) from an asymptotic argument with a numerical step. The values here are the largest certified at
  an explicit finite n.

## Layout

| path | contents |
|---|---|
| `certificates/board/` | the 21 board records, certified by the frozen verifier (table below) |
| `certificates/arb/` | inputs certified only by the 128-bit evaluator (n = 240, 320, 480), the n = 120 cross-check, the three points with margin, and our 45-copy input (`depolarizing_n45_0.06497.json`) |
| `verifier/` | the frozen verifier: `verify.py` and its ball-arithmetic engine `icball.c`, the tests, `PINS.json` and `freeze.sh` (version 2) |
| `data/records.json` | the record snapshot the verifier reads (thresholds, win bars, where Q = 0 is listed) |
| `data/sources/kg2026/`, `seeds/` | the two 45-copy states of Krohn-Grimberghe, arXiv:2608.15870 (CC BY 4.0, copied unmodified from doi:10.5281/zenodo.21962925), as raw payloads and as input files |
| `arb/` | the 128-bit evaluator: `arb_ic.py`, `rigcore.py`, `rigeig.py` |
| `arb/validation/` | its validation suite: `validate_brute.py` (252 brute-force cases, 26 + 9 fifty-digit cases via `hpbrute.py`), `test_eig.py` (spectral unit test), `validate_records.py` (Krohn-Grimberghe states; negative controls), `make_perturbed.py`, `compare_indep.py`, `test_tail.py` (test of the series-truncation bound, added 2026-11-10), and `testdata/` |
| `arb/logs/` | logs and JSON reports of the original runs (October 2026), including the validation logs |
| `indep/` | the independent floating-point evaluator `indep_ic.py` (same master-equation representation as the 128-bit evaluator, separate code), its tests and audit logs; `run_final.py` for the replay |
| `kgcheck/` | the certificate emitter and logs of the check with Krohn-Grimberghe's unmodified verifier (their verifier is not redistributed here; get it from their Zenodo record) |
| `second/` | the second evaluator `kernel_ic.py` (generating-function formula; NumPy and mpmath only), `validate_kernel.py`, `run_headline.py`, `run_mp.py`, and the logs of its first runs |
| `search/` | the search programs as they were run (`*.py`, `qlib/icfloat.py`) and their logs; not needed to check any certificate, and not packaged for reuse (paths refer to the project layout) |
| `zero_capacity/` | the semidefinite program for the antidegradability points (needs cvxpy; not part of the replay) |
| `audits.json` | the audit record of October 2026 |
| `expected/` | expected results of every replay step |
| `reference_out/` | the outputs of the reference replay runs of 10 November 2026 |
| `replay.sh`, `check_outputs.py`, `requirements.txt` | the replay script, its checker, pinned versions |
| `build_release.py` | regenerates the copied part of this directory from the project and rewrites `SHA256SUMS` |
| `SHA256SUMS` | checksums of every file here |

The logs and the audit record name the inputs by their original paths in the project: `candidates/` is
`certificates/board/` here; `work/lead/arbcands/`, `work/lead/cands/` and `work/lead/cand_dep_n45_06497.json` are in
`certificates/arb/`; `results/rigorous/` is `arb/` (scripts other than the three evaluator files are in
`arb/validation/`); `results/indep/` is `indep/`; `results/kgcheck/` is `kgcheck/`; `results/second/` is `second/`.

## Requirements

- bash and a C compiler (`cc`) for the frozen verifier, which otherwise needs only the Python standard library.
- Python with the packages of `requirements.txt` for everything else. Reference versions: Python 3.14.3,
  python-flint 0.9.0 (FLINT 3.6.0), NumPy 2.5.4, mpmath 1.4.1, SciPy 1.18.1; Apple clang 21.0.0 on arm64, macOS 26.5.
  python-flint 0.9.0 did not install on Python 3.9 in an external review; use a recent Python.
- Memory: about 1.2 GB for the verifier at n = 480, about 2.5 GB per process for the 128-bit evaluator at n = 480,
  1.5 GB for the kernel evaluator at n = 480.

## Replay

```bash
PYTHON=/path/to/python ./replay.sh quick      # or: extra | board | headline | full
```

Run it from a copy of this directory. It writes logs to `out/` and a summary to `out/SUMMARY.txt`, compares every
step with `expected/`, and exits with status 0 only if every step passed. `ONLY="step ..."` selects steps; `NPROC`
sets the number of processes of the 128-bit runs.

| mode | steps | time of the reference run |
|---|---|---|
| `quick` | `SHA256SUMS`; versions; build of the verifier engine; verifier tests; board record #20 with the verifier; spectral unit test, 252 + 26 brute-force cases, 30 early-truncation cases, Krohn-Grimberghe states with the 128-bit evaluator; board records #28 and #27 with the 128-bit evaluator; validation of the kernel evaluator; kernel evaluator on the headline inputs and the board records (float64); master-equation float evaluator on three board records | 6 minutes (4 processes) |
| `extra` | verifier on board records #28, #27, #30 and on the depolarizing headline input (UNDECIDED); 9 more fifty-digit cases; early truncation at n = 80; negative controls (p just above the zero crossing; one perturbed amplitude) and the three points with margin with the 128-bit evaluator; board record #20, the n = 120 and n = 240 X-Z inputs with the 128-bit evaluator; kernel evaluator at 30 digits on the n = 45, n = 80 and n = 120 inputs | 3.6 hours (NPROC=2, shared machine); the two n = 240 depolarizing runs take 18 minutes each and the X-Z n = 240 run 1.5 hours with 2 processes |
| `board` | the 128-bit evaluator on the 17 board records that are not final records; each enclosure must lie inside the frozen verifier's | 2.7 hours (NPROC=2, shared machine) |
| `headline` | the three headline values with the 128-bit evaluator, one command each: `headline_dep_n480`, `headline_xz_n320`, `headline_2p_n320` | 3.37 h, 1.28 h and 1.66 h (process counts 3, 4, 3; shared machine, the first run under heavy load) |
| `full` | all of the above | |

How a step is judged (`check_outputs.py`):

- A run of the 128-bit evaluator or of the verifier is `PASS identical` if the sign (or verdict) and the digits of the
  enclosure equal those of the reference run, `PASS consistent` if the sign is the same and the enclosure overlaps the
  reference one, and `FAIL` otherwise. With the pinned versions the digits were identical in the reference runs.
  Another FLINT version, C library or compiler may change the last digits of an enclosure; the certified statement is
  the sign.
- Where the frozen verifier's enclosure of the same input is known (all board records, two points with margin, the
  n = 120 input, the depolarizing headline input), a run of the 128-bit evaluator also fails unless its enclosure lies
  inside the verifier's.
- A run of the 128-bit evaluator also fails if its total block trace does not contain 1. (`arb_ic.py` prints this test
  and does not stop on it; the evaluator was left byte-identical to version 1, so the check is made here.)
- A verifier step fails if the hash of `verify.py`, `icball.c` or `records.json` differs from `verifier/PINS.json`.
  The hash of the compiled engine is machine-specific; a different one is reported, not failed.
- The non-rigorous evaluators are compared with the midpoint of the certified enclosure: within 5e-14 bits in double
  precision, within 1e-27 bits at 30 digits.

Single commands:

```bash
python3 verifier/tests/test_verifier.py
python3 verifier/verify.py certificates/board/20261027-053852_lead_9e078e22e51cd14b.json      # record #20, 30 s
python3 arb/arb_ic.py certificates/board/20261027-124857_lead_2a327ea83c2b10f7.json --nproc 4  # record #28, 30 s
python3 arb/arb_ic.py certificates/arb/depolarizing_n480_0.06688201.json --nproc 3             # headline, about 3 h
python3 arb/arb_ic.py certificates/arb/independent_xz_n320_0.12770564.json --nproc 4           # headline, 1.3 to 2.2 h
python3 arb/arb_ic.py certificates/arb/two_pauli_n320_0.12677472.json --nproc 3                # headline, 1.7 to 2.9 h
python3 second/kernel_ic.py certificates/arb/depolarizing_n480_0.06688201.json                 # 35 s, not rigorous
```

## Reference runs (10 November 2026)

All runs were started from a fresh copy of this directory on an Apple M4 (4 + 6 cores, macOS 26.5) that was shared with
other jobs. Their outputs are in `reference_out/`, with a `README.txt` that states exactly which scripts changed
between runs.

| run | outcome |
|---|---|
| `quick` | all steps PASS; every 128-bit and verifier enclosure identical to the reference digits |
| `extra` | every verifier, 128-bit and 30-digit step PASS identical: verifier on #28, #27, #30; verifier on the depolarizing headline input, UNDECIDED with I_c in [-4.0218e-9, 2.8387e-9]; negative controls NEGATIVE; points with margin; #20, X-Z n = 120 and n = 240 with the 128-bit evaluator; 7 inputs at 30 digits. The run as a whole reported FAIL because of three `make_perturbed` lines, a defect of the replay script of that moment (byte comparison of a free-text field), not of any evaluation; see `reference_out/README.txt` |
| `extra`, steps `make_perturbed`, `arb_test_tail` and `second_mp` only | PASS: the corrected script; the tail-bound test (34 cases); 8 inputs at 30 digits, now including the X-Z input at n = 120 |
| `board` | the 17 other board records with the 128-bit evaluator: all POSITIVE, each enclosure inside the frozen verifier's |
| `headline` | all three PASS with digits identical to the runs of October 2026: depolarizing n = 480, I_c in [485190991292415409011, 485190991292415409190]e-30 bits (3.37 h); X-Z n = 320, [1722045867045543156425, 1722045867045543156528]e-30 (1.28 h); two-Pauli n = 320, [1189800134543454507921, 1189800134543454509050]e-30 (1.66 h) |


## The 128-bit evaluator: identification

| file | SHA-256 |
|---|---|
| `arb/arb_ic.py` | `8a5887956e2b3814b500a3508815351aa6f7fa270d9fdbc8bd7b8f5e151fd5a1` |
| `arb/rigcore.py` | `3a0142b4f95fbf2b8d4c5422ef7cf7e351c5d198f47cd84519c1f5357a7f3f7d` |
| `arb/rigeig.py` | `4e84498eab2803d390eb150fab167e2177a40fc4fc433aae1dfc012fae56a286` |

These three files are byte-identical to those of version 1 and have not changed since their validation on
27 October 2026. Verifier pins (first 16 hexadecimal digits of the SHA-256): `verify.py` 4f007e83056ea3a1, `icball.c`
00b09294b432e40a, `records.json` 3901d412c50f84f7.

## Board records

| post | channel | n | p | I_c enclosure (bits) | file in `certificates/board/` |
|---|---|---|---|---|---|
| 2 | depolarizing | 150 | 0.06629494 | [1.5821e-9, 2.1906e-9] | `20261026-210916_lead_481ce3870c2d4e21.json` |
| 3 | independent X-Z | 48 | 0.12179852 | [9.9207e-9, 1.0197e-8] | `20261026-210926_lead_c9b944ed9a6f8a4a.json` |
| 4 | two-Pauli | 48 | 0.12086867 | [9.8282e-9, 1.0237e-8] | `20261026-210955_lead_ad42e025717074f4.json` |
| 6 | depolarizing | 240 | 0.06646039 | [2.9994e-9, 4.4447e-9] | `20261026-230622_lead_32025a4e90f6200f.json` |
| 9 | independent X-Z | 64 | 0.12331552 | [9.5442e-9, 1.0323e-8] | `20261026-231833_lead_cc6a44c7cf5414bb.json` |
| 10 | two-Pauli | 64 | 0.12235314 | [9.0366e-9, 1.0657e-8] | `20261026-232813_lead_9c563111d96087e5.json` |
| 11 | depolarizing | 240 | 0.06651279 | [1.0198e-11, 1.4597e-9] | `20261027-010137_lead_5a8d8b8895df26cb.json` |
| 12 | independent X-Z | 64 | 0.12332757 | [2.4427e-11, 7.8861e-10] | `20261027-010209_lead_ad2fcf9c940ceb22.json` |
| 13 | two-Pauli | 64 | 0.12236279 | [1.3766e-10, 1.7557e-9] | `20261027-010328_lead_713c6cf1f278172f.json` |
| 15 | depolarizing | 240 | 0.06662756 | [4.0998e-12, 2.1039e-9] | `20261027-033842_lead_38eff8cefeab3c5d.json` |
| 16 | two-Pauli | 80 | 0.1233958 | [5.1234e-11, 7.1850e-9] | `20261027-045237_lead_673b3e61d692d1a4.json` |
| 17 | independent X-Z | 80 | 0.12437928 | [8.7420e-12, 2.4725e-9] | `20261027-045237_lead_4675821f074f5b91.json` |
| 18 | depolarizing | 240 | 0.06667253 | [9.7537e-12, 2.1687e-9] | `20261027-050500_lead_ba11e08853949e86.json` |
| 20 | depolarizing | 240 | 0.0666847 | [9.3916e-12, 2.1094e-9] | `20261027-053852_lead_9e078e22e51cd14b.json` |
| 21 | two-Pauli | 80 | 0.1236419 | [6.1335e-11, 7.6129e-9] | `20261027-074257_lead_d1c9e883b77efe7c.json` |
| 22 | independent X-Z | 80 | 0.12459005 | [2.2419e-11, 2.6302e-9] | `20261027-074257_lead_834edf5841fffd54.json` |
| 24 | independent X-Z | 80 | 0.12463166 | [3.2608e-12, 2.6262e-9] | `20261027-101825_lead_46b0d3127de5b834.json` |
| 25 | two-Pauli | 80 | 0.12367023 | [2.2868e-12, 7.5983e-9] | `20261027-101825_lead_6b2ddb012ce45bff.json` |
| 27 | two-Pauli | 80 | 0.12367528 | [1.3223e-13, 7.5729e-9] | `20261027-124857_lead_585247a2e5ee8280.json` |
| 28 | independent X-Z | 80 | 0.12464123 | [3.5713e-11, 2.6056e-9] | `20261027-124857_lead_2a327ea83c2b10f7.json` |
| 30 | depolarizing | 480 | 0.06679861 | [2.4950e-11, 6.8873e-9] | `20261028-024810_lead_5f09c5df10a5804f.json` |

Records #2 to #28 were certified by version 1 of the verifier and #30 by version 2. Version 2 only raised the size limit
for the depolarizing channel from n = 240 to n = 480. The replay re-runs #20 (`quick`) and #27, #28, #30 (`extra`).

## Values certified by the 128-bit evaluator

| channel | n | p | I_c (bits, lower end) | file in `certificates/arb/` | second evaluation (float64, not rigorous) |
|---|---|---|---|---|---|
| depolarizing | 480 | 0.06688201 (3p = 0.20064603) | 4.85190991292415409011e-10 | `depolarizing_n480_0.06688201.json` | 4.851870977e-10; at p + 1.1e-5: -4.84e-11 |
| independent X-Z | 240 | 0.12737956 | 1.730665272872527453149e-9 | `independent_xz_n240_0.12737956.json` | 1.730666119e-9 |
| independent X-Z | 320 | 0.12770564 | 1.722045867045543156425e-9 | `independent_xz_n320_0.12770564.json` | 1.722044636e-9; at p + 1.1e-5: -1.716e-10 |
| two-Pauli | 320 | 0.12677472 | 1.189800134543454507921e-9 | `two_pauli_n320_0.12677472.json` | 1.189803864e-9; at p + 1.1e-5: -1.175e-10 |
| independent X-Z | 120 | 0.1259 | 2.7601677157813408892899e-8 | `independent_xz_n120_0.1259_board_path1.json` | |
| depolarizing | 240 | 0.06666451 | 2.3132493724068961211e-9 | `depolarizing_n240_0.06666451.json` | |
| independent X-Z | 80 | 0.12463577 | 2.6004754111423794411e-9 | `independent_xz_n80_0.12463577.json` | |
| two-Pauli | 80 | 0.12366746 | 5.6654186637896394375e-9 | `two_pauli_n80_0.12366746.json` | |

At 30 digits the kernel evaluator agrees with the 128-bit evaluator to within 3e-30 bits on the four certified points
of the Krohn-Grimberghe states, on our 45-copy input, on board records #28 and #27 and on the X-Z input at n = 120
(`reference_out/extra_rerun/second_mp.log`). The headline inputs are too large for its 30-digit mode (17 minutes at
n = 120, cost growing like n^4).

## Notes

- **Copied, sanitised and patched files.** `build_release.py` copies every file from the project. Evaluators,
  verifier, certificates, data and test data are byte-identical copies. Logs are copied with the absolute path of the
  project directory removed (the original runs printed absolute input paths). Three scripts that address inputs by
  their project paths have exactly those path strings replaced: `arb/validation/validate_records.py`,
  `arb/validation/make_perturbed.py`, `indep/run_audit.py`. `python3 build_release.py --check` (run inside the
  project) verifies that this directory equals what the script generates.
- **The verifier's "Q = 0" alarm.** The frozen record file lists the antidegradability point for the depolarizing
  channel (1/12) but the entanglement-breaking points for X-Z (0.293) and two-Pauli (1/4). The verifier accepts only
  p <= 1/4, so this alarm is live for the depolarizing channel, cannot fire for X-Z, and can fire for two-Pauli only at
  p = 1/4. The tighter points, (1 - 1/sqrt 2)/2 and 1/6, are in the paper (Appendix A); the frozen files were not changed.
- **`data/records.json`** attributes arXiv:2609.39747 to "Bhalerao, Leditzky et al."; the author list is Ahmed,
  Bhalerao, Lee, Leditzky, Leung, Schaeffer, Smith. The file is hash-pinned and was left as it is.
- **Not replayed here:** the frozen verifier on the other 17 board records (an external referee re-ran all 21 with it
  on 2026-11-09/10 and obtained the listed enclosures; `board` mode runs the 128-bit evaluator on those 17), the 128-bit
  run on board record #30 (n = 480; log of October 2026 in `arb/logs/`), the check with Krohn-Grimberghe's verifier
  (40 to 100 minutes; logs in `kgcheck/logs/`), and the search.

## Changes since capacity-thresholds-v1

Version 1 is `papers/capacity-thresholds/artifacts/` at tag `capacity-thresholds-v1` of the site repository
(67 files).

- **Unchanged:** all 66 data and code files of version 1 are present at the same paths with the same SHA-256
  (both evaluators, the verifier, the records, the 21 board certificates, the 8 inputs of `certificates/arb/`,
  `audits.json`, and the 24 files of `arb/logs/`). No certificate, input or certifying program changed.
- **Rewritten:** `README.md` (this file). It no longer says that all values of `certificates/arb/` are "beyond the
  verifier's size limits": the depolarizing n = 480 input is within the size limit and undecided by the verifier.
- **Added:**
  - `arb/validation/` (the validation suite of the 128-bit evaluator and its 122 test inputs) and the validation and
    remaining run logs in `arb/logs/`; one test in it is new, `test_tail.py` (2026-11-10), which cuts the series early
    and checks that the truncation bound still gives valid enclosures;
  - `indep/`, `kgcheck/`, `seeds/`, `certificates/arb/depolarizing_n45_0.06497.json`;
  - `second/` (a second, non-rigorous evaluator based on a different formula, written on 2026-11-10);
  - `search/` (search programs and logs), `zero_capacity/`;
  - `replay.sh`, `check_outputs.py`, `expected/`, `reference_out/`, `requirements.txt`, `build_release.py`;
  - new evidence: the 128-bit evaluator on the 17 board records it had not been run on (`reference_out/board/`; all
    positive and inside the verifier's enclosures), and the frozen verifier's UNDECIDED outcome on the depolarizing
    headline input (`reference_out/extra/verifier_dep480_headline.json`).
- **Paper:** the manuscript was revised on 2026-11-10 (parameter domains and the definition of the thresholds; the
  comparison with arXiv:2609.39747; the size-limit wording and the UNDECIDED outcome; an appendix with the soundness
  argument of both evaluators; the second evaluation; prior art; no AI byline). The PDF of version 1 is superseded.
- **Location:** version 2 is `papers/capacity-thresholds/artifacts/` at tag `capacity-thresholds-v2` of the same
  repository; the Data availability statement of the manuscript cites that tag.
