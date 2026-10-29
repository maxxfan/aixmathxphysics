# Certificates and checkers: certified thresholds for positive quantum capacity

Release for the paper *Certified thresholds for positive quantum capacity of the depolarizing, independent X-Z and
two-Pauli channels* (`../capacity-thresholds.pdf`). Each certificate is an explicit rank-two permutation-symmetric
input on n qubits: two rows of amplitudes in the Dicke basis, taken as exact double-precision numbers, together with
an exact rational noise level p. A certificate proves Q > 0 by certifying that the coherent information of the input
is positive.

## Layout

| path | contents |
|---|---|
| `certificates/board/` | the 21 board records, certified by the frozen verifier (table below) |
| `certificates/arb/` | inputs certified only by the 128-bit evaluator: beyond the verifier's size limits (n = 240, 320, 480), the n = 120 cross-check, and the three points with margin |
| `verifier/` | the frozen verifier: `verify.py` and its ball-arithmetic engine `icball.c`, the tests, `PINS.json` and `freeze.sh` (version 2) |
| `data/records.json` | the record snapshot the verifier reads (thresholds, win bars, where Q = 0 is proven) |
| `data/sources/kg2026/` | the two 45-copy states of Krohn-Grimberghe, arXiv:2608.15870, used by the tests (CC BY 4.0, copied unmodified from doi:10.5281/zenodo.21962925) |
| `arb/` | the independent 128-bit evaluator (`arb_ic.py`, `rigcore.py`, `rigeig.py`) and `logs/` of its runs quoted in the paper |
| `audits.json` | the audit record: Krohn-Grimberghe's verifier on our 45-copy input, the independent floating-point re-evaluation, and the 128-bit re-certification |
| `SHA256SUMS` | checksums of every file here |

The logs and the audit record name the inputs by their original paths in the project:
`candidates/` is `certificates/board/` here, and `work/lead/arbcands/` and `work/lead/cands/` are `certificates/arb/`.

## Re-running

The frozen verifier needs only Python 3 (standard library) and a C compiler. It builds `verifier/icball` on first use.

```bash
python3 verifier/tests/test_verifier.py
python3 verifier/verify.py certificates/board/20261027-053852_lead_9e078e22e51cd14b.json
```

The tests take about 40 s. The second command re-certifies board record #20 (depolarizing, n = 240, p = 0.0666847)
in about 30 s, with the enclosure I_c in [9.3916e-12, 2.1094e-9] bits. The n = 480 record takes about 10 min and
1.2 GB. The hash of the compiled engine in `PINS.json` is machine-specific, so on another machine `bash
verifier/freeze.sh` rebuilds the engine, re-runs the tests and re-pins it. The hashes of `verify.py`, `icball.c` and
`records.json` should still match the ones shipped.

The 128-bit evaluator needs Python with python-flint 0.9.0 and NumPy:

```bash
python3 arb/arb_ic.py certificates/board/20261027-124857_lead_2a327ea83c2b10f7.json
python3 arb/arb_ic.py certificates/arb/independent_xz_n320_0.12770564.json --nproc 4
```

The first command re-certifies record #28 (independent X-Z, n = 80) in about 15 s, with I_c =
1.605068705609815603826e-9 bits. The n = 320 and n = 480 inputs take 1 to 3 hours with 3 to 4 processes.

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
for the depolarizing channel from n = 240 to n = 480.

## Values certified by the 128-bit evaluator

| channel | n | p | I_c (bits, lower end) | file in `certificates/arb/` |
|---|---|---|---|---|
| depolarizing | 480 | 0.06688201 (3p = 0.20064603) | 4.85190991292415409011e-10 | `depolarizing_n480_0.06688201.json` |
| independent X-Z | 240 | 0.12737956 | 1.730665272872527453149e-9 | `independent_xz_n240_0.12737956.json` |
| independent X-Z | 320 | 0.12770564 | 1.722045867045543156425e-9 | `independent_xz_n320_0.12770564.json` |
| two-Pauli | 320 | 0.12677472 | 1.189800134543454507921e-9 | `two_pauli_n320_0.12677472.json` |
| independent X-Z | 120 | 0.1259 | 2.7601677157813408892899e-8 | `independent_xz_n120_0.1259_board_path1.json` |
| depolarizing | 240 | 0.06666451 | 2.3132493724068961211e-9 | `depolarizing_n240_0.06666451.json` |
| independent X-Z | 80 | 0.12463577 | 2.6004754111423794411e-9 | `independent_xz_n80_0.12463577.json` |
| two-Pauli | 80 | 0.12366746 | 5.6654186637896394375e-9 | `two_pauli_n80_0.12366746.json` |

The depolarizing and independent X-Z channels are semigroups in p, so a certificate at p proves Q > 0 on [0, p]. A
two-Pauli certificate holds only at its own p.
