Outputs of the replay runs of 10 November 2026 (release/replay.sh), one directory per run.  Each run was started from
a fresh copy of release/ in a scratch directory; absolute paths of the scratch copies have been removed from the
logs.  Machine: Apple M4 (4 + 6 cores), macOS 26.5, shared with other jobs (load 6 to 58 during the runs), Python
3.14.3, python-flint 0.9.0 (FLINT 3.6.0), NumPy 2.5.4, mpmath 1.4.1, SciPy 1.18.1, Apple clang 21.0.0.

quick/        ./replay.sh quick, all steps PASS (6 minutes, 4 processes).
extra/        NPROC=2 ./replay.sh extra.  Every verifier, 128-bit and 30-digit step is PASS identical.  The run as a
              whole reported FAIL because of the three make_perturbed lines: that version of replay.sh compared the
              regenerated perturbed inputs with the shipped ones byte by byte, and they differ in the free-text
              "note" field, which names the source file by its path (candidates/... in the project,
              certificates/board/... in the release).  Family, p and all amplitudes are identical.  replay.sh now
              compares those fields and restores the shipped files.  This version also did not yet contain the step
              arb_test_tail.
extra_rerun/  ONLY="make_perturbed arb_test_tail" ./replay.sh extra with the corrected script: PASS.
board/        NPROC=2 ./replay.sh board: the 128-bit evaluator on the 17 board records that are not final records.
              See the comment lines in board/SUMMARY.txt: the reference enclosures of these 17 inputs come from this
              run; what the lines establish is the sign and the containment in the frozen verifier's enclosures.
headline/     ./replay.sh headline: the three headline certificates, one after the other (process counts 3, 4, 3).

How the trees used by these runs differ from the released one (checked by comparing SHA-256 of every file):

headline/     Same three evaluator files (arb/arb_ic.py, rigcore.py, rigeig.py), inputs, verifier, validation suite.
              Earlier versions of replay.sh, check_outputs.py, build_release.py and expected/arb.json (the headline
              entries had the same enclosures); second/, search/, indep/run_final.py, arb/validation/test_tail.py,
              requirements.txt, README.md and the other expected/*.json files were added afterwards.  The three
              headline logs were re-checked with the released check_outputs.py and expected/arb.json: PASS identical.
extra/        Earlier versions of replay.sh, check_outputs.py, build_release.py, expected/arb.json and
              expected/text.json; arb/validation/test_tail.py and README.md added afterwards.
board/        expected/arb.json (the 17 board entries were added from this run) and README.md.
quick/, extra_rerun/ (make_perturbed, arb_test_tail)
              README.md and reference_out/ only.
extra_rerun/ (second_mp)
              README.md and reference_out/ only.

After reference_out/ and README.md were added, ./replay.sh quick was run once more on the final tree; its summary is in
the project at reviews/exam-20261109-evidence/final_quick_SUMMARY.txt (not part of the release, which cannot contain
the output of a run on itself).
