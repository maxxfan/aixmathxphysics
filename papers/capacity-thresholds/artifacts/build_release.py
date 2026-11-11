#!/usr/bin/env python3
"""Regenerate the copied part of this release from the problem directory, and SHA256SUMS.

Run from the repository:   python3 problems/p08-capacity-thresholds/release/build_release.py
(standard library only).  The problem directory is the parent of release/.

What is authored here and NOT regenerated: README.md, replay.sh, check_outputs.py,
requirements.txt, build_release.py, expected/.  Everything else is a copy of a file of the
problem directory, by one of three rules:

  identical   byte-for-byte copy (evaluators, verifier, certificates, data, test data);
  sanitised   text copy with the absolute path of the problem directory removed (logs: the
              original runs printed absolute input paths);
  patched     text copy with the path substitutions of PATCHES below, and nothing else (four
              scripts that address inputs by their path in the problem directory).

`--check` regenerates into a temporary directory and compares it with what is here instead of
writing; exit status 1 on any difference.
"""
import argparse
import hashlib
import os
import re
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
# absolute path of the problem directory on the machine that produced a log
ABS_PREFIX = re.compile(r'/[^\s"\']*/problems/p08-capacity-thresholds/')
AUTHORED = {'README.md', 'replay.sh', 'check_outputs.py', 'requirements.txt', 'build_release.py'}
AUTHORED_DIRS = ('expected',)
NOT_RELEASED_DIRS = ('out', '__pycache__')

# released path of each input the patched scripts name by its problem-directory path
PATCHES = [
    ("'candidates/", "'certificates/board/"),
    ("'work/lead/cand_dep_n45_06497.json'", "'certificates/arb/depolarizing_n45_0.06497.json'"),
]
PATCH_INDEP = PATCHES + [
    ("PROB = os.path.dirname(os.path.dirname(HERE))", "PROB = os.path.dirname(HERE)"),
]

MARGIN_POINTS = ('depolarizing_n240_0.06666451.json', 'independent_xz_n80_0.12463577.json',
                 'two_pauli_n80_0.12366746.json')


def plan():
    """list of (rule, source relative to ROOT, destination relative to release/, patches)."""
    P = []

    def tree(src, dst, rule='identical', patches=None, skip=()):
        base = os.path.join(ROOT, src)
        for dp, dn, fn in os.walk(base):
            dn[:] = sorted(d for d in dn if d not in ('__pycache__',))
            for f in sorted(fn):
                if f.startswith('.') or f.endswith('.pyc') or f in skip:
                    continue
                rel = os.path.relpath(os.path.join(dp, f), base)
                P.append((rule, os.path.join(src, rel), os.path.join(dst, rel), patches))

    tree('candidates', 'certificates/board')
    tree('work/lead/arbcands', 'certificates/arb')
    for f in MARGIN_POINTS:
        P.append(('identical', 'work/lead/cands/' + f, 'certificates/arb/' + f, None))
    P.append(('identical', 'work/lead/cand_dep_n45_06497.json',
              'certificates/arb/depolarizing_n45_0.06497.json', None))
    for f in ('verify.py', 'icball.c', 'PINS.json', 'freeze.sh', 'tests/test_verifier.py'):
        P.append(('identical', 'verifier/' + f, 'verifier/' + f, None))
    P.append(('identical', 'data/records.json', 'data/records.json', None))
    tree('data/sources/kg2026', 'data/sources/kg2026')
    tree('seeds', 'seeds')
    for f in ('arb_ic.py', 'rigcore.py', 'rigeig.py'):
        P.append(('identical', 'results/rigorous/' + f, 'arb/' + f, None))
    for f in ('validate_brute.py', 'hpbrute.py', 'test_eig.py', 'compare_indep.py', 'test_tail.py'):
        P.append(('identical', 'results/rigorous/' + f, 'arb/validation/' + f, None))
    P.append(('patched', 'results/rigorous/validate_records.py', 'arb/validation/validate_records.py', PATCHES))
    P.append(('patched', 'results/rigorous/make_perturbed.py', 'arb/validation/make_perturbed.py', PATCHES[:1]))
    tree('results/rigorous/testdata', 'arb/validation/testdata')
    # queue_*.sh are the shell queues of the original runs (machine-specific paths); not released
    tree('results/rigorous/logs', 'arb/logs', rule='sanitised',
         skip=('queue_final.sh', 'queue_indep.sh', 'queue_pert.sh'))
    for f in sorted(os.listdir(os.path.join(ROOT, 'results/indep'))):
        if f.startswith('.') or f == '__pycache__':
            continue
        if f == 'run_audit.py':
            P.append(('patched', 'results/indep/' + f, 'indep/' + f, PATCH_INDEP))
        elif f.endswith('.txt'):
            P.append(('sanitised', 'results/indep/' + f, 'indep/' + f, None))
        else:
            P.append(('identical', 'results/indep/' + f, 'indep/' + f, None))
    tree('results/kgcheck', 'kgcheck', rule='sanitised')
    if os.path.isdir(os.path.join(ROOT, 'results/second')):
        tree('results/second', 'second', rule='sanitised')
    if os.path.isdir(os.path.join(ROOT, 'results/replay_reference')):
        tree('results/replay_reference', 'reference_out', rule='sanitised')
    # the search programs as run, with their logs (not needed to check any certificate)
    for f in sorted(os.listdir(os.path.join(ROOT, 'work/lead'))):
        if f.endswith('.py'):
            P.append(('sanitised', 'work/lead/' + f, 'search/' + f, None))
    P.append(('sanitised', 'qlib/icfloat.py', 'search/qlib/icfloat.py', None))
    tree('work/lead/logs', 'search/logs', rule='sanitised')
    P.append(('identical', 'results/zero_capacity_sdp.py', 'zero_capacity/zero_capacity_sdp.py', None))
    P.append(('identical', 'results/audits.json', 'audits.json', None))
    return P


def render(rule, src, patches):
    with open(os.path.join(ROOT, src), 'rb') as f:
        data = f.read()
    if rule == 'identical':
        return data
    try:
        text = data.decode('utf-8')
    except UnicodeDecodeError:
        return data
    text = ABS_PREFIX.sub('', text)
    if rule == 'patched':
        for old, new in patches:
            if old not in text:
                raise SystemExit(f'patch {old!r} does not apply to {src}')
            text = text.replace(old, new)
    return text.encode('utf-8')


def write_tree(dest):
    out = []
    for rule, src, dst, patches in plan():
        data = render(rule, src, patches)
        path = os.path.join(dest, dst)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if os.path.exists(path):
            os.chmod(path, 0o644)
        with open(path, 'wb') as f:
            f.write(data)
        mode = os.stat(os.path.join(ROOT, src)).st_mode
        os.chmod(path, 0o755 if (mode & 0o100) else 0o644)
        out.append((rule, src, dst))
    return out


def released_files(base):
    files = []
    for dp, dn, fn in os.walk(base):
        dn[:] = sorted(d for d in dn if d not in NOT_RELEASED_DIRS and not d.startswith('.'))
        for f in sorted(fn):
            rel = os.path.relpath(os.path.join(dp, f), base)
            if f.startswith('.') or f.endswith('.pyc') or rel == 'SHA256SUMS':
                continue
            if rel.startswith('verifier/icball') and not rel.endswith('.c'):
                continue     # the compiled engine is built on the user's machine
            files.append(rel)
    return sorted(files)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    if a.check:
        bad = 0
        with tempfile.TemporaryDirectory() as tmp:
            done = write_tree(tmp)
            for rule, src, dst in done:
                here = os.path.join(HERE, dst)
                if not os.path.exists(here) or sha256(here) != sha256(os.path.join(tmp, dst)):
                    print('DIFFERS', dst)
                    bad += 1
            known = {d for _, _, d in done}
            for rel in released_files(HERE):
                top = rel.split(os.sep)[0]
                if rel not in known and rel not in AUTHORED and top not in AUTHORED_DIRS:
                    print('NOT GENERATED AND NOT AUTHORED', rel)
                    bad += 1
        sums = {}
        with open(os.path.join(HERE, 'SHA256SUMS')) as f:
            for ln in f:
                h, rel = ln.rstrip('\n').split('  ', 1)
                sums[rel] = h
        for rel in released_files(HERE):
            if sums.get(rel) != sha256(os.path.join(HERE, rel)):
                print('SHA256SUMS STALE', rel)
                bad += 1
        for rel in sums:
            if not os.path.exists(os.path.join(HERE, rel)):
                print('SHA256SUMS LISTS A MISSING FILE', rel)
                bad += 1
        print('build_release --check:', 'OK' if not bad else f'{bad} problem(s)')
        sys.exit(1 if bad else 0)
    done = write_tree(HERE)
    files = released_files(HERE)
    with open(os.path.join(HERE, 'SHA256SUMS'), 'w') as f:
        for rel in files:
            f.write(f'{sha256(os.path.join(HERE, rel))}  {rel}\n')
    n = {r: sum(1 for x in done if x[0] == r) for r in ('identical', 'sanitised', 'patched')}
    print(f'release: {len(done)} files regenerated ({n["identical"]} identical, {n["sanitised"]} sanitised, '
          f'{n["patched"]} patched); SHA256SUMS lists {len(files)} files')


if __name__ == '__main__':
    main()
