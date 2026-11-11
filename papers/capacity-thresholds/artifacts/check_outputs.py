#!/usr/bin/env python3
"""Compare one output of replay.sh with the expected result in expected/.  Standard library only.

  check_outputs.py arb      NAME LOG    a log of arb/arb_ic.py           (expected/arb.json)
  check_outputs.py verifier NAME JSON   `verifier/verify.py --json` out  (expected/verifier.json)
  check_outputs.py text     NAME LOG    a validation log                 (expected/text.json)
  check_outputs.py second   NAME LOG    a log of second/kernel_ic.py     (expected/second.json)
  check_outputs.py pert     A B         two perturbed-input files: same family, p, amplitudes
  check_outputs.py env                  versions against requirements.txt (never fails; prints)

Prints one line starting with PASS or FAIL; exit status 0 or 1.

Verdicts for arb and verifier runs:
  PASS identical   same sign/verdict and the same enclosure digits as the reference run
  PASS consistent  same sign/verdict, the enclosure overlaps the reference one, but the digits
                   differ (another FLINT version, C library or compiler can change the last digits
                   of a rigorous enclosure; the certified statement is the sign)
  FAIL             anything else, including a block-trace test that does not contain 1
"""
import json
import os
import re
import sys
from decimal import Decimal

HERE = os.path.dirname(os.path.abspath(__file__))


def expected(kind, name):
    with open(os.path.join(HERE, 'expected', kind + '.json')) as f:
        d = json.load(f)
    if name not in d:
        print(f'FAIL {kind} {name}: no expected entry')
        sys.exit(1)
    return d[name]


def done(ok, msg):
    print(('PASS ' if ok else 'FAIL ') + msg)
    sys.exit(0 if ok else 1)


def check_arb(name, log):
    e = expected('arb', name)
    text = open(log).read()
    m = re.search(r'^family\s+(\S+)\s+n (\d+)\s+rank (\d+)\s+p (\S+)', text, re.M)
    x = re.search(r'^I_c exact \[(-?\d+), (-?\d+)\] x 1e-30', text, re.M)
    s = re.search(r'^sign\s+(\S+)', text, re.M)
    t = re.search(r'\(contains 1: (\w+)\)', text)
    if not (m and x and s and t):
        done(False, f'arb {name}: log incomplete ({log})')
    fam, n, p = m.group(1), int(m.group(2)), m.group(4)
    lo, hi = int(x.group(1)), int(x.group(2))
    sign = s.group(1)
    head = f'arb {name}: {fam} n={n} p={p} sign {sign} I_c in [{lo}, {hi}]e-30 bits'
    if (fam, n, p) != (e['target'], e['n'], e['p']):
        done(False, head + f'; wrong input (expected {e["target"]} n={e["n"]} p={e["p"]})')
    if t.group(1) != 'True':
        done(False, head + '; total block trace does not contain 1')
    if sign != e['sign']:
        done(False, head + f'; expected sign {e["sign"]}')
    if 'verifier_lo' in e:
        # the 128-bit enclosure must lie inside the frozen verifier's enclosure of the same input
        vlo, vhi = Decimal(e['verifier_lo']), Decimal(e['verifier_hi'])
        if not (vlo <= Decimal(lo) * Decimal('1e-30') and Decimal(hi) * Decimal('1e-30') <= vhi):
            done(False, head + f'; not inside the verifier enclosure [{e["verifier_lo"]}, {e["verifier_hi"]}]')
        head += f'; inside the verifier enclosure [{e["verifier_lo"]}, {e["verifier_hi"]}]'
    elo, ehi = int(e['lo_e30']), int(e['hi_e30'])
    if (lo, hi) == (elo, ehi):
        done(True, 'identical   ' + head)
    if lo <= ehi and elo <= hi:
        done(True, 'consistent  ' + head + f'; reference [{elo}, {ehi}]e-30')
    done(False, head + f'; disjoint from the reference [{elo}, {ehi}]e-30')


def check_verifier(name, path):
    e = expected('verifier', name)
    try:
        rep = json.load(open(path))
    except Exception as ex:
        done(False, f'verifier {name}: unreadable report ({ex})')
    with open(os.path.join(HERE, 'verifier', 'PINS.json')) as f:
        pins = json.load(f)
    src = rep.get('sources', {})
    bad_src = [k for k in ('verify.py', 'icball.c', 'records.json') if src.get(k) != pins.get(k)]
    verdict = rep.get('verdict')
    if 'I_c_lo' not in rep:
        done(False, f'verifier {name}: verdict {verdict}: {rep.get("reason", "")}')
    lo, hi = float(rep['I_c_lo']), float(rep['I_c_hi'])
    got = (f'{lo:.4e}', f'{hi:.4e}')
    head = (f'verifier {name}: {rep.get("verifier")} {rep["target"]} n={rep["n"]} p={rep["p"]} '
            f'{verdict} I_c in [{got[0]}, {got[1]}] bits')
    if src.get('icball') == pins.get('icball'):
        binary = ''
    elif 'icball' not in src:
        binary = '; note: the engine was compiled during this run'
    else:
        binary = '; note: the hash of the compiled engine differs from the pinned one (it depends on machine and compiler)'
    if bad_src:
        done(False, head + '; SOURCE PINS DIFFER: ' + ','.join(bad_src))
    if (rep['target'], rep['n'], str(rep['p'])) != (e['target'], e['n'], e['p_fraction']):
        done(False, head + '; wrong input')
    if verdict != e['verdict']:
        done(False, head + f'; expected verdict {e["verdict"]}')
    if got == (e['lo'], e['hi']):
        done(True, 'identical   ' + head + binary)
    if lo <= float(e['hi']) and float(e['lo']) <= hi:
        done(True, 'consistent  ' + head + f'; reference [{e["lo"]}, {e["hi"]}]' + binary)
    done(False, head + f'; disjoint from the reference [{e["lo"]}, {e["hi"]}]')


def check_text(name, log):
    e = expected('text', name)
    text = open(log).read()
    problems = []
    for pat, count in e.get('require', []):
        k = len(re.findall(pat, text, re.M))
        if k != count:
            problems.append(f'/{pat}/ found {k} times, expected {count}')
    for pat in e.get('forbid', []):
        k = len(re.findall(pat, text, re.M))
        if k:
            problems.append(f'/{pat}/ found {k} times, expected 0')
    done(not problems, f'text {name}: ' + ('; '.join(problems) if problems else e.get('says', 'ok')))


def check_second(name, log):
    """second/kernel_ic.py result lines:  RESULT <tag> target n p mode I_c"""
    e = expected('second', name)
    text = open(log).read()
    rows = {}
    for m in re.finditer(r'^RESULT (\S+) (\S+) n=(\d+) p=(\S+) (\S+) I_c=(\S+)', text, re.M):
        rows[m.group(1)] = (m.group(2), int(m.group(3)), m.group(4), m.group(5), Decimal(m.group(6)))
    problems = []
    notes = []
    for tag, want in e['results'].items():
        if tag not in rows:
            problems.append(f'{tag}: missing')
            continue
        fam, n, p, mode, val = rows[tag]
        if (fam, n, p) != (want['target'], want['n'], want['p']):
            problems.append(f'{tag}: wrong input {fam} n={n} p={p}')
            continue
        if (val > 0) != (want['sign'] == 'POSITIVE'):
            problems.append(f'{tag}: I_c = {val:.6e}, expected sign {want["sign"]}')
            continue
        if 'reference' in want:
            ref = Decimal(want['reference'])
            tol = Decimal(want['tolerance'])
            if abs(val - ref) > tol:
                problems.append(f'{tag}: I_c = {val:.10e} differs from the reference {ref:.10e} by '
                                f'{abs(val - ref):.1e} > {tol:.0e}')
                continue
            notes.append(f'{tag} {val:+.6e} (ref {abs(val - ref):.0e} away)')
        else:
            notes.append(f'{tag} {val:+.4e}')
    done(not problems, f'second {name}: ' + ('; '.join(problems) if problems else ', '.join(notes)))


def check_pert(shipped, regenerated):
    a, b = json.load(open(shipped)), json.load(open(regenerated))
    same = all(a.get(k) == b.get(k) for k in ('target', 'p', 'rows'))
    done(same, 'pert ' + os.path.basename(shipped) + (': regenerated input has the same family, p and amplitudes '
         'as the shipped one' if same else ': regenerated input DIFFERS from the shipped one'))


def check_env():
    import platform
    print('python', platform.python_version(), platform.machine(), platform.system())
    for mod in ('numpy', 'flint', 'mpmath', 'scipy'):
        try:
            m = __import__(mod)
            extra = ''
            if mod == 'flint':
                extra = ' (FLINT library %s)' % getattr(m, '__FLINT_VERSION__', '?')
            print(f'{mod} {getattr(m, "__version__", "?")}{extra}')
        except Exception as ex:
            print(f'{mod} MISSING ({type(ex).__name__})')
    req = os.path.join(HERE, 'requirements.txt')
    if os.path.exists(req):
        print('reference versions (requirements.txt): ' +
              ', '.join(ln.strip() for ln in open(req) if ln.strip() and not ln.startswith('#')))


def main():
    a = sys.argv[1:]
    if a[:1] == ['env']:
        check_env()
    elif len(a) == 3 and a[0] == 'arb':
        check_arb(a[1], a[2])
    elif len(a) == 3 and a[0] == 'verifier':
        check_verifier(a[1], a[2])
    elif len(a) == 3 and a[0] == 'text':
        check_text(a[1], a[2])
    elif len(a) == 3 and a[0] == 'pert':
        check_pert(a[1], a[2])
    elif len(a) == 3 and a[0] == 'second':
        check_second(a[1], a[2])
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == '__main__':
    main()
