"""Locate an input file in the release layout (certificates/, seeds/) or in the problem directory."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SEARCH = [
    ('..', 'certificates', 'arb'), ('..', 'certificates', 'board'), ('..', 'seeds'),          # release
    ('..', '..', 'work', 'lead', 'arbcands'), ('..', '..', 'work', 'lead', 'cands'),            # repository
    ('..', '..', 'candidates'), ('..', '..', 'seeds'), ('..', '..', 'work', 'lead'),
]
ALIAS = {'depolarizing_n45_0.06497.json': 'cand_dep_n45_06497.json'}


def find(name):
    for parts in SEARCH:
        for nm in (name, ALIAS.get(name, name)):
            path = os.path.join(HERE, *parts, nm)
            if os.path.exists(path):
                return os.path.normpath(path)
    raise FileNotFoundError(name)
