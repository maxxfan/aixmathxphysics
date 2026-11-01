#!/usr/bin/env python3
"""Negative controls for the replay (run by `replay.py --controls`, or directly).

Each control feeds a deliberately wrong input to the checks of replay.py and must be REJECTED:
tampered SOS identities and lambdas, an SOS certificate whose identity still holds but whose Y is indefinite,
an understated shift for a genuine certificate, tampered strategies, tampered inequality coefficients,
tampered printed windows, a tampered file, and adversarial inputs to the residual kernel rescheck.c
(indefinite matrices, overflow, and a randomized soundness battery with subnormal and large entries whose
every acceptance is re-proved by exact Bareiss elimination).  Ground truth for the PSD controls is
established by exact Bareiss elimination, so a control only counts if the claim it makes is really false.
"""
import copy
import json
import pathlib
import random
import sys
from fractions import Fraction as Q

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "verifier"))
sys.path.insert(0, str(ROOT / "replay"))
import exactcheck as X  # noqa: E402
import replay as RP     # noqa: E402
import verify as V      # noqa: E402


def load(rel):
    return json.loads((ROOT / rel).read_text())


def dense_to_upper(Y):
    return [[str(Y[i][j]) for j in range(i, len(Y))] for i in range(len(Y))]


def verifier_rejects(cand):
    try:
        rep = V.check_candidate(cand, lock=False)
    except Exception as e:  # noqa: BLE001
        return True, f"{type(e).__name__}: {e}"[:160]
    return rep.get("verdict") in ("INVALID", "VERIFIER_ALARM"), f"{rep.get('verdict')}: {rep.get('reason', '')}"[:160]


def exact_rejects(cand, cat, sigma=None, tau=None):
    try:
        out = X.check(cand, cat, sigma, tau)
    except X.Bad as e:
        return True, f"Bad: {e}"[:160]
    if tau is not None and not out.get("tau_ok", True):
        return True, f"tau_exact {out['tau_exact_float']:.3e} > stated {float(tau):.3e}"
    if out.get("bareiss_at_stated_tau") is False:
        return True, "Bareiss: not positive definite at the stated shift"
    return False, "accepted"


def pick_upper(man, need_words=("", "A1"), max_n=120, need_sigma=False):
    for e in man["inequalities"]:
        u = e["upper"]
        if u and u["n_words"] <= max_n and (not need_sigma or u["sigma"] > 0):
            cand = load(u["file"])
            if all(w in cand["words"] for w in need_words):
                return e, cand
    raise RuntimeError("no suitable certificate")


def run_all(man):
    cat = RP.catalog()
    out = []

    def add(name, rejected, detail):
        out.append({"control": name, "rejected": bool(rejected), "detail": detail})

    # ---------------------------------------------------------------- SOS certificates
    e, cand = pick_upper(man)
    D, Y = X.load_Y(cand["Y"], len(cand["words"]))
    c = copy.deepcopy(cand)
    Y1 = [row[:] for row in Y]
    Y1[0][1] += 1
    Y1[1][0] += 1
    c["Y"] = {"denominator": str(D), "upper": dense_to_upper(Y1)}
    r1, d1 = verifier_rejects(c)
    r2, d2 = exact_rejects(c, cat)
    add(f"SOS identity tampered (Y_01 + 1/D), {e['name']} N={e['upper']['n_words']}", r1 and r2, f"verifier {d1}; independent {d2}")

    c = copy.deepcopy(cand)
    c["lambda"] = str(Q(e["upper"]["lambda_exact"]) - Q(1, D))
    r1, d1 = verifier_rejects(c)
    r2, d2 = exact_rejects(c, cat)
    add(f"stated lambda lowered by 1/D, {e['name']}", r1 and r2, f"verifier {d1}; independent {d2}")

    # identity still holds, Y indefinite: add t (E_0a + E_a0) - 2t E_aa with w_0 = 1, w_a = A1 (contributes 2t A1 - 2t A1 = 0)
    i0, ia = cand["words"].index(""), cand["words"].index("A1")
    Y2 = [row[:] for row in Y]
    t = 8 * D
    Y2[i0][ia] += t
    Y2[ia][i0] += t
    Y2[ia][ia] -= 2 * t
    c = copy.deepcopy(cand)
    c["Y"] = {"denominator": str(D), "upper": dense_to_upper(Y2)}
    lam, _, _, _ = X.sos_lambda(c, cat[e["name"]])
    tau = Q(e["upper"]["tau_stated_exact"])
    truth = not X.bareiss_pd(Y2, D, tau)
    r1, d1 = verifier_rejects(c)
    r2, d2 = exact_rejects(c, cat, float.fromhex(e["upper"]["sigma_hex"]), tau)
    add(f"identity intact but Y indefinite (gauge direction, t = 8), {e['name']}",
        truth and r1 and r2 and lam == Q(e["upper"]["lambda_exact"]),
        f"identity still holds, lambda unchanged; Bareiss: not PD at the stated tau = {truth}; verifier {d1}; independent {d2}")

    # genuine certificate, understated shift (only counts if Bareiss proves the understated claim false)
    e, cand = pick_upper(man, need_words=(), need_sigma=True)
    D, Y = X.load_Y(cand["Y"], len(cand["words"]))
    tau = Q(e["upper"]["tau_stated_exact"])
    done = False
    for tp in (tau / 10, tau / 100, tau / 1000, Q(0)):
        if X.bareiss_pd(Y, D, tp):
            continue
        c = copy.deepcopy(cand)
        r, d = exact_rejects(c, cat, float(tp / 2), tp)
        add(f"understated shift tau' = {float(tp):.2e} (stated {float(tau):.2e}; Bareiss: Y + tau' I not PD), {e['name']}", r, d)
        done = True
        break
    if not done:
        add(f"understated shift, {e['name']}", False, "no understated shift with a provably false claim was found")

    # ---------------------------------------------------------------- strategies and windows (J97)
    ej = next(x for x in man["inequalities"] if x["name"] == "J97")
    cj = load(ej["lower"]["file"])
    c = copy.deepcopy(cj)
    c["psi"][0] = str(Q(c["psi"][0]) + Q(1, 1000))
    rep = V.check_candidate(c, lock=False)
    bad = RP.compare_lower(ej, c, rep, X.check(c, cat))
    add("strategy tampered (psi_0 + 1/1000), J97", bool(bad), "; ".join(bad)[:200])

    cat2 = copy.deepcopy(cat)
    cat2["J97"]["bAB"][0][0] += 1
    bad = RP.compare_lower(ej, cj, {"violation": ej["lower"]["violation_down10"]}, X.check(cj, cat2))
    add("inequality coefficient tampered (bAB_11 + 1): strategy value, J97", bool(bad), "; ".join(bad)[:200])
    r, d = exact_rejects(load(ej["upper"]["file"]), cat2)
    add("inequality coefficient tampered (bAB_11 + 1): SOS identity, J97", r, d)

    e2 = copy.deepcopy(ej)
    e2["window"]["lower"] = f"{float(Q(ej['window']['lower']) + Q(1, 10 ** 10)):.10f}"
    bad = RP.compare_lower(e2, cj, V.check_candidate(cj, lock=False), X.check(cj, cat))
    add(f"printed lower raised by 1e-10 ({ej['window']['lower']} -> {e2['window']['lower']}), J97", bool(bad), "; ".join(bad)[:200])

    e3 = copy.deepcopy(ej)
    e3["window"]["upper"] = f"{float(Q(ej['window']['upper']) - Q(1, 10 ** 10)):.10f}"
    cu = load(ej["upper"]["file"])
    ex = X.check(cu, cat, float.fromhex(ej["upper"]["sigma_hex"]), Q(ej["upper"]["tau_stated_exact"]))
    bad, _ = RP.compare_upper(e3, cu, V.check_candidate(cu, lock=False), ex, Q(cat["J97"]["b0"]))
    add(f"printed upper lowered by 1e-10 ({ej['window']['upper']} -> {e3['window']['upper']}), J97", bool(bad), "; ".join(bad)[:200])

    raw = (ROOT / ej["lower"]["file"]).read_bytes()
    k = raw.index(b"1", 200)
    tampered = raw[:k] + b"2" + raw[k + 1:]
    import hashlib
    add("certificate file tampered (one digit), J97 lower", hashlib.sha256(tampered).hexdigest() != ej["lower"]["sha256"],
        "sha256 differs from the manifest")

    # ---------------------------------------------------------------- residual kernel rescheck.c
    V.ensure_rescheck()
    Yi, Di = [[1, 2], [2, 1]], 1                     # eigenvalues 3 and -1
    for sg in (0.0, 0.5, 0.9):
        f = V.rescheck(2, V.packed_lower(Yi, Di), sg)
        g = X.exact_shift(Yi, Di, sg)
        add(f"kernel: indefinite [[1,2],[2,1]] (lambda_min = -1) at sigma = {sg}", f is None and (g is None or g >= 1),
            f"rescheck {'FAIL' if f is None else 'OK %r' % f}; exact {'no factor' if g is None else float(g)}")
    big = 2 ** 670                                       # entries 2^670: the squared residual overflows binary64
    f = V.rescheck(2, V.packed_lower([[3 * big, big], [big, 3 * big]], 1), 0.0)
    add("kernel: residual overflow (entries 2^670) is refused, not returned", f is None, f"rescheck {'FAIL' if f is None else f}")

    rng = random.Random(20261031)
    trials, unsound, accepted = 0, [], 0
    for trial in range(300):
        n = rng.randint(2, 8)
        kind = rng.choice(["gram", "gram_minus", "random"])
        scale = rng.choice([0, 0, -1000, -1050, 500])
        r = rng.randint(1, n)
        if kind == "random":
            M = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1):
                    M[i][j] = M[j][i] = rng.randint(-9, 9)
            den = 1
        else:
            B = [[rng.randint(-5, 5) for _ in range(n)] for _ in range(r)]
            M = [[sum(B[k][i] * B[k][j] for k in range(r)) for j in range(n)] for i in range(n)]
            den = 1
            if kind == "gram_minus":
                den = 2 ** 20
                M = [[M[i][j] * den - (1 if i == j else 0) for j in range(n)] for i in range(n)]
        if scale < 0:
            den *= 2 ** (-scale)
        elif scale > 0:
            M = [[x * 2 ** scale for x in row] for row in M]
        for sg_rel in (0.0, 2.0 ** -40, 2.0 ** -8):
            sigma = sg_rel * 2.0 ** scale
            trials += 1
            f = V.rescheck(n, V.packed_lower(M, den), sigma)
            if f is not None:
                accepted += 1
                if not X.bareiss_pd(M, den, Q(sigma) + Q(f)):
                    unsound.append(f"rescheck trial {trial} ({kind}, n={n}, scale 2^{scale}, sigma {sigma})")
            g = X.exact_shift(M, den, sigma)
            if g is not None and not X.bareiss_pd(M, den, g):
                unsound.append(f"exact_shift trial {trial} ({kind}, n={n}, scale 2^{scale})")
    add(f"kernel soundness battery: {trials} random matrices (singular Gram, indefinite, subnormal 2^-1000/2^-1050 "
        f"and large 2^500 scales); every accepted shift re-proved by exact Bareiss", not unsound,
        f"rescheck accepted {accepted}/{trials}; unsound acceptances: {unsound[:3] or 'none'}")

    rejected = sum(1 for o in out if o["rejected"])
    return {"total": len(out), "rejected": rejected, "controls": out,
            "problems": [o["control"] for o in out if not o["rejected"]]}


if __name__ == "__main__":
    res = run_all(RP.manifest())
    for o in res["controls"]:
        print(("REJECTED  " if o["rejected"] else "ACCEPTED! ") + o["control"] + "\n          " + o["detail"])
    print(f"{res['rejected']}/{res['total']} controls rejected")
    sys.exit(0 if not res["problems"] else 1)
