#!/usr/bin/env python3
"""Emit a Krohn-Grimberghe endpoint certificate for an n = 45 rank-two depolarizing candidate.

Certificate format: schema "physics-proofs/flattened-rank-two-endpoint/v1" of the CC BY 4.0
artifact of arXiv:2608.15870 (artifact/verifier/verify.py).  That verifier rebuilds every
symmetric-subspace block exactly from the payload, so a certificate carries only two kinds of
content:

  * hints that it re-checks: per block and per side (output B, joint RB) a dyadic Hermitian
    centre matrix C ("center_dyadic") and a dyadic approximate eigenbasis V
    ("candidate_basis");
  * replay fields that it recomputes and compares for exact equality: Frobenius bounds, polar
    bounds, Gershgorin discs, merged components, entropy enclosures, scaled sums, traces and the
    POS/NEG/UNRESOLVED classification.

This emitter computes the hints itself (C = the exact block centre rounded to a dyadic grid of
about 2^-kc relative to the largest entry; V = eigenvectors of C from an mpmath eigh at `prec`
bits, rounded to 2^-kv) and fills the replay fields with the verifier's own exact helper
functions, loaded from the verify.py given by --verifier (stdlib-only module, imported as a
library, never modified).  The cubic matrix products are done in exact integer arithmetic
instead of Fractions; exact rationals are unique, so the values are the ones the verifier gets.

usage:
  emit_kg_certificate.py CANDIDATE [P] --verifier .../verify.py --out CERT.json [--jobs 8]

CANDIDATE is a p08 candidate JSON ("rows": 2 rows of 46 [re, im] numbers in the Dicke basis,
or "payload_c128_b64" with "rank": 2) or a raw 1472-byte .c128 payload.  P is an exact
rational "a/b" or a decimal string; it defaults to the candidate's "p".  The payload sha256
and the exact p to pass to verify.py are printed at the end.

Needs numpy and mpmath (the repo .venv).  Run with python -I.

Measured 2026-10-26 (--jobs 7-8): about 8.5 min and 26 MB per certificate; the unmodified
verify.py then takes about 40 min (Python 3.14), dominated by its exact block rebuild.
Calibration: KG payload A at 16239/250000 -> ACCEPT POS, coherent width 1.9e-27 bits (their
own float64-hint certificate: width 3.2e-9 bits).
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import importlib.util
import json
import math
import multiprocessing
import struct
import sys
import time
from fractions import Fraction
from pathlib import Path

N = 45
RANK = 2
SCHEMA = "physics-proofs/flattened-rank-two-endpoint/v1"
CONSTANTS = {
    "sqrt_bits": "224",
    "log_terms": "96",
    "log_grid_digits": "70",
    "accumulation_digits": "50",
}

KG = None  # the KG verify.py module (exact helper functions)


# ----------------------------------------------------------------------------- inputs

def load_kg(path: str):
    spec = importlib.util.spec_from_file_location("kg_verify_lib", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.N != N or module.SCHEMA != SCHEMA:
        raise SystemExit("verify.py is not the n = 45 v1 verifier this emitter targets")
    return module


def num(x) -> float:
    """Same number rules as the p08 verifier: JSON numbers, decimal or hex-float strings."""
    if isinstance(x, bool):
        raise ValueError("booleans are not amplitudes")
    if isinstance(x, (int, float)):
        v = float(x)
    elif isinstance(x, str):
        v = float.fromhex(x) if x.strip().lower().lstrip("+-").startswith("0x") else float(x)
    else:
        raise ValueError("amplitudes must be numbers or numeric strings")
    if not math.isfinite(v):
        raise ValueError("non-finite amplitude")
    return v


def candidate_payload(path: Path) -> tuple[bytes, str | None]:
    """Return (1472-byte little-endian complex128 payload, row-major 2 x 46; candidate p)."""
    raw = path.read_bytes()
    if path.suffix == ".c128":
        payload, p_text = raw, None
    else:
        cand = json.loads(raw)
        p_text = cand.get("p")
        if "payload_c128_b64" in cand:
            if cand.get("rank") != RANK:
                raise SystemExit("payload_c128_b64 candidates need rank 2")
            payload = base64.b64decode(cand["payload_c128_b64"], validate=True)
        else:
            rows = cand["rows"]
            if len(rows) != RANK:
                raise SystemExit("the KG verifier handles rank-two states only")
            values = []
            for row in rows:
                if len(row) != N + 1:
                    raise SystemExit("the KG verifier is hard-wired to n = 45 (46 Dicke amplitudes)")
                for z in row:
                    if isinstance(z, (list, tuple)) and len(z) == 2:
                        values += [num(z[0]), num(z[1])]
                    else:
                        values += [num(z), 0.0]
            payload = struct.pack("<%dd" % len(values), *values)
    if len(payload) != 2 * (N + 1) * 16:
        raise SystemExit("payload must be 2 x 46 complex128 = 1472 bytes")
    return payload, p_text


def parse_p(text: str) -> Fraction:
    p = Fraction(text)
    if not 0 < p <= Fraction(1, 4):
        raise SystemExit("need 0 < p <= 1/4")
    return p


# ----------------------------------------------------------------------------- JSON encoders

def q_json(value: Fraction) -> list[str]:
    value = Fraction(value)
    return [str(value.numerator), str(value.denominator)]


def dyad_json(mantissa: int, exponent: int) -> list[str]:
    """Canonical dyadic mantissa * 2^exponent: odd mantissa, or ["0", "0"]."""
    if mantissa == 0:
        return ["0", "0"]
    while mantissa % 2 == 0:
        mantissa //= 2
        exponent += 1
    return [str(mantissa), str(exponent)]


def gd_json(re_int: int, im_int: int, exponent: int) -> list[str]:
    return dyad_json(re_int, exponent) + dyad_json(im_int, exponent)


# ----------------------------------------------------------------------------- hints

def round_center(center, kc: int):
    """Hermitian dyadic rounding of the exact centre: C = (Cr + i Ci) / 2^ec."""
    d = len(center)
    maxabs = max(max(abs(z[0]), abs(z[1])) for row in center for z in row)
    ec = kc - math.frexp(float(maxabs))[1]
    scale = 1 << ec
    cr = [[0] * d for _ in range(d)]
    ci = [[0] * d for _ in range(d)]
    for i in range(d):
        for j in range(i, d):
            re = round(center[i][j][0] * scale)
            im = round(center[i][j][1] * scale) if i != j else 0
            cr[i][j], ci[i][j] = re, im
            cr[j][i], ci[j][i] = re, -im
    return cr, ci, ec


def eigenbasis(cr, ci, kv: int, prec: int):
    """Eigenvectors (columns) of the integer Hermitian matrix Cr + i Ci, rounded to 2^-kv."""
    import mpmath

    d = len(cr)
    mpmath.mp.prec = prec
    shift = max(abs(x).bit_length() for row in cr + ci for x in row)
    is_complex = any(x != 0 for row in ci for x in row)
    a = mpmath.matrix(d, d)
    for i in range(d):
        for j in range(d):
            re = mpmath.ldexp(mpmath.mpf(cr[i][j]), -shift)
            if is_complex:
                a[i, j] = mpmath.mpc(re, mpmath.ldexp(mpmath.mpf(ci[i][j]), -shift))
            else:
                a[i, j] = re
    _, q = mpmath.eigh(a)
    unit = mpmath.ldexp(mpmath.mpf(1), kv)
    vr = [[int(mpmath.nint(mpmath.re(q[i, j]) * unit)) for j in range(d)] for i in range(d)]
    vi = [[int(mpmath.nint(mpmath.im(q[i, j]) * unit)) for j in range(d)] for i in range(d)]
    return vr, vi


# ----------------------------------------------------------------------------- exact replay

def _obj(matrix):
    import numpy as np

    out = np.empty((len(matrix), len(matrix[0])), dtype=object)
    for i, row in enumerate(matrix):
        for j, x in enumerate(row):
            out[i, j] = int(x)
    return out


def spectral_record(center, radius, cr, ci, ec, vr, vi, kv, dimension, specht):
    """Mirror of KG._verify_spectral that writes the record instead of checking it."""
    d = dimension
    sqrt_bounds = KG._sqrt_bounds
    Cr, Ci, Vr, Vi = _obj(cr), _obj(ci), _obj(vr), _obj(vi)

    # algebraic radius of the exact reference block
    algebraic_square = sum(value * value for row in radius for value in row)
    algebraic = sqrt_bounds(algebraic_square, bits=KG.SQRT_BITS)[1]

    # dyadic centre vs exact centre
    den_c = 1 << ec
    conversion_square = Fraction(0)
    for i in range(d):
        for j in range(d):
            dr = center[i][j][0] - Fraction(cr[i][j], den_c)
            di = center[i][j][1] - Fraction(ci[i][j], den_c)
            conversion_square += dr * dr + di * di
    conversion = sqrt_bounds(conversion_square, bits=KG.FROBENIUS_BITS)[1]
    total = algebraic + conversion

    # Gram defect of V (integer arithmetic, scale 2^(2 kv))
    Gr = Vr.T.dot(Vr) + Vi.T.dot(Vi)
    Gi = Vr.T.dot(Vi) - Vi.T.dot(Vr)
    one = 1 << (2 * kv)
    for i in range(d):
        Gr[i, i] -= one
    gram_square = Fraction(int((Gr * Gr + Gi * Gi).sum()), 1 << (4 * kv))
    gram_defect = sqrt_bounds(gram_square, bits=KG.FROBENIUS_BITS)[1]
    if not gram_defect < 1:
        raise RuntimeError("candidate basis too far from unitary")
    sqrt_minus_lower = sqrt_bounds(1 - gram_defect, bits=KG.SQRT_BITS)[0]
    sqrt_plus_upper = sqrt_bounds(1 + gram_defect, bits=KG.SQRT_BITS)[1]
    polar_distance = KG._qmax(1 - sqrt_minus_lower, sqrt_plus_upper - 1)
    basis_norm = sqrt_plus_upper
    center_square = Fraction(int((Cr * Cr + Ci * Ci).sum()), 1 << (2 * ec))
    center_norm = sqrt_bounds(center_square, bits=KG.FROBENIUS_BITS)[1]
    transform_error = total + polar_distance * center_norm * (1 + basis_norm)

    # transformed = V^H C V, scale 2^(ec + 2 kv)
    Wr = Cr.dot(Vr) - Ci.dot(Vi)
    Wi = Cr.dot(Vi) + Ci.dot(Vr)
    Tr = Vr.T.dot(Wr) + Vi.T.dot(Wi)
    Ti = Vr.T.dot(Wi) - Vi.T.dot(Wr)
    den_t = 1 << (ec + 2 * kv)

    rows = []
    discs = []
    max_radius = Fraction(0)
    for i in range(d):
        if Ti[i, i] != 0:
            raise RuntimeError("transformed diagonal not real")
        diagonal = Fraction(int(Tr[i, i]), den_t)
        offdiagonal = sum(
            KG._abs_upper((Fraction(int(Tr[i, j]), den_t), Fraction(int(Ti[i, j]), den_t)))
            for j in range(d)
            if j != i
        )
        radius_value = offdiagonal + d * transform_error
        max_radius = max(max_radius, radius_value)
        lower, upper = diagonal - radius_value, diagonal + radius_value
        rows.append({"offdiagonal_sum": q_json(offdiagonal), "disc": [q_json(lower), q_json(upper)]})
        discs.append((lower, upper))

    discs.sort(key=KG.cmp_to_key(KG._disc_cmp))
    merged = []
    for lower, upper in discs:
        if not merged or KG._qlt(merged[-1][1], lower):
            merged.append((lower, upper, 1))
        else:
            old_lower, old_upper, count = merged[-1]
            merged[-1] = (old_lower, KG._qmax(old_upper, upper), count + 1)

    components = []
    entropy_lower = Fraction(0)
    entropy_upper = Fraction(0)
    for lower, upper, count in merged:
        (clipped_lower, clipped_upper, term_lower, term_upper,
         at_lower, at_upper, critical) = KG._entropy_range(lower, upper)
        components.append({
            "interval": [q_json(lower), q_json(upper)],
            "multiplicity": str(count),
            "clipped": [q_json(clipped_lower), q_json(clipped_upper)],
            "endpoint_entropy": [
                [q_json(at_lower[0]), q_json(at_lower[1])],
                [q_json(at_upper[0]), q_json(at_upper[1])],
            ],
            "critical_cap_used": bool(critical),
            "entropy_range": [q_json(term_lower), q_json(term_upper)],
        })
        entropy_lower += count * term_lower
        entropy_upper += count * term_upper
    scaled_lower = KG._floor_scaled(specht * entropy_lower)
    scaled_upper = KG._ceil_scaled(specht * entropy_upper)

    record = {
        "dimension": str(d),
        "center_dyadic": {
            "shape": [str(d), str(d)],
            "entries": [gd_json(cr[i][j], ci[i][j], -ec) for i in range(d) for j in range(d)],
        },
        "candidate_basis": {
            "shape": [str(d), str(d)],
            "entries": [gd_json(vr[i][j], vi[i][j], -kv) for i in range(d) for j in range(d)],
        },
        "bounds": {
            "algebraic_frobenius": q_json(algebraic),
            "center_conversion_frobenius": q_json(conversion),
            "total_matrix_frobenius": q_json(total),
            "gram_defect_frobenius": q_json(gram_defect),
            "sqrt_one_minus_lower": q_json(sqrt_minus_lower),
            "sqrt_one_plus_upper": q_json(sqrt_plus_upper),
            "polar_distance": q_json(polar_distance),
            "basis_norm_upper": q_json(basis_norm),
            "center_norm_upper": q_json(center_norm),
            "transform_error": q_json(transform_error),
        },
        "rows": rows,
        "components": components,
        "scaled_entropy": [str(scaled_lower), str(scaled_upper)],
    }
    diag = {
        "gram_defect": float(gram_defect),
        "transform_error_rel": float(transform_error / center_norm),
        "max_disc_radius_rel": float(max_radius / center_norm),
        "components": len(merged),
        "max_multiplicity": max(c for _, _, c in merged),
        "width_bits": float(Fraction(scaled_upper - scaled_lower, KG.SCALE)),
    }
    return record, scaled_lower, scaled_upper, diag


# ----------------------------------------------------------------------------- per block

def _init_worker(verify_path: str) -> None:
    global KG
    KG = load_kg(verify_path)


def block_job(job):
    index, payload, p, kc, kv, prec, self_check = job
    started = time.time()
    state, norm = KG._exact_state(KG._decode_payload(payload))
    first_row = 23 + index
    degree = 2 * first_row - N
    specht = KG._specht_dimension(N, first_row)
    joint_center, joint_radius = KG._reference_block(state, norm, n=N, first_row=first_row, p=p)
    output_center, output_radius = KG._partial_trace(joint_center, joint_radius)
    t_reference = time.time() - started

    trace_center = sum(output_center[i][i][0] for i in range(degree + 1))
    trace_radius = sum(output_radius[i][i] for i in range(degree + 1))

    sides = {}
    diags = {}
    for name, center, radius in (
        ("output", output_center, output_radius),
        ("joint", joint_center, joint_radius),
    ):
        cr, ci, ec = round_center(center, kc)
        vr, vi = eigenbasis(cr, ci, kv, prec)
        record, lo, hi, diag = spectral_record(
            center, radius, cr, ci, ec, vr, vi, kv, len(center), specht)
        if self_check:
            # the verifier's own check of this record (Fraction products; fine for small blocks)
            got = KG._verify_spectral(center, radius, record, dimension=len(center),
                                      specht_dimension=specht)
            if got != (lo, hi):
                raise RuntimeError(f"self-check mismatch in block {index} {name}")
        sides[name] = (record, lo, hi)
        diags[name] = diag

    block = {
        "first_row": str(first_row),
        "degree": str(degree),
        "specht_dimension": str(specht),
        "trace": {"center": q_json(trace_center), "radius": q_json(trace_radius)},
        "output": sides["output"][0],
        "joint": sides["joint"][0],
    }
    return {
        "index": index,
        "block": block,
        "specht": specht,
        "trace": (trace_center, trace_radius),
        "output": sides["output"][1:],
        "joint": sides["joint"][1:],
        "diag": diags,
        "seconds": (round(t_reference, 1), round(time.time() - started, 1)),
        "self_checked": self_check,
    }


# ----------------------------------------------------------------------------- main

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("candidate", type=Path)
    ap.add_argument("p", nargs="?", help="exact p ('a/b' or decimal); default: candidate's p")
    ap.add_argument("--verifier", required=True, help="path to the KG artifact's verifier/verify.py")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--kc", type=int, default=144, help="centre grid: 2^-kc times the largest entry")
    ap.add_argument("--kv", type=int, default=112, help="eigenbasis grid 2^-kv")
    ap.add_argument("--prec", type=int, default=256, help="mpmath eigh precision (bits)")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--self-check-below", type=int, default=6,
                    help="replay the verifier's _verify_spectral on blocks with index below this")
    args = ap.parse_args(argv)

    global KG
    KG = load_kg(args.verifier)
    payload, p_candidate = candidate_payload(args.candidate)
    p_text = args.p if args.p is not None else p_candidate
    if p_text is None:
        raise SystemExit("no p given")
    p = parse_p(p_text)
    sha = hashlib.sha256(payload).hexdigest()
    KG._exact_state(KG._decode_payload(payload))  # rank-two and shape check
    print(f"payload sha256 {sha}  p = {p.numerator}/{p.denominator}", flush=True)

    jobs = [(i, payload, p, args.kc, args.kv, args.prec, i < args.self_check_below)
            for i in reversed(range(23))]
    results = {}
    started = time.time()
    ctx = multiprocessing.get_context("spawn")
    with ctx.Pool(args.jobs, initializer=_init_worker, initargs=(args.verifier,)) as pool:
        for res in pool.imap_unordered(block_job, jobs):
            results[res["index"]] = res
            d = res["diag"]
            print(f"block {res['index']:2d} first_row {23 + res['index']} "
                  f"t_ref {res['seconds'][0]}s t_all {res['seconds'][1]}s "
                  f"out width {d['output']['width_bits']:.3e} joint width {d['joint']['width_bits']:.3e} "
                  f"gram {d['joint']['gram_defect']:.2e} rad_rel {d['joint']['max_disc_radius_rel']:.2e} "
                  f"comps {d['output']['components']}/{d['joint']['components']} "
                  f"maxmult {d['output']['max_multiplicity']}/{d['joint']['max_multiplicity']}"
                  f"{' self-checked' if res['self_checked'] else ''}", flush=True)

    blocks = []
    output_lower = output_upper = joint_lower = joint_upper = 0
    global_center = Fraction(0)
    global_radius = Fraction(0)
    for i in range(23):
        res = results[i]
        blocks.append(res["block"])
        output_lower += res["output"][0]
        output_upper += res["output"][1]
        joint_lower += res["joint"][0]
        joint_upper += res["joint"][1]
        global_center += res["specht"] * res["trace"][0]
        global_radius += res["specht"] * res["trace"][1]
    if not (global_center - global_radius <= 1 <= global_center + global_radius):
        raise SystemExit("global trace enclosure misses 1")
    coherent_lower = output_lower - joint_upper
    coherent_upper = output_upper - joint_lower
    classification = ("POS" if coherent_lower > 0 else "NEG" if coherent_upper < 0
                      else "UNRESOLVED")
    document = {
        "schema": SCHEMA,
        "n": str(N),
        "payload_hex": payload.hex(),
        "payload_sha256": sha,
        "p": q_json(p),
        "constants": dict(CONSTANTS),
        "blocks": blocks,
        "global_trace": {"center": q_json(global_center), "radius": q_json(global_radius)},
        "scaled": {
            "output_lower": str(output_lower),
            "output_upper": str(output_upper),
            "joint_lower": str(joint_lower),
            "joint_upper": str(joint_upper),
            "coherent_lower": str(coherent_lower),
            "coherent_upper": str(coherent_upper),
        },
        "classification": classification,
    }
    KG._validate_document_shape(document)
    data = (json.dumps(document, allow_nan=False, ensure_ascii=True,
                       separators=(",", ":"), sort_keys=False).encode("ascii") + b"\n")
    if len(data) > KG.MAX_CERTIFICATE_BYTES:
        raise SystemExit("certificate too large")
    args.out.write_bytes(data)
    # what the verifier's loader does (rejects JSON numbers, re-checks shape and canonical bytes)
    if KG._load_certificate(args.out) != document:
        raise SystemExit("certificate does not reload identically")

    print(f"wrote {args.out} ({len(data)} bytes) in {time.time() - started:.0f}s")
    print(f"classification {classification}: coherent in "
          f"[{coherent_lower / KG.SCALE:.6e}, {coherent_upper / KG.SCALE:.6e}] bits "
          f"(scaled [{coherent_lower},{coherent_upper}]x1e-50)")
    print(f"verify with: python3 -I verify.py {args.out} {sha} {p.numerator}/{p.denominator}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
