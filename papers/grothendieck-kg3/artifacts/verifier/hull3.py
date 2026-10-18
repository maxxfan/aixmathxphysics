"""Exact 3D convex hull and inradius of conv{points} (p04 verifier, stdlib only).  Version 1.1.

Points are rational: each is (X, Y, Z, W) with integers X, Y, Z and W > 0, meaning (X/W, Y/W, Z/W).

hull(points)        -> list of facets (i, j, k), counter-clockwise seen from outside
                       (randomized incremental construction with a conflict graph, exact predicates)
validate(points, F) -> raises HullError unless all of the following hold, each checked EXACTLY:
                       (1) every directed edge of F occurs once and its reverse occurs once (F is a closed,
                           consistently oriented triangulated surface, i.e. a simplicial 2-cycle),
                       (2) every facet has positive orientation, det(a, b, c) > 0, i.e. the origin lies strictly
                           on the inner side of its plane,
                       (3) a generic ray from the origin (one whose direction avoids the cones over all edges) meets
                           exactly one facet; three such rays are tested,
                       (4) every input point lies on the closed inner side of every facet plane (v1.1: exact for all
                           point/facet pairs; v1.0 skipped pairs that a floating-point screen found more than 1e-9
                           inside),
                       (5) the Euler characteristic of F is 2.
inradius2(planes)   -> min over facets of dist(origin, plane)^2, an exact Fraction.

Why (1)-(3) suffice for the certified bound conv(points) >= eta*Ball with eta^2 = inradius2: the radial projection
r: F -> S^2, x -> x/|x|, is well defined by (2) and maps the 2-cycle F to a 2-cycle on S^2, so it has a degree. For a
regular direction d (one not in the cone over any edge) the degree is the signed number of facets whose cone contains
d; by (2) every facet counts with sign +1, so by (3) the degree is 1. A map of nonzero degree is onto, so every unit
direction u lies in the cone over some facet T = (a, b, c). The ray t*u meets T at t = h_T / (n_T . u) >= h_T, where n_T is
the unit normal and h_T > 0 the distance of the plane of T from the origin; that point lies in T, a subset of
conv(points), and so does the segment from the origin to it. Hence eta*u is in conv(points) for every unit u, with
eta = min_T h_T.  Check (4) is not needed for this lower bound: it shows that the facets are supporting, i.e. that eta
is the exact inradius rather than a lower bound.  Check (5) is a consistency check.
"""
import random
from fractions import Fraction


class HullError(Exception):
    pass


def _sub(p, q):
    """p - q for rational points given as Fraction triples."""
    return (p[0] - q[0], p[1] - q[1], p[2] - q[2])


def _cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def _dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def _lcm(a, b):
    from math import gcd
    return a // gcd(a, b) * b


def _plane(P, i, j, k):
    """Integer plane (a, b, c, d) through rational points P[i], P[j], P[k] (Fraction triples),
    oriented by the right-hand rule on (i, j, k): value(X,Y,Z,W) = aX + bY + cZ - dW."""
    n = _cross(_sub(P[j], P[i]), _sub(P[k], P[i]))
    h = _dot(n, P[i])
    den = 1
    for x in (n[0], n[1], n[2], h):
        den = _lcm(den, x.denominator)
    a, b, c, d = (int(x * den) for x in (n[0], n[1], n[2], h))
    if a == 0 and b == 0 and c == 0:
        raise HullError("degenerate facet")
    return (a, b, c, d)


def _val(pl, q):
    """Sign-carrying value of the plane at homogeneous point q = (X, Y, Z, W), W > 0."""
    return pl[0] * q[0] + pl[1] * q[1] + pl[2] * q[2] - pl[3] * q[3]


def to_fracs(points):
    return [(Fraction(X, W), Fraction(Y, W), Fraction(Z, W)) for (X, Y, Z, W) in points]


def hull(points, seed=12345):
    n = len(points)
    if n < 4:
        raise HullError("need at least 4 points")
    for q in points:
        if len(q) != 4 or not all(isinstance(x, int) for x in q) or q[3] <= 0:
            raise HullError("points must be (X, Y, Z, W) integers with W > 0")
    P = to_fracs(points)
    # initial tetrahedron: four affinely independent points
    i0 = 0
    i1 = next((i for i in range(n) if P[i] != P[i0]), None)
    if i1 is None:
        raise HullError("all points coincide")
    d01 = _sub(P[i1], P[i0])
    i2 = next((i for i in range(n) if _cross(d01, _sub(P[i], P[i0])) != (0, 0, 0)), None)
    if i2 is None:
        raise HullError("points are collinear")
    nrm = _cross(d01, _sub(P[i2], P[i0]))
    i3 = next((i for i in range(n) if _dot(nrm, _sub(P[i], P[i0])) != 0), None)
    if i3 is None:
        raise HullError("points are coplanar (the vectors do not span R^3)")
    if _dot(nrm, _sub(P[i3], P[i0])) > 0:  # make (i0, i1, i2) face away from i3
        i1, i2 = i2, i1
    # interior reference point (homogeneous): centroid of the tetrahedron
    cen = tuple(sum(P[t][c] for t in (i0, i1, i2, i3)) / 4 for c in range(3))
    cden = _lcm(_lcm(cen[0].denominator, cen[1].denominator), cen[2].denominator)
    C = (int(cen[0] * cden), int(cen[1] * cden), int(cen[2] * cden), cden)

    facets = {}     # id -> (i, j, k)
    planes = {}     # id -> integer plane, outward
    edge = {}       # directed edge (a, b) -> facet id that owns it
    conflict_f = {}  # facet id -> set of point indices strictly outside
    conflict_p = {}  # point index -> set of facet ids that see it
    nid = [0]

    def add_facet(i, j, k):
        pl = _plane(P, i, j, k)
        if _val(pl, C) >= 0:
            raise HullError("orientation error (interior point not strictly inside a new facet)")
        f = nid[0]
        nid[0] += 1
        facets[f] = (i, j, k)
        planes[f] = pl
        for e in ((i, j), (j, k), (k, i)):
            if e in edge:
                raise HullError("non-manifold edge during construction")
            edge[e] = f
        conflict_f[f] = set()
        return f

    def remove_facet(f):
        i, j, k = facets.pop(f)
        planes.pop(f)
        for e in ((i, j), (j, k), (k, i)):
            del edge[e]
        for p in conflict_f.pop(f):
            s = conflict_p.get(p)
            if s is not None:
                s.discard(f)

    init = [add_facet(i0, i1, i2), add_facet(i0, i3, i1), add_facet(i1, i3, i2), add_facet(i2, i3, i0)]
    rest = [i for i in range(n) if i not in (i0, i1, i2, i3)]
    rng = random.Random(seed)
    rng.shuffle(rest)
    for p in rest:
        conflict_p[p] = set()
        for f in init:
            if _val(planes[f], points[p]) > 0:
                conflict_f[f].add(p)
                conflict_p[p].add(f)

    for p in rest:
        vis = conflict_p.pop(p)
        if not vis:
            continue  # inside or on the current hull
        # horizon: directed edges of visible facets whose neighbour is not visible
        horizon = []
        for f in vis:
            i, j, k = facets[f]
            for (a, b) in ((i, j), (j, k), (k, i)):
                g = edge.get((b, a))
                if g is None:
                    raise HullError("open surface during construction")
                if g not in vis:
                    horizon.append((a, b, f, g))
        cand_for = []
        for (a, b, f, g) in horizon:
            cand_for.append((a, b, (conflict_f[f] | conflict_f[g]) - {p}))
        for f in list(vis):
            remove_facet(f)
        for (a, b, cands) in cand_for:
            nf = add_facet(a, b, p)
            pl = planes[nf]
            for q in cands:
                if q in conflict_p and _val(pl, points[q]) > 0:
                    conflict_f[nf].add(q)
                    conflict_p[q].add(nf)
    return [facets[f] for f in sorted(facets)]


def _ray_hits(Pf, F, d):
    """Count facets whose cone (from the origin) contains direction d; None if d is non-generic
    (lies on the boundary of some facet cone)."""
    hits = 0
    for (i, j, k) in F:
        a, b, c = Pf[i], Pf[j], Pf[k]
        det = _dot(a, _cross(b, c))
        if det == 0:
            return None  # facet plane through the origin: rejected elsewhere
        x = _dot(d, _cross(b, c)) / det
        y = _dot(a, _cross(d, c)) / det
        z = _dot(a, _cross(b, d)) / det
        if x == 0 or y == 0 or z == 0:
            if x >= 0 and y >= 0 and z >= 0:
                return None
            continue
        if x > 0 and y > 0 and z > 0:
            hits += 1
    return hits


def validate(points, F, seed=777):
    if not F:
        raise HullError("empty facet list")
    n = len(points)
    # (1) closed oriented surface
    dedges = {}
    for t, (i, j, k) in enumerate(F):
        if len({i, j, k}) != 3 or not all(0 <= x < n for x in (i, j, k)):
            raise HullError(f"bad facet {t}")
        for e in ((i, j), (j, k), (k, i)):
            if e in dedges:
                raise HullError(f"directed edge {e} used twice")
            dedges[e] = t
    for (a, b) in dedges:
        if (b, a) not in dedges:
            raise HullError(f"edge {(a, b)} has no reverse (surface not closed)")
    Pf = to_fracs(points)
    planes = [_plane(Pf, i, j, k) for (i, j, k) in F]
    # (2) positive orientation: the plane offset d is a positive multiple of det(a, b, c); d > 0 <=> origin strictly inside
    for pl in planes:
        if pl[3] <= 0:
            raise HullError("origin not strictly inside a facet plane")
    # (4) every point on the closed inner side of every facet plane, exactly, for all pairs (v1.1)
    for (X, Y, Z, W) in points:
        if max(abs(X), abs(Y), abs(Z)) > 8 * W:
            raise HullError("coordinates must have absolute value <= 8")
    for pl in planes:
        a, b, c, d = pl
        for idx, (X, Y, Z, W) in enumerate(points):
            if a * X + b * Y + c * Z - d * W > 0:
                raise HullError(f"point {idx} lies outside a facet plane")
    # (3) degree one: generic rays from the origin meet exactly one facet
    rng = random.Random(seed)
    done = 0
    tries = 0
    while done < 3:
        tries += 1
        if tries > 50:
            raise HullError("could not find a generic ray")
        d = (Fraction(rng.randint(-10 ** 9, 10 ** 9)), Fraction(rng.randint(-10 ** 9, 10 ** 9)),
             Fraction(rng.randint(-10 ** 9, 10 ** 9)))
        if d == (0, 0, 0):
            continue
        h = _ray_hits(Pf, F, d)
        if h is None:
            continue
        if h != 1:
            raise HullError(f"a ray from the origin meets {h} facets (surface does not cover the boundary once)")
        done += 1
    # Euler characteristic of the facet surface
    verts = {x for f in F for x in f}
    if len(verts) - len(dedges) // 2 + len(F) != 2:
        raise HullError("Euler characteristic is not 2")
    return planes


def inradius2(planes):
    """min over facets of d^2 / |n|^2 for planes (a, b, c, d): squared distance from the origin."""
    best = None
    for (a, b, c, d) in planes:
        r = Fraction(d * d, a * a + b * b + c * c)
        if best is None or r < best:
            best = r
    return best


def brute_facet_planes(points):
    """Reference implementation for tests: all supporting planes through point triples (O(n^4))."""
    Pf = to_fracs(points)
    n = len(points)
    out = set()
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                try:
                    pl = _plane(Pf, i, j, k)
                except HullError:
                    continue
                vals = [_val(pl, q) for q in points]
                if all(v <= 0 for v in vals):
                    out.add(_normalize_plane(pl))
                elif all(v >= 0 for v in vals):
                    out.add(_normalize_plane(tuple(-x for x in pl)))
    return out


def _normalize_plane(pl):
    from math import gcd
    g = 0
    for x in pl:
        g = gcd(g, abs(x))
    return tuple(x // g for x in pl)


def facet_planes_normalized(planes):
    return {_normalize_plane(pl) for pl in planes}
