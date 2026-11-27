# -*- coding: utf-8 -*-
"""fq_diag_lib.py -- helpers for genus2_diag.py: the DIAGONAL, COLOUR-DEPENDENT constant scalar frame
R = R' diag(phi_i(z)),  phi_i(z) = exp(a1_i z + a2_i z^3),  ten independent constants (a1_0..a1_4, a2_0..a2_4),
vertex CohFT phi_i . 1 at colour i (no colour mixing).  Imports fq_phi only for phi_inv_coeffs / PhiVertexTable,
called with colour argument 0 (no hard-coded covariant phase) and the colour's own (a1_i, a2_i).

Global variables: index i (0..4) = a1_i (weight 1), index 5+i = a2_i (weight 3); a point is a 10-tuple in F_p.
Local variables of a coloured graph (colour tuple col): vars = [(1,c) for c in cs] + [(2,c) for c in cs], cs = the
distinct colours in order of first appearance; local monomials = exponent tuples over vars of weight <= 3
(5 for one colour: 1, x, x^2, x^3, p;  12 for two: 1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3, p, q).

  diag_frame(rows, pt, name)             gsum_q.Frame for the point pt
  covariant_point(R, a1, a2)             a1_i = zeta^-i a1, a2_i = zeta^-3i a2 (the covariant scalar of fq_phi)
  local_vars(col), local_monos(col), mono_weight, mono_global, mono_name, mono_eval
  solve_coeffs(vals, pts, col, ring)     {mono: NEl} from the first m = len(monos) points (F_p matrix inverse)
  eval_local(T, col, pt, ring)           Sum_mono T[mono] mono(pt)
  fp_inv(M, p), fp_rref(rows, p) -> (rref, pivots), fp_nullspace(rows, n, p)
"""
import itertools, random, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_phi import phi_inv_coeffs, PhiVertexTable

WEIGHT = (1, 1, 1, 1, 1, 3, 3, 3, 3, 3)
GNAMES = ["a1_%d" % i for i in range(5)] + ["a2_%d" % i for i in range(5)]


def diag_frame(rows, pt, name):
    """Frame with rows R phi_i^-1, phi_i = exp(a1_i z + a2_i z^3), a1_i = pt[i], a2_i = pt[5+i]; vertex classes phi_i . 1"""
    from gsum_q import Frame
    R = rows.R
    pt = tuple(R.of(x) for x in pt)
    hfac = {i: phi_inv_coeffs(R, pt[i], pt[5 + i], 0, rows.Nz) for i in range(5)}
    F = Frame(rows, name, hfac, None)
    F.vt = {i: PhiVertexTable(R, pt[i], pt[5 + i], 0) for i in range(5)}
    F.diag_point = pt
    return F


def covariant_point(R, a1, a2):
    a1 = R.of(a1); a2 = R.of(a2)
    return tuple(R.mul(a1, R.z(-i)) for i in range(5)) + tuple(R.mul(a2, R.z(-3 * i)) for i in range(5))


def zero_point():
    return (0,) * 10


def random_point(rng, p):
    return tuple(rng.randrange(1, p) for _ in range(10))


# ------------------------------------------------------------------ local monomials
def colours_of(col):
    cs = []
    for c in col:
        if c not in cs: cs.append(c)
    return cs


def local_vars(col):
    cs = colours_of(col)
    return [(1, c) for c in cs] + [(2, c) for c in cs]


def var_weight(v): return 1 if v[0] == 1 else 3


def local_monos(col):
    """exponent tuples over local_vars(col) with weight <= 3, ordered by (weight, lexicographic)"""
    vs = local_vars(col)
    out = []
    for ex in itertools.product(range(4), repeat=len(vs)):
        w = sum(e * var_weight(v) for e, v in zip(ex, vs))
        if w <= 3: out.append(ex)
    out.sort(key=lambda ex: (sum(e * var_weight(v) for e, v in zip(ex, vs)), [-e for e in ex]))
    return out


def mono_weight(mono, vs):
    return sum(e * var_weight(v) for e, v in zip(mono, vs))


def mono_global(mono, vs):
    """10-tuple of exponents over the global variables"""
    g = [0] * 10
    for e, (k, c) in zip(mono, vs):
        g[c if k == 1 else 5 + c] += e
    return tuple(g)


def gmono_weight(g): return sum(e * w for e, w in zip(g, WEIGHT))


def gmono_name(g):
    parts = []
    for e, n in zip(g, GNAMES):
        if e: parts.append(n if e == 1 else "%s^%d" % (n, e))
    return " ".join(parts) if parts else "1"


def mono_name(mono, vs): return gmono_name(mono_global(mono, vs))


def local_name(mono, vs):
    """short local names: x, y (a1 of the first/second colour), p, q (a2)"""
    nm = ["x", "y"] if len(vs) == 4 else ["x"]
    nm = nm + (["p", "q"] if len(vs) == 4 else ["p"])
    parts = []
    for e, n in zip(mono, nm):
        if e: parts.append(n if e == 1 else "%s^%d" % (n, e))
    return " ".join(parts) if parts else "1"


def mono_eval(R, mono, vs, pt):
    r = R.one
    for e, (k, c) in zip(mono, vs):
        if e: r = R.mul(r, R.pow(pt[c if k == 1 else 5 + c], e))
    return r


def gmono_eval(R, g, pt):
    r = R.one
    for e, x in zip(g, pt):
        if e: r = R.mul(r, R.pow(x, e))
    return r


# ------------------------------------------------------------------ F_p linear algebra
def fp_inv(M, p):
    n = len(M)
    A = [[x % p for x in row] + [int(i == j) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c]), None)
        if piv is None: raise ZeroDivisionError("singular matrix over F_p")
        A[c], A[piv] = A[piv], A[c]
        inv = pow(A[c][c], p - 2, p); A[c] = [x * inv % p for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]; A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def fp_rref(rows, p):
    """reduced row echelon form of a list of rows (lists of equal length) over F_p; returns (rref rows, pivot columns)"""
    A = [[x % p for x in r] for r in rows if any(x % p for x in r)]
    if not A: return [], []
    n = len(A[0]); pivots = []; r0 = 0
    for c in range(n):
        piv = next((r for r in range(r0, len(A)) if A[r][c]), None)
        if piv is None: continue
        A[r0], A[piv] = A[piv], A[r0]
        inv = pow(A[r0][c], p - 2, p); A[r0] = [x * inv % p for x in A[r0]]
        for r in range(len(A)):
            if r != r0 and A[r][c]:
                f = A[r][c]; A[r] = [(x - f * y) % p for x, y in zip(A[r], A[r0])]
        pivots.append(c); r0 += 1
        if r0 == len(A): break
    return A[:r0], pivots


def fp_nullspace(rows, n, p):
    """basis of {v in F_p^n : row . v = 0 for all rows}"""
    A, piv = fp_rref(rows, p)
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = [0] * n; v[f] = 1
        for r, c in zip(A, piv): v[c] = (-r[f]) % p
        basis.append(v)
    return basis


# ------------------------------------------------------------------ the interpolation solve on NEls
def solve_coeffs(vals, pts, col, ring):
    """vals[j] = graph-sum NEl at point pts[j]; uses the first m = #monos points: T = M^-1 vals.  Returns ({mono: NEl}, m)."""
    R = ring.R; p = R.p
    vs = local_vars(col); monos = local_monos(col); m = len(monos)
    M = [[mono_eval(R, mono, vs, pts[j]) for mono in monos] for j in range(m)]
    Minv = fp_inv(M, p)
    T = {}
    for k, mono in enumerate(monos):
        acc = ring.const(0)
        for j in range(m):
            c = Minv[k][j]
            if c: acc = acc.add(vals[j].scale(c), ring)
        T[mono] = ring.reduce(acc)
    return T, m


def eval_local(T, col, pt, ring):
    R = ring.R; vs = local_vars(col)
    acc = ring.const(0)
    for mono, el in T.items():
        c = mono_eval(R, mono, vs, pt)
        if c and not el.iszero(): acc = acc.add(el.scale(c), ring)
    return ring.reduce(acc)


def _selftest():
    import os, sys
    from fq_fp import FpRing, P1
    R = FpRing(P1); p = R.p
    assert len(local_monos((3,))) == 5 and len(local_monos((2, 2))) == 5 and len(local_monos((0, 3))) == 12
    vs = local_vars((0, 3)); ms = local_monos((0, 3))
    assert [local_name(m, vs) for m in ms] == ["1", "x", "y", "x^2", "x y", "y^2", "x^3", "x^2 y", "x y^2", "y^3", "p", "q"], [local_name(m, vs) for m in ms]
    assert mono_global(ms[4], vs) == (1, 0, 0, 1, 0, 0, 0, 0, 0, 0) and mono_global(ms[11], vs) == (0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
    assert all(mono_weight(m, vs) == gmono_weight(mono_global(m, vs)) for m in ms)
    rng = random.Random(5)
    M = [[rng.randrange(p) for _ in range(4)] for _ in range(4)]
    Mi = fp_inv(M, p)
    assert all(sum(M[i][k] * Mi[k][j] for k in range(4)) % p == int(i == j) for i in range(4) for j in range(4))
    rows = [[1, 2, 3, 0], [2, 4, 6, 0], [0, 0, 0, 1]]
    A, piv = fp_rref(rows, p); assert piv == [0, 3] and len(A) == 2
    ns = fp_nullspace(rows, 4, p); assert len(ns) == 2 and all(sum(r[i] * v[i] for i in range(4)) % p == 0 for v in ns for r in rows)
    cp = covariant_point(R, 7, 11)
    assert cp[0] == 7 and cp[5] == 11 and cp[1] == R.mul(7, R.z(-1)) and cp[7] == R.mul(11, R.z(-6))
    print("fq_diag_lib self-tests ok")


if __name__ == "__main__":
    _selftest()
