# -*- coding: utf-8 -*-
"""fq_scalar_lib.py -- helpers for genus2_scalar.py (PLAN_scalar_2026-11-04.md); imports only fq_fp/fq_ring types.

  MONOS, NODES          the monomial basis (1, a1, a1^2, a1^3, a2) and the five interpolation nodes
  vander_inverse()      exact inverse (Fractions) of the 5x5 Vandermonde-type matrix M[j][k] = mono_k(node_j)
  interpolate(vals, ring)   {node: NEl} -> [T0..T4] (NEl), T_k = Sum_j Minv[k][j] v_j
  eval_poly(T, a1, a2, ring) -> NEl  Sum_k T_k mono_k(a1, a2)
  F_p polynomials (lists low->high): p_trim, p_add, p_sub, p_scale, p_mul, p_divmod, p_gcd, p_monic, p_eval,
                   p_powmod, p_roots (roots in F_p of a polynomial of degree <= 3 via gcd with x^p - x and splitting)
  qzeta_coords(R, xs)   Q(zeta)-coordinates (c0..c3, basis 1, z, z^2, z^3) of x = x_1 from its Galois conjugates
                        x_m = sigma_m(x) (m = 1..4) by solving Sum_k c_k zeta^{mk} = x_m over F_p; qzeta_fmt(R, cs)
  ser_coef(S, k), ser_classes(S, R) -> set of exponents mod 5 of the non-zero coefficients
"""
import itertools
from fractions import Fraction as Fr

MONOS = [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1)]                 # a1^e1 a2^e2
NODES = [(Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(-1), Fr(0)), (Fr(2), Fr(0)), (Fr(0), Fr(1))]
MONO_NAMES = ["1", "a1", "a1^2", "a1^3", "a2"]


def _inv_matrix(M):
    n = len(M)
    A = [list(map(Fr, row)) + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def vander_inverse():
    M = [[a1 ** e1 * a2 ** e2 for (e1, e2) in MONOS] for (a1, a2) in NODES]
    Minv = _inv_matrix(M)
    # check
    for i in range(5):
        for j in range(5):
            assert sum(Minv[i][k] * M[k][j] for k in range(5)) == (1 if i == j else 0)
    return Minv


def interpolate(vals, ring, Minv=None):
    """vals: {node (tuple of Fr): NEl}; returns [T0, T1, T2, T3, T4] with v(node) = Sum_k T_k mono_k(node)"""
    if Minv is None: Minv = vander_inverse()
    out = []
    for k in range(5):
        acc = ring.const(0)
        for j, node in enumerate(NODES):
            c = Minv[k][j]
            if c != 0: acc = acc.add(vals[node].scale(c), ring)
        out.append(acc)
    return out


def eval_poly(T, a1, a2, ring):
    R = ring.R
    a1 = R.of(a1); a2 = R.of(a2)
    acc = ring.const(0)
    for k, (e1, e2) in enumerate(MONOS):
        c = R.mul(R.pow(a1, e1), R.pow(a2, e2))
        if c: acc = acc.add(T[k].scale(c), ring)
    return acc


# ------------------------------------------------------------------ polynomials over F_p (lists, low -> high)
def p_trim(a):
    while a and a[-1] == 0: a = a[:-1]
    return a
def p_add(a, b, p):
    n = max(len(a), len(b)); return p_trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(n)])
def p_scale(a, c, p): return p_trim([x * c % p for x in a])
def p_sub(a, b, p): return p_add(a, p_scale(b, p - 1, p), p)
def p_mul(a, b, p):
    if not a or not b: return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): out[i + j] = (out[i + j] + x * y) % p
    return p_trim(out)
def p_divmod(a, b, p):
    a = p_trim(list(a)); b = p_trim(list(b))
    if not b: raise ZeroDivisionError
    inv = pow(b[-1], p - 2, p); q = [0] * max(len(a) - len(b) + 1, 0)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % p; d = len(a) - len(b); q[d] = c
        for i, y in enumerate(b): a[d + i] = (a[d + i] - c * y) % p
        a = p_trim(a)
    return p_trim(q), a
def p_monic(a, p):
    a = p_trim(a)
    return p_scale(a, pow(a[-1], p - 2, p), p) if a else a
def p_gcd(a, b, p):
    a = p_trim(list(a)); b = p_trim(list(b))
    while b: a, b = b, p_divmod(a, b, p)[1]
    return p_monic(a, p)
def p_eval(a, x, p):
    r = 0
    for c in reversed(a): r = (r * x + c) % p
    return r
def p_powmod(base, e, mod, p):
    r = [1]; b = p_divmod(base, mod, p)[1]
    while e:
        if e & 1: r = p_divmod(p_mul(r, b, p), mod, p)[1]
        b = p_divmod(p_mul(b, b, p), mod, p)[1]; e >>= 1
    return r
def p_roots(a, p, R=None):
    """all roots in F_p of a (degree <= 3 expected; general via gcd with x^p - x and random splitting)"""
    import random
    a = p_monic(a, p)
    if len(a) <= 1: return []
    xp = p_powmod([0, 1], p, a, p)
    g = p_gcd(p_sub(xp, [0, 1], p), a, p)           # product of the distinct linear factors
    roots = []
    stack = [g]; rng = random.Random(1)
    while stack:
        f = stack.pop()
        if len(f) <= 1: continue
        if len(f) == 2: roots.append((-f[0]) % p); continue
        while True:
            t = rng.randrange(p)
            h = p_powmod([t, 1], (p - 1) // 2, f, p)      # (x+t)^((p-1)/2) mod f
            d = p_gcd(p_sub(h, [1], p), f, p)
            if 1 < len(d) < len(f):
                stack.append(d); stack.append(p_divmod(f, d, p)[0]); break
    return sorted(roots)
def p_fmt(a, R, var="x"):
    if not a: return "0"
    return " + ".join("(%s)%s" % (R.fmt(c), "" if k == 0 else (" %s" % var if k == 1 else " %s^%d" % (var, k))) for k, c in enumerate(a) if c)


# ------------------------------------------------------------------ Q(zeta) coordinates from Galois conjugates
def qzeta_coords(R, xs):
    """xs[m-1] = sigma_m(x), m = 1..4 (sigma_m: zeta -> zeta^m).  Solve Sum_{k=0}^3 c_k zeta^{mk} = x_m; returns (c0..c3)."""
    p = R.p
    A = [[R.z(m * k) for k in range(4)] + [xs[m - 1] % p] for m in range(1, 5)]
    for c in range(4):
        piv = next(r for r in range(c, 4) if A[r][c])
        A[c], A[piv] = A[piv], A[c]
        inv = pow(A[c][c], p - 2, p); A[c] = [x * inv % p for x in A[c]]
        for r in range(4):
            if r != c and A[r][c]:
                f = A[r][c]; A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
    return tuple(A[k][4] for k in range(4))


def qzeta_ratrec(R, cs):
    """rational reconstruction of the coordinates; None if any fails"""
    out = []
    for c in cs:
        r = R.ratrec(c)
        if r is None: return None
        out.append(r)
    return tuple(out)


def qzeta_fmt(R, cs):
    rs = qzeta_ratrec(R, cs)
    if rs is None: return "[Q(zeta) coords not reconstructible: %s]" % (cs,)
    if all(r == 0 for r in rs[1:]): return str(rs[0])
    parts = []
    for k, r in enumerate(rs):
        if r == 0: continue
        parts.append(("%s" % r) if k == 0 else ("(%s) z%s" % (r, "" if k == 1 else "^%d" % k)))
    return " + ".join(parts) if parts else "0"


def galois_colouring(col, m):
    """the colouring with lambda_i -> lambda_{m i}: colour i -> m i mod 5 (sigma_m on the coloured term)"""
    return tuple((m * c) % 5 for c in col)


# ------------------------------------------------------------------ series helpers
def ser_coef(S, k):
    if k >= S.prec: raise ValueError("coefficient s^%d beyond precision O(s^%d)" % (k, S.prec))
    if k < S.v: return S.R.zero
    return S.c[k - S.v]


def ser_classes(S, R):
    return set((S.v + i) % 5 for i, x in enumerate(S.c) if not R.iszero(x))


def ser_val(S):
    return S.v if S.c else None
