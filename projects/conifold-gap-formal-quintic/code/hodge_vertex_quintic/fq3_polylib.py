# -*- coding: utf-8 -*-
"""fq3_polylib.py -- small exact multivariate polynomial / Macaulay-matrix toolkit over F_p (pure Python) for the genus-3 O(u^3) analysis.

Polynomials in nv variables are dicts {exponent tuple: coefficient mod p}.
  padd, psub, pscale, pmul, pconst, pvar, peval, pdeg, is_hom
  det_poly(M)                      determinant of a small square matrix of polynomials (Laplace expansion along the first row, memoised minors)
  monomials(nv, d)                 list of exponent tuples of total degree d
  Echelon(ncols, p)                incremental row echelon form: .add(row) -> True if the rank grew; .rank; .contains(vec)
  macaulay(polys, nv, D, p)        Echelon of the degree-D part of the homogeneous ideal (rows = poly x monomial of degree D - deg poly)
  binary_gcd_degree(forms, p)      degree of the gcd of binary forms given as coefficient lists
"""
import itertools


def padd(a, b, p):
    r = dict(a)
    for k, v in b.items():
        x = (r.get(k, 0) + v) % p
        if x: r[k] = x
        else: r.pop(k, None)
    return r


def pscale(a, c, p):
    c %= p
    return {k: v * c % p for k, v in a.items()} if c else {}


def psub(a, b, p): return padd(a, pscale(b, p - 1, p), p)


def pmul(a, b, p):
    r = {}
    if len(a) > len(b): a, b = b, a
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb)); r[k] = (r.get(k, 0) + va * vb) % p
    return {k: v for k, v in r.items() if v}


def pconst(c, nv, p): return {(0,) * nv: c % p} if c % p else {}


def pvar(i, nv):
    e = [0] * nv; e[i] = 1
    return {tuple(e): 1}


def peval(a, x, p):
    s = 0
    for k, v in a.items():
        t = v
        for e, xi in zip(k, x):
            if e: t = t * pow(xi, e, p) % p
        s = (s + t) % p
    return s


def pdeg(a): return max((sum(k) for k in a), default=-1)


def is_hom(a): return len(set(sum(k) for k in a)) <= 1


def det_poly(M, p, nv):
    n = len(M); memo = {}
    def minor(r0, colset):
        """determinant of rows r0.. and the columns in colset (tuple)"""
        if r0 == n: return pconst(1, nv, p)
        key = (r0, colset); v = memo.get(key)
        if v is not None: return v
        tot = {}
        for pos, c in enumerate(colset):
            e = M[r0][c]
            if not e: continue
            sub = minor(r0 + 1, colset[:pos] + colset[pos + 1:])
            if not sub: continue
            term = pmul(e, sub, p)
            tot = padd(tot, term, p) if pos % 2 == 0 else psub(tot, term, p)
        memo[key] = tot
        return tot
    return minor(0, tuple(range(n)))


def monomials(nv, d):
    out = []
    def rec(i, left, cur):
        if i == nv - 1: out.append(tuple(cur + [left])); return
        for e in range(left, -1, -1): rec(i + 1, left - e, cur + [e])
    rec(0, d, [])
    return out


class Echelon(object):
    def __init__(self, ncols, p):
        self.n = ncols; self.p = p; self.piv = {}; self.rank = 0

    def _reduce(self, row):
        p = self.p; piv = self.piv; n = self.n; c = 0
        while c < n:
            f = row[c]
            if f:
                pr = piv.get(c)
                if pr is None: return row, c
                row[c:] = [(a - f * b) % p for a, b in zip(row[c:], pr[c:])]
            c += 1
        return row, None

    def add(self, row):
        row, c = self._reduce(list(row))
        if c is None: return False
        inv = pow(row[c], self.p - 2, self.p); row[c:] = [a * inv % self.p for a in row[c:]]
        self.piv[c] = row; self.rank += 1
        return True

    def contains(self, vec):
        return self._reduce(list(vec))[1] is None


def macaulay(polys, nv, D, p, log=None, stop_full=True, extra=None):
    """Echelon form of I_D for homogeneous polys; rows poly x monomial; stops early when the rank equals the number of monomials"""
    mons = monomials(nv, D); index = {m: i for i, m in enumerate(mons)}; E = Echelon(len(mons), p); nrows = 0
    by_deg = {}
    for f in polys:
        if f: by_deg.setdefault(pdeg(f), []).append(f)
    for d in sorted(by_deg):
        if d > D: continue
        mult = monomials(nv, D - d)
        for m in mult:
            for f in by_deg[d]:
                row = [0] * len(mons)
                for k, v in f.items(): row[index[tuple(x + y for x, y in zip(k, m))]] = v
                E.add(row); nrows += 1
                if log and nrows % 500 == 0: log("         degree %d: %d rows processed, rank %d / %d" % (D, nrows, E.rank, len(mons)))
                if stop_full and E.rank == len(mons): return E, mons, index, nrows
    return E, mons, index, nrows


def binary_gcd_degree(forms, p):
    """forms: list of coefficient lists [a_0, .., a_d] of univariate polynomials (dehomogenised binary forms, all of the same formal degree d);
    returns (degree of the gcd of the dehomogenised polynomials, number of common roots at infinity = min over forms of d - deg)"""
    def trim(a):
        a = [x % p for x in a]
        while a and a[-1] == 0: a.pop()
        return a
    def rem(a, b):
        a = list(a); inv = pow(b[-1], p - 2, p)
        while len(a) >= len(b):
            f = a[-1] * inv % p; s = len(a) - len(b)
            for i, x in enumerate(b): a[s + i] = (a[s + i] - f * x) % p
            a = trim(a)
            if not a: break
        return a
    d = max(len(f) for f in forms) - 1; inf = min(d - (len(trim(f)) - 1) for f in forms if trim(f)) if any(trim(f) for f in forms) else d
    g = []
    for f in forms:
        f = trim(f)
        if not f: continue
        if not g: g = f; continue
        a, b = (g, f) if len(g) >= len(f) else (f, g)
        while b: a, b = b, rem(a, b)
        g = a
    return (len(g) - 1 if g else None), inf
