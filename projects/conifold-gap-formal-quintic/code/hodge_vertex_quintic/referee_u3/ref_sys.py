# -*- coding: utf-8 -*-
"""ref_sys.py -- referee ELIMINATION: the circulant system as explicit polynomials / over arbitrary coefficient rings.
Only the DATA of CircSys are used (cubic coefficients S.C[i].coef, quadric matrices S.Q[i], rows S.L2, S.J3, S.m1); the low-part
parametrisation, all evaluations and all eliminations are re-implemented here."""
import os, sys
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ref_lib import *

CUB = [4, 9, 10, 11, 12, 18, 19, 20, 32]
QUA = [3, 6, 7, 8, 17]


class MRing(object):
    """multivariate polynomials in nv variables over F_p (dicts)"""
    def __init__(self, nv, p): self.nv = nv; self.p = p; self.zero = {}; self.one = mconst(1, nv, p)
    def add(self, a, b): return madd(a, b, self.p)
    def sub(self, a, b): return msub(a, b, self.p)
    def mul(self, a, b): return mmul(a, b, self.p)
    def scal(self, c, a): return mscale(a, c, self.p)
    def const(self, c): return mconst(c, self.nv, self.p)
    def var(self, i): return mvar(i, self.nv)
    def iszero(self, a): return not a


class IRing(object):
    def __init__(self, p): self.p = p; self.zero = 0; self.one = 1
    def add(self, a, b): return (a + b) % self.p
    def sub(self, a, b): return (a - b) % self.p
    def mul(self, a, b): return a * b % self.p
    def scal(self, c, a): return a * c % self.p
    def const(self, c): return c % self.p
    def iszero(self, a): return a % self.p == 0
    def inv(self, a): return pow(a, self.p - 2, self.p)


class KRing(object):
    """K = F_p[x]/(f), f monic irreducible (a field); elements = trimmed coefficient lists"""
    def __init__(self, f, p): self.f = umonic(f, p); self.p = p; self.zero = []; self.one = [1]; self.d = len(self.f) - 1
    def add(self, a, b): return uadd(a, b, self.p)
    def sub(self, a, b): return usub(a, b, self.p)
    def mul(self, a, b): return umod(umul(a, b, self.p), self.f, self.p)
    def scal(self, c, a): return uscale(a, c, self.p)
    def const(self, c): return utrim([c], self.p)
    def iszero(self, a): return not a
    def gen(self): return umod([0, 1], self.f, self.p)
    def inv(self, a):
        p = self.p; r0, r1 = self.f, utrim(a, p); s0, s1 = [], [1]          # s_i * a = r_i mod f
        if not r1: raise ZeroDivisionError
        while r1:
            q, r = udivmod(r0, r1, p); r0, r1 = r1, r; s0, s1 = s1, usub(s0, umul(q, s1, p), p)
        if len(r0) != 1: raise ZeroDivisionError("non-unit: f is reducible")
        return umod(uscale(s0, pow(r0[0], p - 2, p), p), self.f, p)


def load(tag):
    from fq3_circsys import CircSys
    S = CircSys(tag)
    return S


class Sys(object):
    """the conditions over a ring R at a 15-vector of ring elements"""
    def __init__(self, S):
        self.S = S; self.p = p = S.p; self.names = S.names
        self.C = {i: dict(S.C[i].coef) for i in CUB}; self.Q = {i: [list(r) for r in S.Q[i]] for i in QUA}
        self.L2 = [list(r) for r in S.L2]; self.J3 = [list(r) for r in S.J3]; self.m1 = list(S.m1)
        # my own solution of m_1 = 0 and the two genus-2 rows for (c2_2; c3_0, c3_1) in terms of (c1_0, c1_1, c1_2, c2_1, c3_2)
        assert all(self.L2[r][j] == 0 for r in range(2) for j in range(8, 15)), "L2 rows involve orders >= 4"
        assert all(self.m1[j] == 0 for j in range(5, 15)) and self.m1[4] % p
        free = [0, 1, 2, 3, 7]; inv4 = pow(self.m1[4], p - 2, p)
        self.e4 = {j: (-self.m1[j]) * inv4 % p for j in free}                      # c2_2 = Sum e4[j] v_j
        a, b, c, d = self.L2[0][5], self.L2[0][6], self.L2[1][5], self.L2[1][6]; det = (a * d - b * c) % p; assert det
        idet = pow(det, p - 2, p); self.e5 = {}; self.e6 = {}
        for j in free:
            r0 = (self.L2[0][j] + self.L2[0][4] * self.e4[j]) % p; r1 = (self.L2[1][j] + self.L2[1][4] * self.e4[j]) % p
            self.e5[j] = (-(d * r0 - b * r1)) * idet % p; self.e6[j] = (-(a * r1 - c * r0)) * idet % p

    def low(self, R, c10, c11, c12, tau, w):
        """the 8 low coordinates as ring elements"""
        fr = {0: c10, 1: c11, 2: c12, 3: tau, 7: w}
        def lin(e):
            s = R.zero
            for j, x in fr.items(): s = R.add(s, R.scal(e[j], x))
            return s
        return [c10, c11, c12, tau, lin(self.e4), lin(self.e5), lin(self.e6), w]

    def cub(self, R, i, v):
        s = R.zero; cache = {}
        for (a, b, c), co in self.C[i].items():
            ab = cache.get((a, b))
            if ab is None: ab = R.mul(v[a], v[b]); cache[(a, b)] = ab
            if R.iszero(ab): continue
            s = R.add(s, R.scal(co, R.mul(ab, v[c])))
        return s

    def quad(self, R, i, v):
        A = self.Q[i]; s = R.zero
        for a in range(15):
            if R.iszero(v[a]): continue
            t = R.zero
            for b in range(15):
                if A[a][b] and not R.iszero(v[b]): t = R.add(t, R.scal(A[a][b], v[b]))
            s = R.add(s, R.mul(v[a], t))
        return s

    def lin(self, R, row, v):
        s = R.zero
        for j in range(15):
            if row[j] and not R.iszero(v[j]): s = R.add(s, R.scal(row[j], v[j]))
        return s

    def residuals(self, R, v):
        out = {}
        for r in range(2): out['L2_%d' % r] = self.lin(R, self.L2[r], v)
        out['m1'] = self.lin(R, self.m1, v)
        for r in range(2): out['J3_%d' % r] = self.lin(R, self.J3[r], v)
        for i in QUA: out['q_%d' % i] = self.quad(R, i, v)
        for i in CUB: out['c_%d' % i] = self.cub(R, i, v)
        return out


def solve_field(R, A, b):
    """Gaussian elimination over a field ring R (needs R.inv): returns (particular solution or None, rank, free columns)"""
    n = len(A[0]); M = [list(r) + [bb] for r, bb in zip(A, b)]; piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, len(M)) if not R.iszero(M[i][c])), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]; inv = R.inv(M[r][c]); M[r] = [R.mul(inv, x) for x in M[r]]
        for i in range(len(M)):
            if i != r and not R.iszero(M[i][c]):
                f = M[i][c]; M[i] = [R.sub(x, R.mul(f, y)) for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == len(M): break
    for i in range(r, len(M)):
        if not R.iszero(M[i][n]): return None, r, [c for c in range(n) if c not in piv]
    x = [R.zero] * n
    for k, c in enumerate(piv): x[c] = M[k][n]
    return x, r, [c for c in range(n) if c not in piv]


def bires(F, G, iv, ie, nv, p, dF=None, dG=None, extra=4):
    """Res with respect to variable ie of two polynomials in the variables (iv, ie) (dicts on nv variables), as a univariate polynomial in
    variable iv, by evaluation / Newton interpolation with the Euclidean resultant at FORMAL degrees (dF, dG) in ie"""
    dF = mdeg(F, ie) if dF is None else dF; dG = mdeg(G, ie) if dG is None else dG
    bound = mdeg(F, iv) * dG + mdeg(G, iv) * dF; n = bound + 1 + extra; xs = list(range(3, n + 3)); ys = []
    for x in xs:
        f = muni(msubs(F, {iv: x}, p), ie, p); g = muni(msubs(G, {iv: x}, p), ie, p); ys.append(ures(f, g, dF, dG, p))
    out = newton_interp(xs, ys, p); assert len(out) - 1 <= bound, (len(out) - 1, bound)
    return out


def fmt_fac(fs): return [len(f) - 1 for f in fs]
