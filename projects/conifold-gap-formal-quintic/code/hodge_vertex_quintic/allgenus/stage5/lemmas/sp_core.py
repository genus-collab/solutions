# -*- coding: utf-8 -*-
"""sp_core.py -- Stage 5 / lemmas: exact formal stationary phase (SP) by Wick contraction over Q (Fractions).

New code, stdlib only; imports no project module.  Notation of FORMAL_SP_LEMMAS.md s.1.
A graded series is a dict {(e, j): c}: e = exponent tuple over the NV variables, j = power of z, c a Fraction.
Doubled weight w2(e, j) = sum_i dw[i] e_i + 2 j, with dw[i] = 1 for an integration variable (weight 1/2) and
dw[i] = 2 for a parameter (weight 1).  Every operation used below is non-decreasing in w2, so truncation at
w2 <= W2 is consistent; all identities are compared on the truncated range.
SP^norm[f; a] := < exp(f_3/z) a >_{zC},  C = -Hess f(0)^-1  (Definition 1.3 of the note): Wick expectation of
x^alpha is z^(|alpha|/2) times the sum over perfect matchings of prod C_ij (computed by the recursion
<x_i x^b> = sum_j C_ij b_j <x^(b - e_j)>, i.e. Gaussian integration by parts).
"""
from fractions import Fraction as Fr
import random


class Ctx:
    def __init__(self, nv, dw, W2):
        self.nv, self.dw, self.W2 = nv, tuple(dw), W2

    def w2(self, k):
        return sum(d * x for d, x in zip(self.dw, k[0])) + 2 * k[1]

    def zero_e(self):
        return (0,) * self.nv

    # ---- ring operations on graded series ---------------------------------------------------------------
    def trunc(self, A):
        return {k: c for k, c in A.items() if c and self.w2(k) <= self.W2}

    def const(self, c):
        return {(self.zero_e(), 0): Fr(c)} if c else {}

    def var(self, i, c=1):
        e = [0] * self.nv; e[i] = 1
        return self.trunc({(tuple(e), 0): Fr(c)})

    def zp(self, A, s):                       # multiply by z^s
        return self.trunc({(e, j + s): c for (e, j), c in A.items()})

    def add(self, A, B, cb=1):
        R = dict(A)
        for k, c in B.items():
            R[k] = R.get(k, 0) + cb * c
        return {k: c for k, c in R.items() if c}

    def scal(self, A, c):
        return {k: v * c for k, v in A.items() if v * c}

    def mul(self, A, B):
        W2, res = self.W2, {}
        Bl = sorted(((self.w2(k), k[0], k[1], c) for k, c in B.items()), key=lambda t: t[0])
        for (ea, ja), ca in A.items():
            wa = self.w2((ea, ja))
            for wb, eb, jb, cb in Bl:
                if wa + wb > W2:
                    break
                k = (tuple(x + y for x, y in zip(ea, eb)), ja + jb)
                res[k] = res.get(k, 0) + ca * cb
        return {k: c for k, c in res.items() if c}

    def power(self, A, n):
        R = self.const(1)
        for _ in range(n):
            R = self.mul(R, A)
        return R

    def exp(self, A):                         # A must have w2 >= 1 on every term
        R, T, m = self.const(1), self.const(1), 0
        while T:
            m += 1; T = self.scal(self.mul(T, A), Fr(1, m)); R = self.add(R, T)
        return R

    def binom_pow(self, A, al):               # (1 + A)^al, A with w2 >= 1
        R, T, m = self.const(1), self.const(1), 0
        while T:
            T = self.scal(self.mul(T, A), Fr(al - m, m + 1)); m += 1; R = self.add(R, T)
        return R

    def deriv(self, A, i):
        R = {}
        for (e, j), c in A.items():
            if e[i]:
                e2 = list(e); e2[i] -= 1; R[(tuple(e2), j)] = c * e[i]
        return R

    def compose(self, A, subs):
        """substitute variable i -> subs[i] (graded series, w2 >= dw[i]) for i in subs (a dict)"""
        cache = {i: [self.const(1)] for i in subs}
        def pw(i, n):
            while len(cache[i]) <= n:
                cache[i].append(self.mul(cache[i][-1], subs[i]))
            return cache[i][n]
        R = {}
        for (e, j), c in A.items():
            keep = tuple(0 if i in subs else x for i, x in enumerate(e))
            T = {(keep, j): c}
            for i in subs:
                if e[i]:
                    T = self.mul(T, pw(i, e[i]))
            R = self.add(R, T)
        return self.trunc(R)

    def split(self, A, xs):
        """split a z-free series by total degree in the variables xs: returns {deg: series}"""
        out = {}
        for (e, j), c in A.items():
            d = sum(e[i] for i in xs)
            out.setdefault(d, {})[(e, j)] = c
        return out

    def quad_matrix(self, Q2, xs):
        """the matrix K of the xs-quadratic part Q2 = (1/2) x^T K x (entries: series in the other variables)"""
        n = len(xs); K = [[{} for _ in range(n)] for _ in range(n)]
        for (e, j), c in Q2.items():
            idx = [a for a, i in enumerate(xs) for _ in range(e[i])]
            rest = tuple(0 if i in xs else x for i, x in enumerate(e))
            a, b = idx
            if a == b:
                K[a][a] = self.add(K[a][a], {(rest, j): 2 * c})
            else:
                K[a][b] = self.add(K[a][b], {(rest, j): c}); K[b][a] = self.add(K[b][a], {(rest, j): c})
        return K

    # ---- matrices of series -----------------------------------------------------------------------------
    def cpart(self, A):
        return A.get((self.zero_e(), 0), Fr(0))

    def inv(self, K):
        """inverse of a matrix of series whose constant part K0 is invertible over Q"""
        n = len(K); K0i = inv_q([[self.cpart(K[a][b]) for b in range(n)] for a in range(n)])
        K0i_s = [[self.const(K0i[a][b]) for b in range(n)] for a in range(n)]
        K1 = [[self.add(K[a][b], self.const(self.cpart(K[a][b])), -1) for b in range(n)] for a in range(n)]
        Mt = mmul(self, K0i_s, K1)                         # K^-1 = sum_r (-K0^-1 K1)^r K0^-1
        T, R = K0i_s, K0i_s
        while any(T[a][b] for a in range(n) for b in range(n)):
            T = [[self.scal(x, -1) for x in row] for row in mmul(self, Mt, T)]
            R = [[self.add(R[a][b], T[a][b]) for b in range(n)] for a in range(n)]
        return R

    def det(self, K):
        n = len(K)
        if n == 1:
            return K[0][0]
        R = {}
        for b in range(n):
            minor = [[K[a][c] for c in range(n) if c != b] for a in range(1, n)]
            R = self.add(R, self.mul(K[0][b], self.det(minor)), (-1) ** b)
        return R

    # ---- Wick expectation and the normalised formal SP --------------------------------------------------
    def moments(self, C, maxdeg):
        """M[alpha] = <x^alpha>_C (series in the other variables) for |alpha| <= maxdeg, C a matrix of series"""
        n = len(C); M = {(0,) * n: self.const(1)}
        def get(al):
            if al in M:
                return M[al]
            if sum(al) % 2:
                return {}
            i = next(k for k in range(n) if al[k]); b = list(al); b[i] -= 1; s = {}
            for jj in range(n):
                if b[jj] and C[i][jj]:
                    bb = list(b); bb[jj] -= 1
                    s = self.add(s, self.mul(C[i][jj], get(tuple(bb))), b[jj])
            M[al] = s
            return s
        return get

    def sp_norm(self, P3, amp, xs, C):
        """< exp(P3/z) amp > over the variables xs with covariance z C (C: matrix of z-free series in the
        other variables).  P3 must have w2 >= 3 on every term, so that P3/z has w2 >= 1."""
        wmin = min([0] + [self.w2(k) for k in amp])          # an amplitude of negative weight needs more phase data
        R = {(e, j - 1): c for (e, j), c in P3.items() if c and self.w2((e, j - 1)) <= self.W2 - wmin}
        E = dict(amp); res = {}; m = 0; mom = self.moments(C, 0)
        assert all(self.w2(k) >= 1 for k in R), "phase part of weight < 1/2"
        assert wmin == 0 or all(not any(k[0]) for row in C for x in row for k in x), \
            "an amplitude of negative weight is supported only with a constant covariance"
        while E:
            for (e, j), c in E.items():
                ex = tuple(e[i] for i in xs)
                if sum(ex) % 2:
                    continue
                rest = tuple(0 if i in xs else x for i, x in enumerate(e))
                T = self.mul({(rest, j + sum(ex) // 2): c}, mom(ex))
                for k, v in T.items():
                    res[k] = res.get(k, 0) + v
            m += 1; E = self.scal(self.mul(E, R), Fr(1, m))
        return {k: c for k, c in res.items() if c}


def mmul(ctx, A, B):
    n, p, q = len(A), len(B), len(B[0])
    return [[_sum(ctx, [ctx.mul(A[a][k], B[k][b]) for k in range(p)]) for b in range(q)] for a in range(n)]


def _sum(ctx, L):
    R = {}
    for x in L:
        R = ctx.add(R, x)
    return R


def inv_q(M):
    n = len(M); A = [[Fr(x) for x in r] + [Fr(int(i == j)) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        r = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[r] = A[r], A[c]
        iv = 1 / A[c][c]; A[c] = [x * iv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]


def rand_poly(ctx, rng, vars_, dmin, dmax, cmax=3, dens=1.0):
    """random polynomial (z-free) in the variables vars_, total degree dmin..dmax, small integer coefficients"""
    out = {}
    def rec(i, e, d):
        if i == len(vars_):
            if dmin <= d <= dmax and rng.random() < dens:
                c = rng.randint(-cmax, cmax)
                if c:
                    ee = [0] * ctx.nv
                    for v, x in zip(vars_, e):
                        ee[v] = x
                    out[(tuple(ee), 0)] = Fr(c)
            return
        for x in range(dmax - d + 1):
            rec(i + 1, e + [x], d + x)
    rec(0, [], 0)
    return ctx.trunc(out)


def zseries(ctx, A):
    """the z-coefficients of a series in z only: {j: c}"""
    out = {}
    for (e, j), c in A.items():
        assert not any(e), "unexpected variable left"
        out[j] = out.get(j, 0) + c
    return {j: c for j, c in out.items() if c}
