"""cr5_cyc.py -- Stage 5 critic: (CYC) by special geometry, exact over Q (stdlib only; imports no project or agent module).

(CYC) (PB_AUDIT s.2.3): with D = q d/dq = -(1-u) d/du, PF4 = D^4 - (1-u) prod_(k=1..4) (D + k/5) (= u D^4 - (1-u)(2D^3 + ...)), and a
conifold frame (varpi0 analytic at u = 0, g the vanishing period, t = g/varpi0),
    PF4 = M := (1/varpi0) D (1/Dt) D (u varpi0^2 Dt^2) D (1/Dt) D (1/varpi0).
The critic's proof (STAGE5_SYNTHESIS_2 s.9) uses: (a) PF4 is formally self-adjoint for D (measure dq/q); (b) the Lagrange concomitant
B(f, g) = Sum_k Sum_(j<k) (-1)^j D^j(A_k g) D^(k-1-j) f is a constant, alternating, non-degenerate form on the solutions; (c) its
jet matrix is anti-triangular with anti-diagonal +-u, hence Sum_ij Binv_ij y_i D^k y_j = 0 (k = 0, 1, 2) and = 1/u (k = 3);
(d) B(varpi0, g) = 0 for every analytic varpi0 (Picard-Lefschetz); (e) special geometry gives M y = 0 for all four solutions.
This module holds the series and operator machinery; the checks (a)-(e) and the operator identity M u^n = PF4 u^n, with failing
controls, are in cr5_cyc_run.py.  Usage: py -I -B cr5_cyc_run.py N   (series to u^N; log logs/cr5_cyc_run_N.log)
"""
import os, sys, time
from fractions import Fraction as Fr
from math import comb

INF = 10 ** 9


class LS:  # Laurent series in u, absolute precision p (exponents >= p unknown)
    __slots__ = ('c', 'p')

    def __init__(self, c=None, p=INF):
        self.p = p
        self.c = {e: v for e, v in (c or {}).items() if v != 0 and e < p}

    def val(self):
        return min(self.c) if self.c else self.p

    def __add__(a, b):
        d = dict(a.c)
        for e, v in b.c.items():
            d[e] = d.get(e, 0) + v
        return LS(d, min(a.p, b.p))

    def __neg__(a):
        return LS({e: -v for e, v in a.c.items()}, a.p)

    def __sub__(a, b):
        return a + (-b)

    def sc(a, s):
        return LS({e: s * v for e, v in a.c.items()}, a.p)

    def __mul__(a, b):
        va, vb = a.val(), b.val()
        p = min(a.p + vb if a.p < INF else INF, b.p + va if b.p < INF else INF)
        d = {}
        for e1, v1 in a.c.items():
            for e2, v2 in b.c.items():
                e = e1 + e2
                if e < p:
                    d[e] = d.get(e, 0) + v1 * v2
        return LS(d, p)

    def inv(a):  # a = c0 + O(u), c0 != 0, finite precision p
        assert a.val() == 0 and a.p < INF
        p, c0 = a.p, a.c[0]
        r = [Fr(0)] * p
        r[0] = 1 / Fr(c0)
        for n in range(1, p):
            s = sum(a.c.get(k, 0) * r[n - k] for k in range(1, n + 1))
            r[n] = -s / c0
        return LS({i: r[i] for i in range(p)}, p)

    def trunc(a, p):
        return LS(a.c, min(p, a.p))


def ls_th(a):  # D = -(1-u) d/du on Laurent series
    d = {}
    for e, v in a.c.items():
        if e:
            d[e - 1] = d.get(e - 1, 0) - e * v
            d[e] = d.get(e, 0) + e * v
    return LS(d, a.p - 1 if a.p < INF else INF)


ELL = LS({0: Fr(1), -1: Fr(-1)})  # D(log u) = -(1-u)/u


class LF:  # Sum_k F[k] (log u)^k
    def __init__(self, F):
        self.F = {k: v for k, v in F.items()}

    def __add__(a, b):
        d = dict(a.F)
        for k, v in b.F.items():
            d[k] = d[k] + v if k in d else v
        return LF(d)

    def sc(a, s):
        return LF({k: v.sc(s) for k, v in a.F.items()})

    def mul(a, b):  # b an LF
        d = {}
        for k1, v1 in a.F.items():
            for k2, v2 in b.F.items():
                d[k1 + k2] = d[k1 + k2] + v1 * v2 if (k1 + k2) in d else v1 * v2
        return LF(d)

    def mls(a, s):  # times a log-free series
        return LF({k: v * s for k, v in a.F.items()})


def lf_th(f):
    d = {}
    for k, v in f.F.items():
        d[k] = d[k] + ls_th(v) if k in d else ls_th(v)
        if k:
            t = (v * ELL).sc(k)
            d[k - 1] = d[k - 1] + t if (k - 1) in d else t
    return LF(d)


def lf(s):
    return LF({0: s})


def is_zero(f, upto):  # every coefficient below u^upto vanishes (and is known)
    for k, v in f.F.items():
        assert v.p >= upto, ('precision', k, v.p, upto)
        if any(e < upto for e in v.c):
            return False
    return True


# ---------------------------------------------------------------- PF4
def poly(d):
    return LS({e: Fr(v) for e, v in d.items()})


A = {4: poly({1: 1}), 3: poly({0: -2, 1: 2}), 2: poly({0: Fr(-7, 5), 1: Fr(7, 5)}),
     1: poly({0: Fr(-2, 5), 1: Fr(2, 5)}), 0: poly({0: Fr(-24, 625), 1: Fr(24, 625)})}


def apply_op(Aop, f, th=lf_th):  # Sum_k A_k D^k f
    out, cur = None, f
    for k in range(0, 5):
        if k in Aop:
            t = cur.mls(Aop[k])
            out = t if out is None else out + t
        cur = th(cur)
    return out


def adjoint(Aop):  # Sum_k (-D)^k o A_k = Sum_k (-1)^k Sum_j C(k,j) D^j(A_k) D^(k-j)
    out = {}
    for k, a in Aop.items():
        dj = a
        for j in range(0, k + 1):
            t = dj.sc((-1) ** k * comb(k, j))
            out[k - j] = out[k - j] + t if (k - j) in out else t
            dj = ls_th(dj)
    return out


def concomitant(Aop, f, g):  # B(f, g) = Sum_k Sum_(j<k) (-1)^j D^j(A_k g) D^(k-1-j) f
    tf = [f]
    for _ in range(4):
        tf.append(lf_th(tf[-1]))
    out = None
    for k, a in Aop.items():
        h = g.mls(a)
        for j in range(0, k):
            t = h.mul(tf[k - 1 - j]).sc((-1) ** j)
            out = t if out is None else out + t
            h = lf_th(h)
    return out
