# -*- coding: utf-8 -*-
"""cy_frob.py -- referee cyc_adversary (Stage 5): PF4 (optionally perturbed by polynomials in u), the closed form
PF4 u^x = Sum_j P_j(x) u^(x+j) (P_j polynomials in the symbol x), and from it the Frobenius basis at u = 0 by a recurrence;
d/dx u^x = u^x log u gives the log solution and its obstruction (a method different from the critic's and pb's codes,
which apply the operator to truncated series).  Split off cy_core.py (same code) to keep every file under 250 lines."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cy_core import Ser


def pf4_poly(pert=None):
    """A_k as dicts {power of u: coefficient}; pert: {k: {m: x}} is added"""
    A = {4: {1: Fr(1)}, 3: {0: Fr(-2), 1: Fr(2)}, 2: {0: Fr(-7, 5), 1: Fr(7, 5)}, 1: {0: Fr(-2, 5), 1: Fr(2, 5)},
         0: {0: Fr(-24, 625), 1: Fr(24, 625)}}
    for k, d in (pert or {}).items():
        for m, x in d.items(): A[k][m] = A[k].get(m, 0) + Fr(x)
    return A


def pf4_ser(A): return [Ser.poly(A[k]) for k in range(5)]


def _pn(p):
    p = list(p)
    while p and p[-1] == 0: p.pop()
    return p


def padd(a, b): return _pn([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(max(len(a), len(b)))])
def pscale(a, s): return _pn([s * x for x in a])
def pder(a): return _pn([i * a[i] for i in range(1, len(a))])


def pmul(a, b):
    if not a or not b: return []
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i + j] += x * y
    return _pn(out)


def peval(a, x):
    r = Fr(0)
    for c in reversed(a): r = r * x + c
    return r


def closed_form(A):
    """{j: P_j}: PF4 u^x = Sum_j P_j(x) u^(x+j), from D u^y = y u^y - y u^(y-1)"""
    tot = {}; T = {0: [Fr(1)]}
    for k in range(5):
        if k:
            new = {}
            for j, c in T.items():
                t = pmul([Fr(j), Fr(1)], c)
                new[j] = padd(new.get(j, []), t); new[j - 1] = padd(new.get(j - 1, []), pscale(t, -1))
            T = new
        for m, a in A[k].items():
            for j, c in T.items(): tot[j + m] = padd(tot.get(j + m, []), pscale(c, a))
    return {j: p for j, p in tot.items() if p}


class Frob(object):
    """Frobenius data at u = 0 of PF4 (or a perturbation A), coefficients through u^(N-1)"""
    def __init__(self, A, N):
        self.A = A; self.N = N; self.P = closed_form(A); self.dP = {j: pder(p) for j, p in self.P.items()}
        assert min(self.P) == -3
        self.jmax = max(self.P)

    def piv(self, n): return peval(self.P[-3], n)

    def _rhs(self, c, m, extra=None):
        s = Fr(0)
        for j, p in self.P.items():
            n = m - j
            if j > -3 and 0 <= n < len(c): s += peval(p, n) * c[n]
        if extra is not None:
            for j, p in self.dP.items():
                n = m - j
                if 0 <= n < len(extra): s += peval(p, n) * extra[n]
        return s

    def series(self, init):
        """power-series solution with (b_0, b_1, b_2) = init (list of N coefficients)"""
        b = [Fr(x) for x in init]
        for m in range(0, self.N - 3): b.append(-self._rhs(b, m) / self.piv(m + 3))
        return b

    def low_orders(self, c):
        """coefficients of u^-3, u^-2, u^-1 of PF4(Sum c_n u^n) (must vanish for every power series)"""
        return [self._rhs(c, m) + (peval(self.P[-3], m + 3) * c[m + 3] if 0 <= m + 3 < len(c) else 0) for m in (-3, -2, -1)]

    def obstruction(self, b):
        """u^-3, u^-2, u^-1 coefficients of PF4(g log u) for g = Sum b_n u^n"""
        return [sum((peval(p, m - j) * b[m - j] for j, p in self.dP.items() if 0 <= m - j < len(b)), Fr(0)) for m in (-3, -2, -1)]

    def logpart(self, b):
        """h (h_0 = h_1 = h_2 = 0) with PF4(g log u + h) = 0, assuming the obstruction of g vanishes"""
        c = [Fr(0)] * 3
        for m in range(0, self.N - 3): c.append(-self._rhs(c, m, extra=b) / self.piv(m + 3))
        return c


def frobenius(A, N):
    """returns dict: F (Frob), f0, f1, f2 (coefficient lists), gamma, g, h, obstructions"""
    F = Frob(A, N)
    f = [F.series([1 if i == k else 0 for i in range(3)]) for k in range(3)]
    ob = [F.obstruction(x) for x in f]
    gamma = -ob[1][2] / ob[2][2]
    g = [x + gamma * y for x, y in zip(f[1], f[2])]
    return dict(F=F, f=f, ob=ob, gamma=gamma, g=g, h=F.logpart(g), obg=F.obstruction(g))
