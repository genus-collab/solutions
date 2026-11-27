# -*- coding: utf-8 -*-
"""lm_core.py -- referee lmi_adversary (Stage 5): core arithmetic, written from scratch (stdlib only).

* Field: Q (Fraction) or F_p (int mod p), chosen by setmod(p) (None = Q).
* LS: Laurent series in the variable e = L - s* with ABSOLUTE precision: all coefficients of e^k with k < prec are known
  exactly (prec = INF: exact Laurent polynomial).  No rescaling of sigma or z is ever made: the Wick sum below is the
  definition of CP s.1.2 evaluated at the critical point L = s* + e, with every Taylor datum an element of LS.
* The Bernoulli numbers, the graded exponential gexp and the Wick engine wick1 are in lm_wick.py.
"""
from fractions import Fraction as Fr

INF = 1 << 40
_M = [None]


def setmod(p):
    _M[0] = p


def getmod():
    return _M[0]


def F(x):
    """rational (int or Fraction) -> field element"""
    p = _M[0]
    if p is None:
        return Fr(x)
    x = Fr(x)
    return x.numerator % p * pow(x.denominator % p, p - 2, p) % p


def finv(x):
    p = _M[0]
    if p is None:
        if x == 0:
            raise ZeroDivisionError('finv(0)')
        return 1 / Fr(x)
    x %= p
    if x == 0:
        raise ZeroDivisionError('finv(0) mod p')
    return pow(x, p - 2, p)


def fpow(x, n):
    p = _M[0]
    if p is None:
        return Fr(x) ** n
    return pow(x % p, n, p) if n >= 0 else pow(finv(x), -n, p)


def fz(x):
    p = _M[0]
    return x % p == 0 if p is not None else x == 0


def fr(x):
    p = _M[0]
    return x % p if p is not None else x


class LS:
    """Sum_i c[i] e^(v+i); coefficients of e^k known for k < prec."""
    __slots__ = ('v', 'c', 'prec')

    def __init__(self, v, c, prec=INF):
        if prec < INF and len(c) > prec - v:
            c = c[:max(0, prec - v)]
        i = 0
        while i < len(c) and fz(c[i]):
            i += 1
        if i:
            c = c[i:]
            v += i
        if not c:
            v = prec if prec < INF else 0
        self.v, self.c, self.prec = v, list(c), prec

    def ezero(self):
        return not self.c and self.prec >= INF

    def coef(self, k):
        if k >= self.prec:
            raise ValueError('coefficient e^%d not known (prec %d)' % (k, self.prec))
        i = k - self.v
        return self.c[i] if 0 <= i < len(self.c) else F(0)

    def __add__(a, b):
        if not isinstance(b, LS):
            b = LS(0, [b])
        prec = min(a.prec, b.prec)
        parts = [x for x in (a, b) if x.c]
        if not parts:
            return LS(0, [], prec)
        lo = min(x.v for x in parts)
        hi = max(x.v + len(x.c) for x in parts)
        if prec < INF:
            hi = min(hi, prec)
        if hi <= lo:
            return LS(0, [], prec)
        c = [0] * (hi - lo)
        for x in parts:
            o = x.v - lo
            for i in range(min(len(x.c), len(c) - o)):
                c[o + i] += x.c[i]
        return LS(lo, [fr(y) for y in c], prec)

    __radd__ = __add__

    def __neg__(a):
        return LS(a.v, [fr(-y) for y in a.c], a.prec)

    def __sub__(a, b):
        return a + (-b if isinstance(b, LS) else LS(0, [fr(-b)]))

    def scal(a, s):
        if fz(s):
            return LS(0, [])
        return LS(a.v, [fr(s * y) for y in a.c], a.prec)

    def __mul__(a, b):
        if not isinstance(b, LS):
            return a.scal(b)
        if a.ezero() or b.ezero():
            return LS(0, [])
        prec = INF
        if b.prec < INF:
            prec = min(prec, a.v + b.prec)
        if a.prec < INF:
            prec = min(prec, b.v + a.prec)
        if not a.c or not b.c:
            return LS(0, [], prec)
        v = a.v + b.v
        n = len(a.c) + len(b.c) - 1
        if prec < INF:
            n = min(n, prec - v)
        if n <= 0:
            return LS(0, [], prec)
        c = [0] * n
        bc = b.c
        lb = len(bc)
        for i, x in enumerate(a.c):
            if i >= n:
                break
            m = min(lb, n - i)
            for j in range(m):
                c[i + j] += x * bc[j]
        return LS(v, [fr(y) for y in c], prec)

    __rmul__ = __mul__

    def inv(a, rel=None):
        """1/a; relative precision = that of a, or rel if a is exact."""
        if not a.c:
            raise ZeroDivisionError('LS.inv of zero')
        rp = (a.prec - a.v) if a.prec < INF else rel
        c0 = finv(a.c[0])
        ac, la = a.c, len(a.c)
        b = [0] * rp
        b[0] = c0
        for k in range(1, rp):
            s = 0
            for i in range(1, min(k, la - 1) + 1):
                s += ac[i] * b[k - i]
            b[k] = fr(-s * c0)
        return LS(-a.v, b, -a.v + rp)

    def trunc(a, prec):
        return LS(a.v, a.c, min(prec, a.prec))


def lsc(x, prec=INF):
    """constant (rational or field element given as rational)"""
    return LS(0, [F(x)], prec)


def lsm(x, k, prec=INF):
    """x e^k"""
    return LS(k, [F(x)], prec)


def lsexp(a, prec=None):
    """exp(a) for a of valuation >= 1; precision a.prec (or prec if a is exact)"""
    P_ = a.prec if a.prec < INF else prec
    if a.c and a.v < 1:
        raise ValueError('lsexp needs valuation >= 1')
    A = [0] * P_
    for i, y in enumerate(a.c):
        if a.v + i < P_:
            A[a.v + i] = y
    E = [0] * P_
    E[0] = F(1)
    for k in range(1, P_):
        s = 0
        for i in range(1, k + 1):
            if A[i]:
                s += i * A[i] * E[k - i]
        E[k] = fr(s * finv(F(k)))
    return LS(0, E, P_)
