# -*- coding: utf-8 -*-
"""rf_rings.py -- Stage 5 / chart / referee_formal: two coefficient rings over Q for the formal SP code of rf_core.

LS: truncated Laurent series in eps, value = Sum_i c[i] eps^(v+i), coefficients known for exponents < p (absolute precision).
    Products, sums and inverses propagate p rigorously (p = INF for exact elements), so a final coefficient is reported only if
    it is known; this is used to expand c_n(L) at L = 1 + eps WITHOUT assuming the shape of c_n (poles must cancel).
TS: Q[[u]]/(u^n), n = TS.n (global), for the hatted u-series and for Taylor expansions in a parameter.
"""
from fractions import Fraction as Fr

INF = 10 ** 9


class LS:
    __slots__ = ('v', 'c', 'p')
    CAP = 40                                       # relative precision of inverses of exact non-monomials

    def __init__(s, v, c, p=INF):
        n = max(0, min(len(c), p - v)); i = 0
        while i < n and c[i] == 0: i += 1
        s.c = [Fr(x) for x in c[i:n]]; s.p = min(p, INF)
        s.v = v + i if s.c else s.p

    @staticmethod
    def of(x):
        x = Fr(x); return LS(0, [x]) if x else LS(INF, [], INF)

    @staticmethod
    def mono(k, a=1): return LS(k, [Fr(a)])            # exact a eps^k

    def __add__(s, o):
        if not isinstance(o, LS): o = LS.of(o)
        p = min(s.p, o.p)
        if not s.c: return LS(o.v, o.c, p)
        if not o.c: return LS(s.v, s.c, p)
        v = min(s.v, o.v); n = min(p, max(s.v + len(s.c), o.v + len(o.c))) - v
        r = [Fr(0)] * max(n, 0)
        for x in (s, o):
            off = x.v - v
            for i, a in enumerate(x.c):
                if off + i < n: r[off + i] += a
        return LS(v, r, p)
    __radd__ = __add__

    def __neg__(s): return LS(s.v, [-a for a in s.c], s.p)
    def __sub__(s, o): return s + (-o if isinstance(o, LS) else LS.of(-Fr(o)))
    def __rsub__(s, o): return (-s) + o

    def __mul__(s, o):
        if not isinstance(o, LS):
            o = Fr(o)
            return LS(s.v, [a * o for a in s.c], s.p) if o else LS(INF, [], INF)
        p = min(s.v + o.p, o.v + s.p, INF)
        if not s.c or not o.c: return LS(p, [], p)
        v = s.v + o.v; n = min(p - v, len(s.c) + len(o.c) - 1)
        r = [Fr(0)] * n
        for i, a in enumerate(s.c):
            if i >= n: break
            for j in range(min(len(o.c), n - i)): r[i + j] += a * o.c[j]
        return LS(v, r, p)
    __rmul__ = __mul__

    def inv(s):
        assert s.c, 'LS: inverse of a (possibly inexact) zero'
        rel = s.p - s.v
        if len(s.c) == 1 and s.p >= INF: return LS(-s.v, [1 / s.c[0]])
        n = min(rel, LS.CAP); a = s.c; b = [1 / a[0]]
        for k in range(1, n):
            t = Fr(0)
            for i in range(1, min(k, len(a) - 1) + 1): t += a[i] * b[k - i]
            b.append(-t * b[0])
        return LS(-s.v, b, -s.v + n)

    def __truediv__(s, o): return s * (o.inv() if isinstance(o, LS) else 1 / Fr(o))
    def __rtruediv__(s, o): return LS.of(o) * s.inv()

    def coef(s, k):
        assert k < s.p, 'LS: coefficient eps^%d not known (precision %d)' % (k, s.p)
        return s.c[k - s.v] if s.v <= k < s.v + len(s.c) else Fr(0)

    def __repr__(s): return 'LS(v=%d,p=%d,%s)' % (s.v, s.p, [str(x) for x in s.c[:4]])


class TS:
    __slots__ = ('c',)
    n = 10

    def __init__(s, c):
        c = list(c)[:TS.n]; s.c = [Fr(x) for x in c] + [Fr(0)] * (TS.n - len(c))

    @staticmethod
    def of(x): return x if isinstance(x, TS) else TS([Fr(x)])

    @staticmethod
    def var(): return TS([0, 1])

    def __add__(s, o):
        if not isinstance(o, TS): o = TS.of(o)
        return TS([a + b for a, b in zip(s.c, o.c)])
    __radd__ = __add__
    def __neg__(s): return TS([-a for a in s.c])
    def __sub__(s, o): return s + (-TS.of(o))
    def __rsub__(s, o): return TS.of(o) + (-s)

    def __mul__(s, o):
        if not isinstance(o, TS):
            o = Fr(o); return TS([a * o for a in s.c])
        n = TS.n; r = [Fr(0)] * n; oc = o.c
        for i, a in enumerate(s.c):
            if a:
                for j in range(n - i):
                    if oc[j]: r[i + j] += a * oc[j]
        return TS(r)
    __rmul__ = __mul__

    def inv(s):
        a = s.c; assert a[0] != 0, 'TS: not invertible'
        b = [1 / a[0]]
        for k in range(1, TS.n):
            b.append(-sum(a[i] * b[k - i] for i in range(1, k + 1)) * b[0])
        return TS(b)

    def __truediv__(s, o): return s * (o.inv() if isinstance(o, TS) else 1 / Fr(o))
    def __rtruediv__(s, o): return TS.of(o) * s.inv()
    def __eq__(s, o): return s.c == TS.of(o).c
    def __hash__(s): return hash(tuple(s.c))
    def __repr__(s): return 'TS' + str([str(x) for x in s.c[:6]])
