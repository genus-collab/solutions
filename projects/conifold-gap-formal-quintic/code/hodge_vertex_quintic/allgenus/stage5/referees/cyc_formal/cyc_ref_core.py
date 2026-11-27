# -*- coding: utf-8 -*-
"""cyc_ref_core.py -- Stage-5 referee 'cyc_formal': exact machinery (Python stdlib only).

Written from scratch for refereeing critic/CYC_PROOF.md; imports nothing from the project or from any agent
(critic/cr5_cyc.py and pb/pb_cyc.py were read, never imported or copied).
  Ser   truncated Laurent series in u: dense list c for the exponents v, v+1, ..; absolute precision N (exponents >= N
        unknown; N = EX means exact).  Coefficients: Fraction, or P2 (polynomials in the frame symbols alpha, beta).
  P2    polynomial in (alpha, beta) over Q, for the symbolic-frame check.
  ONE, U, DL = 1, u, D log u = 1 - 1/u.  REl and Op are in cyc_ref_ops.py; jets, PN and PF4 in cyc_ref_jet.py.
D = -(1-u) d/du throughout:  D u^e = -e u^(e-1) + e u^e.
"""
from fractions import Fraction as Fr
from math import comb

EX = 1 << 40


def _z(x):
    return x == 0


def _inv(c):
    return c.inv() if isinstance(c, P2) else 1 / Fr(c)


class P2:
    """polynomial in (alpha, beta) over Q: dict (i, j) -> Fraction"""
    __slots__ = ('d',)
    __hash__ = None

    def __init__(self, d=None):
        self.d = {k: Fr(v) for k, v in (d or {}).items() if v != 0}

    @staticmethod
    def lift(x):
        return x if isinstance(x, P2) else P2({(0, 0): x})

    def __add__(a, b):
        b = P2.lift(b); d = dict(a.d)
        for k, v in b.d.items():
            d[k] = d.get(k, 0) + v
        return P2(d)
    __radd__ = __add__

    def __neg__(a):
        return P2({k: -v for k, v in a.d.items()})

    def __sub__(a, b):
        return a + (-P2.lift(b))

    def __rsub__(a, b):
        return P2.lift(b) + (-a)

    def __mul__(a, b):
        b = P2.lift(b); d = {}
        for (i, j), x in a.d.items():
            for (k, m), y in b.d.items():
                d[(i + k, j + m)] = d.get((i + k, j + m), 0) + x * y
        return P2(d)
    __rmul__ = __mul__

    def __eq__(a, b):
        return a.d == P2.lift(b).d

    def inv(a):
        assert set(a.d) == {(0, 0)}, 'P2: only non-zero constants are inverted'
        return P2({(0, 0): 1 / a.d[(0, 0)]})


class Ser:
    __slots__ = ('v', 'c', 'N')

    def __init__(self, v, c, N=EX):
        c = list(c)
        if N < EX:
            c = c[:max(0, N - v)]
        i = 0
        while i < len(c) and _z(c[i]):
            i += 1
        c = c[i:]; v += i
        while c and _z(c[-1]):
            c.pop()
        self.v, self.c, self.N = (v if c else 0), c, N

    @staticmethod
    def poly(d, N=EX):                       # {exponent: coefficient}
        d = {e: x for e, x in d.items() if not _z(x)}
        if not d:
            return Ser(0, [], N)
        lo, hi = min(d), max(d)
        return Ser(lo, [d.get(e, Fr(0)) for e in range(lo, hi + 1)], N)

    def val(self):
        return self.v if self.c else self.N

    def coef(self, e):
        assert e < self.N, ('unknown coefficient', e, self.N)
        k = e - self.v
        return self.c[k] if 0 <= k < len(self.c) else Fr(0)

    def __add__(a, b):
        N = min(a.N, b.N)
        if not a.c or not b.c:
            s = b if not a.c else a
            return Ser(s.v, s.c, N)
        v = min(a.v, b.v); hi = max(a.v + len(a.c), b.v + len(b.c))
        if N < EX:
            hi = min(hi, N)
        out = [Fr(0)] * max(0, hi - v)
        for s in (a, b):
            for k, x in enumerate(s.c):
                if s.v + k < hi:
                    out[s.v + k - v] = out[s.v + k - v] + x
        return Ser(v, out, N)

    def __neg__(a):
        return Ser(a.v, [-x for x in a.c], a.N)

    def __sub__(a, b):
        return a + (-b)

    def sc(a, s):
        return Ser(a.v, [s * x for x in a.c], a.N)

    def __mul__(a, b):
        N = min(a.N + b.val(), b.N + a.val())
        if N > EX // 2:
            N = EX
        if not a.c or not b.c:
            return Ser(0, [], N)
        v = a.v + b.v; L = len(a.c) + len(b.c) - 1
        if N < EX:
            L = min(L, N - v)
        out = [Fr(0)] * max(0, L)
        for i, x in enumerate(a.c):
            if _z(x):
                continue
            for j in range(min(len(b.c), L - i)):
                out[i + j] = out[i + j] + x * b.c[j]
        return Ser(v, out, N)

    def D(a):
        N = a.N - 1 if a.N < EX else EX
        out = [Fr(0)] * (len(a.c) + 1)
        for k, x in enumerate(a.c):
            e = a.v + k
            if e:
                out[k] = out[k] - e * x
                out[k + 1] = out[k + 1] + e * x
        return Ser(a.v - 1, out, N)

    def inv(a, N=None):
        assert a.c and a.v == 0, 'not a unit'
        N = a.N if N is None else min(N, a.N)
        assert N < EX
        i0 = _inv(a.c[0]); r = [i0]
        for n in range(1, N):
            s = Fr(0)
            for k in range(1, min(n, len(a.c) - 1) + 1):
                s = s + a.c[k] * r[n - k]
            r.append(-(s * i0))
        return Ser(0, r, N)

    def trunc(a, N):
        return Ser(a.v, a.c, min(N, a.N))

    def zero_below(a, M):                    # every coefficient of u^e, e < M, is known and zero
        assert a.N >= M, ('precision', a.N, M)
        return all(_z(x) for k, x in enumerate(a.c) if a.v + k < M)


ONE = Ser(0, [Fr(1)])
U = Ser(1, [Fr(1)])
DL = Ser(-1, [Fr(-1), Fr(1)])                # D log u = 1 - 1/u
