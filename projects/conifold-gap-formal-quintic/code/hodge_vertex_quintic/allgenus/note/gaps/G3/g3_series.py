# -*- coding: utf-8 -*-
"""g3_series.py -- truncated power / Laurent series in u over Q (exact), and the model coordinate v(u) = u tau^(-5/4).

A series mod u^N is a list of N Fractions.  A Laurent series is a pair (o, list) meaning u^o * Sum list[i] u^i, known mod
u^(o + len(list)).  tau(u) = 1/2F1(1,4/5;9/5;u), 2F1(1,4/5;9/5;u) = Sum_n 4 u^n/(5n+4), so v(u) = u * F(u)^(5/4).
Imported by g3_decomp.py (gap cluster G3).  Nothing here writes a file.
"""
from fractions import Fraction as Fr


def zero(N):
    return [Fr(0)] * N


def one(N):
    s = zero(N)
    s[0] = Fr(1)
    return s


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def scal(c, a):
    return [c * x for x in a]


def mul(a, b):
    N = len(a)
    out = zero(N)
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                if b[j]:
                    out[i + j] += x * b[j]
    return out


def ordu(a):
    for i, x in enumerate(a):
        if x:
            return i
    return None


def binom_pow(f, alpha):
    """f^alpha for f(0) = 1 (binomial series), alpha in Q."""
    N = len(f)
    x = f[:]
    x[0] -= 1
    out, term, c = one(N), one(N), Fr(1)
    for m in range(1, N):
        c = c * (alpha - m + 1) / m
        term = mul(term, x)
        out = add(out, scal(c, term))
    return out


def hyp_F(N):
    return [Fr(4, 5 * n + 4) for n in range(N)]


def v_unit(N, kind="model"):
    """the unit w(u) with v(u) = u * w(u) ("model": w = F^(5/4)); kind "shift" returns the control v(u) = 1 + u (no u factor)."""
    if kind == "model":
        return binom_pow(hyp_F(N), Fr(5, 4))
    raise ValueError(kind)


class VPowers:
    """v(u)^e for integers e >= 0 (and e < 0 for the model coordinate), as Laurent series (o, list) of length N."""

    def __init__(self, N, kind="model"):
        self.N, self.kind = N, kind
        if kind == "model":
            self.w = v_unit(N + 8)[:N]
        elif kind == "shift":                      # control: v(u) = 1 + u, v(0) != 0
            self.w = [Fr(1), Fr(1)] + [Fr(0)] * (N - 2)
        self.cache = {}

    def get(self, e):
        if e in self.cache:
            return self.cache[e]
        N = self.N
        if e > 0:
            o, prev = self.get(e - 1)
            s = mul(prev, self.w)
        elif e == 0:
            s = one(N)
        else:
            s = binom_pow(self.w, Fr(e))
        out = (e, s) if self.kind == "model" else (0, s)
        self.cache[e] = out
        return out

    def trunc(self, e, N):
        """v(u)^e as an ordinary series mod u^N (requires e >= 0)."""
        o, s = self.get(e)
        out = zero(N)
        for i, x in enumerate(s):
            if o + i < N:
                out[o + i] += x
        return out


def lmul(A, B):
    """product of Laurent series of the same length."""
    (oa, a), (ob, b) = A, B
    return (oa + ob, mul(a, b))


def ladd(A, B, N):
    """sum of Laurent series, truncated to absolute order < N (returns ordinary list indexed from the min offset)."""
    (oa, a), (ob, b) = A, B
    o = min(oa, ob)
    out = [Fr(0)] * (N - o)
    for (oo, s) in ((oa, a), (ob, b)):
        for i, x in enumerate(s):
            if oo + i < N:
                out[oo + i - o] += x
    return (o, out)
