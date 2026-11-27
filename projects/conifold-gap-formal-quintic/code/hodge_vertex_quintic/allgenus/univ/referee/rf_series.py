# -*- coding: utf-8 -*-
"""rf_series.py -- REFEREE (univ): exact truncated Laurent series in u over Q with explicit absolute precision.

LS(o, c, prec): coefficients c[i] of u^(o+i) for o+i < prec (len(c) = prec - o); every coefficient below prec is KNOWN.
Operations propagate the absolute precision exactly (mul: min(a.v + b.prec, b.v + a.prec) with true valuations; inv: prec - 2v;
deriv: prec - 1).  No module of the worker is imported; nothing of the project either.
"""
from fractions import Fraction as Fr


class LS(object):
    __slots__ = ("o", "c", "prec")

    def __init__(self, o, c, prec):
        c = [Fr(x) for x in c][:max(0, prec - o)]
        c += [Fr(0)] * (prec - o - len(c))
        self.o = o; self.c = c; self.prec = prec
        self._norm()

    def _norm(self):
        i = 0
        while i < len(self.c) and self.c[i] == 0: i += 1
        if i: self.o += i; self.c = self.c[i:]

    @staticmethod
    def poly(coefs, prec, o=0): return LS(o, list(coefs), prec)

    @staticmethod
    def const(x, prec): return LS(0, [x], prec)

    @staticmethod
    def zero(prec): return LS(prec, [], prec)

    def iszero(self): return not self.c

    def v(self): return self.o if self.c else self.prec

    def coef(self, n):
        assert n < self.prec, "coefficient u^%d beyond precision %d" % (n, self.prec)
        return self.c[n - self.o] if n >= self.o else Fr(0)

    def trunc(self, prec): return LS(self.o, self.c, min(prec, self.prec))

    def add(self, b):
        pr = min(self.prec, b.prec); o = min(self.v(), b.v(), pr)
        return LS(o, [self.coef(n) + b.coef(n) if n < pr else 0 for n in range(o, pr)], pr)

    def neg(self): return LS(self.o, [-x for x in self.c], self.prec)

    def sub(self, b): return self.add(b.neg())

    def scale(self, x):
        x = Fr(x)
        if x == 0: return LS.zero(self.prec)
        return LS(self.o, [y * x for y in self.c], self.prec)

    def shift(self, k): return LS(self.o + k, self.c, self.prec + k)          # times u^k

    def mul(self, b):
        va, vb = self.v(), b.v(); pr = min(va + b.prec, vb + self.prec)
        if self.iszero() or b.iszero(): return LS.zero(pr)
        o = va + vb; n = pr - o; out = [Fr(0)] * max(n, 0); A = self.c; B = b.c
        for i in range(min(n, len(A))):
            x = A[i]
            if x:
                for j in range(min(n - i, len(B))):
                    y = B[j]
                    if y: out[i + j] += x * y
        return LS(o, out, pr)

    def inv(self):
        assert not self.iszero(), "inverse of zero"
        v = self.v(); r = self.prec - v; A = self.c; a0 = A[0]; out = [1 / a0]
        for n in range(1, r):
            out.append(-sum(A[j] * out[n - j] for j in range(1, min(n, len(A) - 1) + 1)) / a0)
        return LS(-v, out, -v + r)

    def div(self, b): return self.mul(b.inv())

    def deriv(self):
        if self.iszero(): return LS.zero(self.prec - 1)
        return LS(self.o - 1, [(self.o + i) * x for i, x in enumerate(self.c)], self.prec - 1)

    def integ(self):
        """antiderivative with zero constant term; asserts that there is no u^-1 term"""
        if self.prec > -1: assert self.coef(-1) == 0, "logarithm in integration"
        if self.iszero(): return LS.zero(self.prec + 1)
        o = self.o
        return LS(o + 1, [(x / (o + i + 1) if o + i != -1 else 0) for i, x in enumerate(self.c)], self.prec + 1)

    def log1(self):
        """log of a series 1 + O(u)"""
        assert self.v() == 0 and self.c[0] == 1
        return self.deriv().mul(self.inv()).integ()

    def exp0(self):
        """exp of a series O(u)"""
        assert self.v() >= 1
        N = self.prec; d = self.deriv(); out = [Fr(1)]
        for n in range(1, N):                       # E' = d E  =>  n e_n = Sum_j (j a_j) e_(n-j)
            out.append(sum(d.coef(j - 1) * out[n - j] for j in range(1, n + 1)) / n)
        return LS(0, out, N)

    def powr(self, r):
        """self^r for a series 1 + O(u), r rational"""
        return self.log1().scale(r).exp0()

    def compose(self, w):
        """self(w(u)) for a regular self and w = O(u)"""
        assert self.v() >= 0 and w.v() >= 1
        pr = min(self.prec * w.v(), w.prec)
        acc = LS.zero(pr); pw = LS.const(1, pr)
        for n in range(0, pr):
            if n >= self.prec: break
            if self.coef(n): acc = acc.add(pw.scale(self.coef(n)))
            pw = pw.mul(w).trunc(pr)
            if pw.v() >= pr: break
        return acc.trunc(pr)

    def fmt(self, n=6):
        return " + ".join("(%s) u^%d" % (self.coef(k), k) for k in range(self.v(), min(self.prec, self.v() + n)))


def tau_series(N):
    """tau = 1/2F1(1, 4/5; 9/5; u) = 1/Sum_n 4 u^n/(5n+4)"""
    return LS(0, [Fr(4, 5 * n + 4) for n in range(N)], N).inv()


def model_v(N):
    """v(u) = u tau^(-5/4)"""
    return tau_series(N).powr(Fr(-5, 4)).shift(1)
