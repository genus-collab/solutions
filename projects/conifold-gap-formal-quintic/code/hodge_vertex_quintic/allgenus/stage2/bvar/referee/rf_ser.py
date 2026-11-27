# -*- coding: utf-8 -*-
"""rf_ser.py -- referee's own truncated Laurent series over Q (Fractions), independent of fq_series and of the worker's bv_* modules.

L(v, c, P): sum_{i} c[i] u^(v+i) + O(u^P)  (P = absolute precision; c may have leading zeros, they are stripped).
Also: residues mod p (fr_mod), power-series helpers on plain lists mod p (pmul, pinv, ppow).
"""
from fractions import Fraction as Fr


class L(object):
    __slots__ = ("v", "c", "P")

    def __init__(self, v, c, P):
        c = [Fr(x) for x in c[:max(P - v, 0)]]
        i = 0
        while i < len(c) and c[i] == 0: i += 1
        self.v = v + i; self.c = c[i:]; self.P = P
        if not self.c: self.v = P

    @staticmethod
    def poly(lst, P): return L(0, list(lst) + [0] * max(0, P - len(lst)), P)

    @staticmethod
    def const(a, P): return L(0, [a] + [0] * (P - 1), P)

    def zero(self): return not self.c
    def co(self, k):
        assert k < self.P, "coefficient u^%d beyond O(u^%d)" % (k, self.P)
        return self.c[k - self.v] if self.v <= k < self.v + len(self.c) else Fr(0)

    def add(self, o):
        P = min(self.P, o.P); v = min(self.v, o.v)
        return L(v, [self.co(k) + o.co(k) for k in range(v, P)], P)

    def neg(self): return L(self.v, [-x for x in self.c], self.P)
    def sub(self, o): return self.add(o.neg())
    def sc(self, a): return L(self.v, [a * x for x in self.c], self.P) if a != 0 else L(self.P, [], self.P)
    def sh(self, k): return L(self.v + k, self.c, self.P + k)

    def mul(self, o):
        if not self.c or not o.c:
            P = min(self.P + (o.v if o.c else o.P), o.P + (self.v if self.c else self.P))
            return L(P, [], P)
        P = min(self.P + o.v, o.P + self.v); v = self.v + o.v; n = P - v; out = [Fr(0)] * max(n, 0)
        for i, x in enumerate(self.c[:n]):
            if x:
                for j, y in enumerate(o.c[:n - i]): out[i + j] += x * y
        return L(v, out, P)

    def inv(self):
        assert self.c, "inverse of zero"
        n = len(self.c); a = self.c; out = [Fr(0)] * n; out[0] = 1 / a[0]
        for k in range(1, n):
            out[k] = -sum(a[j] * out[k - j] for j in range(1, k + 1)) * out[0]
        return L(-self.v, out, -self.v + n)

    def der(self):
        return L(self.v - 1, [(self.v + i) * x for i, x in enumerate(self.c)], self.P - 1)

    def integ(self):
        """antiderivative with zero constant (asserts no u^-1 term)"""
        assert self.co(-1) == 0 if self.P > -1 else True, "log term in an integral"
        return L(self.v + 1, [x / (self.v + i + 1) for i, x in enumerate(self.c)], self.P + 1)

    def powfrac(self, e, lead):
        """self^e for a power series with self.v = 0, given lead = (c0)^e (rational): binomial recursion y = c0^e (1 + x)^e"""
        assert self.v == 0 and self.c
        c0 = self.c[0]; x = [a / c0 for a in self.c]; n = len(x); y = [Fr(0)] * n; y[0] = Fr(1)
        # y' (1 + x) = e x' y  (coefficient recursion)
        for k in range(1, n):
            s = Fr(0)
            for j in range(1, k + 1):
                s += (e * j - (k - j)) * x[j] * y[k - j]
            y[k] = s / k
        return L(0, [lead * t for t in y], self.P)

    def compose(self, g):
        """self(g), self a power series, g with g.v >= 1"""
        assert self.v >= 0 and g.v >= 1
        P = min(g.P, self.P * g.v if self.P > 0 else g.P); res = L(P, [], P); one = L.const(1, P)
        for k in range(self.P - 1, -1, -1):
            res = res.mul(g).add(one.sc(self.co(k)))
            res = L(res.v, res.c, min(res.P, P))
        return res

    def revert(self):
        """compositional inverse of a u + O(u^2), a != 0"""
        a = self.co(1); assert self.v == 1 and a != 0
        P = self.P; g = L(1, [1 / a] + [0] * (P - 2), P)
        for n in range(2, P):
            d = self.compose(g).co(n)
            g = L(1, [g.co(k) for k in range(1, n)] + [-d / a] + [0] * (P - n - 1), P)
        return g

    def __repr__(self):
        t = ["(%s)u^%d" % (x, self.v + i) for i, x in enumerate(self.c[:5]) if x]
        return " + ".join(t) + " + O(u^%d)" % self.P


def fr_mod(x, p):
    x = Fr(x); return x.numerator % p * pow(x.denominator % p, p - 2, p) % p


def pmul(a, b, K, p):
    out = [0] * K
    for i in range(K):
        if a[i]:
            for j in range(K - i): out[i + j] = (out[i + j] + a[i] * b[j]) % p
    return out


def pinv(a, K, p):
    out = [0] * K; i0 = pow(a[0], p - 2, p); out[0] = i0
    for k in range(1, K):
        out[k] = -sum(a[j] * out[k - j] for j in range(1, k + 1)) * i0 % p
    return out


def ppow(a, n, K, p):
    if n < 0: a = pinv(a, K, p); n = -n
    out = [1] + [0] * (K - 1)
    for _ in range(n): out = pmul(out, a, K, p)
    return out
