"""Truncated Laurent series in one variable with Fraction coefficients.
A series is (v, c) meaning sum_{k>=0} c[k] x^(v+k), known modulo x^(v+len(c)).
"""
from fractions import Fraction as Fr

class S:
    __slots__ = ("v", "c")
    def __init__(self, v, c):
        self.v = v
        self.c = [Fr(x) for x in c]
    @staticmethod
    def const(a, N):
        return S(0, [Fr(a)] + [Fr(0)] * (N - 1))
    @staticmethod
    def mono(n, N, a=1):
        # a*x^n known to absolute precision x^N
        return S(n, [Fr(a)] + [Fr(0)] * (N - n - 1))
    def prec(self):
        return self.v + len(self.c)
    def coef(self, n):
        k = n - self.v
        if k < 0:
            return Fr(0)
        if k >= len(self.c):
            raise ValueError("beyond precision: want %d, prec %d" % (n, self.prec()))
        return self.c[k]
    def __add__(self, o):
        if not isinstance(o, S):
            o = S(0, [Fr(o)] + [Fr(0)] * (max(self.prec(), 1) - 1))
        v = min(self.v, o.v)
        p = min(self.prec(), o.prec())
        return S(v, [self.coef(n) + o.coef(n) for n in range(v, p)])
    __radd__ = __add__
    def __neg__(self):
        return S(self.v, [-x for x in self.c])
    def __sub__(self, o):
        if not isinstance(o, S):
            return self + (-Fr(o))
        return self + (-o)
    def __rsub__(self, o):
        return (-self) + o
    def __mul__(self, o):
        if not isinstance(o, S):
            return S(self.v, [x * Fr(o) for x in self.c])
        n = min(len(self.c), len(o.c))
        c = [Fr(0)] * n
        for i in range(n):
            a = self.c[i]
            if a == 0:
                continue
            for j in range(n - i):
                c[i + j] += a * o.c[j]
        return S(self.v + o.v, c)
    __rmul__ = __mul__
    def norm(self):
        # strip leading zeros
        k = 0
        while k < len(self.c) and self.c[k] == 0:
            k += 1
        return S(self.v + k, self.c[k:])
    def inv(self):
        s = self.norm()
        if not s.c:
            raise ZeroDivisionError
        n = len(s.c)
        a0 = s.c[0]
        b = [Fr(0)] * n
        b[0] = 1 / a0
        for k in range(1, n):
            b[k] = -sum(s.c[j] * b[k - j] for j in range(1, k + 1)) / a0
        return S(-s.v, b)
    def __truediv__(self, o):
        if not isinstance(o, S):
            return self * (Fr(1) / Fr(o))
        return self * o.inv()
    def d(self):
        # d/dx
        return S(self.v - 1, [(self.v + k) * self.c[k] for k in range(len(self.c))])
    def __pow__(self, n):
        r = S(0, [Fr(1)] + [Fr(0)] * (len(self.c) - 1))
        base = self
        if n < 0:
            base = self.inv()
            n = -n
        for _ in range(n):
            r = r * base
        return r
    def iszero(self):
        return all(x == 0 for x in self.c)
    def lowest(self):
        s = self.norm()
        return (s.v, s.c[0]) if s.c else (self.prec(), None)
    def __repr__(self):
        t = []
        for k, x in enumerate(self.c):
            if x != 0:
                t.append("%s*u^%d" % (x, self.v + k))
        return " + ".join(t[:12]) + " + O(u^%d)" % self.prec()
