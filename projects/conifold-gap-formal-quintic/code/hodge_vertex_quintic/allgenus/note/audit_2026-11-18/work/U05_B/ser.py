# Truncated Laurent series in u with exact rational coefficients.
# S(c, prec): sum_e c[e] u^e + O(u^prec).  Independent sanity-check code (referee U05_B).
from fractions import Fraction as Fr


class S:
    __slots__ = ("c", "prec")

    def __init__(self, c, prec):
        self.c = {e: Fr(x) for e, x in c.items() if x != 0 and e < prec}
        self.prec = prec

    @staticmethod
    def const(x, prec):
        return S({0: Fr(x)}, prec)

    @staticmethod
    def u(prec):
        return S({1: Fr(1)}, prec)

    def val(self):
        return min(self.c) if self.c else self.prec

    def iszero(self):
        return not self.c

    def __add__(self, o):
        if not isinstance(o, S):
            o = S.const(o, self.prec)
        p = min(self.prec, o.prec)
        d = dict(self.c)
        for e, x in o.c.items():
            d[e] = d.get(e, 0) + x
        return S(d, p)

    __radd__ = __add__

    def __neg__(self):
        return S({e: -x for e, x in self.c.items()}, self.prec)

    def __sub__(self, o):
        if not isinstance(o, S):
            o = S.const(o, self.prec)
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, S):
            o = Fr(o)
            return S({e: x * o for e, x in self.c.items()}, self.prec)
        p = min(self.prec + o.val(), o.prec + self.val())
        d = {}
        for e1, x1 in self.c.items():
            for e2, x2 in o.c.items():
                e = e1 + e2
                if e < p:
                    d[e] = d.get(e, 0) + x1 * x2
        return S(d, p)

    __rmul__ = __mul__

    def inv(self):
        assert self.c, "inverse of zero series"
        v = self.val()
        r = self.prec - v  # relative precision
        a0 = self.c[v]
        # unit part a = a0*(1 + sum_{i>=1} a_i u^i)
        a = [self.c.get(v + i, Fr(0)) / a0 for i in range(r)]
        b = [Fr(0)] * r
        b[0] = Fr(1)
        for n in range(1, r):
            b[n] = -sum(a[i] * b[n - i] for i in range(1, n + 1))
        return S({-v + i: b[i] / a0 for i in range(r)}, -v + r)

    def __truediv__(self, o):
        if not isinstance(o, S):
            return self * (Fr(1) / Fr(o))
        return self * o.inv()

    def __rtruediv__(self, o):
        return self.inv() * o

    def d(self):
        return S({e - 1: e * x for e, x in self.c.items() if e != 0}, self.prec - 1)

    def integ(self):
        assert self.c.get(-1, 0) == 0, "logarithm would arise"
        assert self.val() >= 0 or True
        return S({e + 1: x / (e + 1) for e, x in self.c.items()}, self.prec + 1)

    def coef(self, e):
        assert e < self.prec, "coefficient beyond precision"
        return self.c.get(e, Fr(0))

    def powfr(self, a):
        """f^a for f = 1 + O(u), a rational (binomial series)."""
        assert self.val() == 0 and self.c[0] == 1
        N = self.prec
        f = [self.c.get(i, Fr(0)) for i in range(N)]
        g = [Fr(0)] * N
        g[0] = Fr(1)
        a = Fr(a)
        # f g' = a f' g  ->  n g_n = sum_{i=1}^{n} (a*i - (n-i)) f_i g_{n-i}
        for n in range(1, N):
            g[n] = sum((a * i - (n - i)) * f[i] * g[n - i] for i in range(1, n + 1)) / n
        return S({i: g[i] for i in range(N)}, N)

    def trunc(self, p):
        return S(self.c, min(p, self.prec))

    def __repr__(self):
        ks = sorted(self.c)
        return " + ".join(f"({self.c[k]})u^{k}" for k in ks[:8]) + f" + O(u^{self.prec})"
