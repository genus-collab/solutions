# -*- coding: utf-8 -*-
"""rc_lser.py -- minimal exact Laurent series in one variable over Fraction.

Written for the independent re-check of the genus-5 cluster results; deliberately NOT built on fq_series.py.
A series is stored as a flat list of Fractions for the exponents LO .. LO + len - 1 (LO fixed, negative allowed),
so multiplication by u^-1 and the Laurent output of the operator D = -(1-u) d/du need no special casing.
"""
from fractions import Fraction as Fr

LO = -12          # lowest exponent carried


class LS(object):
    __slots__ = ("c", "n")

    def __init__(self, c):
        self.c = list(c)
        self.n = len(c)

    # ---------- constructors -------------------------------------------------
    @staticmethod
    def zero(n):
        return LS([Fr(0)] * n)

    @staticmethod
    def mono(e, coef, n):
        """coef * u^e"""
        a = [Fr(0)] * n
        a[e - LO] = Fr(coef)
        return LS(a)

    @staticmethod
    def from_pairs(pairs, n):
        a = [Fr(0)] * n
        for e, v in pairs:
            a[e - LO] += Fr(v)
        return LS(a)

    @staticmethod
    def from_u0(coeffs, n):
        """power series given by its coefficients from u^0 upwards"""
        return LS.from_pairs(list(enumerate(coeffs)), n)

    # ---------- access -------------------------------------------------------
    def coef(self, e):
        i = e - LO
        return self.c[i] if 0 <= i < self.n else Fr(0)

    def val(self):
        for i, x in enumerate(self.c):
            if x:
                return i + LO
        return None

    def u0list(self, k):
        """coefficients of u^0 .. u^(k-1)"""
        return [self.coef(e) for e in range(k)]

    # ---------- arithmetic ---------------------------------------------------
    def add(self, o):
        return LS([a + b for a, b in zip(self.c, o.c)])

    def sub(self, o):
        return LS([a - b for a, b in zip(self.c, o.c)])

    def neg(self):
        return LS([-a for a in self.c])

    def scale(self, k):
        k = Fr(k)
        return LS([a * k for a in self.c])

    def shift(self, d):
        """multiply by u^d"""
        a = [Fr(0)] * self.n
        for i, x in enumerate(self.c):
            if x and 0 <= i + d < self.n:
                a[i + d] = x
        return LS(a)

    def mul(self, o):
        n = self.n
        a = [Fr(0)] * n
        sc = self.c
        oc = o.c
        for i, x in enumerate(sc):
            if not x:
                continue
            # exponent of sc[i] is i + LO, of oc[j] is j + LO; the product sits at index i + j + LO
            jlo = max(0, -LO - i)             # i + j + LO >= 0
            jhi = min(o.n, n - i - LO)        # i + j + LO <= n - 1
            for j in range(jlo, jhi):
                y = oc[j]
                if y:
                    a[i + j + LO] += x * y
        return LS(a)

    def deriv(self):
        """d/du"""
        a = [Fr(0)] * self.n
        for i, x in enumerate(self.c):
            e = i + LO
            if x and 0 <= i - 1 < self.n:
                a[i - 1] = x * e
        return LS(a)

    def inv(self):
        """inverse of a series of valuation v; result has valuation -v"""
        v = self.val()
        assert v is not None, "cannot invert zero"
        f = self.shift(-v)                  # valuation 0
        a0 = f.coef(0)
        out = [Fr(0)] * self.n
        top = self.n + LO                   # exponents 0 .. top-1 of the valuation-0 inverse
        g = [Fr(0)] * max(top, 1)
        g[0] = 1 / a0
        for m in range(1, max(top, 1)):
            s = Fr(0)
            for j in range(1, m + 1):
                fj = f.coef(j)
                if fj:
                    s += fj * g[m - j]
            g[m] = -s / a0
        for m in range(max(top, 1)):
            if 0 <= m - LO < self.n:
                out[m - LO] = g[m]
        return LS(out).shift(-v)

    def powi(self, k):
        assert k >= 0
        r = LS.mono(0, 1, self.n)
        b = self
        while k:
            if k & 1:
                r = r.mul(b)
            b = b.mul(b)
            k >>= 1
        return r


def compose_u0(f_u0, g, n):
    """f(g) for f given by coefficients from u^0 (list) and g a series of valuation >= 1"""
    out = LS.zero(n)
    pw = LS.mono(0, 1, n)
    for a in f_u0:
        if a:
            out = out.add(pw.scale(a))
        pw = pw.mul(g)
    return out


def reverse_u1(f, k):
    """compositional inverse of f = c1 u + c2 u^2 + ... (c1 != 0), returned as a list of k coefficients
    of t^1 .. t^k (index 0 holds [t^1]).  Lagrange-free: solve f(w(t)) = t order by order."""
    c = [f.coef(e) for e in range(0, k + 2)]
    assert c[0] == 0 and c[1] != 0, "reverse_u1 needs valuation exactly 1"
    w = [Fr(0)] * (k + 1)                # w[m] = [t^m] w
    w[1] = 1 / c[1]
    for m in range(2, k + 1):
        # [t^m] f(w) = c1 w[m] + [t^m] sum_{j>=2} c_j w^j ; the latter uses only w[1..m-1]
        wt = [Fr(0)] * (m + 1)
        for i in range(1, m):
            wt[i] = w[i]
        acc = Fr(0)
        pw = [Fr(0)] * (m + 1)
        pw[0] = Fr(1)
        for j in range(1, m + 1):
            nw = [Fr(0)] * (m + 1)
            for i in range(m + 1):
                if pw[i]:
                    for l in range(1, m + 1 - i):
                        if wt[l]:
                            nw[i + l] += pw[i] * wt[l]
            pw = nw
            if j >= 2 and j < len(c) and c[j]:
                acc += c[j] * pw[m]
        w[m] = -acc / c[1]
    return w[1:]
