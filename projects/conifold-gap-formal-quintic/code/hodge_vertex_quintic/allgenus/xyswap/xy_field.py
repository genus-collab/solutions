# -*- coding: utf-8 -*-
"""xy_field.py -- coefficient fields and series helpers for the direct topological recursion (stdlib only, exact).

Two interchangeable coefficient fields containing a primitive fifth root of unity zeta:
  QZ : Q(zeta), elements = 4-tuples of Fractions in the basis 1, zeta, zeta^2, zeta^3   (exact, slow)
  FP : Z/p, p = 2^89-1 or 2^61-1 (both are 1 mod 5)                                    (fast; evidence mod p only)
Both expose: zero(), one(), frac(Fraction), zeta(j), inv1mz(j) = 1/(1 - zeta^j) (j != 0 mod 5),
             inv_a_minus_zeta(a, j) = 1/(a - zeta^j) for a rational a with a^5 != 1, and el.rational().
Laurent series in one variable t are dicts {exponent: coefficient}; rational power series are lists of Fractions.
Nothing here writes a file.
"""
from fractions import Fraction as Fr

P89 = 2 ** 89 - 1
P61 = 2 ** 61 - 1


class QZ:
    __slots__ = ("c",)
    name = "Q(zeta_5)"

    def __init__(self, c): self.c = c

    @staticmethod
    def zero(): return QZ((Fr(0), Fr(0), Fr(0), Fr(0)))

    @staticmethod
    def one(): return QZ((Fr(1), Fr(0), Fr(0), Fr(0)))

    @staticmethod
    def frac(a): return QZ((Fr(a), Fr(0), Fr(0), Fr(0)))

    @staticmethod
    def zeta(j):
        j %= 5
        if j == 4: return QZ((Fr(-1), Fr(-1), Fr(-1), Fr(-1)))
        c = [Fr(0)] * 4; c[j] = Fr(1)
        return QZ(tuple(c))

    @staticmethod
    def inv1mz(j):
        # (1 - zeta) Sum_{k=1..4} k zeta^k = -5
        out = QZ.zero()
        for k in range(1, 5):
            out = out + QZ.zeta(j * k).scale(Fr(-k, 5))
        return out

    @staticmethod
    def inv_a_minus_zeta(a, j):
        a = Fr(a); out = QZ.zero()
        for k in range(5):
            out = out + QZ.zeta(j * k).scale(a ** (4 - k))
        return out.scale(1 / (a ** 5 - 1))

    def __add__(self, o):
        a, b = self.c, o.c
        return QZ((a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]))

    def __sub__(self, o):
        a, b = self.c, o.c
        return QZ((a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]))

    def __neg__(self):
        a = self.c
        return QZ((-a[0], -a[1], -a[2], -a[3]))

    def __mul__(self, o):
        a, b = self.c, o.c
        d = [Fr(0)] * 7
        for i in range(4):
            ai = a[i]
            if ai:
                for k in range(4):
                    if b[k]: d[i + k] += ai * b[k]
        d4 = d[4]
        return QZ((d[0] - d4 + d[5], d[1] - d4 + d[6], d[2] - d4, d[3] - d4))

    def scale(self, s):
        a = self.c
        return QZ((a[0] * s, a[1] * s, a[2] * s, a[3] * s))

    def is_zero(self):
        a = self.c
        return not (a[0] or a[1] or a[2] or a[3])

    def rational(self):
        a = self.c
        return a[0] if not (a[1] or a[2] or a[3]) else None

    def __repr__(self): return "QZ(%s)" % (", ".join(str(x) for x in self.c),)


def make_fp(p):
    z = 1; g = 2
    while z == 1:
        z = pow(g, (p - 1) // 5, p); g += 1
    zp = [pow(z, j, p) for j in range(5)]

    class FP:
        __slots__ = ("v",)
        name = "Z/(%d)" % p
        P = p

        def __init__(self, v): self.v = v % p

        @staticmethod
        def zero(): return FP(0)

        @staticmethod
        def one(): return FP(1)

        @staticmethod
        def frac(a):
            a = Fr(a); return FP(a.numerator * pow(a.denominator, -1, p))

        @staticmethod
        def zeta(j): return FP(zp[j % 5])

        @staticmethod
        def inv1mz(j): return FP(pow(1 - zp[j % 5], -1, p))

        @staticmethod
        def inv_a_minus_zeta(a, j): return FP(pow(FP.frac(a).v - zp[j % 5], -1, p))

        def __add__(self, o): return FP(self.v + o.v)

        def __sub__(self, o): return FP(self.v - o.v)

        def __neg__(self): return FP(-self.v)

        def __mul__(self, o): return FP(self.v * o.v)

        def scale(self, s):
            s = Fr(s); return FP(self.v * s.numerator * pow(s.denominator, -1, p))

        def is_zero(self): return self.v == 0

        def rational(self, bound=None):
            # rational reconstruction (Wang); None if no small fraction exists
            bound = bound or int((p // 2) ** 0.5)
            r0, r1, s0, s1 = p, self.v, 0, 1
            while r1 > bound:
                q = r0 // r1
                r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
            if s1 == 0 or abs(s1) > bound: return None
            return Fr(r1, s1)

        def __repr__(self): return "FP(%d)" % self.v
    return FP


# ---------------- rational power series (lists of Fractions, index = power) ----------------
def ps_mul(a, b, n):
    out = [Fr(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai:
            for k, bk in enumerate(b[: n + 1 - i]):
                if bk: out[i + k] += ai * bk
    return out


def ps_inv(a, n):
    out = [Fr(0)] * (n + 1); out[0] = 1 / a[0]
    for k in range(1, n + 1):
        s = Fr(0)
        for i in range(1, k + 1):
            if i < len(a) and a[i]: s += a[i] * out[k - i]
        out[k] = -s / a[0]
    return out


def ps_pow(a, m, n):
    """a^m for an integer m (a[0] != 0 if m < 0)"""
    if m < 0: a = ps_inv(a, n); m = -m
    out = [Fr(1)] + [Fr(0)] * n
    for _ in range(m): out = ps_mul(out, a, n)
    return out


def ps_compose(f, g, n):
    """f(g(t)) mod t^(n+1), g[0] = 0"""
    out = [Fr(0)] * (n + 1); pw = [Fr(1)] + [Fr(0)] * n
    for k, fk in enumerate(f[: n + 1]):
        if fk:
            for i in range(n + 1): out[i] += fk * pw[i]
        pw = ps_mul(pw, g, n)
    return out


# ---------------- Laurent series with field coefficients: dict {e: coeff} ----------------
def l_mul(a, b, top):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = ea + eb
            if e <= top:
                p = ca * cb
                if e in out: out[e] = out[e] + p
                else: out[e] = p
    return out


def l_add_into(acc, a, coeff=None):
    for e, c in a.items():
        if coeff is not None: c = c * coeff
        if e in acc: acc[e] = acc[e] + c
        else: acc[e] = c


def l_scale(a, coeff):
    return {e: c * coeff for e, c in a.items()}
