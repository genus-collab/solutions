# -*- coding: utf-8 -*-
"""rc_core.py -- critic of the Stage-5 referees (referees/critic): exact arithmetic for the (CYC) re-check.

New code (Python stdlib only); shares nothing with the critic's cr5_cyc*, pb's pb_cyc, the referees' cyc_ref_* / cy_* or the
project's fq_* modules.
  LS : truncated Laurent series in u over Q with ABSOLUTE precision (coefficients known for exponents < prec).
  RL : elements of R = K[l], K = Q((u)), l an indeterminate with D l = 1 - 1/u (l = log u); a list of LS (coefficient of l^n).
  D  : D = q d/dq = -(1 - u) d/du on LS and on RL;  Dx : the twisted derivation g -> x (1 - 1/u) g + D g, i.e.
       u^-x D (u^x g), used to test operator identities on u^x for several non-integer rational x (an operator of order <= 4 with
       coefficients in K is determined by its values on u^x at five distinct x).
"""
from fractions import Fraction as Fr

BIG = 10 ** 9
INV_CAP = 160          # relative precision used when inverting an exact series


class LS:
    __slots__ = ('c', 'prec')

    def __init__(self, c, prec):
        self.prec = prec
        self.c = {k: Fr(v) for k, v in c.items() if v != 0 and k < prec}

    @staticmethod
    def mono(k, a=1, prec=BIG):
        return LS({k: a}, prec)

    def val(self):
        return min(self.c) if self.c else self.prec

    def coef(self, k):
        if k >= self.prec:
            raise ValueError('coefficient u^%d beyond precision %d' % (k, self.prec))
        return self.c.get(k, Fr(0))

    def __add__(self, o):
        if not isinstance(o, LS):
            o = LS({0: o}, BIG)
        p = min(self.prec, o.prec)
        d = dict(self.c)
        for k, v in o.c.items():
            d[k] = d.get(k, 0) + v
        return LS(d, p)

    __radd__ = __add__

    def __neg__(self):
        return LS({k: -v for k, v in self.c.items()}, self.prec)

    def __sub__(self, o):
        return self + (-o if isinstance(o, LS) else LS({0: -Fr(o)}, BIG))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, LS):
            return LS({k: v * o for k, v in self.c.items()}, self.prec)
        p = min(self.prec + o.val(), o.prec + self.val())
        d = {}
        for k1, v1 in self.c.items():
            for k2, v2 in o.c.items():
                k = k1 + k2
                if k < p:
                    d[k] = d.get(k, 0) + v1 * v2
        return LS(d, p)

    __rmul__ = __mul__

    def inv(self):
        """1/self, coefficient recursion b_n = -(1/a_0) sum_(k=1..n) a_k b_(n-k); an exact (infinite-precision) series is
        inverted to relative precision INV_CAP."""
        if not self.c:
            raise ZeroDivisionError('series zero to its precision')
        v = self.val()
        rel = min(self.prec - v, INV_CAP)
        a = [self.c.get(v + k, Fr(0)) for k in range(rel)]
        b = [Fr(1) / a[0]]
        for n in range(1, rel):
            b.append(-b[0] * sum(a[k] * b[n - k] for k in range(1, n + 1) if a[k]))
        return LS({k - v: b[k] for k in range(rel)}, rel - v)

    def __truediv__(self, o):
        if not isinstance(o, LS):
            return self * Fr(1, 1) * (Fr(1) / Fr(o))
        return self * o.inv()

    def der(self):
        return LS({k - 1: k * v for k, v in self.c.items() if k != 0}, self.prec - 1)

    def D(self):
        return LS({0: -1, 1: 1}, BIG) * self.der()     # -(1-u) d/du

    def Dx(self, x):
        return LS({0: x, -1: -x}, BIG) * self + self.D()

    def zero_below(self, n):
        return all(self.coef(k) == 0 for k in range(min(self.val(), n), n))


ONE = LS({0: 1}, BIG)
U = LS({1: 1}, BIG)
DL = LS({0: 1, -1: -1}, BIG)        # D l = 1 - 1/u


class RL:
    """a_0 + a_1 l + ... (a_n in LS)."""
    __slots__ = ('a',)

    def __init__(self, a):
        self.a = list(a) if a else [LS({}, BIG)]

    @staticmethod
    def of(s):
        return RL([s])

    def __add__(self, o):
        n = max(len(self.a), len(o.a))
        z = LS({}, BIG)
        return RL([(self.a[i] if i < len(self.a) else z) + (o.a[i] if i < len(o.a) else z) for i in range(n)])

    def __neg__(self):
        return RL([-x for x in self.a])

    def __sub__(self, o):
        return self + (-o)

    def scale(self, s):           # multiplication by an element of K (LS) or a rational
        return RL([x * s for x in self.a])

    def __mul__(self, o):
        z = LS({}, BIG)
        out = [z] * (len(self.a) + len(o.a) - 1)
        for i, x in enumerate(self.a):
            for j, y in enumerate(o.a):
                out[i + j] = out[i + j] + x * y
        return RL(out)

    def D(self):
        n = len(self.a)
        out = [self.a[i].D() for i in range(n)]
        for i in range(1, n):
            out[i - 1] = out[i - 1] + self.a[i] * DL * i
        return RL(out)

    def logfree(self, upto):
        return all(x.zero_below(upto) for x in self.a[1:])


def apply_op(coeffs, f):
    """sum_k coeffs[k] D^k f for f in RL (coeffs: list of LS)."""
    out, g = None, f
    for k, ck in enumerate(coeffs):
        t = g.scale(ck)
        out = t if out is None else out + t
        g = g.D()
    return out


def twisted(chain, x, start=None):
    """u^-x (chain) u^x . start, where chain is a list applied right to left of ('m', LS) multiplications and ('D',) steps."""
    g = start if start is not None else ONE
    for step in reversed(chain):
        g = g * step[1] if step[0] == 'm' else g.Dx(x)
    return g
