# -*- coding: utf-8 -*-
"""rf_ser.py -- referee code (allgenus/tr/referee): fields and truncated Laurent series.  Written independently of tr_*.py.

Fields (uniform interface F.zero, F.one, F.of(x), F.iszero(a), F.name; elements support + - * / and unary -):
  QField()        exact rationals (fractions.Fraction)
  Q5Field()       exact Q(zeta), zeta^5 = 1 (basis 1, zeta, zeta^2, zeta^3), element class Q5
  FpField(p, r)   integers mod p with a fixed primitive r-th root of unity F.w (needs p = 1 mod r); F.rat(a) = rational reconstruction
Laurent series: class Ls(v, c, hi): Sum_k c[k] q^(v+k), coefficients known exactly for degrees v..hi (len(c) = hi - v + 1).
"""
from fractions import Fraction as Fr
import math

P61 = (1 << 61) - 1
P89 = (1 << 89) - 1


class QField(object):
    name = "Q"
    def __init__(self): self.zero = Fr(0); self.one = Fr(1)
    def of(self, x): return Fr(x)
    def iszero(self, a): return a == 0


class Q5(object):
    __slots__ = ('a',)
    def __init__(self, a): self.a = a                          # tuple of 4 Fractions
    def _c(self, o):
        if isinstance(o, Q5): return o
        return Q5((Fr(o), Fr(0), Fr(0), Fr(0)))
    def __add__(self, o): o = self._c(o); return Q5(tuple(x + y for x, y in zip(self.a, o.a)))
    __radd__ = __add__
    def __sub__(self, o): o = self._c(o); return Q5(tuple(x - y for x, y in zip(self.a, o.a)))
    def __rsub__(self, o): return self._c(o) - self
    def __neg__(self): return Q5(tuple(-x for x in self.a))
    def __mul__(self, o):
        if not isinstance(o, Q5): o = Fr(o); return Q5(tuple(x * o for x in self.a))
        r = [Fr(0)] * 7
        for i, x in enumerate(self.a):
            if x:
                for j, y in enumerate(o.a):
                    if y: r[i + j] += x * y
        r[0] += r[5]; r[1] += r[6]                              # zeta^5 = 1
        return Q5((r[0] - r[4], r[1] - r[4], r[2] - r[4], r[3] - r[4]))   # zeta^4 = -1 - zeta - zeta^2 - zeta^3
    __rmul__ = __mul__
    def gal(self, m):
        v = [Fr(0)] * 5
        for k, x in enumerate(self.a): v[(m * k) % 5] += x
        return Q5((v[0] - v[4], v[1] - v[4], v[2] - v[4], v[3] - v[4]))
    def inv(self):
        num = self.gal(2) * self.gal(3) * self.gal(4); n = (self * num).a
        assert n[1] == n[2] == n[3] == 0 and n[0] != 0
        return num * (1 / n[0])
    def __truediv__(self, o):
        if not isinstance(o, Q5): return self * (1 / Fr(o))
        return self * o.inv()
    def __rtruediv__(self, o): return self.inv() * o
    def __eq__(self, o): return self.a == self._c(o).a
    def __hash__(self): return hash(self.a)
    def iszero(self): return not any(self.a)
    def __repr__(self): return "Q5%s" % (tuple(str(x) for x in self.a),)


class Q5Field(object):
    name = "Q(zeta5)"
    def __init__(self):
        self.zero = Q5((Fr(0),) * 4); self.one = Q5((Fr(1), Fr(0), Fr(0), Fr(0)))
        self.w = Q5((Fr(0), Fr(1), Fr(0), Fr(0)))
    def of(self, x): return x if isinstance(x, Q5) else Q5((Fr(x), Fr(0), Fr(0), Fr(0)))
    def iszero(self, a): return a.iszero()
    def pw(self, k):
        r = self.one
        for _ in range(k % 5): r = r * self.w
        return r


def FpField(p, r=5):
    class E(object):
        __slots__ = ('v',)
        def __init__(self, v): self.v = v % p
        def __add__(self, o): return E(self.v + (o.v if type(o) is E else F.of(o).v))
        __radd__ = __add__
        def __sub__(self, o): return E(self.v - (o.v if type(o) is E else F.of(o).v))
        def __rsub__(self, o): return E(F.of(o).v - self.v)
        def __neg__(self): return E(-self.v)
        def __mul__(self, o): return E(self.v * (o.v if type(o) is E else F.of(o).v))
        __rmul__ = __mul__
        def __truediv__(self, o):
            ov = o.v if type(o) is E else F.of(o).v
            assert ov % p, "division by 0 mod p"
            return E(self.v * pow(ov, p - 2, p))
        def __rtruediv__(self, o): return F.of(o) / self
        def __eq__(self, o): return self.v == (o.v if type(o) is E else F.of(o).v)
        def __hash__(self): return hash(self.v)
        def iszero(self): return self.v == 0
        def __repr__(self): return "Fp(%d)" % self.v

    class Fld(object):
        pass
    F = Fld(); F.name = "F_p, p = %d" % p; F.p = p; F.E = E; F.zero = E(0); F.one = E(1)
    def of(x):
        if type(x) is E: return x
        if isinstance(x, Fr): return E(x.numerator * pow(x.denominator % p, p - 2, p))
        return E(x)
    F.of = of; F.iszero = lambda a: a.v == 0
    assert (p - 1) % r == 0
    g = 2
    while True:
        w = pow(g, (p - 1) // r, p)
        if all(pow(w, r // q, p) != 1 for q in range(2, r + 1) if r % q == 0 and all(q % s for s in range(2, q))): break
        g += 1
    F.w = E(w); F.wv = w
    F.pw = lambda k: E(pow(w, k % r, p))
    def rat(a):
        a = a.v if type(a) is E else a % p
        if a == 0: return Fr(0)
        B = math.isqrt(p // 2); r0, r1, s0, s1 = p, a, 0, 1
        while r1 > B:
            q = r0 // r1; r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
        if s1 == 0 or abs(s1) > B: return None
        f = Fr(r1, s1)
        return f if (f.numerator - a * f.denominator) % p == 0 else None
    F.rat = rat
    return F


class Ls(object):
    """truncated Laurent series: degrees v..hi known"""
    __slots__ = ('v', 'c', 'hi')
    def __init__(self, v, c, hi):
        self.v = v; self.hi = hi; n = hi - v + 1
        self.c = list(c[:n]) if n > 0 else []

    def coef(self, d, F):
        assert d <= self.hi, ("coefficient beyond precision", d, self.hi)
        k = d - self.v
        return self.c[k] if 0 <= k < len(self.c) else F.zero

    def mul(self, o, F):
        v = self.v + o.v; hi = min(self.hi + o.v, o.hi + self.v); n = hi - v + 1
        out = [F.zero] * max(n, 0)
        for i, x in enumerate(self.c[:n]):
            if F.iszero(x): continue
            for j, y in enumerate(o.c[:n - i]):
                out[i + j] = out[i + j] + x * y
        return Ls(v, out, hi)

    def add(self, o, F, s=None):
        v = min(self.v, o.v); hi = min(self.hi, o.hi); out = []
        for d in range(v, hi + 1):
            a = self.coef(d, F) if d >= self.v else F.zero
            b = o.coef(d, F) if d >= o.v else F.zero
            out.append(a + (b if s is None else b * s))
        return Ls(v, out, hi)

    def scale(self, s): return Ls(self.v, [x * s for x in self.c], self.hi)

    def flip(self):
        """f(-q)"""
        return Ls(self.v, [x if (self.v + k) % 2 == 0 else -x for k, x in enumerate(self.c)], self.hi)

    def normal(self, F):
        k = 0
        while k < len(self.c) and F.iszero(self.c[k]): k += 1
        assert k < len(self.c), "series is zero to its precision"
        return Ls(self.v + k, self.c[k:], self.hi)

    def inv(self, F):
        s = self.normal(F); n = s.hi - s.v; a0i = F.one / s.c[0]; out = [a0i] + [F.zero] * n
        for m in range(1, n + 1):
            acc = F.zero
            for k in range(1, m + 1): acc = acc + s.c[k] * out[m - k]
            out[m] = -(acc * a0i)
        return Ls(-s.v, out, -s.v + n)

    def deriv(self, F):
        return Ls(self.v - 1, [x * F.of(self.v + k) for k, x in enumerate(self.c)], self.hi - 1)


def ps_power(a, alpha, N, F):
    """(a)^alpha for a power series a = 1 + a1 q + .. (list), alpha a Fraction: coefficients 0..N  (n b_n = Sum (alpha k - n + k) a_k b_{n-k})"""
    assert a[0] == F.one
    b = [F.one] + [F.zero] * N; al = Fr(alpha)
    for n in range(1, N + 1):
        acc = F.zero
        for k in range(1, min(n, len(a) - 1) + 1):
            acc = acc + a[k] * b[n - k] * F.of(al * k - (n - k))
        b[n] = acc / F.of(n)
    return b


def bernoulli(n):
    B = [Fr(1)]
    for m in range(1, n + 1):
        B.append(-sum(math.comb(m + 1, k) * B[k] for k in range(m)) / Fr(m + 1))
    return B[n]


if __name__ == "__main__":
    Q5F = Q5Field(); z = Q5F.w
    assert z * z * z * z * z == Q5F.one and (z * (z + 3)).inv() * (z * (z + 3)) == Q5F.one
    for p in (P61, P89):
        F = FpField(p); x = F.of(Fr(-7, 11)); assert F.rat(x) == Fr(-7, 11) and F.w * F.w * F.w * F.w * F.w == F.one
    F = QField(); s = Ls(-1, [F.one, F.of(2), F.of(3)], 1); t = s.mul(s.inv(F), F)
    assert t.v == 0 and t.c[0] == 1 and all(c == 0 for c in t.c[1:])
    assert bernoulli(4) == Fr(-1, 30) and bernoulli(12) == Fr(-691, 2730)
    print("rf_ser self-test ok")
