# g4_field.py -- exact arithmetic for the G4 gap checks (stdlib only).
# Qz: the cyclotomic field Q(zeta), zeta = exp(2 pi i/5), elements a0 + a1 z + a2 z^2 + a3 z^3 (Fractions).
# Power series helpers on lists of field elements (index = power).
from fractions import Fraction as Fr


class Qz(object):
    __slots__ = ('c',)

    def __init__(self, c):
        self.c = tuple(Fr(x) for x in c)

    @staticmethod
    def of(q):
        return Qz((q, 0, 0, 0))

    @staticmethod
    def zeta(k):
        k %= 5
        if k == 4:
            return Qz((-1, -1, -1, -1))
        v = [0, 0, 0, 0]; v[k] = 1
        return Qz(v)

    def __add__(self, o):
        o = _co(o); return Qz(tuple(a + b for a, b in zip(self.c, o.c)))
    __radd__ = __add__

    def __neg__(self):
        return Qz(tuple(-a for a in self.c))

    def __sub__(self, o):
        return self + (-_co(o))

    def __rsub__(self, o):
        return _co(o) - self

    def __mul__(self, o):
        o = _co(o); r = [Fr(0)] * 7
        for i, a in enumerate(self.c):
            if a == 0: continue
            for j, b in enumerate(o.c):
                if b: r[i + j] += a * b
        for k in (6, 5, 4):                      # z^4 = -1 - z - z^2 - z^3
            t = r[k]
            if t:
                r[k] = Fr(0)
                for m in range(4): r[k - 4 + m] -= t
        return Qz(r[:4])
    __rmul__ = __mul__

    def galois(self, k):
        out = Qz.of(0)
        for i, a in enumerate(self.c):
            if a: out = out + Qz.zeta(i * k) * a
        return out

    def inv(self):
        num = self.galois(2) * self.galois(3) * self.galois(4)
        n = (self * num).c
        assert n[1] == n[2] == n[3] == 0 and n[0] != 0
        return num * (Fr(1) / n[0])

    def __truediv__(self, o):
        return self * _co(o).inv()

    def __rtruediv__(self, o):
        return _co(o) * self.inv()

    def iszero(self):
        return all(a == 0 for a in self.c)

    def __eq__(self, o):
        return (self - _co(o)).iszero()

    def __hash__(self):
        return hash(self.c)

    def rational(self):
        assert self.c[1] == self.c[2] == self.c[3] == 0
        return self.c[0]

    def __repr__(self):
        return 'Qz(%s)' % ', '.join(str(a) for a in self.c)


def _co(o):
    return o if isinstance(o, Qz) else Qz.of(o)


ZERO = Qz.of(0); ONE = Qz.of(1)


def smul(a, b, N):
    """product of series a, b truncated at order N (inclusive)"""
    r = [ZERO] * (N + 1)
    for i, x in enumerate(a[:N + 1]):
        if x.iszero(): continue
        for j, y in enumerate(b[:N + 1 - i]):
            if not y.iszero(): r[i + j] = r[i + j] + x * y
    return r


def sinv(a, N):
    """1/a for a series with a[0] != 0"""
    r = [ZERO] * (N + 1); r[0] = ONE / a[0]
    for n in range(1, N + 1):
        s = ZERO
        for k in range(1, min(n, len(a) - 1) + 1):
            s = s + a[k] * r[n - k]
        r[n] = -(s * r[0])
    return r


def ssqrt1(a, N):
    """sqrt of a series with a[0] == 1"""
    assert a[0] == ONE
    r = [ZERO] * (N + 1); r[0] = ONE
    for n in range(1, N + 1):
        s = a[n] if n < len(a) else ZERO
        for k in range(1, n):
            s = s - r[k] * r[n - k]
        r[n] = s * Fr(1, 2)
    return r


def srevert(a, N):
    """compositional inverse of a series a = a1 t + a2 t^2 + ... (a[0] = 0, a[1] != 0)"""
    assert a[0].iszero()
    b = [ZERO] * (N + 1); b[1] = ONE / a[1]
    for n in range(2, N + 1):
        # coefficient of s^n in a(b(s)) must vanish
        pw = [ZERO] * (N + 1); pw[0] = ONE
        tot = ZERO
        for k in range(1, n + 1):
            pw = smul(pw, b, n)
            if k < len(a): tot = tot + a[k] * pw[n]
        b[n] = -(tot / a[1])
    return b


def scompose(f, g, N):
    """f(g(s)) with g[0] = 0"""
    r = [ZERO] * (N + 1); pw = [ZERO] * (N + 1); pw[0] = ONE
    for k, fk in enumerate(f[:N + 1]):
        if k > 0: pw = smul(pw, g, N)
        if not fk.iszero():
            r = [x + fk * y for x, y in zip(r, pw)]
    return r


def dfact(n):
    r = 1
    while n > 1:
        r *= n; n -= 2
    return r


def bernoulli(n):
    B = [Fr(0)] * (n + 1); B[0] = Fr(1)
    from math import comb
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B
