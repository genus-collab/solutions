# rf4_field.py -- referee's own exact arithmetic for the G4 referee checks (stdlib only; independent of g4_field.py).
# K = Q(zeta), zeta = exp(2 pi i/5), stored as a dict-free 5-vector over Q in the redundant basis 1, z, .., z^4 and
# normalised by subtracting the coefficient of z^4 (1 + z + .. + z^4 = 0). Truncated power series = lists of K elements.
from fractions import Fraction as F


class K(object):
    __slots__ = ('a',)

    def __init__(self, a):
        a = [F(x) for x in a]
        while len(a) < 5:
            a.append(F(0))
        c4 = a[4]
        self.a = tuple(a[i] - c4 for i in range(4))   # canonical: coefficients of 1, z, z^2, z^3

    @staticmethod
    def q(r):
        return K([r])

    @staticmethod
    def z(k):
        v = [0] * 5
        v[k % 5] = 1
        return K(v)

    def _full(self):
        return list(self.a) + [F(0)]

    def __add__(s, o):
        o = kk(o); return K([x + y for x, y in zip(s.a, o.a)])
    __radd__ = __add__

    def __neg__(s):
        return K([-x for x in s.a])

    def __sub__(s, o):
        return s + (-kk(o))

    def __rsub__(s, o):
        return kk(o) - s

    def __mul__(s, o):
        o = kk(o)
        r = [F(0)] * 5
        for i, x in enumerate(s.a):
            if x:
                for j, y in enumerate(o.a):
                    if y:
                        r[(i + j) % 5] += x * y
        return K(r)
    __rmul__ = __mul__

    def conj(s, k):          # Galois automorphism z -> z^k
        r = [F(0)] * 5
        for i, x in enumerate(s.a):
            r[(i * k) % 5] += x
        return K(r)

    def inv(s):
        n = s.conj(2) * s.conj(3) * s.conj(4)
        d = s * n
        assert all(x == 0 for x in d.a[1:]) and d.a[0] != 0, 'norm not rational'
        return n * (F(1) / d.a[0])

    def __truediv__(s, o):
        return s * kk(o).inv()

    def __rtruediv__(s, o):
        return kk(o) * s.inv()

    def __eq__(s, o):
        return s.a == kk(o).a

    def __hash__(s):
        return hash(s.a)

    def is0(s):
        return all(x == 0 for x in s.a)

    def rat(s):
        assert all(x == 0 for x in s.a[1:]), 'not rational: %r' % (s,)
        return s.a[0]

    def __pow__(s, n):
        if n < 0:
            return s.inv() ** (-n)
        r = K.q(1); b = s
        while n:
            if n & 1: r = r * b
            b = b * b; n >>= 1
        return r

    def __repr__(s):
        return 'K(%s)' % ','.join(str(x) for x in s.a)


def kk(o):
    return o if isinstance(o, K) else K.q(o)


ZERO = K.q(0); ONE = K.q(1)


# ---- truncated power series in one variable: list c[0..N] ----
def smul(a, b, N):
    r = [ZERO] * (N + 1)
    for i in range(min(N, len(a) - 1) + 1):
        if a[i].is0(): continue
        for j in range(min(N - i, len(b) - 1) + 1):
            if not b[j].is0():
                r[i + j] = r[i + j] + a[i] * b[j]
    return r


def sinv(a, N):
    r = [ZERO] * (N + 1)
    i0 = a[0].inv(); r[0] = i0
    for n in range(1, N + 1):
        acc = ZERO
        for k in range(1, min(n, len(a) - 1) + 1):
            acc = acc + a[k] * r[n - k]
        r[n] = -(acc * i0)
    return r


def spow(a, e, N):
    """a^e for integer e (a[0] invertible if e < 0)"""
    if e < 0:
        a = sinv(a, N); e = -e
    r = [ONE] + [ZERO] * N
    for _ in range(e):
        r = smul(r, a, N)
    return r


def ssqrt(a, N):
    """square root of a series with a[0] = 1"""
    assert a[0] == ONE
    r = [ZERO] * (N + 1); r[0] = ONE
    for n in range(1, N + 1):
        acc = a[n] if n < len(a) else ZERO
        for k in range(1, n):
            acc = acc - r[k] * r[n - k]
        r[n] = acc * F(1, 2)
    return r


def sderiv(a):
    return [a[k] * k for k in range(1, len(a))] + [ZERO]


def bernoulli(n):
    from math import comb
    B = [F(1)] + [F(0)] * n
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


def dfact(n):
    r = 1
    while n > 1:
        r *= n; n -= 2
    return r
