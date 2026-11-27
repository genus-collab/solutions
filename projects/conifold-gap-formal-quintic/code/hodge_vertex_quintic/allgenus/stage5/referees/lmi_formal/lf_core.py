"""lf_core: exact arithmetic for the lmi_formal referee (Python stdlib only).

Contents: Bernoulli numbers and polynomials; b_k = B_2k/(2k(2k-1)); generalised binomials;
truncated power series in one variable (lists of Fractions); the ring TS = Q[eps]/(eps^K);
the field Q(zeta_5) (class Cyc5, used only to compute h from its definition); the power sums
P_r(x) = Sum_(j=0..4) (zeta^j - x)^-r and p_r(x) = Sum_(j=1..4) (zeta^j - x)^-r at rational x,
from the generating functions 5 y^4/(1 - y^5) and -Phi5'(y)/Phi5(y), y = x + t (no zeta needed).
Written for this referee report; imports no project, worker, referee or critic module.
"""
from fractions import Fraction as Fr
from math import comb, factorial


def bernoulli_numbers(N):
    """B_0..B_N, convention B_1 = -1/2."""
    B = [Fr(0)] * (N + 1)
    B[0] = Fr(1)
    for n in range(1, N + 1):
        B[n] = -sum(comb(n + 1, k) * B[k] for k in range(n)) / (n + 1)
    return B


_BN = bernoulli_numbers(260)


def bern(n):
    return _BN[n]


def bern_poly(n, x):
    """Bernoulli polynomial B_n(x) = Sum_k binom(n, k) B_k x^(n-k) (B_1(x) = x - 1/2)."""
    x = Fr(x)
    return sum(comb(n, k) * _BN[k] * x ** (n - k) for k in range(n + 1))


def bk(k, shift=0):
    """b_k = B_2k/(2k(2k-1)); `shift` is added (used only by failing controls)."""
    return _BN[2 * k] / (2 * k * (2 * k - 1)) + shift


def gbinom(a, r):
    """Generalised binomial a(a-1)..(a-r+1)/r! for a in Q, r >= 0."""
    a = Fr(a)
    out = Fr(1)
    for i in range(r):
        out *= (a - i)
    return out / factorial(r)


def dfact(n):
    """Double factorial for odd n >= -1 (dfact(-1) = 1)."""
    out = 1
    while n > 1:
        out *= n
        n -= 2
    return out


# ---------------- truncated power series in one variable (lists of Fractions) ----------------

def ser_mul(a, b, N):
    out = [Fr(0)] * N
    for i, ai in enumerate(a[:N]):
        if ai:
            for j, bj in enumerate(b[:N - i]):
                if bj:
                    out[i + j] += ai * bj
    return out


def ser_inv(a, N):
    """1/a for a power series with a[0] != 0."""
    out = [Fr(0)] * N
    out[0] = Fr(1) / a[0]
    for n in range(1, N):
        s = sum(a[k] * out[n - k] for k in range(1, min(n, len(a) - 1) + 1))
        out[n] = -s / a[0]
    return out


def ser_exp(a, N):
    """exp(a) for a power series with a[0] = 0: n f_n = Sum_k k a_k f_(n-k)."""
    a = list(a) + [Fr(0)] * max(0, N - len(a))
    assert a[0] == 0
    f = [Fr(0)] * N
    f[0] = Fr(1)
    for n in range(1, N):
        f[n] = sum(k * a[k] * f[n - k] for k in range(1, n + 1)) / n
    return f


# ---------------- the ring TS = Q[eps]/(eps^K) ----------------

class TS:
    __slots__ = ("c",)
    K = 8

    def __init__(self, c=None):
        c = list(c) if c is not None else []
        c = [Fr(x) for x in c[:TS.K]]
        self.c = c + [Fr(0)] * (TS.K - len(c))

    @staticmethod
    def const(x):
        return TS([x])

    def __bool__(self):
        return any(self.c)

    def __add__(self, o):
        if not isinstance(o, TS):
            o = TS.const(o)
        return TS([x + y for x, y in zip(self.c, o.c)])

    __radd__ = __add__

    def __neg__(self):
        return TS([-x for x in self.c])

    def __sub__(self, o):
        return self + (-o if isinstance(o, TS) else TS.const(-Fr(o)))

    def __mul__(self, o):
        if not isinstance(o, TS):
            o = Fr(o)
            return TS([x * o for x in self.c])
        K = TS.K
        out = [Fr(0)] * K
        for i, x in enumerate(self.c):
            if x:
                for j in range(K - i):
                    y = o.c[j]
                    if y:
                        out[i + j] += x * y
        return TS(out)

    __rmul__ = __mul__

    def inv(self):
        return TS(ser_inv(self.c, TS.K))

    def __eq__(self, o):
        if not isinstance(o, TS):
            o = TS.const(o)
        return self.c == o.c

    def __repr__(self):
        return "TS(%s)" % ", ".join(str(x) for x in self.c)


def eps_ts():
    return TS([0, 1])


# ---------------- Q(zeta_5), zeta^4 = -1 - zeta - zeta^2 - zeta^3 ----------------

class Cyc5:
    __slots__ = ("c",)

    def __init__(self, c):
        c = [Fr(x) for x in c] + [Fr(0)] * (5 - len(c))
        c4 = c[4]
        self.c = tuple(c[i] - c4 for i in range(4))   # subtract c4 (1 + z + .. + z^4) = 0

    def full(self):
        return list(self.c) + [Fr(0)]

    def __add__(self, o):
        return Cyc5([x + y for x, y in zip(self.full(), o.full())])

    def __mul__(self, o):
        if not isinstance(o, Cyc5):
            return Cyc5([x * Fr(o) for x in self.full()])
        out = [Fr(0)] * 5
        for i, x in enumerate(self.full()):
            for j, y in enumerate(o.full()):
                out[(i + j) % 5] += x * y
        return Cyc5(out)

    def conj(self, k):
        out = [Fr(0)] * 5
        for i, x in enumerate(self.full()):
            out[(i * k) % 5] += x
        return Cyc5(out)

    def inv(self):
        y = self.conj(2) * self.conj(3) * self.conj(4)
        n = (self * y).c
        assert n[1] == n[2] == n[3] == 0 and n[0] != 0
        return y * (Fr(1) / n[0])

    def rational(self):
        assert self.c[1] == self.c[2] == self.c[3] == 0, self.c
        return self.c[0]


def zeta_pow(j):
    c = [0] * 5
    c[j % 5] = 1
    return Cyc5(c)


# ---------------- power sums at rational x ----------------

def _poly_shift(coeffs, x, N):
    """coefficients in t of Sum_i coeffs[i] (x + t)^i, truncated at t^N."""
    out = [Fr(0)] * N
    for i, a in enumerate(coeffs):
        for r in range(min(i, N - 1) + 1):
            out[r] += a * comb(i, r) * Fr(x) ** (i - r)
    return out


def power_sums_all(x, R):
    """P[r] = Sum_(j=0..4) (zeta^j - x)^-r for r = 1..R (P[0] unused), x rational, x^5 != 1."""
    num = _poly_shift([0, 0, 0, 0, 5], x, R)
    den = _poly_shift([1, 0, 0, 0, 0, -1], x, R)
    s = ser_mul(num, ser_inv(den, R), R)
    return [None] + [s[r - 1] for r in range(1, R + 1)]


def power_sums_reg(x, R):
    """p[r] = Sum_(j=1..4) (zeta^j - x)^-r for r = 1..R, from -Phi5'(y)/Phi5(y), y = x + t."""
    num = _poly_shift([-1, -2, -3, -4], x, R)
    den = _poly_shift([1, 1, 1, 1, 1], x, R)
    s = ser_mul(num, ser_inv(den, R), R)
    return [None] + [s[r - 1] for r in range(1, R + 1)]
