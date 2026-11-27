# -*- coding: utf-8 -*-
"""rf_ring.py -- REFEREE library for the Stage-4 gamma task (stdlib only; imports NO worker module and NO project module).

Fields: FQ (exact rationals) and FP(p) (prime field).  Truncated power series are Python lists of field elements, index = exponent.
Mod p, long products use Kronecker substitution (pack the coefficients into one integer, multiply, unpack).
Also: Bernoulli numbers (recurrence Sum_{k<=n} C(n+1,k) B_k = 0), Bernoulli polynomials, and exact arithmetic in Q(zeta_5)
(elements = 5 coefficients on 1, zeta, .., zeta^4 modulo 1 + zeta + .. + zeta^4; trace and Galois action by zeta -> zeta^a).
"""
from fractions import Fraction as Fr
from math import comb

P61 = (1 << 61) - 1
P89 = (1 << 89) - 1
P127 = (1 << 127) - 1


class FQ:
    name = "Q"; p = None; zero = Fr(0); one = Fr(1)
    def of(self, x): return Fr(x)
    def inv(self, x): return 1 / x
    def red(self, x): return x


class FP:
    def __init__(self, p):
        self.p = p; self.zero = 0; self.one = 1
        self.name = "p=2^%d-1" % p.bit_length() if p == (1 << p.bit_length()) - 1 else "p=%d" % p
    def of(self, x):
        x = Fr(x); return x.numerator * pow(x.denominator, -1, self.p) % self.p
    def inv(self, x): return pow(x, -1, self.p)
    def red(self, x): return x % self.p


# ------------------------------------------------------------------ products
def _nbytes(p, n, k=1):
    return (2 * p.bit_length() + n.bit_length() + k.bit_length() + 8) // 8 + 1


def _pack(a, n, B):
    return int.from_bytes(b"".join(x.to_bytes(B, "little") for x in a[:n]), "little")


def _unpack(A, n, B, p):
    m = max(n * B, (A.bit_length() + 7) // 8); bs = A.to_bytes(m, "little")
    return [int.from_bytes(bs[i * B:(i + 1) * B], "little") % p for i in range(n)]


def mul(F, a, b, n):
    """a*b mod t^n"""
    if F.p is not None and n > 20 and len(a) > 8 and len(b) > 8:
        B = _nbytes(F.p, n); return _unpack(_pack(a, n, B) * _pack(b, n, B), n, B, F.p)
    out = [F.zero] * n
    for i, x in enumerate(a[:n]):
        if x:
            for j, y in enumerate(b[:n - i]):
                out[i + j] += x * y
    return [F.red(x) for x in out]


def mulsum(F, pairs, n):
    """Sum of a*b over (a, b) in pairs, mod t^n (one unpack mod p)"""
    if F.p is not None and n > 20:
        B = _nbytes(F.p, n, len(pairs) + 1); acc = 0
        for a, b in pairs: acc += _pack(a, n, B) * _pack(b, n, B)
        return _unpack(acc & ((1 << (8 * B * n)) - 1), n, B, F.p)
    out = [F.zero] * n
    for a, b in pairs:
        for i, x in enumerate(a[:n]):
            if x:
                for j, y in enumerate(b[:n - i]):
                    out[i + j] += x * y
    return [F.red(x) for x in out]


def pad(a, n, F):
    a = list(a[:n]); return a + [F.zero] * (n - len(a))


def add(F, a, b, c=1):
    n = max(len(a), len(b)); a = pad(a, n, F); b = pad(b, n, F)
    return [F.red(x + c * y) for x, y in zip(a, b)]


def scal(F, a, c): return [F.red(c * x) for x in a]


def inv(F, a, n):
    """1/a mod t^n (a[0] invertible)"""
    a = pad(a, n, F); i0 = F.inv(a[0]); out = [F.zero] * n; out[0] = i0
    for k in range(1, n):
        out[k] = F.red(-i0 * sum(a[j] * out[k - j] for j in range(1, k + 1)))
    return out


def deriv(F, a): return [F.red(k * a[k]) for k in range(1, len(a))]


def integ(F, a):
    return [F.zero] + [F.red(x * F.inv(F.of(k + 1))) for k, x in enumerate(a)]


def log1(F, a, n):
    """log(a), a[0] = 1, mod t^n"""
    a = pad(a, n, F); assert a[0] == F.one
    return integ(F, mul(F, deriv(F, a), inv(F, a, n), n - 1))[:n]


def exp0(F, f, n):
    """exp(f), f[0] = 0, mod t^n (k e_k = Sum_j j f_j e_(k-j))"""
    f = pad(f, n, F); assert f[0] == F.zero
    jf = [F.red(j * f[j]) for j in range(n)]; out = [F.zero] * n; out[0] = F.one
    for k in range(1, n):
        out[k] = F.red(sum(jf[j] * out[k - j] for j in range(1, k + 1)) * F.inv(F.of(k)))
    return out


def powr(F, a, r, n):
    """a^r (a[0] = 1, r rational) mod t^n"""
    return exp0(F, scal(F, log1(F, a, n), F.of(r)), n)


def reversion(F, phi, n):
    """t(y) with y = t/phi(t) (phi[0] = 1), mod y^n, by Lagrange: [y^k] t = (1/k) [t^(k-1)] phi^k"""
    out = [F.zero] * n; pw = [F.one] + [F.zero] * (n - 1)
    for k in range(1, n):
        pw = mul(F, pw, phi, n); out[k] = F.red(pw[k - 1] * F.inv(F.of(k)))
    return out


def compose(F, a, t, n):
    """a(t(y)) mod y^n, t[0] = 0 (Horner; used only at small order)"""
    out = [F.zero] * n
    for c in reversed(pad(a, n, F)):
        out = mul(F, out, t, n); out[0] = F.red(out[0] + c)
    return out


# ------------------------------------------------------------------ Bernoulli
_BERN = [Fr(1)]


def bernoulli(n):
    while len(_BERN) <= n:
        m = len(_BERN); _BERN.append(-sum(comb(m + 1, k) * _BERN[k] for k in range(m)) / (m + 1))
    return _BERN[n]


def bernoulli_poly(n, x):
    return sum(comb(n, k) * bernoulli(k) * Fr(x) ** (n - k) for k in range(n + 1))


def stirling_b(k):
    """b_k = B_2k/(2k(2k-1)): log Gamma(y) = (y-1/2) log y - y + log(2 pi)/2 + Sum_k b_k y^(1-2k)"""
    return bernoulli(2 * k) / (2 * k * (2 * k - 1))


# ------------------------------------------------------------------ Q(zeta_5)
def cz_red(c):
    """5 coefficients (exponents mod 5) -> canonical 4 coefficients (zeta^4 = -1 - zeta - zeta^2 - zeta^3)"""
    return tuple(Fr(c[i]) - Fr(c[4]) for i in range(4))


def cz_full(a): return list(a) + [Fr(0)]


def cz_mul(a, b):
    out = [Fr(0)] * 5
    for i, x in enumerate(cz_full(a)):
        if x:
            for j, y in enumerate(cz_full(b)):
                out[(i + j) % 5] += x * y
    return cz_red(out)


def cz_gal(a, g):
    out = [Fr(0)] * 5
    for i, x in enumerate(cz_full(a)): out[(g * i) % 5] += x
    return cz_red(out)


def cz_trace(a):
    """Tr_{Q(zeta)/Q}: Tr(1) = 4, Tr(zeta^i) = -1"""
    return 4 * a[0] - a[1] - a[2] - a[3]


def cz_inv(a):
    num = cz_mul(cz_mul(cz_gal(a, 2), cz_gal(a, 3)), cz_gal(a, 4)); nrm = cz_mul(a, num)
    assert nrm[1] == nrm[2] == nrm[3] == 0 and nrm[0] != 0
    return tuple(x / nrm[0] for x in num)
