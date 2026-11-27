# -*- coding: utf-8 -*-
"""gm_series.py -- truncated power series over F_p or Q (stdlib only), used by the Gamma-route scripts.

A series is a list of coefficients (index = power), all of the same length n (truncation O(x^n)).
Ring R: R.of(Fraction) -> element; R.add/mul; R.inv(a).  Rq = exact Fractions; Rp(p) = residues mod p.
Functions: smul, sinv, sder, sint, slog1 (log of a series with constant term 1), sexp0 (exp of a series with zero constant
term), spow (rational power of a series with constant term 1), reversion_powers (Lagrange inversion: powers of the inverse series).
Bernoulli numbers B_n (B_1 = -1/2) and Bernoulli polynomials B_n(x) are exact (Fractions).
"""
from fractions import Fraction as Fr
from math import comb

P89 = (1 << 89) - 1
P61 = (1 << 61) - 1


class Rq(object):
    name = "Q"
    def of(self, x): return Fr(x)
    def inv(self, a): return 1 / Fr(a)
    def zero(self): return Fr(0)
    def one(self): return Fr(1)
    def norm(self, a): return a


class Rp(object):
    def __init__(self, p): self.p = p; self.name = "p=2^%d-1" % (p.bit_length())
    def of(self, x):
        x = Fr(x); return x.numerator % self.p * pow(x.denominator % self.p, self.p - 2, self.p) % self.p
    def inv(self, a): return pow(a % self.p, self.p - 2, self.p)
    def zero(self): return 0
    def one(self): return 1
    def norm(self, a): return a % self.p


def bernoulli_list(N):
    B = [Fr(0)] * (N + 1); B[0] = Fr(1)
    for m in range(1, N + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


def bernoulli_poly(n, x, B=None):
    B = B or bernoulli_list(n)
    return sum(comb(n, k) * B[k] * Fr(x) ** (n - k) for k in range(n + 1))


def smul(R, a, b, n=None):
    n = n or len(a); out = [R.zero()] * n
    for i, x in enumerate(a[:n]):
        if x:
            for j in range(min(len(b), n - i)):
                out[i + j] += x * b[j]
    return [R.norm(c) for c in out]


def sinv(R, a, n=None):
    n = n or len(a); i0 = R.inv(a[0]); out = [R.zero()] * n; out[0] = i0
    for k in range(1, n):
        s = sum(a[j] * out[k - j] for j in range(1, min(k, len(a) - 1) + 1))
        out[k] = R.norm(-s * i0)
    return out


def sder(R, a): return [R.norm(k * a[k]) for k in range(1, len(a))] + [R.zero()]


def sint(R, a):
    """integral from 0 (same length; the top coefficient is dropped)"""
    return [R.zero()] + [R.norm(a[k] * R.inv(k + 1)) for k in range(len(a) - 1)]


def slog1(R, a):
    """log(a), a[0] = 1"""
    n = len(a); da = sder(R, a); q = smul(R, da, sinv(R, a, n), n)
    return sint(R, q)


def sexp0(R, a):
    """exp(a), a[0] = 0"""
    n = len(a); out = [R.zero()] * n; out[0] = R.one()
    for k in range(1, n):
        s = sum(j * a[j] * out[k - j] for j in range(1, k + 1))
        out[k] = R.norm(s * R.inv(k))
    return out


def spow(R, a, e):
    """a^e for a[0] = 1, e rational"""
    l = slog1(R, a); return sexp0(R, [R.norm(R.of(e) * x) for x in l])


def reversion_powers(R, w, n, rmax):
    """Y(t) = t w(t), w[0] = 1.  Returns T with T[r][k] = [Y^k] t(Y)^r for 0 <= r <= rmax, k < n (Lagrange:
    [Y^k] t^r = (r/k) [t^(k-r)] w^(-k))."""
    winv = sinv(R, w, n); T = [[R.zero()] * n for _ in range(rmax + 1)]; T[0][0] = R.one()
    Wk = [R.one()] + [R.zero()] * (n - 1)
    for k in range(1, n):
        Wk = smul(R, Wk, winv, n)
        for r in range(1, min(k, rmax) + 1):
            T[r][k] = R.norm(R.of(Fr(r, k)) * Wk[k - r])
    return T


def compose(R, f, T, n):
    """f(t(Y)) mod Y^n, f a series in t with f[0] allowed; T from reversion_powers"""
    out = [R.zero()] * n
    for r, fr in enumerate(f):
        if not fr or r >= len(T): continue
        Tr = T[r]
        for k in range(r, n):
            if Tr[k]: out[k] += fr * Tr[k]
    return [R.norm(c) for c in out]
