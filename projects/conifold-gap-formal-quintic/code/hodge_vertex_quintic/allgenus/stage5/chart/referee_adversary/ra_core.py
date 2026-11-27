# -*- coding: utf-8 -*-
"""ra_core.py -- stage5/chart/referee_adversary: fields, Bernoulli numbers, the rings R used by ra_sp.py.

Written by the adversarial referee; imports no project module.  stdlib only.
Fields:  FQ (exact rationals, Fraction)   and   FP(p) (integers mod a prime p = 1 mod 10, so zeta_5 lies in F_p).
Rings:   RSym(F, K) = F[eps]/(eps^K)  (elements: lists of length K)   -- the Laurent expansion at the colour point;
         RNum(F, v) = F with eps := v (elements: field elements)      -- evaluation at one point L = sigma* + v.
Both rings expose the same small interface, so ra_sp.py runs unchanged on either.
"""
from fractions import Fraction as Fr
from math import comb


def is_prime(n):
    """deterministic Miller-Rabin for n < 3.3e24 (bases = first 13 primes)"""
    if n < 2:
        return False
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
    for p in bases:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in bases:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_1mod10_below(start):
    p = start - ((start - 1) % 10)
    while not is_prime(p):
        p -= 10
    return p


# Referee primes: chosen below 10^18 and 10^21, different from the worker's (4611686018427388081, 1208925819614629174706261)
# and from the Stage-4 primes (2^61-1, 2^89-1, 2^107-1, 2^127-1).
PR1 = prime_1mod10_below(10 ** 18)
PR2 = prime_1mod10_below(10 ** 21)


class FQ:
    name = 'Q'
    zero = Fr(0)
    one = Fr(1)

    def c(self, x):
        return Fr(x)

    def red(self, x):
        return x

    def inv(self, x):
        return 1 / Fr(x)

    def iszero(self, x):
        return x == 0


class FP:
    def __init__(self, p):
        assert is_prime(p) and p % 10 == 1
        self.p = p; self.name = 'F_p, p = %d' % p; self.zero = 0; self.one = 1

    def c(self, x):
        x = Fr(x)
        return x.numerator * pow(x.denominator, -1, self.p) % self.p

    def red(self, x):
        return x % self.p

    def inv(self, x):
        return pow(x % self.p, -1, self.p)

    def iszero(self, x):
        return x % self.p == 0

    def zeta5(self):
        """a primitive 5th root of unity in F_p"""
        for a in range(2, 1000):
            z = pow(a, (self.p - 1) // 5, self.p)
            if z != 1:
                return z
        raise ValueError('no zeta')


def bernoulli(nmax):
    """B_0..B_nmax (B_1 = -1/2), exact"""
    B = [Fr(0)] * (nmax + 1); B[0] = Fr(1)
    for n in range(1, nmax + 1):
        B[n] = -sum(comb(n + 1, k) * B[k] for k in range(n)) / (n + 1)
    return B


def stirling_b(kmax):
    """b_k = B_2k/(2k(2k-1)), k = 1..kmax (index 0 unused)"""
    B = bernoulli(2 * kmax)
    return [Fr(0)] + [B[2 * k] / (2 * k * (2 * k - 1)) for k in range(1, kmax + 1)]


def gbinom(a, r):
    """generalised binomial coefficient binom(a, r), a rational, r >= 0"""
    out = Fr(1)
    for i in range(r):
        out = out * (Fr(a) - i) / (i + 1)
    return out


class RSym:
    """F[eps]/(eps^K)"""

    def __init__(self, F, K):
        self.F = F; self.K = K; self.name = '%s[eps]/(eps^%d)' % (F.name, K)

    def const(self, x):
        return [self.F.c(x)] + [self.F.zero] * (self.K - 1)

    def zero(self):
        return [self.F.zero] * self.K

    def epsp(self, r, x=1):
        """x * eps^r"""
        out = [self.F.zero] * self.K
        if r < self.K:
            out[r] = self.F.c(x)
        return out

    def add(self, a, b):
        red = self.F.red
        return [red(x + y) for x, y in zip(a, b)]

    def sub(self, a, b):
        red = self.F.red
        return [red(x - y) for x, y in zip(a, b)]

    def scal(self, x, a):
        x = self.F.c(x); red = self.F.red
        return [red(x * y) for y in a]

    def mul(self, a, b):
        K = self.K; out = [0] * K
        for i, ai in enumerate(a):
            if ai == 0:
                continue
            for j in range(K - i):
                bj = b[j]
                if bj != 0:
                    out[i + j] += ai * bj
        red = self.F.red
        return [red(x) for x in out]

    def inv(self, a):
        F = self.F; K = self.K; i0 = F.inv(a[0]); out = [F.zero] * K; out[0] = i0
        for n in range(1, K):
            s = 0
            for i in range(1, n + 1):
                s += a[i] * out[n - i]
            out[n] = F.red(-s * i0)
        return out

    def iszero(self, a):
        return all(self.F.iszero(x) for x in a)


class RNum:
    """the field F itself, with eps := v (evaluation at one point)"""

    def __init__(self, F, v):
        self.F = F; self.v = F.c(v); self.K = 10 ** 9; self.name = '%s at eps = %s' % (F.name, v)

    def const(self, x):
        return self.F.c(x)

    def zero(self):
        return self.F.zero

    def epsp(self, r, x=1):
        return self.F.red(self.F.c(x) * pow(self.v, r, self.F.p) if isinstance(self.F, FP) else Fr(x) * self.v ** r)

    def add(self, a, b):
        return self.F.red(a + b)

    def sub(self, a, b):
        return self.F.red(a - b)

    def scal(self, x, a):
        return self.F.red(self.F.c(x) * a)

    def mul(self, a, b):
        return self.F.red(a * b)

    def inv(self, a):
        return self.F.inv(a)

    def iszero(self, a):
        return self.F.iszero(a)
