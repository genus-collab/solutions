# -*- coding: utf-8 -*-
"""ch_core.py -- Stage 5 / chart: shared exact machinery (stdlib only; imports no project or agent module).

Fields / rings with one interface (zero, one, of, add, sub, neg, mul, smul, inv, iszero, eq):
  QZ  exact Q(zeta_5), basis (1, zeta, zeta^2, zeta^3), zeta^4 = -1 - zeta - zeta^2 - zeta^3;
  FP  F_p with p = 1 mod 5 and zeta a primitive 5th root of unity in F_p;
  RS  truncated power series ring Q[[w]]/(w^D) (only units are inverted) -- for the Fubini tests.
Graded multivariate polynomials P = {deg: {mono: coeff}}; formal stationary phase N[Phi, A; p] by Wick's theorem (CHART_PROOF s.1.2):
  N_n = Sum_m Sum_V < [Phi_3^V A_m / V!]_(2(n-m+V)) >_C ,  C = -H^-1.
Also: Bernoulli numbers, b_k = B_2k/(2k(2k-1)), and the QRR series h_(w_0) by Newton sums (no zeta arithmetic).
"""
from fractions import Fraction as Fr
from math import comb


class QZ:
    name = 'Q(zeta5)'
    zero = (Fr(0),) * 4
    one = (Fr(1), Fr(0), Fr(0), Fr(0))

    def of(self, x): return (Fr(x), Fr(0), Fr(0), Fr(0))
    def add(self, a, b): return (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])
    def sub(self, a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3])
    def neg(self, a): return (-a[0], -a[1], -a[2], -a[3])
    def smul(self, s, a): s = Fr(s); return (s * a[0], s * a[1], s * a[2], s * a[3])
    def iszero(self, a): return not (a[0] or a[1] or a[2] or a[3])
    def isunit(self, a): return not self.iszero(a)
    def eq(self, a, b): return self.iszero(self.sub(a, b))

    def mul(self, a, b):
        a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
        c4 = a1 * b3 + a2 * b2 + a3 * b1
        return (a0 * b0 - c4 + a2 * b3 + a3 * b2, a0 * b1 + a1 * b0 - c4 + a3 * b3,
                a0 * b2 + a1 * b1 + a2 * b0 - c4, a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0 - c4)

    def zeta(self, j):
        j %= 5
        if j == 4: return (Fr(-1),) * 4
        return tuple(Fr(1) if i == j else Fr(0) for i in range(4))

    def conj(self, a, k):                       # Galois action zeta -> zeta^k
        r = self.zero
        for i in range(4):
            if a[i]: r = self.add(r, self.smul(a[i], self.zeta(i * k)))
        return r

    def inv(self, a):
        b = self.mul(self.mul(self.conj(a, 2), self.conj(a, 3)), self.conj(a, 4))
        n = self.mul(a, b); assert n[1] == n[2] == n[3] == 0 and n[0] != 0, 'norm'
        return self.smul(1 / n[0], b)

    def rat(self, a):
        assert a[1] == a[2] == a[3] == 0, 'not rational: %s' % (a,)
        return a[0]


class FP:
    def __init__(self, p):
        self.p = p; self.name = 'F_p, p = %d' % p; self.zero = 0; self.one = 1
        assert p % 5 == 1
        c = 2
        while pow(c, (p - 1) // 5, p) == 1: c += 1
        self.z = pow(c, (p - 1) // 5, p)

    def of(self, x): x = Fr(x); return x.numerator * pow(x.denominator, -1, self.p) % self.p
    def add(self, a, b): return (a + b) % self.p
    def sub(self, a, b): return (a - b) % self.p
    def neg(self, a): return (-a) % self.p
    def smul(self, s, a): return self.of(s) * a % self.p
    def mul(self, a, b): return a * b % self.p
    def iszero(self, a): return a % self.p == 0
    def isunit(self, a): return a % self.p != 0
    def eq(self, a, b): return (a - b) % self.p == 0
    def zeta(self, j): return pow(self.z, j % 5, self.p)
    def inv(self, a): return pow(a, -1, self.p)


class QF:
    """the field Q (Fractions) with the same interface"""
    name = 'Q'; zero = Fr(0); one = Fr(1)
    def of(self, x): return Fr(x)
    def add(self, a, b): return a + b
    def sub(self, a, b): return a - b
    def neg(self, a): return -a
    def mul(self, a, b): return a * b
    def smul(self, s, a): return Fr(s) * a
    def iszero(self, a): return a == 0
    def isunit(self, a): return a != 0
    def eq(self, a, b): return a == b
    def inv(self, a): return 1 / a


class RS:
    """Q[[w]]/(w^D); elements are tuples of D Fractions."""
    def __init__(self, D):
        self.D = D; self.name = 'Q[[w]]/(w^%d)' % D
        self.zero = (Fr(0),) * D; self.one = (Fr(1),) + (Fr(0),) * (D - 1)
    def of(self, x): return (Fr(x),) + (Fr(0),) * (self.D - 1)
    def add(self, a, b): return tuple(x + y for x, y in zip(a, b))
    def sub(self, a, b): return tuple(x - y for x, y in zip(a, b))
    def neg(self, a): return tuple(-x for x in a)
    def smul(self, s, a): s = Fr(s); return tuple(s * x for x in a)
    def iszero(self, a): return not any(a)
    def isunit(self, a): return a[0] != 0
    def eq(self, a, b): return self.iszero(self.sub(a, b))
    def mul(self, a, b):
        D = self.D; c = [Fr(0)] * D
        for i, x in enumerate(a):
            if x:
                for j in range(D - i):
                    if b[j]: c[i + j] += x * b[j]
        return tuple(c)
    def inv(self, a):
        D = self.D; assert a[0] != 0, 'not a unit'; c = [Fr(0)] * D; c[0] = 1 / a[0]
        for m in range(1, D):
            c[m] = -sum(a[i] * c[m - i] for i in range(1, m + 1)) / a[0]
        return tuple(c)


# two primes = 1 mod 10, NOT those of the Stage-4 critic (2^61-1, 2^89-1); both < 3.3e24, so Miller-Rabin with the 13 prime
# bases 2..41 is a primality proof for them
PRIMES = (4611686018427388081, 1208925819614629174706261)


def bern(N):
    B = [Fr(0)] * (N + 1); B[0] = Fr(1)
    for m in range(1, N + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


def bk_list(K):
    """b_k = B_2k/(2k(2k-1)), k = 1..K (index 0 unused)"""
    B = bern(2 * K)
    return [None] + [B[2 * k] / (2 * k * (2 * k - 1)) for k in range(1, K + 1)]


def h_newton(N):
    """[h_0..h_N], h = h_(w_0)(z) = exp(Sum_k b_k z^(2k-1) (5^(1-2k) - p_(2k-1))), p_r = Sum_(j=1..4) (1-zeta^j)^(-r),
    by Newton's identities for the roots of y^4 - 2y^3 + 2y^2 - y + 1/5 (no zeta arithmetic)."""
    e = [Fr(2), Fr(2), Fr(1), Fr(1, 5)]
    p = [Fr(4)]
    for r in range(1, N + 2):
        s = sum((-1) ** (i - 1) * e[i - 1] * p[r - i] for i in range(1, min(r - 1, 4) + 1))
        if r <= 4: s += (-1) ** (r - 1) * r * e[r - 1]
        p.append(s)
    bk = bk_list(N // 2 + 2); H = [Fr(0)] * (N + 1)
    for k in range(1, N // 2 + 2):
        if 2 * k - 1 <= N: H[2 * k - 1] = bk[k] * (Fr(1, 5 ** (2 * k - 1)) - p[2 * k - 1])
    return sexp_q(H, N)


def sexp_q(H, N):
    """exp of a Q-series with H[0] = 0, through z^N"""
    h = [Fr(1)] + [Fr(0)] * N
    for m in range(1, N + 1):
        h[m] = sum(i * H[i] * h[m - i] for i in range(1, m + 1)) / m
    return h


# ---------------- univariate truncated series over a field F (lists, length K+1) ----------------
def smulf(F, a, b, K):
    c = [F.zero] * (K + 1)
    for i, x in enumerate(a[:K + 1]):
        if F.iszero(x): continue
        for j in range(min(len(b), K + 1 - i)):
            if not F.iszero(b[j]): c[i + j] = F.add(c[i + j], F.mul(x, b[j]))
    return c


def log1p_lin(F, c, K):
    """log(1 + t/c) = Sum_k (-1)^(k+1) t^k/(k c^k), k = 1..K (c a field element)"""
    ic = F.inv(c); out = [F.zero] * (K + 1); pw = F.one
    for k in range(1, K + 1):
        pw = F.mul(pw, ic); out[k] = F.smul(Fr((-1) ** (k + 1), k), pw)
    return out


def geom(F, c, K, e=-1):
    """(1 + t/c)^e for an integer e, through t^K"""
    ic = F.inv(c); out = [F.zero] * (K + 1); pw = F.one
    for k in range(K + 1):
        b = Fr(1)
        for i in range(k): b = b * (e - i) / (i + 1)
        out[k] = F.smul(b, pw); pw = F.mul(pw, ic)
    return out
