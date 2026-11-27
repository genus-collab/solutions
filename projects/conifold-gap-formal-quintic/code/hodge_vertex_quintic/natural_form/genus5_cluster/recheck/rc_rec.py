# -*- coding: utf-8 -*-
"""rc_rec.py -- CRT over the two Mersenne primes and rational reconstruction (own implementation).

p89 = 2^89 - 1, p61 = 2^61 - 1.  crt(r89, r61) -> residue mod M = p89 p61.
ratrec(r, M) -> (n, d, certified) with n/d == r mod M, gcd(n, d) = 1, d > 0; 'certified' is True when
2 |n| d < M, which makes n/d the unique such rational (standard uniqueness bound for rational reconstruction).
"""
from fractions import Fraction as Fr
from math import gcd

P89 = (1 << 89) - 1
P61 = (1 << 61) - 1
M = P89 * P61


def crt(r89, r61, p=P89, q=P61):
    """the residue mod p*q congruent to r89 mod p and r61 mod q"""
    m = p * q
    inv = pow(q % p, -1, p)
    t = (r89 - r61) % p * inv % p
    return (r61 + q * t) % m


def ratrec(r, mod=M):
    """rational reconstruction by the extended Euclidean algorithm on (mod, r)"""
    r = r % mod
    a0, b0 = mod, 0
    a1, b1 = r, 1
    bound = isqrt(mod // 2)
    while a1 > bound:
        q = a0 // a1
        a0, a1 = a1, a0 - q * a1
        b0, b1 = b1, b0 - q * b1
    n, d = a1, b1
    if d < 0:
        n, d = -n, -d
    gg = gcd(abs(n), d) or 1
    n //= gg
    d //= gg
    certified = d > 0 and 2 * abs(n) * d < mod and (n - r * d) % mod == 0
    return n, d, certified


def isqrt(x):
    if x < 0:
        raise ValueError
    if x == 0:
        return 0
    k = (x.bit_length() + 1) // 2
    g = 1 << k
    while True:
        h = (g + x // g) // 2
        if h >= g:
            break
        g = h
    return g


def margin(n, d, mod=M):
    """how far inside the uniqueness bound the reconstruction sits: mod / (2 |n| d)"""
    if n == 0:
        return float('inf')
    return mod / (2 * abs(n) * d)


def recon_pair(r89, r61):
    """-> (Fraction or None, n, d, certified, margin)"""
    r = crt(r89, r61)
    n, d, cert = ratrec(r)
    return (Fr(n, d) if cert else None), n, d, cert, margin(n, d)


def tofp(x, p):
    """Fraction -> F_p"""
    x = Fr(x)
    return x.numerator % p * pow(x.denominator % p, -1, p) % p


def bern(n):
    """Bernoulli number B_n as a Fraction (Akiyama-Tanigawa, independent of fq_fp)"""
    A = [Fr(1, m + 1) for m in range(n + 1)]
    for j in range(n):
        for m in range(n - j):
            A[m] = (m + 1) * (A[m] - A[m + 1])
    return A[0]
