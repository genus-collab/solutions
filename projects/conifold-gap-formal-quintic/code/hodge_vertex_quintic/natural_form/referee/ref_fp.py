# -*- coding: utf-8 -*-
"""ref_fp.py -- REFEREE code (independent of the fq_* / fq3_* modules): arithmetic helpers over F_p.

  primes: P89 = 2^89-1, P61 = 2^61-1, P521 = 2^521-1 (all Mersenne primes with exponent = 1 mod 4, hence = 1 mod 5).
  find_zeta(p)            a primitive 5th root of unity mod p (smallest base g >= 2 with g^((p-1)/5) != 1).
  u-series: plain lists of K residues (index = power of u).  smul, sinv, stheta, sone_minus (multiplication by 1-u).
  polynomials in u: plain lists (no truncation).
  5x5 matrices: lists of lists.  mmul, madd, mT, mid, mzero, mscale, miszero.
  bivariate matrix series (zhat^n u^k): dict {(n, k): 5x5 matrix}, truncated at n <= NZ, k < K:  bmul, bclean, blog, bexp, bsub.
  ratrec(a, p)            rational reconstruction with the symmetric bound sqrt(p/2) (None if it fails);  fmt(a, p).
  bernoulli(n)            Fraction (B_1 = -1/2).
"""
from fractions import Fraction as Fr
import math

P89 = (1 << 89) - 1
P61 = (1 << 61) - 1
P521 = (1 << 521) - 1


def inv(a, p): return pow(a % p, p - 2, p)


def of(x, p):
    x = Fr(x)
    return x.numerator % p * inv(x.denominator, p) % p


def find_zeta(p):
    assert p % 5 == 1
    e = (p - 1) // 5
    for g in range(2, 1000):
        z = pow(g, e, p)
        if z != 1:
            assert pow(z, 5, p) == 1 and (1 + z + z * z + pow(z, 3, p) + pow(z, 4, p)) % p == 0
            return z
    raise RuntimeError("no zeta")


# ------------------------------------------------------------------ u-series (lists of length K)
def smul(a, b, p, K=None):
    if K is None: K = min(len(a), len(b))
    out = [0] * K
    for i, x in enumerate(a[:K]):
        if not x: continue
        for j in range(min(len(b), K - i)):
            y = b[j]
            if y: out[i + j] = (out[i + j] + x * y) % p
    return out


def sinv(a, p):
    K = len(a); out = [0] * K; a0i = inv(a[0], p); out[0] = a0i
    for k in range(1, K):
        s = sum(a[j] * out[k - j] for j in range(1, k + 1)) % p
        out[k] = (-s * a0i) % p
    return out


def stheta(a, p): return [k * x % p for k, x in enumerate(a)]
def sscale(a, c, p): return [x * c % p for x in a]
def sadd(a, b, p): return [(x + y) % p for x, y in zip(a, b)]
def sone_minus(a, p):
    """(1 - u) a, same length (truncating)"""
    return [(a[k] - (a[k - 1] if k else 0)) % p for k in range(len(a))]
def sderiv(a, p):
    """d/du; the result has one coefficient fewer"""
    return [(k + 1) * a[k + 1] % p for k in range(len(a) - 1)]


# ------------------------------------------------------------------ matrices
def mid(): return [[1 if i == j else 0 for j in range(5)] for i in range(5)]
def mzero(): return [[0] * 5 for _ in range(5)]
def mmul(A, B, p): return [[sum(A[i][k] * B[k][j] for k in range(5)) % p for j in range(5)] for i in range(5)]
def madd(A, B, p, c=1): return [[(A[i][j] + c * B[i][j]) % p for j in range(5)] for i in range(5)]
def mscale(A, c, p): return [[x * c % p for x in r] for r in A]
def mT(A): return [list(r) for r in zip(*A)]
def miszero(A): return not any(x for r in A for x in r)


# ------------------------------------------------------------------ bivariate matrix series
def bclean(X): return {key: A for key, A in X.items() if not miszero(A)}


def bmul(X, Y, p, NZ, K):
    out = {}
    for (n1, k1), A in X.items():
        for (n2, k2), B in Y.items():
            n = n1 + n2; k = k1 + k2
            if n > NZ or k >= K: continue
            C = mmul(A, B, p); D = out.get((n, k))
            out[(n, k)] = C if D is None else madd(D, C, p)
    return bclean(out)


def bsub(X, Y, p):
    out = dict(X)
    for key, B in Y.items():
        A = out.get(key, mzero()); out[key] = madd(A, B, p, p - 1)
    return bclean(out)


def bscale(X, c, p): return bclean({key: mscale(A, c, p) for key, A in X.items()})


def blog(X, p, NZ, K, nil):
    """log of X = 1 + Y where Y^nil = 0 within the truncation (caller supplies nil)"""
    Y = bclean({key: A for key, A in X.items() if key != (0, 0)})
    assert X.get((0, 0)) == mid()
    out = {}; power = {(0, 0): mid()}
    for m in range(1, nil + 1):
        power = bmul(power, Y, p, NZ, K)
        if not power: break
        c = inv(m, p) * (1 if m % 2 else p - 1) % p
        for key, A in power.items():
            B = out.get(key, mzero()); out[key] = madd(B, A, p, c)
    return bclean(out)


def bexp(Y, p, NZ, K, nil):
    out = {(0, 0): mid()}; power = {(0, 0): mid()}; f = 1
    for m in range(1, nil + 1):
        power = bmul(power, Y, p, NZ, K)
        if not power: break
        f = f * m % p; c = inv(f, p)
        for key, A in power.items():
            B = out.get(key, mzero()); out[key] = madd(B, A, p, c)
    return bclean(out)


def bflip_T(X, p):
    """X(-zhat)^T"""
    return {(n, k): mscale(mT(A), 1 if n % 2 == 0 else p - 1, p) for (n, k), A in X.items()}


# ------------------------------------------------------------------ rational reconstruction
def ratrec(a, p, B=None):
    a %= p
    if a == 0: return Fr(0)
    if B is None: B = math.isqrt(p // 2)
    r0, r1 = p, a; s0, s1 = 0, 1
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B: return None
    return Fr(r1, s1)


def fmt(a, p, zeta=None, small=None):
    """rational if reconstructible with numerator/denominator below 'small' (default p^(1/3)), else r*zeta^k, else residue"""
    if small is None: small = int(round(p ** (1.0 / 3)))
    r = ratrec(a, p)
    if r is not None and abs(r.numerator) <= small and r.denominator <= small: return str(r)
    if zeta is not None:
        for k in range(1, 5):
            r = ratrec(a * pow(zeta, 5 - k, p) % p, p)
            if r is not None and abs(r.numerator) <= small and r.denominator <= small: return "%s*zeta^%d" % (r, k)
    return "[%d mod p]" % (a % p) if p < (1 << 100) else "[non-rational residue]"


_bern = {}
def bernoulli(n):
    if n in _bern: return _bern[n]
    B = [Fr(1)]
    for m in range(1, n + 1):
        B.append(-sum(Fr(math.comb(m + 1, j)) * B[j] for j in range(m)) / (m + 1))
    for m, b in enumerate(B): _bern[m] = b
    return _bern[n]


if __name__ == "__main__":
    for p in (P61, P89, P521):
        z = find_zeta(p); assert ratrec(of(Fr(-7, 13), p), p) == Fr(-7, 13)
        print("p bits", p.bit_length(), "zeta ok")
    assert bernoulli(1) == Fr(-1, 2) and bernoulli(6) == Fr(1, 42) and bernoulli(12) == Fr(-691, 2730)
    a = [1, 3, 5, 7]; assert smul(a, sinv(a, P61), P61) == [1, 0, 0, 0]
    print("ref_fp ok")
