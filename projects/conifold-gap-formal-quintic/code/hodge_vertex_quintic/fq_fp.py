# -*- coding: utf-8 -*-
"""fq_fp.py -- coefficient rings for the formal-quintic code (PLAN_2026-11-03.md, step 1-2).

Two rings with the same small interface (used by fq_series.Ser and the drivers):
  FpRing(p)  : integers mod a prime p with p = 1 mod 5; carries a primitive fifth
               root of unity  zeta  (zeta**5 = 1, zeta != 1) and rational reconstruction.
  QRing()    : fractions.Fraction (exact, no zeta).
Interface:  R.of(x) (int/Fraction -> element), R.add, R.sub, R.mul, R.neg, R.inv,
            R.zero, R.one, R.iszero(a), R.eq(a,b), R.pow(a,n), R.fmt(a) (printable form).
Conventions:  zeta^(1/2) := zeta^3  (so (zeta^(3i))^(1/2) = zeta^(4i)), used for Delta_i^(1/2).
Primes used:  P1 = 2^61-1, P2 = 2^89-1 (both prime and = 1 mod 5); default P2.
Pure Python; no third-party modules.
"""
from fractions import Fraction as Fr
import math

P1 = (1 << 61) - 1
P2 = (1 << 89) - 1


class FpRing(object):
    def __init__(self, p=P2):
        assert p % 5 == 1, "need p = 1 mod 5"
        self.p = p
        self.zero = 0
        self.one = 1
        self._inv = {}
        self.zeta = self._find_zeta()
        self.zpow = [pow(self.zeta, k, p) for k in range(5)]
        self.zhalf = self.zpow[3]          # zeta^(1/2) := zeta^3
        self.is_fp = True

    # --- basic ops
    def of(self, x):
        if isinstance(x, Fr):
            return x.numerator % self.p * self.inv(x.denominator % self.p) % self.p
        return x % self.p

    def add(self, a, b): return (a + b) % self.p
    def sub(self, a, b): return (a - b) % self.p
    def mul(self, a, b): return (a * b) % self.p
    def neg(self, a): return (-a) % self.p
    def pow(self, a, n):
        if n < 0: return pow(self.inv(a), -n, self.p)
        return pow(a, n, self.p)

    def inv(self, a):
        a %= self.p
        r = self._inv.get(a)
        if r is None:
            if a == 0: raise ZeroDivisionError("inverse of 0 mod p")
            r = pow(a, self.p - 2, self.p)
            if len(self._inv) < 200000: self._inv[a] = r
        return r

    def div(self, a, b): return a * self.inv(b) % self.p
    def iszero(self, a): return a % self.p == 0
    def eq(self, a, b): return (a - b) % self.p == 0

    # --- zeta
    def _find_zeta(self):
        p = self.p
        e = (p - 1) // 5
        for g in range(2, 1000):
            z = pow(g, e, p)
            if z != 1:
                assert pow(z, 5, p) == 1
                return z
        raise RuntimeError("no zeta found")

    def z(self, k):
        """zeta^k"""
        return self.zpow[k % 5]

    # --- rational reconstruction
    RATREC_BOUND = 10 ** 10

    def ratrec(self, a, B=None):
        """rational a/b with |a|,|b| <= B (default min(sqrt(p/2)/100, 10^10)) congruent to a mod p, or None.
        The cap 10^10 keeps spurious reconstructions rare (rate ~ B^2/p): use p = 2^89-1 for reporting."""
        p = self.p
        a %= p
        if a == 0: return Fr(0)
        r0, r1 = p, a; s0, s1 = 0, 1
        if B is None: B = min(math.isqrt(p // 2) // 100, self.RATREC_BOUND)
        while r1 > B:
            q = r0 // r1
            r0, r1 = r1, r0 - q * r1
            s0, s1 = s1, s0 - q * s1
        if s1 == 0 or abs(s1) > B: return None
        return Fr(r1, s1)

    def fmt(self, a):
        """printable: rational if reconstructible, else 'r*zeta^k' if a*zeta^-k is, else residue"""
        a %= self.p
        r = self.ratrec(a)
        if r is not None: return str(r)
        for k in range(1, 5):
            r = self.ratrec(a * self.zpow[(5 - k) % 5] % self.p)
            if r is not None: return "%s*zeta^%d" % (r, k)
        # try a + b zeta ... (sum of two): skip, print residue
        return "[%d mod p]" % a

    def fmt_zeta_poly(self, a):
        """try to write a as c0 + c1 zeta + c2 zeta^2 + c3 zeta^3 with small rationals (basis 1,z,z^2,z^3)."""
        # brute-force small search is not feasible; use fmt
        return self.fmt(a)


class QRing(object):
    """Fractions; same interface (no zeta)."""
    def __init__(self):
        self.zero = Fr(0)
        self.one = Fr(1)
        self.is_fp = False

    def of(self, x): return Fr(x)
    def add(self, a, b): return a + b
    def sub(self, a, b): return a - b
    def mul(self, a, b): return a * b
    def neg(self, a): return -a
    def pow(self, a, n): return a ** n
    def inv(self, a):
        if a == 0: raise ZeroDivisionError("inverse of 0")
        return 1 / a
    def div(self, a, b): return a / b
    def iszero(self, a): return a == 0
    def eq(self, a, b): return a == b
    def ratrec(self, a): return Fr(a)
    def fmt(self, a): return str(a)


_bern = {}
def bernoulli(n):
    """B_n (B_1 = -1/2), Fraction; Akiyama-Tanigawa."""
    if n in _bern: return _bern[n]
    A = [Fr(0)] * (n + 1)
    for m in range(n + 1):
        A[m] = Fr(1, m + 1)
        for j in range(m, 0, -1):
            A[j - 1] = j * (A[j - 1] - A[j])
        _bern.setdefault(m, A[0])
    if n >= 1: _bern[1] = Fr(-1, 2)
    return _bern[n]


if __name__ == "__main__":
    for p in (P1, P2):
        R = FpRing(p)
        assert pow(R.zeta, 5, p) == 1 and R.zeta != 1
        # 1 + z + z^2 + z^3 + z^4 = 0
        assert sum(R.zpow) % p == 0
        # prod_{k=1}^4 (1 - zeta^k) = 5
        pr = 1
        for k in range(1, 5): pr = pr * (1 - R.z(k)) % p
        assert pr == 5
        assert R.ratrec(R.of(Fr(-7, 13))) == Fr(-7, 13)
        assert R.fmt(R.of(Fr(3, 4)) * R.z(2) % p) == "3/4*zeta^2"
        assert R.mul(R.zhalf, R.zhalf) == R.zeta
        print("FpRing(%d) ok, zeta=%d" % (p, R.zeta))
    assert bernoulli(2) == Fr(1, 6) and bernoulli(4) == Fr(-1, 30) and bernoulli(12) == Fr(-691, 2730)
    print("bernoulli ok")
