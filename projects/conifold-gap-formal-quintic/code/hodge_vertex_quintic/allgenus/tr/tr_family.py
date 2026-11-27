# -*- coding: utf-8 -*-
"""tr_family.py -- the family x = kappa (p^r/r + 1/p), y = p (the line of H_{0;1,r} through the special point) and a prime search
for p = 1 mod (r+1) (Miller-Rabin).  Split off tr_core.py (file-size rule); used by tr_generic.py."""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from tr_core import ModelCurve, binom


class FamilyCurve(ModelCurve):
    """x = kappa (p^r/r + 1/p), y = p (the line of H_{0;1,r} through the special point); N = r + 1 branch points
    c_j = zeta^-j, zeta of order N (the field must provide F.z of order N); r = 4 is ModelCurve."""
    def __init__(self, F, N, r, kappa=1):
        self.F = F; self.N = N; self.r = r; self.kappa = Fr(kappa); self.kap = F.of(self.kappa)
        self.name = "x = %s (p^%d/%d + 1/p), y = p" % (self.kappa, r, r)
        self.nbp = r + 1
        self.c = [F.z(-j) for j in range(r + 1)]
        self._loc = {}

    def x_taylor(self, a, N):
        F = self.F; ai = F.one / a; out = []; r = self.r
        for k in range(N + 1):
            v = F.zero
            if k <= r: v = v + F.of(binom(r, k) / r) * self._pw(a, r - k)
            v = v + F.of((-1) ** k) * self._pw(ai, k + 1)
            out.append(v * self.kap)
        return out

    def ydx_taylor(self, j, N):
        """[t^k] of y x' = kappa (p^r - 1/p) at c_j + t"""
        F = self.F; a = self.c[j]; ai = F.one / a; out = []; r = self.r
        for k in range(N + 1):
            v = F.zero
            if k <= r: v = v + F.of(binom(r, k)) * self._pw(a, r - k)
            v = v - F.of((-1) ** k) * self._pw(ai, k + 1)
            out.append(v * self.kap)
        return out


def find_prime(N, start=(1 << 61)):
    """smallest prime p >= start with p = 1 mod N (Miller-Rabin with fixed bases, deterministic below 3.3e24)"""
    def isprime(n):
        if n < 2: return False
        for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
            if n % q == 0: return n == q
        d = n - 1; s = 0
        while d % 2 == 0: d //= 2; s += 1
        for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
            x = pow(a, d, n)
            if x in (1, n - 1): continue
            for _ in range(s - 1):
                x = x * x % n
                if x == n - 1: break
            else: return False
        return True
    p = start - (start % N) + 1
    while not isprime(p): p += N
    return p
