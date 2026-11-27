# -*- coding: utf-8 -*-
"""lm_harness.py -- referee lmi_adversary: logging, PASS/FAIL accounting, comparisons, primes (stdlib only).
Logs are written to ./logs/<name>.log next to the scripts; they print relative names only."""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


class Run:
    def __init__(self, name):
        os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
        self.name = name
        self.fh = open(os.path.join(HERE, 'logs', name + '.log'), 'w', encoding='utf-8')
        self.npass = self.nfail = 0
        self.t0 = time.time()
        self.say('log logs/%s.log  (python %s)' % (name, sys.version.split()[0]))

    def say(self, s):
        print(s)
        self.fh.write(s + '\n')
        self.fh.flush()

    def check(self, tag, ok, info=''):
        if ok:
            self.npass += 1
        else:
            self.nfail += 1
        self.say('%s  %s  %s' % ('PASS' if ok else 'FAIL', tag, info))
        return ok

    def done(self):
        self.say('%d checks: %d PASS, %d FAIL;  %.1f s' % (self.npass + self.nfail, self.npass, self.nfail,
                                                           time.time() - self.t0))
        self.say('ALL CHECKS PASSED' if self.nfail == 0 else 'SOME CHECKS FAILED')
        self.fh.close()


def ls_equal(x, y):
    """compare two LS on every coefficient known for both (from the lower valuation); returns (equal, count)"""
    lo = min(x.v if x.c else x.prec, y.v if y.c else y.prec)
    hi = min(x.prec, y.prec)
    if hi >= (1 << 39):                       # both exact: compare the stored coefficients only
        lo = min(x.v, y.v)
        hi = max(x.v + len(x.c), y.v + len(y.c))
    n = 0
    for k in range(lo, hi):
        if x.coef(k) != y.coef(k):
            return False, n
        n += 1
    return True, n


def regular(cs, upto=None):
    """True iff no c_n (n <= upto) has a known non-zero coefficient of a negative power of e, and prec >= 1"""
    for n, x in enumerate(cs):
        if upto is not None and n > upto:
            break
        if x.prec < 1 or (x.c and x.v < 0):
            return False
    return True


def poles(cs):
    return [(-x.v if (x.c and x.v < 0) else 0) for x in cs]


def values(cs):
    return [x.coef(0) for x in cs]


def first_pole(cs):
    for n, x in enumerate(cs):
        if x.c and x.v < 0:
            return n, -x.v, x.c[0]
    return None


def is_prime(n):
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:                      # deterministic for n < 3.3e24
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


def prime_below(b, mod=20, res=1):
    p = b - 1
    while p % mod != res:
        p -= 1
    while not is_prime(p):
        p -= mod
    return p


def root_of_unity(p, k):
    """a primitive k-th root of unity mod p (k prime, k | p - 1)"""
    for g in range(2, 1000):
        z = pow(g, (p - 1) // k, p)
        if z != 1:
            return z
    raise ValueError


def sqrt_m1(p):
    for g in range(2, 1000):
        z = pow(g, (p - 1) // 4, p)
        if z * z % p == p - 1:
            return z
    raise ValueError
