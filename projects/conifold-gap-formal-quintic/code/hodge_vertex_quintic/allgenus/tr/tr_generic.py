# -*- coding: utf-8 -*-
"""tr_generic.py -- TR at generic points of the Hurwitz space H_{0;1,4} and along the family H_{0;1,r} (TR_MODEL.md s.4.2).

GenericCurve(F, N, rng): x = P(z) + v/z, P = z^4/4 + al z^3 + be z^2 + ga z, y = z, with the five branch points chosen in F_p
(c_0..c_3 random, c_4 fixed by e_4(c) = 0, which is the vanishing of the z^1-coefficient of z^2 x'(z) = z^5 + 3 al z^4 + 2 be z^3 + ga z^2 - v);
TR is run at all five branch points (no symmetry used).  The pole position a_1 is not varied: z -> z - a_1, y -> y - a_1 does not
change the omega_{g,n} (y enters only through y(q) - y(sigma(q))).
Run: py -I -B tr_generic.py  -> tr_generic.log: F_2, F_3 at 2 random points x 2 primes, and on x = p^r/r + 1/p, r = 1..6.
"""
import os, sys, time, random
from fractions import Fraction as Fr
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from tr_field import make_fp, P89, P61
from tr_core import Local, binom
from tr_family import FamilyCurve, find_prime
from tr_rec import TR
from tr_compare import bernoulli


class GenericCurve(object):
    symmetric = False; nbp = 5; name = "generic"

    def __init__(self, F, N, rng):
        self.F = F; self.N = N; self._loc = {}
        def esym(xs, k):
            tot = F.zero
            for cm in combinations(xs, k):
                pr = F.one
                for x in cm: pr = pr * x
                tot = tot + pr
            return tot
        while True:
            c = [F.of(rng.randrange(1, F.p)) for _ in range(4)]
            e3 = esym(c, 3); e4 = esym(c, 4)
            if e3.iszero(): continue
            c = c + [-(e4 / e3)]
            if len(set(x.v for x in c)) == 5 and not any(x.iszero() for x in c): break
        self.c = c; E = [esym(c, k) for k in range(6)]; assert E[4].iszero()
        self.al = -(E[1]) / 3; self.be = E[2] / 2; self.ga = -(E[3]); self.v = E[5]
        self.co = {4: F.of(Fr(1, 4)), 3: self.al, 2: self.be, 1: self.ga}

    def _poly_taylor(self, co, a, N):
        F = self.F; out = []
        for k in range(N + 1):
            s = F.zero
            for e, ce in co.items():
                if k <= e:
                    pw = F.one
                    for _ in range(e - k): pw = pw * a
                    s = s + ce * pw * F.of(binom(e, k))
            out.append(s)
        return out

    def x_taylor(self, a, N):
        F = self.F; out = self._poly_taylor(self.co, a, N); ai = F.one / a; pw = ai
        for k in range(N + 1):
            out[k] = out[k] + self.v * pw * F.of((-1) ** k); pw = pw * ai
        return out

    def local(self, j):
        L = self._loc.get(j)
        if L is None:
            L = Local(self.F, self.x_taylor(self.c[j], self.N + 2), self.c[j], self.N); self._loc[j] = L
        return L

    def ydx_taylor(self, j, N):
        """[t^k] of y x' = z P'(z) - v/z at c_j + t"""
        F = self.F; a = self.c[j]; ai = F.one / a
        out = self._poly_taylor({e: ce * e for e, ce in self.co.items()}, a, N); pw = ai
        for k in range(N + 1):
            out[k] = out[k] - self.v * pw * F.of((-1) ** k); pw = pw * ai
        return out


if __name__ == "__main__":
    f = open(os.path.join(HERE, "tr_generic.log"), "w"); t0 = time.time(); fails = 0
    def log(s): print(s, flush=True); f.write(s + "\n")
    log("tr_generic.py: F_g at generic points of H_{0;1,4} and on the family x = p^r/r + 1/p (TR mod p, rational reconstruction)")
    for seed in (1, 2):
        for p in (P89, P61):
            F = make_fp(p); C = GenericCurve(F, 2 * 15 + 8, random.Random(seed)); tr = TR(C, full=True)
            for g in (2, 3):
                val = tr.Fg(g); vv = F.one
                for _ in range(2 * g - 2): vv = vv * C.v
                ok = (val * vv - F.of(bernoulli(2 * g) / (2 * g * (2 * g - 2)))).iszero(); fails += not ok
                log("  [%s] seed %d, p = 2^%d-1, (al, be, ga, v) = %s: F_%d v^(2g-2) = %s = B_%d/(2g(2g-2)); symmetry failures of omega_{1,2}, omega_{0,4}: %d, %d" % (
                    "PASS" if ok else "FAIL", seed, p.bit_length(), [F.fmt(x) for x in (C.al, C.be, C.ga, C.v)], g, F.fmt(val * vv), 2 * g,
                    tr.symmetric_failures(1, 2), tr.symmetric_failures(0, 4)))
    for r in range(1, 7):
        p = find_prime(r + 1, 1 << 80); F = make_fp(p, order=r + 1); C = FamilyCurve(F, 2 * 15 + 8, r); tr = TR(C)
        for g in (2, 3):
            val = tr.Fg(g); ok = (val - F.of(bernoulli(2 * g) / (2 * g * (2 * g - 2)))).iszero(); fails += not ok
            log("  [%s] r = %d (p = %d = 1 mod %d): F_%d(x = p^r/r + 1/p) = %s" % ("PASS" if ok else "FAIL", r, p, r + 1, g, F.fmt(val)))
    log("== %s (%.1f s, %d failures)" % ("ALL CHECKS PASSED" if not fails else "SOME CHECKS FAILED", time.time() - t0, fails))
    f.close()
