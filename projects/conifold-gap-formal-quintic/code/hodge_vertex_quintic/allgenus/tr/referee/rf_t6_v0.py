# -*- coding: utf-8 -*-
"""rf_t6_v0.py -- referee check of TR_MODEL.md s.5 (u -> 0 at fixed p) with the referee engine, both primes.

Homogeneity (re-derived in REFEREE.md): omega^(u)_{g,n}(p) = Sum_m c_m u^((2-2g-n) + |m|/5) prod dp_l/p_l^(m_l+1), where
omega^(1) = Sum_m c_m prod dp_l/p_l^(m_l+1) is the expansion at p = infinity of TR on x = p^4/4 + 1/p (|p| > 1).  Regular at u = 0 iff
c_m = 0 for |m| < 5(2g-2+n).  Expansion: dp/(p - c_j)^M = Sum_{m >= M-1} binom(m, M-1) c_j^(m-M+1) dp/p^(m+1).
Control: the same expansion restricted to the terms with first-slot branch point c_0 only must NOT vanish at low |m| (the vanishing is
a cancellation between the five branch points).  Writes rf_t6_v0.log.  Run: py -I -B rf_t6_v0.py
"""
import os, sys, time
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from rf_ser import FpField, P61, P89
from rf_tr import TR, Mmax
from rf_curves import special

CASES = ((0, 3), (1, 1), (0, 4), (1, 2), (2, 1), (1, 3), (2, 2), (3, 1))
CLAIM = {(1, 1): {(5,): Fr(-5, 8)}, (2, 1): {(15,): Fr(2079, 128)}, (3, 1): {(25,): Fr(3132675, 1024)},
         (0, 3): {(1, 1, 3): Fr(-3), (1, 2, 2): Fr(-4), (1, 3, 1): Fr(-3), (2, 1, 2): Fr(-4), (2, 2, 1): Fr(-4), (3, 1, 1): Fr(-3)},
         (1, 2): {(1, 9): Fr(45, 8), (3, 7): Fr(-21, 8), (5, 5): Fr(25, 8), (7, 3): Fr(-21, 8), (9, 1): Fr(45, 8)}}


def expand(F, c, W, n, mtot):
    cur = dict(W)
    for slot in range(n):
        new = {}
        for key, v in cur.items():
            j, M = key[slot]; fixed = sum((e if isinstance(e, int) else e[1] - 1) for l, e in enumerate(key) if l != slot)
            cj = c[j]
            for m in range(M - 1, mtot - fixed + 1):
                w = v * F.of(comb(m, M - 1)); e = m - M + 1; pw = F.one
                for _ in range(e): pw = pw * cj
                nk = key[:slot] + (m,) + key[slot + 1:]
                new[nk] = new[nk] + w * pw if nk in new else w * pw
        cur = {k: x for k, x in new.items() if not F.iszero(x)}
    return cur


def main():
    log = open(os.path.join(HERE, "rf_t6_v0.log"), "w"); t0 = time.time(); fails = [0]
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    def check(ok, s):
        say("  [%s] %s" % ("PASS" if ok else "FAIL", s)); fails[0] += (not ok)
    say("rf_t6_v0.py  %s" % time.strftime("%Y-%m-%d %H:%M"))
    for p in (P61, P89):
        F = FpField(p); tag = "2^%d-1" % p.bit_length(); tr = TR(special(F, 1), 2 * Mmax(2, 2) + 6)
        for g, n in CASES:
            reg = 2 * g - 2 + n; mt = 5 * reg
            E = expand(F, tr.C.c, tr.omega(g, n), n, mt)
            low = {m: x for m, x in E.items() if sum(m) < mt}
            off5 = [m for m in E if sum(m) % 5]
            top = {m: x for m, x in E.items() if sum(m) == mt}
            check(not low and not off5 and top, "%s omega_{%d,%d}: all c_m with |m| < %d vanish (%d non-zero c_m with |m| <= the bound, all at the bound), none with |m| != 0 mod 5, %d non-zero at |m| = %d" % (
                tag, g, n, mt, sum(1 for _ in E), len(top), mt))
            if (g, n) in CLAIM:
                ok = all(F.iszero(top.get(m, F.zero) - F.of(v)) for m, v in CLAIM[(g, n)].items()) and set(top) == set(CLAIM[(g, n)])
                check(ok, "%s omega_{%d,%d}: limit u -> 0 = %s" % (tag, g, n, ", ".join("%s: %s" % (m, F.rat(top[m])) for m in sorted(top))))
            W0 = {k: v for k, v in tr.omega(g, n).items() if k[0][0] == 0}
            E0 = expand(F, tr.C.c, W0, n, mt)
            low0 = [m for m, x in E0.items() if sum(m) < mt]
            check(len(low0) > 0, "%s control: the c_0-part of the first slot alone has %d non-zero c_m with |m| < %d" % (tag, len(low0), mt))
    say("== %s (%.1f s, %d failures)" % ("ALL CHECKS PASSED" if not fails[0] else "SOME CHECKS FAILED", time.time() - t0, fails[0]))
    log.close()


if __name__ == "__main__":
    main()
