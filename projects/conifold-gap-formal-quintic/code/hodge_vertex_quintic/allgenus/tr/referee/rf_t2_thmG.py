# -*- coding: utf-8 -*-
"""rf_t2_thmG.py -- referee test of TR_MODEL.md s.4 (Theorem G: F_g = B_2g/(2g(2g-2)) v^(2-2g) for x = P(p) + v/(p - b), y = p).

(a) EXACT over Q at random points with RATIONAL branch points (rf_curves.rational_point), deg P = r = 1..5, F_2 and F_3;
(b) the pole moved to b != 0 (equivalently y -> y + const relative to the pole), r = 4;
(c) the x-moment Sum Res x omega_{g,1} (zero <=> F_g invariant under y -> y + const);
(d) F_4 (and F_5 for r = 2) mod 2^61-1 and 2^89-1 at rational points (rational data reduced mod p), r = 2, 3, 4.
Writes rf_t2_thmG.log.  Run: py -I -B rf_t2_thmG.py
"""
import os, sys, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from rf_ser import QField, FpField, P61, P89, bernoulli
from rf_tr import TR, Mmax
from rf_curves import rational_point


def pred(g, v): return bernoulli(2 * g) / (2 * g * (2 * g - 2)) * Fr(v) ** (2 - 2 * g)


def main():
    log = open(os.path.join(HERE, "rf_t2_thmG.log"), "w"); t0 = time.time(); fails = [0]
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    def check(ok, s):
        say("  [%s] %s" % ("PASS" if ok else "FAIL", s)); fails[0] += (not ok)
    say("rf_t2_thmG.py  %s" % time.strftime("%Y-%m-%d %H:%M"))
    F = QField()
    for r, seed, gmax in ((1, 11, 3), (2, 12, 3), (3, 13, 3), (4, 14, 3), (4, 15, 2), (5, 16, 2)):
        C, Pc, v, c = rational_point(F, r, random.Random(seed))
        say("  curve r = %d (seed %d): P coeffs %s, v = %s, branch points %s" % (r, seed, [str(x) for x in Pc], v, [str(x) for x in c]))
        tr = TR(C, 2 * Mmax(gmax - 1, 2) + 6)
        for g in range(2, gmax + 1):
            t = time.time(); val = tr.Fg(g)
            check(val == pred(g, v), "exact Q, r = %d: F_%d = %s, B_2g/(2g(2g-2)) v^(2-2g) = %s (%.1f s)" % (r, g, val, pred(g, v), time.time() - t))
        for g in range(2, gmax + 1):
            say("     x-moment Sum Res x omega_{%d,1} = %s" % (g, tr.xmoment(tr.omega(g, 1))))
    for b in (Fr(3, 2), Fr(-2, 7)):
        C, Pc, v, c = rational_point(F, 4, random.Random(14), b)
        tr = TR(C, 2 * Mmax(2, 2) + 6)
        for g in (2, 3):
            val = tr.Fg(g)
            check(val == pred(g, v), "exact Q, r = 4, pole at b = %s (y shifted relative to the pole): F_%d = %s (expected %s)" % (b, g, val, pred(g, v)))
    for p in (P61, P89):
        Fp = FpField(p)
        for r, seed, gmax in ((4, 21, 4), (3, 22, 4), (2, 23, 5)):
            C, Pc, v, c = rational_point(Fp, r, random.Random(seed))
            tr = TR(C, 2 * Mmax(gmax - 1, 2) + 6)
            for g in range(4, gmax + 1):
                t = time.time(); val = tr.Fg(g)
                check(Fp.iszero(val - Fp.of(pred(g, v))), "mod 2^%d-1, r = %d, v = %s: F_%d v^(2g-2) = %s (expected %s) (%.1f s)" % (
                    p.bit_length(), r, v, g, Fp.rat(val * Fp.of(Fr(v) ** (2 * g - 2))), bernoulli(2 * g) / (2 * g * (2 * g - 2)), time.time() - t))
    say("== %s (%.1f s, %d failures)" % ("ALL CHECKS PASSED" if not fails[0] else "SOME CHECKS FAILED", time.time() - t0, fails[0]))
    log.close()


if __name__ == "__main__":
    main()
