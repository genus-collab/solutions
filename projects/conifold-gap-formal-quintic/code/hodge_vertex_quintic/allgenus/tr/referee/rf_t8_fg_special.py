# -*- coding: utf-8 -*-
"""rf_t8_fg_special.py -- referee: F_4 (both primes) and F_5 (2^61-1) of TR on the special curve x = W = p^4/4 + 1/p, y = p, with the
referee engine (all five branch points, no symmetry), against B_2g/(2g(2g-2)) and the frozen vertex values c_g = (-5)^(1-g) F_g
(c_4 = 1/180000, c_5 = 1/660000 on record).  Writes rf_t8_fg_special.log.  Run: py -I -B rf_t8_fg_special.py
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from rf_ser import FpField, P61, P89, bernoulli
from rf_tr import TR, Mmax
from rf_curves import special

RECORD = {4: Fr(1, 180000), 5: Fr(1, 660000)}


def main():
    log = open(os.path.join(HERE, "rf_t8_fg_special.log"), "w"); t0 = time.time(); fails = [0]
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    def check(ok, s):
        say("  [%s] %s" % ("PASS" if ok else "FAIL", s)); fails[0] += (not ok)
    say("rf_t8_fg_special.py  %s" % time.strftime("%Y-%m-%d %H:%M"))
    for p, gmax in ((P89, 4), (P61, 5)):
        F = FpField(p); tr = TR(special(F, 1), 2 * Mmax(gmax - 1, 2) + 6)
        for g in range(4, gmax + 1):
            t = time.time(); val = tr.Fg(g); pred = bernoulli(2 * g) / (2 * g * (2 * g - 2)); fr = F.rat(val * F.of(Fr(-5) ** (1 - g)))
            check(F.iszero(val - F.of(pred)) and fr == RECORD[g], "mod 2^%d-1: F_%d(x = W) = %s (B_2g/(2g(2g-2)) = %s); (-5)^(1-g) F_g = %s, record c_%d = %s (%.1f s)" % (
                p.bit_length(), g, F.rat(val), pred, fr, g, RECORD[g], time.time() - t))
    say("== %s (%.1f s, %d failures)" % ("ALL CHECKS PASSED" if not fails[0] else "SOME CHECKS FAILED", time.time() - t0, fails[0]))
    log.close()


if __name__ == "__main__":
    main()
