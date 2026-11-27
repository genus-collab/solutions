# -*- coding: utf-8 -*-
"""rf_t1_regress.py -- regression of the referee's TR engine (rf_tr.py): Airy curve against intersection numbers of psi classes
(table below, standard values), and the Gaussian curve x = p + 1/p, y = p (two rational branch points, exact over Q) against
F_g = B_2g/(2g(2g-2)) (the classical Euler-characteristic values, the expected EO-convention result).  Writes rf_t1_regress.log.
Run: py -I -B rf_t1_regress.py
"""
import os, sys, time, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from rf_ser import QField, FpField, P61, bernoulli
from rf_tr import TR
from rf_curves import airy, pole_curve

PSI = {(0, (0, 0, 0)): Fr(1), (1, (1,)): Fr(1, 24), (0, (0, 0, 0, 1)): Fr(1), (1, (1, 1)): Fr(1, 24), (1, (0, 2)): Fr(1, 24),
       (2, (4,)): Fr(1, 1152), (1, (1, 1, 1)): Fr(1, 12), (1, (0, 1, 2)): Fr(1, 12), (1, (0, 0, 3)): Fr(1, 24),
       (2, (1, 4)): Fr(1, 384), (2, (2, 3)): Fr(29, 5760), (2, (0, 5)): Fr(1, 1152), (3, (7,)): Fr(1, 82944)}


def dfact(n):
    r = 1
    while n > 1: r *= n; n -= 2
    return r


def predicted(g, n):
    out = {}
    for (gg, ds), val in PSI.items():
        if gg != g or len(ds) != n: continue
        for perm in set(itertools.permutations(ds)):
            c = val
            for d in perm: c *= dfact(2 * d + 1)
            out[tuple((0, 2 * d + 2) for d in perm)] = c
    return out


def main():
    log = open(os.path.join(HERE, "rf_t1_regress.log"), "w"); t0 = time.time(); fails = [0]
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    def check(ok, s):
        say("  [%s] %s" % ("PASS" if ok else "FAIL", s)); fails[0] += (not ok)
    say("rf_t1_regress.py (referee engine rf_tr.py)  %s" % time.strftime("%Y-%m-%d %H:%M"))
    for F in (QField(), FpField(P61)):
        tr = TR(airy(F), 34)
        for g, n in ((0, 3), (1, 1), (0, 4), (1, 2), (2, 1), (1, 3), (2, 2), (3, 1)):
            W = tr.omega(g, n); P = predicted(g, n); keys = set(W) | set(P)
            bad = [k for k in keys if not F.iszero(W.get(k, F.zero) - F.of(P.get(k, 0)))]
            check(not bad and tr.symmetric(g, n) == 0, "%s Airy omega_{%d,%d} = psi-class intersection numbers x prod (2d+1)!! dp/p^(2d+2): %d keys, %d mismatches" % (("Q" if F.name == "Q" else "2^61-1"), g, n, len(keys), len(bad)))
        for g, n in ((0, 3), (1, 1), (1, 2), (0, 4), (2, 1), (1, 3)):
            D = tr.dilaton(tr.omega(g, n + 1)); low = tr.omega(g, n)
            ok = set(D) == set(k for k, v in low.items() if not F.iszero(v)) and all(F.iszero(D[k] - low[k] * F.of(2 - 2 * g - n)) for k in D)
            check(ok, "%s Airy dilaton: Sum Res Phi omega_{%d,%d} = (2-2g-n) omega_{%d,%d}" % (("Q" if F.name == "Q" else "2^61-1"), g, n + 1, g, n))
    # Gaussian curve: x = p + 1/p, y = p; branch points +1, -1
    F = QField()
    C = pole_curve(F, [0, 1], 1, 0, [F.one, -F.one], "x = p + 1/p")
    tr = TR(C, 40)
    for g in (2, 3, 4):
        t = time.time(); val = tr.Fg(g); pred = bernoulli(2 * g) / (2 * g * (2 * g - 2))
        check(val == pred, "Gaussian curve x = p + 1/p (exact Q): F_%d = %s, B_2g/(2g(2g-2)) = %s (%.1f s)" % (g, val, pred, time.time() - t))
    check(tr.symmetric(1, 2) == 0 and tr.symmetric(0, 4) == 0 and tr.symmetric(2, 2) == 0, "Gaussian curve: omega_{1,2}, omega_{0,4}, omega_{2,2} symmetric")
    for g, n in ((1, 1), (2, 1), (0, 3), (1, 2)):
        D = tr.dilaton(tr.omega(g, n + 1)); low = tr.omega(g, n)
        ok = all(F.iszero(D.get(k, F.zero) - low.get(k, F.zero) * F.of(2 - 2 * g - n)) for k in set(D) | set(low))
        check(ok, "Gaussian curve dilaton equation omega_{%d,%d} -> omega_{%d,%d}" % (g, n + 1, g, n))
    say("== %s (%.1f s, %d failures)" % ("ALL CHECKS PASSED" if not fails[0] else "SOME CHECKS FAILED", time.time() - t0, fails[0]))
    log.close()


if __name__ == "__main__":
    main()
