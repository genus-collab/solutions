# -*- coding: utf-8 -*-
"""rf_t7_extra.py -- referee extras (own code), both primes:
(a) Phi-pairing  P_{i,d} := Sum_j Res_{c_j} Phi d((d/dx)^d xit^i) = -delta_{d,1}/c_i for d = 0..10, kappa = 1, 4, -1
    (this is what makes the TR free energy F_g = (1/(2-2g)) Sum Res Phi omega_{g,1} equal the CohFT value <tau_1(1)>_g/(2g-2) in ALL
    genera; the worker checked d <= 3; REFEREE.md s.3 gives the residue proof);
(b) primary sector at the special point: all flat primary correlators <t0(H^a) t0(H^b) t0(H^c)>_2 (genus 2, n = 3) and
    <t0(H^a) t0(H^b)>_3 (genus 3, n = 2), and, for information, <t0(H^a) t0(H^b) t0(H^c)>_1;  prediction of the primary-sector lemma
    (g >= 2): only a = b = c = 1 survives, with <H,H,H>_2 = d^3/du^3 (c_2 u^-2) = -24 c_2 = -1/50 and <H,H>_3 = 20 c_3 = 1/1260.
Writes rf_t7_extra.log.  Run: py -I -B rf_t7_extra.py
"""
import os, sys, time, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from rf_ser import FpField, P61, P89
from rf_tr import TR, Mmax
from rf_curves import special
from rf_xi import XiSpecial


def main():
    log = open(os.path.join(HERE, "rf_t7_extra.log"), "w"); t0 = time.time(); fails = [0]
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    def check(ok, s):
        say("  [%s] %s" % ("PASS" if ok else "FAIL", s)); fails[0] += (not ok)
    say("rf_t7_extra.py  %s" % time.strftime("%Y-%m-%d %H:%M"))
    res = {}
    for p in (P61, P89):
        F = FpField(p); tag = "2^%d-1" % p.bit_length()
        for kappa in (1, 4, -1):
            C = special(F, kappa); tr = TR(C, 8); XB = XiSpecial(F, kappa, 10); ok = True
            for i in range(5):
                for d in range(11):
                    P = sum((v * tr.phi_res(j, M) for (j, M), v in XB.dxi[(i, d)].items()), F.zero)
                    target = -(F.one / XB.c[i]) if d == 1 else F.zero
                    ok = ok and F.iszero(P - target)
            check(ok, "%s kappa = %s: Sum_j Res_{c_j} Phi d xit^i_d = -delta_{d,1}/c_i for i = 0..4, d = 0..10" % (tag, kappa))
        tr = TR(special(F, 1), 2 * Mmax(2, 3) + 6); XB = XiSpecial(F, 1, 8)
        for g, n in ((1, 3), (2, 3), (3, 2)):
            t = time.time(); col, nbad = XB.frozen(tr.omega(g, n), g, n)
            nz = {}
            for avec in itertools.combinations_with_replacement(range(5), n):
                v = XB.flat(col, avec, (0,) * n)
                if not F.iszero(v): nz[avec] = F.rat(v)
            res.setdefault((g, n), []).append(nz)
            say("  %s (%d,%d): non-zero flat primary correlators <prod t0(H^a)> (a-tuples): %s (%.1f s, %d remainders)" % (
                tag, g, n, ", ".join("%s: %s" % (k, v) for k, v in sorted(nz.items())), time.time() - t, nbad))
    for (g, n), pred in (((2, 3), {(1, 1, 1): Fr(-1, 50)}), ((3, 2), {(1, 1): Fr(1, 1260)})):
        check(all(r == pred for r in res[(g, n)]), "primary-sector lemma at (%d,%d), both primes: only H^1 insertions survive, value %s" % (
            g, n, ", ".join("%s: %s" % (k, v) for k, v in pred.items())))
    say("== %s (%.1f s, %d failures)" % ("ALL CHECKS PASSED" if not fails[0] else "SOME CHECKS FAILED", time.time() - t0, fails[0]))
    log.close()


if __name__ == "__main__":
    main()
