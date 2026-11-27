# -*- coding: utf-8 -*-
"""ms_onept.py -- task 2(a): the one-point law  <tau_k(H^b)>_h = PiC(tau_k H^b) d_v F_h  (h >= 2), and what happens at h = 1.

For every admissible single insertion (b + k = 1 mod 5, k <= 3h - 2) prints: v-exponent e, the correlator at v = 1, PiC (v = 1),
the prediction PiC F_h'(1), the remainder, and the alternative linear functional dW (W = w^H/w^1_x, the ancestor genus-0 point:
dW = 1 on H, -v on tau_1(1), 0 elsewhere).  Run: py -I -B ms_onept.py [hmax]
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_corr as mc


def dW(ins):
    a, k = ins
    if (a, k) == (1, 0): return Fr(1)
    if (a, k) == (0, 1): return Fr(-1)
    return Fr(0)


def run(hmax):
    allok = True; sing_ok = True; nsing = 0; ntot = 0
    for h in range(1, hmax + 1):
        print("genus %d:" % h)
        for k in range(0, 3 * h - 1):
            for b in range(5):
                if (b + k) % 5 != 1: continue
                ins = ((b, k),)
                x = mc.corr(h, ins); e = mc.vexp(h, ins); p = mc.pic((b, k))
                pred = p * mc.Fder(h, 1); rem = x - pred
                ntot += 1
                if h >= 2:
                    if rem != 0: allok = False
                    if e < 0:
                        nsing += 1
                        if rem != 0: sing_ok = False
                print("   <t%d(H%d)>_%d: e = %2d  value %-28s PiC %-10s PiC F_h' %-28s remainder %-24s dW-pred %s"
                      % (k, b, h, e, x, p, pred, rem, dW((b, k)) * mc.Fder(h, 1)))
    print("ONE-POINT LAW for h >= 2 on all %d admissible insertions (h <= %d): %s" % (ntot, hmax, "HOLDS" if allok else "FAILS somewhere"))
    print("ONE-POINT LAW on the %d singular ones (e < 0), h >= 2: %s" % (nsing, "HOLDS" if sing_ok else "FAILS"))


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
