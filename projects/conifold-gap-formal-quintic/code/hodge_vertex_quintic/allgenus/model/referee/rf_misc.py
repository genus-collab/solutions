# -*- coding: utf-8 -*-
"""rf_misc.py -- spot values quoted in the worker's notes: F_4, <kappa(H^2,4)>_1 (claimed v^3/32), the H^2-tower exponents e = V >= 1
for all pure kappa(H^2, k) correlators in the referee tables (h <= 3, n <= 3), and the regular-part example <tau_3(H^3)>_2.
Run: py -I -B rf_misc.py"""
import itertools, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_flat as rf
import rf_tables as rt


def kcorr(h, gl):
    """exact <kappa(H^b1,k1) ...>_h via the multilinear expansion"""
    tot = Fr(0)
    for combo in itertools.product(*[rf.kappa(b, k) for b, k in gl]):
        ins = tuple(x[1] for x in combo)
        if sum(k for a, k in ins) > 3 * h - 3 + len(ins): continue
        c = Fr(1)
        for cf, _ in combo: c *= cf
        tot += c * rt.exact_val(h, ins)
    return tot


if __name__ == "__main__":
    print("F_4 =", rt.exact_val(4, ()), "(record 1/180000)")
    print("<kappa(H^2,4)>_1 =", kcorr(1, ((2, 4),)), "at v = 1, v-exponent", rf.vexp(1, ((2, 4),)), "(worker: v^3/32)")
    print("<tau_3(H^3)>_2 =", rt.exact_val(2, ((3, 3),)), "(worker: 7/5760, e = 0)")
    bad = 0; seen = 0
    for h, n in ((1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2)):
        for ks in itertools.combinations_with_replacement(range(0, 3 * h + n + 4), n):
            if sum(k + 1 for k in ks) % 5: continue
            B = sum(k + 1 for k in ks) // 5; e = rf.vexp(h, tuple((2, k) for k in ks)); Vh = 2 - 2 * h + 4 * B - n
            if e != Vh: bad += 1
            val = kcorr(h, tuple((2, k) for k in ks))
            if val != 0:
                seen += 1
                if e < 1: bad += 1; print("   H^2 tower: non-zero with e = %s < 1: h=%d ks=%s value %s" % (e, h, ks, val))
    print("H^2 tower: e = 2 - 2h + 4B - n for all tuples; %d non-zero pure kappa(H^2) correlators, all with e >= 1: %s" % (seen, bad == 0))
