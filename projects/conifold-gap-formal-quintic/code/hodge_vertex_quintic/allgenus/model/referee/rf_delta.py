# -*- coding: utf-8 -*-
"""rf_delta.py -- the worker's Delta_g(m) (first-order effect of a scalar factor exp(c z^m), m odd, on the genus-g free energy of the
frozen theory) recomputed from the referee tables:
   Delta_g(m) = <tau_(m+1)(1)>_g + (1/2) Sum_(j=0)^(m-1) (-1)^j [ (1/5) Sum_a <tau_(m-1-j)(H^a) tau_j(H^(3-a))>_(g-1)
                                                             + (1/5) Sum_a Sum_(h1+h2=g, h_i>=1) <tau_(m-1-j)(H^a)>_h1 <tau_j(H^(3-a))>_h2 ],
(leaf c z^(m+1) 1; edge c (w^m + z^m)/(w + z) eta^-1, eta^-1 = (1/5) Sum_a N_a N_(3-a) H^a (x) H^(3-a), N_a N_(3-a) = 1).
Derivation of the formula: REFEREE.md s.3 item 10.  Run: py -I -B rf_delta.py
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_tables as rt
import rf_flat as rf


def corr(h, ins):
    n = len(ins)
    if h < 0 or 2 * h - 2 + n <= 0: return Fr(0)
    if sum(k for a, k in ins) > 3 * h - 3 + n: return Fr(0)
    return rt.exact_val(h, tuple(ins))


def delta(g, m):
    leaf = corr(g, ((0, m + 1),)); loop = Fr(0); br = Fr(0)
    for j in range(m):
        s = (-1) ** j
        for a in range(5):
            b = (3 - a) % 5
            if g - 1 >= 1 or (g - 1 == 0 and False): loop += s * corr(g - 1, ((a, m - 1 - j), (b, j))) / 5
            for h1 in range(1, g):
                br += s * corr(h1, ((a, m - 1 - j),)) * corr(g - h1, ((b, j),)) / 5
    return leaf, loop / 2, br / 2


if __name__ == "__main__":
    for g, m in ((2, 1), (2, 3), (2, 5), (3, 1), (3, 3), (3, 5), (4, 1), (4, 3), (4, 5), (4, 7)):
        lf, lp, br = delta(g, m)
        print("g = %d, m = %d: leaf %s  loop %s  bridges %s  ->  Delta = %s   (v-exponent %s)"
              % (g, m, lf, lp, br, lf + lp + br, rf.vexp(g, ((0, m + 1),))))
