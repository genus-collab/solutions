# -*- coding: utf-8 -*-
"""rf_run_e.py -- referee run E (exact over Q, no TR):
 (1) F_g, g = 2..6, from Corollary 1.2 (eq. (1.3)) of arXiv:2502.18115 (rf_cor12) at: two GUE points (normalisation control),
     the line x = z^4 + 4/z, a generic point of the family with rational branch points, a random rational point (c0 != 0,
     irrational branch points); residue at z = 0 and at infinity printed separately;
 (2) sensitivity control: a double pole x = z^4 + 2 z^2 + 3/z + W/z^2 (the formula must NOT give c_2 V^-2);
 (3) the Borel-sum identity of XYSWAP s.2.3: Sum_r (r-1)! [eps^2g][t^(r+1)] exp(E_0) = -2^(2g) B_2g/(2g),
     E_0 = -Sum_k eps^(2k) t^(2k+1)/(2k+1), g = 1..10.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import factorial
from collections import Counter
from rf_cor12 import F_cor12
from rf_ser import bernoulli, partitions
from rf_fam import family_x

t0 = time.time()
B = bernoulli(30)
cases = [("GUE x = z + 1/z", {1: Fr(1), -1: Fr(1)}, Fr(1)), ("GUE x = z + 3/z", {1: Fr(1), -1: Fr(3)}, Fr(3)),
         ("line x = z^4 + 4/z", {4: Fr(1), -1: Fr(4)}, Fr(4))]
X, par = family_x([1, -2, 3, -6, Fr(-3, 2)])
cases.append(("family r=(1,-2,3,-6,-3/2)", X, par[3]))
cases.append(("random x = z^4+(2/7)z^3-(5/3)z^2+(11/2)z+3-(9/5)/z",
              {4: Fr(1), 3: Fr(2, 7), 2: Fr(-5, 3), 1: Fr(11, 2), 0: Fr(3), -1: Fr(-9, 5)}, Fr(-9, 5)))
print("(1) F_g from arXiv:2502.18115 Cor. 1.2")
for name, X, V in cases:
    for g in range(2, 7):
        f0, fi, F = F_cor12(X, g)
        ex = B[2 * g] / (2 * g * (2 * g - 2)) * V ** (2 - 2 * g)
        print("  %-52s g=%d: from z=0 %s, from z=inf %s; = B_2g/(2g(2g-2)) V^(2-2g): %s" % (name, g, f0, fi, F == ex))
print("(2) control, double pole: x = z^4 + 2z^2 + 3/z + W/z^2 (c_2 V^-2 = %s)" % (Fr(-1, 240) / 9))
for W in (Fr(1), Fr(1, 10), Fr(-3, 7)):
    f0, fi, F = F_cor12({4: Fr(1), 2: Fr(2), -1: Fr(3), -2: W}, 2)
    print("  W = %s: F_2 = %s (z=0: %s, z=inf: %s); equals c_2 V^-2: %s" % (W, F, f0, fi, F == Fr(-1, 2160)))
print("(3) Borel-sum identity")
for g in range(1, 11):
    tot = Fr(0)
    for lam in partitions(g):
        c = Fr(1)
        for k in lam:
            c *= Fr(-1, 2 * k + 1)
        for v in Counter(lam).values():
            c /= factorial(v)
        r = 2 * g + len(lam) - 1
        tot += factorial(r - 1) * c
    ex = -Fr(2) ** (2 * g) * B[2 * g] / (2 * g)
    print("  g=%d: %s ; -2^2g B_2g/(2g) = %s ; equal %s" % (g, tot, ex, tot == ex))
print("done %.1fs" % (time.time() - t0))
