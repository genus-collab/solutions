# -*- coding: utf-8 -*-
"""rf_run_b.py -- referee run B: the model lemma ML1 OFF the symmetric line, exact over Q, own TR.
Path in the Hurwitz family: branch points r(s) = (1, 2, -2/3, s, -s), so (S = s^2)
   c3 = -28/9, c2 = -2 S, c1 = 16/3 + (28/3) S, V = (16/3) S:   V -> 0 as s -> 0 while c3, c1 stay non-zero; two branch points
   collide with the pole z = 0 and the curve degenerates to x0 = z^4 - (28/9) z^3 + (16/3) z (branch points 1, 2, -2/3).
ML1 predicts: every correlator with all legs in L (polynomials in z) is a polynomial in (c, V), hence a polynomial in S of
degree <= floor(W/2), W = 5(2-2g-n) + Sum deg (weighted homogeneity), and = 0 if W < 0.  Test: fit through deg+1 values of S,
check at the remaining values; compare the S -> 0 value with direct TR on the limiting curve x0.  Controls: legs log z,
z^-5, 1/z, x (function of x: must vanish), mixed pairs.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from rf_tr import TR
from rf_fam import family_x, leg_lp, leg_log
from rf_ser import lp_mul

t0 = time.time()
SV = [Fr(1, 2), Fr(1, 3), Fr(1, 4), Fr(1, 5), Fr(5, 7), Fr(3, 4), Fr(2, 5), Fr(3, 5), Fr(1, 6), Fr(5, 6), Fr(4, 7), Fr(3, 7)]
ORDER = [(0, 3), (1, 1), (0, 4), (1, 2), (2, 1)]
L1 = [("z^%d" % d, {d: Fr(1)}, d) for d in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15)]
L2 = [((a, b), a + b) for (a, b) in ((1, 4), (1, 9), (2, 3), (5, 5), (5, 10), (3, 8), (4, 6), (1, 14))]
L3 = [((a, b, c), a + b + c) for (a, b, c) in ((1, 1, 3), (1, 2, 2), (1, 1, 8), (2, 3, 5))]
BAD = [("-log z", leg_log({0: Fr(-1)})), ("z^-5", leg_lp({-5: Fr(1)})), ("z^-1", leg_lp({-1: Fr(1)}))]


def interp_check(xs, ys, deg):
    """fit a polynomial of degree deg through the first deg+1 points; return (ok on the rest, value at 0)"""
    n = deg + 1
    if n > len(xs) - 2:
        return None, None
    def P(t):
        tot = Fr(0)
        for i in range(n):
            li = Fr(1)
            for j in range(n):
                if j != i:
                    li *= (t - xs[j]) / (xs[i] - xs[j])
            tot += ys[i] * li
        return tot
    return all(P(xs[k]) == ys[k] for k in range(n, len(xs))), P(Fr(0))


data = {}
for s in SV:
    X, par = family_x([1, 2, Fr(-2, 3), s, -s])
    tr = TR(X, [1, 2, Fr(-2, 3), s, -s], 9)
    for gn in ORDER:
        tr.step(*gn)
    for g in (1, 2):
        for name, L, d in L1:
            data.setdefault(("g%d" % g, name), []).append(tr.pair(g, 1, [leg_lp(L)]))
        for name, lg in BAD:
            data.setdefault(("g%d" % g, name), []).append(tr.pair(g, 1, [lg]))
        xl = lp_mul(X, X)
        data.setdefault(("g%d" % g, "x^2 (function of x)"), []).append(tr.pair(g, 1, [leg_lp(xl)]))
    for (a, b), d in L2:
        data.setdefault(("(1,2)", "z^%d,z^%d" % (a, b)), []).append(tr.pair(1, 2, [leg_lp({a: Fr(1)}), leg_lp({b: Fr(1)})]))
    data.setdefault(("(1,2)", "z^5,-log z (mixed)"), []).append(tr.pair(1, 2, [leg_lp({5: Fr(1)}), BAD[0][1]]))
    data.setdefault(("(1,2)", "z^5,z^-5 (mixed)"), []).append(tr.pair(1, 2, [leg_lp({5: Fr(1)}), BAD[1][1]]))
    for (a, b, c), d in L3:
        data.setdefault(("(0,3)", "z^%d,z^%d,z^%d" % (a, b, c)), []).append(tr.pair(0, 3, [leg_lp({a: Fr(1)}), leg_lp({b: Fr(1)}), leg_lp({c: Fr(1)})]))
    print("s = %s done (%.1fs)" % (s, time.time() - t0), flush=True)

# limiting curve (V = 0): three branch points
X0 = {4: Fr(1), 3: Fr(-28, 9), 1: Fr(16, 3)}
tr0 = TR(X0, [1, 2, Fr(-2, 3)], 9)
for gn in ORDER:
    tr0.step(*gn)
lim = {}
for g in (1, 2):
    for name, L, d in L1:
        lim[("g%d" % g, name)] = tr0.pair(g, 1, [leg_lp(L)])
for (a, b), d in L2:
    lim[("(1,2)", "z^%d,z^%d" % (a, b))] = tr0.pair(1, 2, [leg_lp({a: Fr(1)}), leg_lp({b: Fr(1)})])
for (a, b, c), d in L3:
    lim[("(0,3)", "z^%d,z^%d,z^%d" % (a, b, c))] = tr0.pair(0, 3, [leg_lp({a: Fr(1)}), leg_lp({b: Fr(1)}), leg_lp({c: Fr(1)})])

xs = [s * s for s in SV]
weights = {}
for g in (1, 2):
    for name, L, d in L1:
        weights[("g%d" % g, name)] = 5 * (2 - 2 * g - 1) + d
for (a, b), d in L2:
    weights[("(1,2)", "z^%d,z^%d" % (a, b))] = 5 * (2 - 2 - 2) + d
for (a, b, c), d in L3:
    weights[("(0,3)", "z^%d,z^%d,z^%d" % (a, b, c))] = 5 * (2 - 0 - 3) + d
print("--- legs in L (polynomials in z): W = weight; ML1 predicts a polynomial in S = s^2 of degree <= floor(W/2), 0 if W < 0")
for key, ys in data.items():
    if key not in weights:
        continue
    W = weights[key]
    if W < 0:
        print("  %s %-14s W=%3d  all values zero: %s" % (key[0], key[1], W, all(y == 0 for y in ys)))
        continue
    ok, v0 = interp_check(xs, ys, W // 2)
    print("  %s %-14s W=%3d  polynomial of degree <= %d in S: %s;  value at S=0 %s;  = TR of the limiting curve x0: %s"
          % (key[0], key[1], W, W // 2, ok, v0, v0 == lim.get(key)))
print("--- controls (bad legs / function of x): S^k * value for the smallest k making it a polynomial of degree <= 12")
for key, ys in data.items():
    if key in weights:
        continue
    if all(y == 0 for y in ys):
        print("  %s %-22s identically zero on the path" % key)
        continue
    for k in range(0, 8):
        ok, v0 = interp_check(xs, [y * x ** k for y, x in zip(ys, xs)], 8)
        if ok:
            print("  %s %-22s S^%d * value is polynomial (deg <= 8); value at S=0 of S^k*value: %s  => pole order %s in S"
                  % (key[0], key[1], k, v0, k if v0 != 0 else "< %d" % k))
            break
    else:
        print("  %s %-22s no pole order <= 7 found" % key)
print("done %.1fs" % (time.time() - t0))
