# -*- coding: utf-8 -*-
"""rf_run_a.py -- referee run A (exact over Q, own code only):
 (1) direct TR at GENERIC points of the family x = z^4 + c3 z^3 + c2 z^2 + c1 z + V/z with rational branch points, and at
     the Gaussian control x = z + V/z;  (2) F_g by the dilaton formula, split into polynomial and log parts of Phi;
 (3) the vanishing period Int_inf^0 omega_{g,1} and the contour lemma Sum_a Res_a log z omega = Phi_g(0) - Phi_g(inf);
 (4) the closed formula of arXiv:2212.00320 Thm 7.1 (own code rf_closed) against direct TR, n = 1, 2, 3;
 (5) the closed-formula primitive T_r at z = 0 and infinity.
Usage: py -I -B rf_run_a.py gmax
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from rf_tr import TR
from rf_fam import family_x, phi_parts, leg_lp, leg_log
from rf_ser import bernoulli
import rf_closed as CL

gmax = int(sys.argv[1]) if len(sys.argv) > 1 else 2
B = bernoulli(2 * gmax + 4)
t0 = time.time()
cases = []
for r in ([1, -2, 3, -6, Fr(-3, 2)], [3, 6, -2, 4, -4]):
    X, par = family_x(r)
    cases.append(("family r=%s (c3,c2,c1,V)=%s" % ([str(a) for a in r], tuple(str(p) for p in par)), X, r, par[3]))
cases.append(("Gaussian control x = z + 9/z", {1: Fr(1), -1: Fr(9)}, [3, -3], Fr(9)))
order = [(0, 3), (1, 1), (0, 4), (1, 2), (2, 1)] + ([(1, 3), (2, 2), (3, 1)] if gmax >= 3 else [])
pts1 = [Fr(7, 3), Fr(-5, 11), Fr(13, 2)]
pts2 = [(Fr(7, 3), Fr(-5, 11)), (Fr(13, 2), Fr(1, 7))]
pts3 = [(Fr(7, 3), Fr(-5, 11), Fr(13, 2)), (Fr(1, 7), Fr(-9, 4), Fr(5, 3))]
for name, X, r, V in cases:
    print("=== %s" % name, flush=True)
    assert all(c for e, c in X.items()), "a parameter vanishes"
    tr = TR(X, r, 6 * gmax - 3)
    for gn in order:
        tr.step(*gn)
        print("  TR omega_%s: %d coefficients, symmetry violations %d (%.1fs)" % (gn, len(tr.om[gn]), tr.asym(*gn), time.time() - t0), flush=True)
    pol, lg = phi_parts(X)
    for g in range(2, gmax + 1):
        Fp = tr.pair(g, 1, [leg_lp(pol)]) / (2 - 2 * g)
        Fl = tr.pair(g, 1, [leg_log({0: lg})]) / (2 - 2 * g)
        ex = B[2 * g] / (2 * g * (2 * g - 2)) * V ** (2 - 2 * g)
        print("  F_%d: polynomial part of Phi gives %s; log part gives %s; B_2g/(2g(2g-2)) V^(2-2g) = %s; EQUAL %s"
              % (g, Fp, Fl, ex, Fp + Fl == ex))
    for g in range(1, gmax + 1):
        per = sum((c * Fr(-1, k) / (-tr.a[j]) ** k for ((j, k),), c in tr.om[(g, 1)].items()), Fr(0))
        lz = tr.pair(g, 1, [leg_log({0: Fr(1)})])
        ex = B[2 * g] / (2 * g) * V ** (1 - 2 * g)
        print("  g=%d: Int_inf^0 omega = %s, (B_2g/2g) V^(1-2g) = %s, EQUAL %s; Sum Res log z omega = %s = Phi(0)-Phi(inf): %s"
              % (g, per, ex, per == ex, lz, lz == per))
    if name.startswith("Gaussian"):
        continue
    C = CL.setup(X)
    for g in range(1, gmax + 1):
        ok = all(CL.omega_g1_value(C, g, z) == tr.value(g, 1, (z,)) for z in pts1)
        terms = CL.omega_g1_terms(C, g)
        reg0 = all(T is None or min(e for e, c in T[0].items() if c) >= T[1] for _, T in terms)
        dinf = max(CL.deg_inf(C, T) for _, T in terms if T is not None)
        phi0 = sum((CL.at_zero(C, T) for _, T in terms if T is not None), Fr(0))
        print("  closed n=1 g=%d: = TR at %d points: %s; each T_r regular at 0: %s; max deg_inf T_r = %d; Phi^A_g(0) = %s,"
              " -(B_2g/2g) V^(1-2g) = %s: %s" % (g, len(pts1), ok, reg0, dinf, phi0, -B[2 * g] / (2 * g) * V ** (1 - 2 * g),
                                                 phi0 == -B[2 * g] / (2 * g) * V ** (1 - 2 * g)), flush=True)
    for g in [g for g in range(1, 3) if (g, 2) in tr.om]:
        ok = all(CL.omega_g2_value(C, g, *b) == tr.value(g, 2, b) for b in pts2)
        print("  closed n=2 g=%d (single edge, edge expanded in hbar) = TR at %d points: %s (%.1fs)" % (g, len(pts2), ok, time.time() - t0), flush=True)
    ok = all(CL.omega_03_value(C, b) == tr.value(0, 3, b) for b in pts3)
    print("  closed (0,3) = Sum_i D_xi[prod_j (zi - zj)^-2] = TR at %d points: %s" % (len(pts3), ok))
    bad = CL.omega_g1_value(C, 1, pts1[0])
    print("  control: closed omega_{1,1} with the opposite sign differs from TR: %s" % (-bad != tr.value(1, 1, (pts1[0],))))
print("done %.1fs" % (time.time() - t0))
