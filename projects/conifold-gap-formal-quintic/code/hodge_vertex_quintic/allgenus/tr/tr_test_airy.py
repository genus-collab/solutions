# -*- coding: utf-8 -*-
"""tr_test_airy.py -- regression of tr_rec on the Airy curve x = -p^2/2, y = p (DNOPS conventions (2-26)-(2-35)):
the local coordinate w = p (x = -w^2/2), Delta = 1, xi = -1/p, xi_d = (d/dx)^d xi = -(2d-1)!!/p^(2d+1), d xi_d = (2d+1)!! e_{0,2d+1},
so omega_{g,n} = Sum <tau_d1 .. tau_dn>_g prod (2d_j+1)!! e_{0,2d_j+1}(p_j)  (Witten-Kontsevich numbers, exact over Q(zeta) and mod p).
Also: dilaton equation Sum Res Phi omega_{g,n+1} = (2-2g-n) omega_{g,n}, string equation, symmetry.  Run: py -I -B tr_test_airy.py
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from tr_field import QzField, make_fp, P89
from tr_core import AiryCurve
from tr_rec import TR

WK = {(1, (1,)): Fr(1, 24), (0, (0, 0, 0)): Fr(1), (0, (1, 0, 0, 0)): Fr(1), (1, (1, 1)): Fr(1, 24), (1, (2, 0)): Fr(1, 24),
      (2, (4,)): Fr(1, 1152), (2, (4, 1)): Fr(1, 384), (2, (3, 2)): Fr(29, 5760), (2, (5, 0)): Fr(1, 1152), (3, (7,)): Fr(1, 82944),
      (0, (2, 0, 0, 0, 0)): Fr(1), (0, (1, 1, 0, 0, 0)): Fr(2), (1, (1, 1, 1)): Fr(1, 12), (1, (2, 1, 0)): Fr(1, 12),
      (1, (3, 0, 0)): Fr(1, 24), (1, (2, 2, 0, 0)) : Fr(1, 6)}


def dfact(n):
    r = 1
    while n > 1: r *= n; n -= 2
    return r


def predicted(g, n):
    """{key: Fraction} for omega_{g,n} of the Airy curve from the WK numbers in the table"""
    import itertools
    out = {}
    for (gg, ds), v in WK.items():
        if gg != g or len(ds) != n: continue
        for perm in set(itertools.permutations(ds)):
            key = tuple((0, 2 * d + 1) for d in perm)
            c = v
            for d in perm: c *= dfact(2 * d + 1)
            out[key] = c
    return out


def run(F, log):
    t0 = time.time(); tr = TR(AiryCurve(F, 40)); ok = True
    for g, n in ((0, 3), (1, 1), (0, 4), (1, 2), (2, 1), (0, 5), (1, 3), (2, 2), (3, 1)):
        W = tr.omega(g, n); P = predicted(g, n)
        if (g, n) == (1, 3): P = {k: v for k, v in P.items()}
        keys = set(W) | set(P)
        bad = [k for k in keys if not (W.get(k, F.zero) - F.of(P.get(k, 0))).iszero()]
        sf = tr.symmetric_failures(g, n)
        log("  Airy omega_{%d,%d}: %d coefficients, %d mismatches with Witten-Kontsevich, symmetry failures %d" % (g, n, len(W), len(bad), sf))
        ok = ok and not bad and sf == 0
    # dilaton and string equations
    for g, n in ((0, 3), (1, 1), (1, 2), (2, 1), (0, 4)):
        Wn1 = tr.omega(g, n + 1); D = tr.dilaton(Wn1); Wn = tr.omega(g, n)
        okd = set(D) == set(Wn) and all((D[k] - Wn[k] * (2 - 2 * g - n)).iszero() for k in Wn)
        log("  dilaton  Sum Res Phi omega_{%d,%d} = (2-2g-n) omega_{%d,%d}: %s" % (g, n + 1, g, n, okd)); ok = ok and okd
    log("  Airy regression over %s: %s (%.1f s)" % (F.name, "PASS" if ok else "FAIL", time.time() - t0))
    return ok


if __name__ == "__main__":
    ok = run(QzField(), print) and run(make_fp(P89), print)
    print("ALL AIRY CHECKS PASSED" if ok else "AIRY CHECK FAILED")
