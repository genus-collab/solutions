# -*- coding: utf-8 -*-
"""rf_run_c.py -- referee run C (exact over Q, own TR): leg dictionary and sign checks at GENERIC points of the family.
 (a) dilaton leg -Phi:  <omega_{1,1}, -Phi> = 5/24 (= rank/24) and <omega_{2,1}, -Phi> = 2 F_2;
 (b) the sign of the TR string equation: <omega_{g,2}, -z, Psi> against <omega_{g,1}, +- dPsi/dx>, Psi = Int z^m dx;
 (c) the string reconciliation with tau_1(H) = -(x log z - Int x dz/z) and the H leg -log z;
 (d) Euler-dilaton leg on the whole family, e = -Phi_pol (log-free part of -Phi): <e>_1, <e, H>_1, <e, H>_2;
 (e) the leg 1/z in multi-point correlators along the degenerating path of run B (regular or singular in S = s^2?);
 (f) the descendant chain Int^k (-(4/5) z^5) dx on the line: coefficients and the first log z.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from rf_tr import TR
from rf_fam import family_x, phi_parts, leg_lp, leg_log, leg_sum, leg_scale
from rf_ser import lp_mul, lp_deriv

t0 = time.time()
ORDER = [(0, 3), (1, 1), (0, 4), (1, 2), (2, 1), (1, 3), (2, 2)]


def int_dx(L, X):
    """Int L(z) dx for a Laurent polynomial L: returns (Laurent polynomial, coefficient of log z)"""
    P = lp_mul(L, lp_deriv(X))
    out, lg = {}, Fr(0)
    for e, c in P.items():
        if e == -1:
            lg += c
        else:
            out[e + 1] = out.get(e + 1, Fr(0)) + c / (e + 1)
    return out, lg


for r in ([1, -2, 3, -6, Fr(-3, 2)], [1, 2, Fr(-2, 3), Fr(1, 3), Fr(-1, 3)]):
    X, (c3, c2, c1, V) = family_x(r)
    print("=== point r = %s, (c3, c2, c1, V) = (%s, %s, %s, %s)" % ([str(a) for a in r], c3, c2, c1, V), flush=True)
    tr = TR(X, r, 11)
    for gn in ORDER:
        tr.step(*gn)
    pol, lg = phi_parts(X)
    mphi = leg_sum(leg_lp({e: -c for e, c in pol.items()}), leg_log({0: -lg}))
    F2 = tr.pair(2, 1, [leg_sum(leg_lp(pol), leg_log({0: lg}))]) / (2 - 4)
    print("  (a) <omega_{1,1}, -Phi> = %s (5/24: %s);  <omega_{2,1}, -Phi> = %s, 2 F_2 = %s: %s"
          % (tr.pair(1, 1, [mphi]), tr.pair(1, 1, [mphi]) == Fr(5, 24), tr.pair(2, 1, [mphi]), 2 * F2, tr.pair(2, 1, [mphi]) == 2 * F2))
    unit = leg_lp({1: Fr(-1)})
    for g, m in ((1, 5), (1, 7), (1, 10), (2, 15), (2, 16)):
        Psi, l = int_dx({m: Fr(1)}, X)
        assert l == 0
        lhs = tr.pair(g, 2, [unit, leg_lp(Psi)])
        rhs = tr.pair(g, 1, [leg_lp({m: Fr(1)})])
        print("  (b) g=%d Psi = Int z^%d dx: <omega_{g,2}, -z, Psi> = %s; <omega_{g,1}, z^%d> = %s; unit acts by +d/dx: %s, by -d/dx: %s"
              % (g, m, lhs, m, rhs, lhs == rhs, lhs == -rhs))
    xdz = {e - 1: c for e, c in X.items()}                  # x/z
    Ixz = {}
    lgx = Fr(0)
    for e, c in xdz.items():
        if e == -1:
            lgx += c
        else:
            Ixz[e + 1] = c / (e + 1)
    assert lgx == 0                                          # Int x dz/z has no log (x has no z^0 term... c0 = 0 here)
    tau1H = leg_sum(leg_log({e: -c for e, c in X.items()}), leg_lp(Ixz))   # -(x log z - Int x dz/z)
    Hleg = leg_log({0: Fr(-1)})
    for g in (1, 2):
        a = tr.pair(g, 2, [unit, tau1H])
        b = tr.pair(g, 1, [Hleg])
        print("  (c) g=%d: <omega_{g,2}, tau_0(1), tau_1(H)> = %s; <omega_{g,1}, tau_0(H)> = %s; equal %s" % (g, a, b, a == b))
    e = leg_lp({k: -c for k, c in pol.items()})
    print("  (d) e = -Phi_pol: <omega_{1,1}, e> = %s (1/8: %s); <omega_{1,2}, e, H> = %s; <omega_{2,2}, e, H> = %s"
          % (tr.pair(1, 1, [e]), tr.pair(1, 1, [e]) == Fr(1, 8), tr.pair(1, 2, [e, Hleg]), tr.pair(2, 2, [e, Hleg])))
    print("  (%.1fs)" % (time.time() - t0), flush=True)

print("=== (e) the leg 1/z in multi-point correlators along r(s) = (1, 2, -2/3, s, -s)")
vals = {}
SV = [Fr(1, 2), Fr(1, 3), Fr(1, 4), Fr(1, 5), Fr(3, 4), Fr(2, 5), Fr(3, 5), Fr(1, 6), Fr(5, 6), Fr(4, 7)]
for s in SV:
    r = [1, 2, Fr(-2, 3), s, -s]
    X, par = family_x(r)
    tr = TR(X, r, 9)
    for gn in ORDER[:5]:
        tr.step(*gn)
    iz = leg_lp({-1: Fr(1)})
    tests = {"(0,3) 1/z,z,z^4": [iz, leg_lp({1: Fr(1)}), leg_lp({4: Fr(1)})],
             "(0,3) 1/z,z^2,z^3": [iz, leg_lp({2: Fr(1)}), leg_lp({3: Fr(1)})],
             "(0,3) 1/z,1/z,z^6": [iz, iz, leg_lp({6: Fr(1)})],
             "(0,3) 1/z,1/z,1/z": [iz, iz, iz],
             "(0,3) v/z,v/z,z^6": [leg_lp({-1: par[3]}), leg_lp({-1: par[3]}), leg_lp({6: Fr(1)})],
             "(1,2) 1/z,z^6": [iz, leg_lp({6: Fr(1)})], "(1,2) 1/z,1/z": [iz, iz], "(1,2) 1/z,z^11": [iz, leg_lp({11: Fr(1)})]}
    for k, legs in tests.items():
        g, n = (0, 3) if k.startswith("(0,3)") else (1, 2)
        vals.setdefault(k, []).append((s * s, tr.pair(g, n, legs)))
for k, lst in vals.items():
    S0, v0 = lst[0]
    pole = None
    for kk in range(0, 6):
        ys = [v * S ** kk for S, v in lst]
        xs = [S for S, v in lst]
        n = 5
        def P(t):
            tot = Fr(0)
            for i in range(n):
                li = Fr(1)
                for j in range(n):
                    if j != i:
                        li *= (t - xs[j]) / (xs[i] - xs[j])
                tot += ys[i] * li
            return tot
        if all(P(xs[i]) == ys[i] for i in range(n, len(xs))):
            pole = (kk, P(Fr(0)))
            break
    print("  %-22s values zero: %s; smallest k with S^k*value polynomial (deg<=4): %s" % (k, all(v == 0 for _, v in lst), pole))

print("=== (f) descendant chain Int^k(-(4/5) z^5) dx on the line x = z^4 + v/z")
for v in (Fr(4), Fr(1, 3)):
    X = {4: Fr(1), -1: v}
    L = {5: Fr(-4, 5)}
    for k in range(1, 6):
        L, l = int_dx(L, X)
        print("  v=%s k=%d: %s ; log z coefficient %s" % (v, k, {e: str(c) for e, c in sorted(L.items())}, l))
        if l:
            break
print("done %.1fs" % (time.time() - t0))
