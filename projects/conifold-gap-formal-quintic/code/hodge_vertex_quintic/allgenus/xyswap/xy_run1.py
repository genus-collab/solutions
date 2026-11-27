# -*- coding: utf-8 -*-
"""xy_run1.py -- first checks: direct TR (xy_tr) against the closed formula (xy_closed) for n = 1; F_g by the dilaton
equation from direct TR; the vanishing-period integral Phi_g(0) - Phi_g(inf) from the closed formula against -B_2g/(2g).
Usage: py -I -B xy_run1.py [qz|fp89|fp61] [gmax_tr] [gmax_closed]
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from xy_field import QZ, make_fp, P89, P61
from xy_tr import TR
from xy_closed import omega_g1, omega_g1_terms, ev_rat, bernoulli


def main():
    fld = sys.argv[1] if len(sys.argv) > 1 else "qz"
    gtr = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    gcl = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    F = QZ if fld == "qz" else make_fp(P89 if fld == "fp89" else P61)
    print("field", F.name, "gmax(TR)", gtr); t0 = time.time()
    order = [(0, 3), (1, 1)]
    if gtr >= 2: order += [(0, 4), (1, 2), (2, 1)]
    if gtr >= 3: order += [(0, 5), (1, 3), (2, 2), (3, 1)]
    tr = TR(F, gtr, max(n for (_, n) in order))
    for (g, n) in order:
        tr.step(g, n)
        print("  TR omega_{%d,%d}: %d coefficients, symmetric failures %d, %.1fs" % (
            g, n, len(tr.om[(g, n)]), tr.symmetric_failures(g, n) if n <= 3 else -1, time.time() - t0))
    v = 4; pts = [Fr(2), Fr(-3, 7), Fr(5, 3), Fr(1, 2)]
    for g in range(1, gtr + 1):
        N, K = omega_g1(g); ok = True
        for z in pts:
            a = ev_rat(N, K, z, v); b = tr.value(g, 1, (z,))
            same = (b - F.frac(a)).is_zero()
            ok = ok and same
            if fld == "qz": print("  omega_{%d,1} at p = %s: closed %s, TR %s, rational %s" % (g, z, a, b.c[0], b.rational() is not None))
        print("CHECK omega_{%d,1}: closed formula (Thm 7.1, EO sign) = direct TR at %d points: %s" % (g, len(pts), ok))
    B = bernoulli(2 * gtr + 2)
    for g in range(2, gtr + 1):
        Fg = tr.free_energy(g)
        expect = B[2 * g] / (2 * g * (2 * g - 2)) * Fr(v) ** (2 - 2 * g)
        print("CHECK F_%d (dilaton, direct TR, v=4) = B_2g/(2g(2g-2)) v^(2-2g) = %s: %s" % (g, expect, (Fg - F.frac(expect)).is_zero()))
    print("--- vanishing-period integral from the closed formula (ABDKS sign), symbolic v ---")
    B = bernoulli(2 * gcl + 2)
    for g in range(1, gcl + 1):
        terms = omega_g1_terms(g)
        reg0 = True; inf_ok = True; val0 = {}
        for (r, N, k, prim) in terms:
            if prim is None: continue
            Np, kp = prim
            if any(ez < 0 for (ez, ev) in Np): reg0 = False
            if max(ez for (ez, ev) in Np) - 5 * kp >= 0: inf_ok = False
            for (ez, ev), c in Np.items():
                if ez == 0: val0[ev - kp] = val0.get(ev - kp, 0) + c * Fr(-1) ** kp
        val0 = {e: c for e, c in val0.items() if c}
        expect = {1 - 2 * g: -B[2 * g] / (2 * g)}
        print("g=%d: each T_r regular at z=0: %s; each T_r -> 0 at z=inf: %s; Phi_g(0) = %s (v-exponent: coeff); "
              "-B_2g/(2g) v^(1-2g): %s; MATCH %s; #terms %d; %.1fs" % (
                  g, reg0, inf_ok, val0, expect, val0 == expect, len(terms), time.time() - t0))


if __name__ == "__main__":
    main()
