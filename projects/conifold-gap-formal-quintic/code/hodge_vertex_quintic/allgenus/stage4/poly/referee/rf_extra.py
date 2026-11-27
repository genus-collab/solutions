# -*- coding: utf-8 -*-
"""rf_extra.py -- REFEREE: side checks with literature and closed forms (own code).
(E1) h_n exactly over Q (Q(zeta_5) traces) against the mod-p implementation; h_5 = -257071/8064000000 (POLY s.1).
(E2) arXiv:0710.0889 (Zagier-Zinger), Thm 3 and the formulas for Phi_2, Phi_3 after Thm 4, n = 5: their F(w,x) is the colour-0 series
     with w = 1/z, x = q, and their expansion (5) is coefficientwise in x with Phi_s(0) = delta_(s,0).  Hence Phi^ZZ = Phi^QRR / h(z):
     A^ZZ_n = Sum_k [z^k](1/h) A^QRR_(n-k).  Checked as q-series mod p for s = 0..3 (a published check of rf_qexp's large-radius WKB).
(E3) closed form of the free constants (observed by the parallel Stage-4 'gamma' worker; re-checked here from rf_qexp's polynomials):
     F_inf(z) = Sum_m f_m z^(5m) = exp( Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l) / (l (5l+1)) ),  f_1 = B_6(1/5)/6 = 4537/3937500.
Run: py -I -B rf_extra.py M      (E3 for m <= M; needs rf_qexp to n = 5M)
"""
import os, sys, time
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_core import P61, P89, md, Ser, bernoulli, hcoef_mod, hcoef_exact
from rf_qexp import wkb_q, hatted_polys


def bern_poly(n, x, B):
    Bm = list(B); Bm[1] = Fr(-1, 2)                       # standard convention for Bernoulli polynomials
    return sum(comb(n, k) * Bm[k] * x ** (n - k) for k in range(n + 1))


def zz_phi(n=5):
    """Zagier-Zinger Phi_0..Phi_3 as dicts {power of L: coefficient} (arXiv:0710.0889 (6) and the display after Thm 4)"""
    P0 = {1: Fr(1)}
    c1 = Fr((n - 2) * (n + 1), 24 * n); P1 = {1: c1, n: -c1}
    c2 = Fr((n + 1) ** 2 * (n - 2) ** 2, 2 * (24 * n) ** 2); P2 = {1: c2, n: -2 * c2, 2 * n - 1: c2}
    c3 = Fr((n + 1) * (n - 2), 30 * (24 * n) ** 3); P3 = {}
    def add(k, v): P3[k] = P3.get(k, 0) + c3 * v
    add(3 * n - 2, 1003 * n ** 4 - 2366 * n ** 3 + 3759 * n ** 2 - 1676 * n - 164)
    add(2 * n - 2, -72 * (n - 1) * (3 * n - 1) * (7 * n ** 2 - 9 * n + 14))
    add(2 * n - 1, 15 * (n + 1) ** 2 * (n - 2) ** 2); add(n, -15 * (n + 1) ** 2 * (n - 2) ** 2)
    add(n - 2, 72 * (n - 1) * (7 * n ** 3 - 17 * n ** 2 + 22 * n - 24))
    add(1, 5 * n ** 4 + 134 * n ** 3 - 447 * n ** 2 + 308 * n - 556)
    return [P0, P1, P2, P3]


def main(M):
    t0 = time.time(); ok = True
    he = hcoef_exact(12)
    e1 = he[5] == Fr(-257071, 8064000000) and all([md(x, p) for x in he] == hcoef_mod(p, 12) for p in (P61, P89))
    print("(E1) h_n exact (Q(zeta_5) traces) = h_n mod p for n <= 12, both primes, and h_5 = -257071/8064000000: %s;  h_1 = %s" % (e1, he[1]))
    ok = ok and e1
    for p in (P61, P89):
        S, A, L = wkb_q(p, 3); h = hcoef_mod(p, 3); K = S.K
        g = S.inv(h + [0] * (K - 4))                                     # 1/h(z) as a series in z (only 4 terms used)
        Lp = {}
        def Lpow(k):
            if k not in Lp: Lp[k] = S.powr([1, (-3125) % p] + [0] * (K - 2), Fr(-k, 5))
            return Lp[k]
        zz = zz_phi(5); res = []
        for s in range(4):
            Azz = S.z()
            for k in range(s + 1): Azz = S.add(Azz, A[s - k], g[k])
            ref = S.z()
            for e, c in zz[s].items(): ref = S.add(ref, Lpow(e), md(c, p))
            res.append(Azz == ref)
        print("(E2) p=%d: Zagier-Zinger Phi_0..Phi_3 (n = 5) = (1/h) x QRR WKB, as q-series mod q^%d: %s" % (p, K, res))
        ok = ok and all(res)
    B = bernoulli(5 * M + 2)
    logG = [Fr(0)] * (M + 1)
    for l in range(1, M + 1): logG[l] = (-1) ** (l + 1) * bern_poly(5 * l + 1, Fr(1, 5), B) / (l * (5 * l + 1))
    print("     B_6(1/5)/6 = %s" % (bern_poly(6, Fr(1, 5), B) / 6))
    for p in (P61, P89):
        polys, _ = hatted_polys(p, 5 * M)
        f = [1] + [polys[5 * m][4 * m] if len(polys[5 * m]) > 4 * m else 0 for m in range(1, M + 1)]
        S = Ser(p, M + 1); Gi = S.exp0([0] + [md(x, p) for x in logG[1:]])
        print("(E3) p=%d: F_inf = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1))) for m <= %d: %s" % (p, M, f == Gi))
        ok = ok and f == Gi
    print("RESULT rf_extra: %s  [%.1f s]" % ("PASS" if ok else "FAIL", time.time() - t0))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 10)
