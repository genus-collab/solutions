# -*- coding: utf-8 -*-
"""rf_run_d.py -- referee run D: the explicit formulas of XYSWAP_MODEL s.2.2 on the line x = z^4 + v/z (v = 4 and v = 1/3)
against the referee's own closed-formula code (rf_closed, itself checked against own TR at generic points in run A):
   omega_{1,1} = -D_x^2(h dz) = -5 z^4 (3 v^2 + 36 v z^5 + 8 z^10)/(4 z^5 - v)^4 dz,  h = z^2/2 + v/(12 z^3);
   omega_{2,1} = -D_x^4(k dz) + (1/2) D_x^5(h^2 dz) = [266112 z^34 - ... - 8 v^6 z^4]/(4 z^5 - v)^10 dz,  k = 1/80 + v/(80 z^5);
and the one-point correlators of the worker's leg table on the line by the residue at infinity (no TR needed):
   <omega_{g,1}, z^d> = -Res_inf z^d omega_{g,1} (omega regular at z = 0).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import rf_closed as CL

PRINTED21 = {34: 266112, 29: -10384448, 24: -20583920, 19: -7658960, 14: -690770, 9: -12053, 4: -8}   # z-exponent: coeff (v^((34-e)/5))
for v in (Fr(4), Fr(1, 3)):
    C = CL.setup({4: Fr(1), -1: v})
    h = ({5: Fr(1, 2), 0: v / 12}, 3, 0)                    # h = (z^5/2 + v/12)/z^3
    kk = ({5: Fr(1, 80), 0: v / 80}, 5, 0)
    hh = ({10: Fr(1, 4), 5: v / 12, 0: v * v / 144}, 6, 0)
    def Dn(F, n):
        for _ in range(n):
            F = CL.D1(C, F)
        return F
    A11 = Dn(h, 2)
    A4, A5 = Dn(kk, 4), Dn(hh, 5)
    ok11 = ok21 = okp11 = okp21 = True
    for z in (Fr(2), Fr(-3, 7), Fr(5, 3), Fr(1, 2)):
        cl11 = CL.omega_g1_value(C, 1, z)
        cl21 = CL.omega_g1_value(C, 2, z)
        ok11 &= (-CL.eval1(C, A11, z) == cl11)
        ok21 &= (-CL.eval1(C, A4, z) + CL.eval1(C, A5, z) / 2 == cl21)
        okp11 &= (-5 * z ** 4 * (3 * v * v + 36 * v * z ** 5 + 8 * z ** 10) / (4 * z ** 5 - v) ** 4 == cl11)
        num = sum((c * v ** ((34 - e) // 5) * z ** e for e, c in PRINTED21.items()), Fr(0))
        okp21 &= (num / (4 * z ** 5 - v) ** 10 == cl21)
    print("v = %s: omega_{1,1} = -D^2(h dz): %s; printed rational omega_{1,1}: %s; omega_{2,1} = -D^4(k dz) + D^5(h^2 dz)/2: %s;"
          " printed rational omega_{2,1}: %s" % (v, ok11, okp11, ok21, okp21))
# residue at infinity on the line, v symbolic through homogeneity: use v = 1 and the weight
v = Fr(1)
C = CL.setup({4: Fr(1), -1: v})
for g in (1, 2):
    terms = CL.omega_g1_terms(C, g)
    # expand omega^EO_{g,1} = -Sum D^r f_r at infinity: N/(z^A q^K) with q = 4 z^5 - 1
    out = {}
    for F, _ in terms:
        N, A, K = F
        # 1/q^K = (4 z^5)^-K (1 - 1/(4 z^5))^-K = Sum_j binom(K+j-1, j) 4^(-K-j) z^(-5K-5j)
        from math import comb
        for j in range(0, 12):
            cf = Fr(comb(K + j - 1, j), 4 ** (K + j))
            for e, c in N.items():
                ee = e - A - 5 * K - 5 * j
                out[ee] = out.get(ee, Fr(0)) - c * cf
    for d in (5, 10, 15, 20):
        val = out.get(-d - 1, Fr(0))                          # Sum_a Res_a z^d omega = -Res_inf = [z^(-d-1)] coefficient
        print("  line v=1: <omega_{%d,1}, z^%d> = %s" % (g, d, val))
