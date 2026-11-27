# -*- coding: utf-8 -*-
"""rc_t0_pf.py -- smoke test of the independent Picard-Fuchs / t_c implementation."""
import os
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rc_pf
import rc_rec
from fractions import Fraction as Fr

tcc, diag = rc_pf.tc_coeffs(0, 0, kmax=12)
print("gamma  =", diag['gamma'])
for k in ('f0', 'f1', 'f2', 'g', 'h', 'varpi0'):
    print("%-7s" % k, diag[k])
print("tc  u^1..u^12:", tcc)
print()
tcc2, diag2 = rc_pf.tc_coeffs(Fr(1, 3), Fr(-2, 7), kmax=10)
print("second frame (alpha,beta) = (1/3,-2/7)")
print("varpi0 ", diag2['varpi0'])
print("tc  u^1..u^10:", tcc2)
print()
print("B_10 =", rc_rec.bern(10), "  |B_10|/(10*8*5^4) =", abs(rc_rec.bern(10)) / (10 * 8 * 5 ** 4))
print("1/660000 =", Fr(1, 660000))
print()
# sanity: the identity change of variable (t_c = u) must be the identity on the polar part
print("(u/t_c)^8 as a u-series, frame (0,0):", rc_pf.u_over_tc_pow(8, 0, 0, 8))
