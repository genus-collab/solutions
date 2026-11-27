# -*- coding: utf-8 -*-
"""rf_s15_rat.py -- display only: CRT (2^61-1, 2^89-1) + rational reconstruction of the own genus-3 tables (frame (0,0)); the shifts
caused by the stored normalisation (control X5); the count of non-zero regular coefficients.  Run: py -I -B rf_s15_rat.py"""
import os, sys, json, math
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_s15_base import P61, P89
from rf_s15_cmp import MINE, WORK
M = P61 * P89


def crt(a, b): return (a * P89 * pow(P89, -1, P61) + b * P61 * pow(P61, -1, P89)) % M


def rr(a):
    r0, r1, s0, s1 = M, a % M, 0, 1; B = math.isqrt(M // 2)
    while r1 > B: q = r0 // r1; r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B: return "?"
    return str(Fr(r1, s1))


def L(f): return {int(k): v for k, v in json.load(open(os.path.join(MINE, f)))['A'].items()}


c61, c89, s61, s89 = L('rf_g3_61_can.json'), L('rf_g3_89_can.json'), L('rf_g3_61_stored.json'), L('rf_g3_89_stored.json')
print("own flat, frame (0,0), CRT display (k = coefficient of u^(k-4)):")
print("  vertex k=0..5:", [rr(crt(c61[1][k], c89[1][k])) for k in range(6)])
for i in (2, 3, 4, 5, 9, 39):
    print("  graph %2d k=4..7:" % i, [rr(crt(c61[i][k], c89[i][k])) for k in range(4, 8)])
print("X5 shifts (stored normalisation minus canonical), non-zero graphs:")
for i in range(1, 43):
    d = [crt(s61[i][k] - c61[i][k], s89[i][k] - c89[i][k]) for k in range(8)]
    if any(d): print("  graph %2d:" % i, [rr(x) if x else "0" for x in d])
nz_e = sum(1 for i in c89 if i != 1 for k in range(4, 8) if c89[i][k])
nz_a = sum(1 for i in c89 for k in range(4, 8) if c89[i][k])
print("non-zero coefficients at k = 4..7: edge graphs %d, all 42 graphs %d" % (nz_e, nz_a))
H = json.load(open(os.path.join(WORK, "hat_g3_89_nat_K8.json")))['A']
print("worker hatted table, same count: edge graphs %d, all graphs %d"
      % (sum(1 for i, v in H.items() if i != '1' for k in range(4, 8) if v[k]), sum(1 for i, v in H.items() for k in range(4, 8) if v[k])))
