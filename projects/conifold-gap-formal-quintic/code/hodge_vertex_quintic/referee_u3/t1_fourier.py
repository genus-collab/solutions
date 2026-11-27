# -*- coding: utf-8 -*-
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
from rlib import *
from fq3_rec import Rec
tag = sys.argv[1] if len(sys.argv) > 1 else "61"
S, d, cs = load(tag); p = S.p; rec = Rec(S.R)
M = d['M']; Mi = d['Minv']; mu0 = d['mu0']
A = cs['Pi_A']; B = cs['Pi_B']
print("residuals A, B, A+3B all zero:", [not any(S.residuals(v).values()) for v in (A, B, [(a + 3 * b) % p for a, b in zip(A, B)])])
NA = ser_from_circ(S, A); NB = ser_from_circ(S, B)
show(S, "A (right-placed Delta)", NA, rec); show(S, "B (right-placed Delta)", NB, rec)
DA = Ad(M, Mi, NA, p); DB = Ad(M, Mi, NB, p)
show(S, "D_A = M0 A M0^-1 (left-placed)", DA, rec); show(S, "D_B = M0 B M0^-1", DB, rec)
show(S, "M0^-1 A M0", Ad(Mi, M, NA, p), rec); show(S, "M0^-1 B M0", Ad(Mi, M, NB, p), rec)
Nt = sp_mats(S.nt, p); show(S, "N1taut", Nt, rec); show(S, "mu0", mu0, rec)
show(S, "Ad_M0 N1taut", Ad(M, Mi, Nt, p), rec)
