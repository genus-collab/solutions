# -*- coding: utf-8 -*-
"""ref12_n2rank.py -- rank of the N_2 system (2 genus-3 one-edge rows J3 and the 5 polar forms grad Q_Gamma(N_1)) on the 15 circulant
coordinates at Delta*, at points of Pi and at the isolated solution on c1_0 = 0.  Rank 7 = always solvable whatever the right-hand side."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_loc import *
for tag in ("61", "89"):
    S = load(tag); p = S.p; I = IRing(p)
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag))); A = [x % p for x in cs['Pi_A']]; B = [x % p for x in cs['Pi_B']]
    d = json.load(open(os.path.join(HERE, "ref9_newsol_%s.json" % tag)))
    pts = [("Delta*", d['delta']), ("Pi: A + 2 B", [(A[j] + 2 * B[j]) % p for j in range(15)]), ("Pi: A + 12345 B", [(A[j] + 12345 * B[j]) % p for j in range(15)]),
           ("isolated solution on c1_0 = 0", [x % p for x in cs['circ_b_solutions'][0]])]
    for name, v in pts:
        M = [list(S.J3[0]), list(S.J3[1])] + [S.quad_grad(i, v) for i in (3, 6, 7, 8, 17)]
        print("p = 2^%s-1  %-32s rank of the 7 x 15 N_2 system: %d;  rank of the five quadric gradients alone: %d; residuals all zero: %s" % (tag, name, rank_field(I, M), rank_field(I, M[2:]), not any(S.residuals(v).values())))
