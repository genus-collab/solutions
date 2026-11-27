# -*- coding: utf-8 -*-
"""ref_alg2.py -- CircSys audit: the restricted quadrics are homogeneous in Delta (c0 = 0, b = 0, dropped silently by CircSys);
S.quad agrees with the 75-variable QuadForm at N_1 = to_sp(v); rank of the high fibre at Pi_B and at a generic point of Pi."""
import os, sys, json, random
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); sys.path.insert(0, PROJ)
from fq3_circsys import CircSys
from fq3_quad import QuadForm, IDX_2EDGE, rank_fp
for tag in ("61", "89"):
    S = CircSys(tag); p = S.p; rng = random.Random(5)
    qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag))); cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag)))
    ok_h = True; ok_e = True
    for i in IDX_2EDGE:
        Q = QuadForm(p, qs['Q'][str(i)]['c0'], qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']); r = Q.restrict(S.nt, S.cols)
        ok_h = ok_h and r.c0 == 0 and not any(r.b)
        v = [rng.randrange(p) for _ in range(15)]; n = S.to_sp(v)
        # own evaluation of the 75-variable quadratic polynomial
        own = (Q.c0 + sum(Q.b[k] * n[k] for k in range(75)) + sum(n[k] * Q.A[k][l] * n[l] for k in range(75) for l in range(75))) % p
        ok_e = ok_e and own == S.quad(i, v)
    print(tag, "restricted quadrics homogeneous (c0 = 0, b = 0):", ok_h, "| S.quad = 75-variable Q at to_sp(v):", ok_e)
    A = cs['Pi_A']; B = cs['Pi_B']
    def fibre_rank(v):
        """rank of the linear system for the high part (c4, c5, c6) at fixed low part: J3 rows, quadric gradients, cubic c4-derivatives"""
        low = list(v[:8]); v0 = low + [0] * 7; rows = []
        for r in range(2): rows.append([S.J3[r][j] for j in range(8, 15)])
        for i in IDX_2EDGE:
            g = S.quad_grad(i, v0); rows.append([g[j] for j in range(8, 15)])
        for g_ in (4, 9, 10, 11, 12, 18, 19, 20, 32):
            b0 = S.cub(g_, v0); rows.append([(S.cub(g_, low + [1, 0]) - b0) % p, (S.cub(g_, low + [0, 1]) - b0) % p] + [0] * 5)
        return rank_fp(rows, p)
    print(tag, "rank of the high-part system (7 unknowns): at Pi_A", fibre_rank(A), "at Pi_B", fibre_rank(B), "at A+2B", fibre_rank([(A[j] + 2 * B[j]) % p for j in range(15)]),
          "at circ_b_sol0", fibre_rank(cs['circ_b_solutions'][0]), "at extras", [fibre_rank(e) for e in cs['circ_b_extra']])
