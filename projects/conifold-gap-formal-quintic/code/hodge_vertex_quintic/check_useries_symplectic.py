# -*- coding: utf-8 -*-
"""Check (report stage): every u-series coefficient (s^0, s^5, s^10) of every colour-summed genus-2 graph is independent
of the constant symplectic M (referee's strengthening of the s^-10, s^-5 statement), at random colour-mixing M."""
import os, sys, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_rows import Rows
from fq_conifold import Conifold
from fq_sser import SRows, SympM
from fq_msum import GraphSummer, outer_graphs

for p in (P1, P2):
    R = FpRing(p); rng = random.Random(12345 + p % 1000)
    t0 = time.time()
    rows = Rows(R, Nz=6); C = Conifold(R, Nu=10); S = SRows(rows, C, prec_abs=24); graphs = outer_graphs()
    res1 = GraphSummer(S, SympM.identity(R)).all_graphs(graphs)
    ok = True; nz_other = 0
    for trial in range(3):
        M = SympM.random(R, rng); resM = GraphSummer(S, M).all_graphs(graphs)
        for idx in resM:
            prec = min(resM[idx].prec, res1[idx].prec)
            for k in range(-10, prec):
                a, b = resM[idx].coef(k), res1[idx].coef(k)
                if k % 5 == 0:
                    if a != b: ok = False; print("DIFF graph %d s^%d" % (idx, k))
                elif a != b: nz_other += 1
    print("prime %d: u-series coefficients (s^-10, s^-5, s^0, s^5, s^10) M-independent for all 7 graphs, 3 random M: %s;"
          " non-u-series coefficients differing from the point frame: %d; %.1f s" % (p, ok, nz_other, time.time() - t0))
