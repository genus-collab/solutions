# -*- coding: utf-8 -*-
"""p8_g2.py -- referee PASSES: the genus-2 statement of genus3_u3.py is only in the logs; re-evaluate the seven genus-2 graphs from the
STORED (n1, n2, n3) of each finished class: edge graphs 2..7 must be O(u^2) (a genus-2 term is u^-2 G), vertex [u^0] = 1/1200, total =
1/1200 - 7/6000 u + 97/72000 u^2 - 6809/18000000 u^3.  Also with N_1 alone (what do N_2, N_3 change at genus 2?).  Run: py -I p8_g2.py"""
import os, sys, json
from fractions import Fraction as Fr
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
from fq_fp import FpRing, P1, P2
from fq3_steps import frame_of
from fq3_par import setup
from fq3_m import sp_mats, MU
from fq3_msum import GraphSummer6, UVertex
from fq_msum import outer_graphs
TOT = (Fr(1, 1200), Fr(-7, 6000), Fr(97, 72000), Fr(-6809, 18000000))
for tag in ("61", "89"):
    p = P2 if tag == "89" else P1; R = FpRing(p); st = setup(frame_of(p, (0, 0), 34, 6))
    for cls in ("circ", "circpi", "circb", "circstar"):
        d = json.load(open(os.path.join(PROJ, "g3_u3_%s_%s.json" % (tag, cls))))
        have = [k for k in ('n1', 'n2', 'n3') if k in d]
        for use in sorted(set([1, len(have)])):
            Ns = [sp_mats(d[k], p) for k in have[:use]]
            MUo = MU(R, st['M0'], Ns, 6, 'right'); G2 = GraphSummer6(st['H'], MUo, UVertex(R, MUo, 6, 21)).all_graphs(outer_graphs())
            edge_ok = all(G2[i].coef(k) == 0 for i in range(2, 8) for k in range(2)); tot = None
            for i in G2: tot = G2[i] if tot is None else tot.add(G2[i])
            print("p=2^%s-1 %-8s with %s: %d graphs; edge graphs O(u^2): %s; vertex [u^0] = 1/1200: %s; vertex [u^1] = -7/6000: %s; total u^0..u^3 as stated: %s; edge graphs non-zero at u^2: %s"
                  % (tag, cls, "+".join(have[:use]), len(G2), edge_ok, G2[1].coef(0) == R.of(Fr(1, 1200)), G2[1].coef(1) == R.of(Fr(-7, 6000)),
                     [tot.coef(k) == R.of(TOT[k]) for k in range(4)], [i for i in range(2, 8) if G2[i].coef(2)]), flush=True)
