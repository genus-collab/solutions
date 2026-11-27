# -*- coding: utf-8 -*-
"""p5_rows.py -- referee PASSES: the 7 x 15 matrix of the N_2 system (J3 rows and polar forms of the five two-edge quadrics) at stored
N_1-solutions, from the stored forms only (no graph evaluation): rank, identically vanishing rows, and the value of the linear form m_1.
Run: py -I p5_rows.py"""
import os, sys, json, random
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fq3_circsys import CircSys
from p3_n2cons import rank
TWO = (3, 6, 7, 8, 17)
for tag in ("61", "89"):
    S = CircSys(tag); p = S.p; cj = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag)))
    A = [x % p for x in cj['Pi_A']]; B = [x % p for x in cj['Pi_B']]; rng = random.Random(7 + int(tag))
    pts = {'pi_A': A, 'pi_B': B, 'pi_A+2B': [(a + 2 * b) % p for a, b in zip(A, B)], 'star': cj['circ_star'], 'circb': cj['circ_b_solutions'][0]}
    for k in range(3):
        x, y = rng.randrange(p), rng.randrange(p); pts['pi_rand%d' % k] = [(x * a + y * b) % p for a, b in zip(A, B)]
    for k, v in enumerate(cj['circ_b_extra']): pts['extra%d' % k] = v
    for nm, v in pts.items():
        v = [t % p for t in v]; bad = [k for k, r in S.residuals(v).items() if r]
        rows = [S.J3[0], S.J3[1]] + [S.quad_grad(i, v) for i in TWO]
        zero = [i for i, r in zip(TWO, rows[2:]) if not any(r)]
        m1 = sum(S.m1[j] * v[j] for j in range(15)) % p
        print("p=2^%s-1 %-9s: N_1-conditions violated %s; rank of the 7 x 15 N_2 matrix %d; rank of the 5 polar rows %d; zero rows: graphs %s; m_1 = %s; c1_0 = %s"
              % (tag, nm, bad, rank(rows, p), rank(rows[2:], p), zero, "0" if m1 == 0 else "non-zero", v[0] if v[0] < 10 else "non-zero"))
