# Single-vertex graph at genus g: full polar part against (-3)^(g-1) beta_{2g-2}/t_c^(2g-2). Evidence only.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import time
from fractions import Fraction as Fr
from gs import *

e1, e2 = Fr(sys.argv[1]), Fr(sys.argv[2])
for g in range(int(sys.argv[3]), int(sys.argv[4]) + 1):
    t0 = time.time()
    N = 2 * g + 4
    zmax = 3 * g + 3
    R = Rows(e1, e2, N, zmax)
    V = Vert(e1, e2)
    Ev = Evaluator(R, V)
    b = beta(2 * g - 2, e1, e2)
    gap = fp((-3) ** (g - 1) * b)
    tcm = R.b.tc.pow(2 - 2 * g)
    for col in (0, 1, 2):
        T = Ev.vertex(g, col, ())
        if col == 0:
            d = T - tcm * gap
            print('g', g, 'colour 0: val T', T.val, ' polar(T-gap/t^(2g-2))', d.polar(), ' val diff', d.val, 'aprec', d.aprec)
        else:
            print('g', g, 'colour', col, 'val', T.val, 'polar', T.polar())
    print('  time', round(time.time() - t0, 1), flush=True)
