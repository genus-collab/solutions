# Fresh graph-by-graph, colour-by-colour polar table at g = 2, 3. Evidence only.
import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sys, time
from itertools import product
from fractions import Fraction as Fr
from gs import *

e1, e2 = Fr(sys.argv[1]), Fr(sys.argv[2])
g = int(sys.argv[3])
N = int(sys.argv[4]) if len(sys.argv) > 4 else 2 * g + 8
zmax = 3 * g + 6
t0 = time.time()
R = Rows(e1, e2, N, zmax)
V = Vert(e1, e2)
Ev = Evaluator(R, V)
gr = stable_graphs(g)
print('genus', g, 'graphs', len(gr), 'eps', e1, e2)
b = beta(2 * g - 2, e1, e2)
gap = fp((-3) ** (g - 1) * b)
tcm = R.b.tc.pow(2 - 2 * g)
total = None
vertexpolar = None
for genera, edges, aut in gr:
    gt = None
    rows = []
    for col in product(range(3), repeat=len(genera)):
        T = Ev.graph(genera, edges, col)
        pol = T.polar()
        rows.append((col, T.val if T.c else None, len(pol)))
        gt = T if gt is None else gt + T
    if len(edges) == 0:
        vertexpolar = gt
    tot = gt * inv(aut)
    total = tot if total is None else total + tot
    negs = [r for r in rows if r[1] is not None and r[1] < 0]
    print(genera, edges, 'aut', aut, 'val(total)', tot.val if tot.c else None,
          'neg colourings', negs if len(edges) else [r for r in rows], flush=True)
diff = total - tcm * gap
print('G_g - gap/t_c^(2g-2): val', diff.val, 'polar', diff.polar(), 'aprec', diff.aprec)
vd = vertexpolar - tcm * gap
print('vertex - gap/t^(2g-2): val', vd.val, 'polar', vd.polar())
# t_c^0 coefficient of the difference = value at xi = 0
c0 = diff.coef(0)
print('t_c^0 coefficient of the difference: mod p', c0, ' rat', ratrec(c0))
print('leading coefficient of G_g in t_c^(2-2g):', ratrec(total.coef(2 - 2 * g), 2 ** 62), ' gap:', (-3) ** (g - 1) * b)
print('time', time.time() - t0)
