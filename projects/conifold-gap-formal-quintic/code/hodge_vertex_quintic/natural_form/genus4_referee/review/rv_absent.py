# -*- coding: utf-8 -*-
"""rv_absent.py -- why do 81 outer graphs with 3..5 edges receive no non-zero piece at all?  (read-only probe)
Lists the absent outer graphs, and for the absent 3-edge graphs evaluates the piece (inner graph = the graph itself, S = all edges)
with (i) the true E_out, (ii) a generic (random) replacement of E_out with the same valuation: if (ii) is non-zero the vanishing is a
property of E_out (a selection rule), not a dead code path."""
import os, sys, json, random, itertools
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_common import PROJ, OUT
from fractions import Fraction as Fr
from fq_fp import P2
from fq4_split import Split
sp = Split(P2, 4, (Fr(0), Fr(0)), 6, 'nat'); p = sp.p; K = sp.K
Eof = [G['E'] for G in sp.inner_graphs]
rows = set(int(k) for k in json.load(open(os.path.join(PROJ, "g4_split_g4_89_nat.json")))['table'])
absent = [i for i in range(len(Eof)) if Eof[i] <= 5 and i not in rows]
pres = [i for i in range(len(Eof)) if 1 <= Eof[i] <= 5 and i in rows]
def desc(i):
    G = sp.inner_graphs[i]; return "%d: gv=%s ed=%s" % (i, G['gv'], G['ed'])
print("absent outer graphs with E <= 5: %d; by E %s" % (len(absent), {E: sum(1 for i in absent if Eof[i] == E) for E in range(1, 6)}))
for i in absent:
    if Eof[i] <= 4: print("   absent  ", desc(i))
for i in pres:
    if Eof[i] == 3: print("   present ", desc(i))
# classification attempt: number of genus-0 vertices / vertices
def feat(i):
    G = sp.inner_graphs[i]; return (Eof[i], sum(1 for g in G['gv'] if g == 0), tuple(sorted((G['gv'][v], len(G['hsAt'][v])) for v in range(G['V']))))
print("features (E, #genus-0 vertices) absent :", sorted(set(feat(i)[:2] for i in absent)))
print("features (E, #genus-0 vertices) present:", sorted(set(feat(i)[:2] for i in pres)))
# [u^1] E_out(a, b) in the frozen flat frame: rank
def rank_mod(Mx):
    Mx = [r[:] for r in Mx]; rk = 0
    for c in range(5):
        piv = next((r for r in range(rk, 5) if Mx[r][c] % p), None)
        if piv is None: continue
        Mx[rk], Mx[piv] = Mx[piv], Mx[rk]; inv = pow(Mx[rk][c], p - 2, p)
        for r in range(5):
            if r != rk and Mx[r][c] % p:
                f = Mx[r][c] * inv % p; Mx[r] = [(x - f * y) % p for x, y in zip(Mx[r], Mx[rk])]
        rk += 1
    return rk
print("rank of [u^1] E_out(a,b) as a 5x5 colour matrix:", {key: rank_mod([[(sp.EO[key][i][j] or [0] * K)[1] for j in range(5)] for i in range(5)]) for key in [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0)]})
print("rank of [u^2] E_out(0,0):", rank_mod([[(sp.EO[(0, 0)][i][j] or [0] * K)[2] for j in range(5)] for i in range(5)]))
# generic replacement
rng = random.Random(5); EOsave = sp.EO
def randser(): return [0] + [rng.randrange(p) for _ in range(K - 1)]
EOrand = {key: [[randser() for j in range(5)] for i in range(5)] for key in sp.EO}
for i in [j for j in absent if Eof[j] == 3][:3]:
    sp.EO = EOsave; sp.EM = {k: [[None] * 5 for _ in range(5)] for k in sp.ER}
    # true E_out on all edges: use variant-id path (S = all edges) by swapping tables
    sp.variant = 'id'; ER = sp.ER; sp.EO = EOsave; r1 = sp.graph(i)
    sp.EO = EOrand; r2 = sp.graph(i); sp.EO = EOsave
    print("   outer graph %s: piece with the true E_out: %s; with a generic E_out of valuation 1: %s" % (desc(i), r1.get(i), "non-zero %s" % [k for k in range(K) if r2[i][k]] if i in r2 else None))
