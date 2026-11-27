# -*- coding: utf-8 -*-
"""rv_mass.py -- numerical test of the regrouping identity (REVIEW.md (a)) with the weights t^(number of inner edges).

Claim tested, for every leg-free stable graph Gamma of genus g (outer graph):
   Sum over pairs (Gt, S), Gt a leg-free stable graph of genus g, S a labelled subset of its edges with Gt/(E(Gt) - S) isomorphic to Gamma,
        of  t^(|E(Gt)| - |S|) / |Aut Gt|
   =  (1/|Aut Gamma|)  prod_{v in V(Gamma)}  N_{g_v, n_v}(t),      N_{g,n}(t) := Sum_{gamma} t^|E(gamma)| / |Aut gamma|,
where gamma runs over the stable graphs of type (g, n) with n LABELLED legs (fq3_graphs.legged_graphs) and Aut gamma fixes the legs.
This is the identity (a) with every vertex/edge/leaf weight replaced by 1 (inner edges by t): it tests the weight 1/|Aut Gt| per labelled
subset, the contraction, the automorphism numbers of the generator and the completeness of the lists.  Exact rational arithmetic.
The contraction used here is the independent one of rv_parts.py (equal to Split.contract on all pairs, see out/parts.log).
The vertex types (g_v, n_v) are generated in order of size; generation stops when the time budget is exhausted, and only the outer graphs
all of whose vertex types are available are tested (the coverage is reported).
Usage: py -I rv_mass.py [genus (default 3)] [budget minutes (default 10)] [types to skip, e.g. 0:8,0:7]   -> out/mass_g<genus>[b].log
"""
import os, sys, time
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_common import PROJ, OUT
from fractions import Fraction as Fr
from fq3_graphs import canon, legged_graphs

g = int(sys.argv[1]) if len(sys.argv) > 1 else 3; budget = 60 * float(sys.argv[2]) if len(sys.argv) > 2 else 600.0
EXCL = set(tuple(int(y) for y in x.split(":")) for x in sys.argv[3].split(",")) if len(sys.argv) > 3 else set()     # optional: types to skip, e.g. 0:8,0:7
LOG = open(os.path.join(OUT, "mass_g%d%s.log" % (g, "b" if EXCL else "")), "w", encoding="utf-8"); FAILED = []
def log(s):
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)


def contract_indep(G, S):
    gen = {v: G['gv'][v] for v in range(G['V'])}
    edges = [[u, w, bool((S >> e) & 1)] for e, (u, w) in enumerate(G['ed'])]
    while True:
        k = next((i for i, e in enumerate(edges) if not e[2]), None)
        if k is None: break
        u, w, _ = edges.pop(k)
        if u == w: gen[u] += 1; continue
        gen[u] += gen.pop(w)
        for e in edges:
            if e[0] == w: e[0] = u
            if e[1] == w: e[1] = u
    names = sorted(gen); idx = {v: i for i, v in enumerate(names)}
    return canon(tuple(gen[v] for v in names), [(idx[u], idx[w]) for u, w, _ in edges], ())


def padd(a, b, c=1):
    out = dict(a)
    for k, v in b.items(): out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v}
def pmul(a, b):
    out = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items(): out[k1 + k2] = out.get(k1 + k2, 0) + v1 * v2
    return {k: v for k, v in out.items() if v}


t0 = time.time()
graphs = legged_graphs(g, 0, 3 * g - 3); pos = {G['key']: i for i, G in enumerate(graphs)}; n = len(graphs)
log("rv_mass: genus %d, %d leg-free stable graphs" % (g, n))
lhs = [dict() for _ in range(n)]; npairs = 0
for G in graphs:
    E = G['E']; w = Fr(1, G['aut'])
    for S in range(1 << E):
        i = pos[contract_indep(G, S)]; k = E - bin(S).count("1"); lhs[i][k] = lhs[i].get(k, 0) + w; npairs += 1
log("   left-hand sides: %d pairs (Gt, S) in %.0f s" % (npairs, time.time() - t0))
types = sorted(set((G['gv'][v], len(G['hsAt'][v])) for G in graphs for v in range(G['V'])), key=lambda x: (3 * x[0] - 3 + x[1], x[1], x[0]))
N = {}; skipped = []
for (gv, nv) in types:
    if time.time() - t0 > budget or (gv, nv) in EXCL: skipped.append((gv, nv)); continue
    t1 = time.time(); L = legged_graphs(gv, nv, 3 * gv - 3 + nv); poly = {}
    for H in L: poly[H['E']] = poly.get(H['E'], 0) + Fr(1, H['aut'])
    N[(gv, nv)] = poly
    log("   N_{%d,%d}(t): %d legged graphs, %.1f s: %s" % (gv, nv, len(L), time.time() - t1, " + ".join("%s t^%d" % (poly[k], k) for k in sorted(poly))))
if skipped: log("   vertex types NOT generated (time budget): %s" % skipped)
bad = []; tested = 0; untested = []
for i, G in enumerate(graphs):
    ty = [(G['gv'][v], len(G['hsAt'][v])) for v in range(G['V'])]
    if any(t not in N for t in ty): untested.append(i); continue
    rhs = {0: Fr(1, G['aut'])}
    for t in ty: rhs = pmul(rhs, N[t])
    tested += 1
    if rhs != {k: v for k, v in lhs[i].items() if v}: bad.append(i)
check("genus %d: regrouping identity with weights t^(inner edges) for %d of %d outer graphs" % (g, tested, n), not bad and tested > 0, bad[:10])
log("   untested outer graphs (vertex type too large for the budget): %d; by E: %s" % (len(untested), {E: sum(1 for i in untested if graphs[i]['E'] == E) for E in range(3 * g - 2) if any(graphs[i]['E'] == E for i in untested)}))
# total: Sum over all (Gt, S) of 1/|Aut Gt| at t = 1 is Sum_Gt 2^E/|Aut Gt|
tot = sum(sum(l.values()) for l in lhs)
check("bookkeeping: Sum of all left-hand sides at t = 1 equals Sum_Gt 2^E/|Aut Gt|", tot == sum(Fr(2 ** G['E'], G['aut']) for G in graphs), str(tot))
log("RV_MASS genus %d %s  (%.0f s)" % (g, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0))
