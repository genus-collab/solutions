# -*- coding: utf-8 -*-
"""rv_parts.py -- read-only analysis of the genus-4 part files of record and an independent test of Split.contract.

(1) For every inner graph t and both primes: Sum over the outer graphs of the nat pieces of t (parts/g4_<p>_nat_run.json) equals the
    point-frame term of t (parts/g4_<p>_id_run.json), all K = 6 coefficients.  (Independent of E_M: see REVIEW.md (d).)
(2) Statistics of the table of record g4_split_g4_<p>_nat.json: which outer graphs received a non-zero term at all, by number of edges.
(3) Split.contract(G, S) against an independent contraction (inner edges contracted one at a time on an explicit multigraph), for all
    379 inner graphs and all 2^E subsets S; the set of keys reached = keys of the 379 stable graphs; S = {} -> vertex graph only;
    S = all edges -> the graph itself; number of outer edges of the key = |S|; total genus of the key = 4.
(4) Cross-prime: the vertex-graph rows of the two primes reconstruct to the same rationals where reconstruction succeeds (u^0..u^3).
Usage: py -I rv_parts.py    -> out/parts.log
"""
import os, sys, json, time
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_common import PROJ, OUT
from fractions import Fraction as Fr
from fq_fp import FpRing, P1, P2
from fq3_graphs import canon, legged_graphs

LOG = open(os.path.join(OUT, "parts.log"), "w", encoding="utf-8"); FAILED = []
def log(s):
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)

t0 = time.time(); K = 6
graphs = legged_graphs(4, 0, 9); nin = len(graphs); Eof = [G['E'] for G in graphs]
check("379 leg-free stable graphs of genus 4; by E: %s" % {E: Eof.count(E) for E in range(10)}, nin == 379)
log("   Sum 1/|Aut| by E (reported only): %s" % {E: str(sum(Fr(1, G['aut']) for G in graphs if G['E'] == E)) for E in range(10)})

# ---------------------------------------------------------------- (1), (2)
for prime, p in (("89", P2), ("61", P1)):
    nat = json.load(open(os.path.join(PROJ, "parts", "g4_%s_nat_run.json" % prime)))['res']
    idr = json.load(open(os.path.join(PROJ, "parts", "g4_%s_id_run.json" % prime)))['res']
    check("2^%s-1: both part files cover the 379 inner graphs" % prime, set(nat) == set(idr) == set(str(t) for t in range(nin)))
    bad = []; nterms = 0
    for t in range(nin):
        s = [0] * K
        for pos, c in nat[str(t)].items():
            nterms += 1
            for k in range(K): s[k] = (s[k] + c[k]) % p
        ref = idr[str(t)]
        assert set(ref) <= {str(t)}, "the id run must only fill the graph itself"
        if s != [x % p for x in ref.get(str(t), [0] * K)]: bad.append(t)
    check("2^%s-1 (1): for each of the 379 inner graphs the nat pieces add up to the point-frame term, u^0..u^5 (%d pieces)" % (prime, nterms), not bad, bad[:10])
    tabd = json.load(open(os.path.join(PROJ, "g4_split_g4_%s_nat.json" % prime)))
    rows = sorted(int(k) for k in tabd['table']); nz = [i for i in rows if any(x % p for x in tabd['table'][str(i)])]
    log("   2^%s-1 (2): rows present in the table of record by E: %s; non-zero rows: %s; failed: %s" % (
        prime, {E: sum(1 for i in rows if Eof[i] == E) for E in range(10)}, nz, tabd['failed']))
    # pieces of the outer graphs with edges: how many non-zero pieces cancel
    cnt = {}
    for t in range(nin):
        for pos, c in nat[str(t)].items():
            if int(pos) != 0 and any(x % p for x in c): cnt[int(pos)] = cnt.get(int(pos), 0) + 1
    log("      outer graphs with edges that receive non-zero pieces (which then cancel): %d; by E: %s; total non-zero pieces %d" % (
        len(cnt), {E: sum(1 for i in cnt if Eof[i] == E) for E in range(1, 10)}, sum(cnt.values())))
    lowest = {}
    for t in range(nin):
        for pos, c in nat[str(t)].items():
            i = int(pos)
            if i == 0: continue
            v = next((k for k in range(K) if c[k] % p), K); lowest[i] = min(lowest.get(i, K), v)
    log("      lowest u-order of a non-zero piece minus E, histogram over these outer graphs: %s" % (
        {d: sum(1 for i in lowest if lowest[i] < K and lowest[i] - Eof[i] == d) for d in range(0, 6)}))
    check("2^%s-1: no piece of an outer graph with E edges is non-zero below u^E" % prime, all(lowest[i] >= min(Eof[i], K) for i in lowest))

# ---------------------------------------------------------------- (4) cross-prime rationals of the vertex graph
R2 = FpRing(P2); R1 = FpRing(P1)
v2 = json.load(open(os.path.join(PROJ, "g4_split_g4_89_nat.json")))['table']['0']; v1 = json.load(open(os.path.join(PROJ, "g4_split_g4_61_nat.json")))['table']['0']
rec = [R2.ratrec(x) for x in v2]
check("(4) vertex graph: rationals reconstructed at 2^89-1 reduce to the residues at 2^61-1 (orders %s)" % [k for k in range(K) if rec[k] is not None],
      all(rec[k] is None or R1.of(rec[k]) == v1[k] % P1 for k in range(K)), [str(r) for r in rec])

# ---------------------------------------------------------------- (3) independent contraction
def contract_indep(G, S):
    """contract the edges not in S one at a time; returns the canonical key"""
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

from fq4_split import Split
sp = Split(P2, 4, (Fr(0), Fr(0)), K, 'id')
keyset = set(sp.outer_keys); reached = set(); mism = 0; nsub = 0; badS = 0
for t, G in enumerate(sp.inner_graphs):
    E = G['E']
    for S in range(1 << E):
        k1 = sp.contract(G, S); k2 = contract_indep(G, S); nsub += 1
        if k1 != k2: mism += 1
        reached.add(k1)
        gtot = sum(k1[0]) + len(k1[1]) - len(k1[0]) + 1
        if k1 not in keyset or len(k1[1]) != bin(S).count("1") or gtot != 4: badS += 1
        if S == 0 and sp.outer_pos[k1] != 0: badS += 1
        if S == (1 << E) - 1 and sp.outer_pos[k1] != t: badS += 1
check("(3) Split.contract = independent contraction on all %d pairs (inner graph, subset)" % nsub, mism == 0, mism)
check("(3) every contraction is one of the 379 keys, has |S| edges and genus 4; S = {} is the vertex graph; S = all is the graph itself", badS == 0, badS)
check("(3) all 379 outer graphs are reached", reached == keyset)
log("RV_PARTS %s  (%.0f s)" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0))
