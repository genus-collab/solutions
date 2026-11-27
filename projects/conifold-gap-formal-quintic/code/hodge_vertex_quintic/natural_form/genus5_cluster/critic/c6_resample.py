# -*- coding: utf-8 -*-
"""c6_resample.py -- the critic re-runs the genus-5 engine on a sample of inner graphs and re-does, locally, the
two checks that CANNOT be redone from the copied results/ files (parts/ was not copied back):
  (A) ADDITIVITY: for each sampled inner graph t, the sum over outer graphs of its split pieces must equal its
      point-frame term, i.e. results/g5_89_id.json table[t].  (This is the collect step's additivity check; the
      copied JSONs record neither its outcome nor whether it ran.)
  (B) graph_dp (the fast one-elimination engine used on the cluster) vs fq4_split.Split.graph (the slow
      subset-by-subset path) on the same graphs.
  (C) valuation structure: every piece going to an outer graph with |S| edges must be O(u^|S|) (this is what makes
      the engine's cap |S| <= K-1 = 7 sound, and what makes the 2891 outer graphs with E >= 8 regular for free).
Run:  py -I -B c6_resample.py [ncheap] [nslow]
"""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)


def main():
    ncheap = int(sys.argv[1]) if len(sys.argv) > 1 else 250
    nslow = int(sys.argv[2]) if len(sys.argv) > 2 else 80
    from g5_engine import make_split, graph_dp
    from fq_fp import P2
    idtab = json.load(open(os.path.join(PKG, "results", "g5_89_id.json")))["table"]
    cost = json.load(open(os.path.join(PKG, "data", "cost_table_g5.json")))
    order = sorted(range(4555), key=lambda t: cost[str(t)]["est_dp"])
    sample = order[:ncheap]
    t0 = time.time()
    sp = make_split(P2, 5, ("0", "0"), 8, "nat")
    K = sp.K
    p = sp.p
    print("Split built in %.0f s" % (time.time() - t0))
    badA = []
    badC = []
    pieces_nz = 0
    t0 = time.time()
    store = {}
    for t in sample:
        res = graph_dp(sp, t, None)
        store[t] = res
        s = [0] * K
        for q, v in res.items():
            if any(v):
                pieces_nz += 1
            Eq = sp.inner_graphs[q]["E"]
            val = next((j for j, x in enumerate(v) if x), K)
            if val < Eq:
                badC.append((t, q, Eq, val))
            for k in range(K):
                s[k] = (s[k] + v[k]) % p
        ref = idtab.get(str(t), [0] * K)
        if s != ref:
            badA.append(t)
    print("(A) additivity on %d inner graphs (the %d cheapest), %.0f s: %d failures %s"
          % (len(sample), ncheap, time.time() - t0, len(badA), badA[:10]))
    print("    non-zero (inner graph, outer graph) pieces seen in the sample: %d" % pieces_nz)
    print("(C) pieces with u-valuation below the number of outer edges: %d %s" % (len(badC), badC[:5]))
    t0 = time.time()
    badB = []
    for t in order[:nslow]:
        a = {q: v for q, v in graph_dp(sp, t, None).items() if any(v)}
        b = {q: v for q, v in sp.graph(t).items() if any(v)}
        if a != b:
            badB.append(t)
    print("(B) graph_dp vs Split.graph on %d inner graphs, %.0f s: %d disagreements %s"
          % (nslow, time.time() - t0, len(badB), badB[:10]))
    # how much of the regularity claim does the sample touch?
    touched = set()
    for t, res in store.items():
        for q, v in res.items():
            if any(v):
                touched.add(q)
    print("    outer graphs with edges receiving a non-zero piece from the sample: %d (of the 627 that the full run"
          " says cancel)" % len([q for q in touched if q != 0]))


main()
