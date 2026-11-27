# -*- coding: utf-8 -*-
"""ig_canon.py -- own canonical form and own |Aut| for all 4555 genus-5 stable graphs, from the INDEPENDENT code
(natural_form/genus4_referee/indep: ind_graphs.canon = minimum over vertex permutations, ind_graphs.n_aut = depth-first
count of half-edge bijections).  The genera and the edge lists are taken from the package list
cluster_genus5_quintic/data/graphs_5_0_12.json; nothing else is taken from it.

Checks performed here:
  * every graph is connected, stable and of genus 5 (own check, not the package's);
  * the 4555 own canonical forms are pairwise distinct (the list is non-redundant as seen by an independent canon);
  * own |Aut| equals the package's |Aut| for every graph;
  * contracting all edges gives the vertex graph and contracting none gives the graph itself.
Writes out/canon_map.json: {"keys": [canonical form per index], "map": {canonical form as a string -> index}, "aut": [...]}.
  py -I -B ig_canon.py
"""
import os, sys, json, time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
IND = os.path.normpath(os.path.join(HERE, "..", "..", "genus4_referee", "indep"))
PKG = os.environ.get("G5PKG", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic"))
sys.path.insert(0, IND)


def keystr(key):
    gv, ed = key
    return "%s|%s" % (",".join(str(x) for x in gv), ",".join("%d-%d" % e for e in ed))


def main():
    from ind_graphs import canon, n_aut, contract, make
    graphs = json.load(open(os.path.join(PKG, "data", "graphs_5_0_12.json")))
    print("package list: %d graphs" % len(graphs))
    t0 = time.time()
    keys = []
    auts = []
    bad = []
    for t, (gk, aut) in enumerate(graphs):
        gv, ed, legs = gk
        assert not legs, t
        gv = tuple(gv)
        ed = [tuple(e) for e in ed]
        V = len(gv)
        E = len(ed)
        val = [0] * V
        for u, w in ed:
            val[u] += 1
            val[w] += 1
        if sum(gv) + E - V + 1 != 5:
            bad.append((t, "genus"))
        if not all(2 * gv[v] - 2 + val[v] > 0 for v in range(V)):
            bad.append((t, "stability"))
        comp = list(range(V))
        for _ in range(V):
            for u, w in ed:
                m = min(comp[u], comp[w])
                comp[u] = comp[w] = m
        if len(set(comp)) != 1:
            bad.append((t, "connected"))
        keys.append(canon(gv, ed))
        a = n_aut(gv, ed)
        auts.append(a)
        if a != aut:
            bad.append((t, "aut own %d package %d" % (a, aut)))
        if t % 500 == 0:
            print("   %d/%d (%.0f s)" % (t, len(graphs), time.time() - t0), flush=True)
    print("own canon + own |Aut| for all %d graphs in %.0f s" % (len(graphs), time.time() - t0))
    ks = [keystr(k) for k in keys]
    print("distinct own canonical forms: %d of %d" % (len(set(ks)), len(ks)))
    assert len(set(ks)) == len(ks), "the package list is redundant as seen by the independent canon"
    assert not bad, bad[:20]
    print("own |Aut| equals the package |Aut| for all %d graphs; max |Aut| = %d" % (len(graphs), max(auts)))
    # contraction sanity with the independent contract()
    cm = {k: i for i, k in enumerate(ks)}
    nchk = 0
    for t, (gk, aut) in enumerate(graphs):
        gv, ed, _ = gk
        G = make(tuple(gv), [tuple(e) for e in ed])
        assert cm[keystr(contract(G, 0))] == 0, t
        assert cm[keystr(contract(G, (1 << G["E"]) - 1))] == t, t
        nchk += 2
    print("contraction sanity: %d checks (contract everything -> vertex graph, contract nothing -> itself) all pass" % nchk)
    json.dump(dict(keys=ks, map=cm, aut=auts), open(os.path.join(OUT, "canon_map.json"), "w"))
    print("wrote %s" % os.path.join(OUT, "canon_map.json"))


if __name__ == "__main__":
    main()
