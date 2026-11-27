# -*- coding: utf-8 -*-
"""c5_contract_audit.py -- independent audit of the map (inner graph, outer-edge set S) -> outer graph index,
i.e. of fq4_split.Split.contract composed with Split.outer_pos.  This map decides WHICH outer graph each piece
of the split sum is added to; nothing in the genus-5 checks constrains it (the total and the per-inner-graph
additivity both sum over S), and it is validated against truly independent code only at genus 2-4.
Here it is recomputed with the critic's own contraction and canonical form (c2_graphlib) and compared.
Run:  py -I -B c5_contract_audit.py [nlarge]
"""
import json, os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import c2_graphlib as L
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)


def my_contract(gv, ed, S):
    """contract the edges NOT in the bitmask S; return (gv2, ed2) of the outer graph (critic's own code)"""
    V = len(gv)
    par = list(range(V))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for e, (u, w) in enumerate(ed):
        if not (S >> e) & 1:
            a, b = find(u), find(w)
            if a != b:
                par[a] = b
    roots = sorted(set(find(v) for v in range(V)))
    cid = {r: i for i, r in enumerate(roots)}
    n = len(roots)
    g2 = [0] * n
    nv = [0] * n
    ne = [0] * n
    for v in range(V):
        c = cid[find(v)]
        g2[c] += gv[v]
        nv[c] += 1
    out = []
    for e, (u, w) in enumerate(ed):
        if (S >> e) & 1:
            out.append((cid[find(u)], cid[find(w)]))
        else:
            ne[cid[find(u)]] += 1
    g2 = [g2[c] + ne[c] - nv[c] + 1 for c in range(n)]
    return tuple(g2), out


def main():
    nlarge = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    from g5_engine import make_split
    from fq_fp import P2
    t0 = time.time()
    sp = make_split(P2, 5, ("0", "0"), 8, "nat")
    print("Split built in %.0f s; %d inner graphs, K = %d" % (time.time() - t0, len(sp.inner_graphs), sp.K))
    # critic's own index: brute-force canonical key -> position in the shipped list
    mine = {}
    for i, G in enumerate(sp.inner_graphs):
        mine[L.canon_brute(G['V'], G['gv'], L.mat_of(G['V'], G['ed']))] = i
    assert len(mine) == len(sp.inner_graphs)
    cap = sp.K - 1
    small = [t for t, G in enumerate(sp.inner_graphs) if G['E'] <= 7]
    random.seed(1)
    big = random.sample([t for t, G in enumerate(sp.inner_graphs) if G['E'] >= 8], nlarge)
    bad = []
    npairs = 0
    t0 = time.time()
    for t in small + big:
        G = sp.inner_graphs[t]
        E = G['E']
        for S in range(1 << E):
            if bin(S).count("1") > cap:
                continue
            npairs += 1
            pkg = sp.outer_pos[sp.contract(G, S)]
            gv2, ed2 = my_contract(G['gv'], G['ed'], S)
            key = L.canon_brute(len(gv2), gv2, L.mat_of(len(gv2), ed2))
            m = mine.get(key)
            if m != pkg:
                bad.append((t, S, pkg, m))
            # the outer graph's edge count must equal |S|
            if sp.inner_graphs[pkg]['E'] != bin(S).count("1"):
                bad.append((t, S, "E mismatch", sp.inner_graphs[pkg]['E']))
    print("checked %d (inner graph, outer-edge set) pairs on %d graphs with E <= 7 and %d with E >= 8 (%.0f s)"
          % (npairs, len(small), len(big), time.time() - t0))
    print("disagreements with the critic's own contraction + canonical form: %d %s" % (len(bad), bad[:5]))
    # the cap itself
    print("cap on |S| used by the engine: %d (= K-1); outer graphs with E >= %d therefore receive no piece at all"
          % (cap, cap + 1))
    n8 = sum(1 for G in sp.inner_graphs if G['E'] >= cap + 1)
    print("   that is %d of the 4554 outer graphs with edges (%.0f%%)" % (n8, 100.0 * n8 / 4554))
    print("   sound only because every entry of E_out is O(u): min valuation of E_out =", min(sp.valO.values()))


main()
