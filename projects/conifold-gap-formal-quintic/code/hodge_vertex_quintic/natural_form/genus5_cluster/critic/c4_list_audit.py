# -*- coding: utf-8 -*-
"""c4_list_audit.py -- audit the SHIPPED genus-5 graph list (data/graphs_5_0_12.json) with the critic's own code:
validity (connected, genus 5, stable, leg-free), |Aut| recomputed by backtracking, pairwise non-isomorphism
(brute-force canonical key), and closure under one-edge contraction.
Run:  py -I -B c4_list_audit.py
"""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import c2_graphlib as L

PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")


def load(g, Emax):
    raw = json.load(open(os.path.join(PKG, "data", "graphs_%d_0_%d.json" % (g, Emax))))
    return [((tuple(k[0]), [tuple(e) for e in k[1]], tuple(k[2])), a) for k, a in raw]


def contract(gv, ed, e):
    u, w = ed[e]
    rest = [x for i, x in enumerate(ed) if i != e]
    if u == w:
        gv2 = list(gv)
        gv2[u] += 1
        return tuple(gv2), rest
    a, b = min(u, w), max(u, w)
    ren = {}
    nv = 0
    for v in range(len(gv)):
        if v == b:
            continue
        ren[v] = nv
        nv += 1
    ren[b] = ren[a]
    gv2 = [0] * nv
    for v in range(len(gv)):
        if v != b:
            gv2[ren[v]] = gv[v]
    gv2[ren[a]] = gv[a] + gv[b]
    return tuple(gv2), [(ren[x], ren[y]) for x, y in rest]


def main():
    g = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    Emax = 3 * g - 3
    flat = load(g, Emax)
    print("shipped list genus %d: %d entries" % (g, len(flat)))
    bad = []
    autbad = []
    t0 = time.time()
    keys = {}
    for n, (k, aut) in enumerate(flat):
        gv, ed, legs = k
        V = len(gv)
        M = L.mat_of(V, ed)
        if legs or not L.connected(V, M) or L.genus(V, gv, M) != g or not L.stable(V, gv, M):
            bad.append((n, k))
        a2 = L.aut_total(V, gv, M)
        if a2 != aut:
            autbad.append((n, k, aut, a2))
        ck = L.canon_brute(V, gv, M)
        keys.setdefault(ck, []).append(n)
    print("  invalid graphs (not connected / wrong genus / unstable / with legs): %d %s" % (len(bad), bad[:5]))
    print("  |Aut| disagreements with the critic's backtracking count: %d %s" % (len(autbad), autbad[:5]))
    dup = {k: v for k, v in keys.items() if len(v) > 1}
    print("  distinct up to isomorphism (critic's brute-force canonical key): %d of %d; duplicates: %d %s"
          % (len(keys), len(flat), len(dup), list(dup.values())[:5]))
    # closure under one-edge contraction
    ckey = {}
    for n, (k, aut) in enumerate(flat):
        gv, ed, legs = k
        ckey[L.canon_brute(len(gv), gv, L.mat_of(len(gv), ed))] = n
    miss = 0
    for n, (k, aut) in enumerate(flat):
        gv, ed, legs = k
        for e in range(len(ed)):
            gv2, ed2 = contract(gv, ed, e)
            k2 = L.canon_brute(len(gv2), gv2, L.mat_of(len(gv2), ed2))
            if k2 not in ckey:
                miss += 1
    print("  one-edge contractions landing outside the list: %d (must be 0)" % miss)
    byE = {}
    for k, a in flat:
        byE[len(k[1])] = byE.get(len(k[1]), 0) + 1
    print("  by number of edges:", sorted(byE.items()))
    print("  %.0f s" % (time.time() - t0))


main()
