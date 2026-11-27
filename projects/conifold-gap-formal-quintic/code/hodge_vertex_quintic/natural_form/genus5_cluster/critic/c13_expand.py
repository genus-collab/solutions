# -*- coding: utf-8 -*-
"""c13_expand.py -- COMPLETENESS of the shipped genus-g graph list, by an independent induction on the number of edges.
Every stable graph with E >= 1 edges becomes a stable graph with E-1 edges when any edge is contracted, so a list that
(i) is right at E = 0 (the single genus-g vertex) and (ii) contains, at every level E, every one-edge EXPANSION of every
graph at level E-1, is complete.  The expansions are the inverse of contraction: split a vertex (distribute its
half-edges and its genus between the two halves, join them by a new edge), or lower a vertex's genus by one and add a
loop.  This script generates them with the critic's own code, canonicalises with c2_graphlib.canon_brute, and compares
level by level with data/graphs_<g>_0_<Emax>.json.  Independent of fq3_graphs.degenerations and of g5_graphs.canon_count.
Run:  py -I -B c13_expand.py [g]
"""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import c2_graphlib as L
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")


def key_of(gv, ed):
    return L.canon_brute(len(gv), tuple(gv), L.mat_of(len(gv), ed))


def expansions(gv, ed):
    """all one-edge expansions (gv2, ed2) of the leg-free stable graph (gv, ed)"""
    V = len(gv)
    out = []
    inc = [[] for _ in range(V)]
    for e, (u, w) in enumerate(ed):
        inc[u].append((e, 0))
        if u == w:
            inc[u].append((e, 1))
        else:
            inc[w].append((e, 1))
    for v in range(V):
        # (b) lower the genus by one and add a loop
        if gv[v] >= 1:
            gv2 = list(gv)
            gv2[v] -= 1
            out.append((tuple(gv2), list(ed) + [(v, v)]))
        # (a) split the vertex
        h = inc[v]
        n = len(h)
        for mask in range(1 << n):
            keep = [h[i] for i in range(n) if (mask >> i) & 1]
            move = [h[i] for i in range(n) if not (mask >> i) & 1]
            for g1 in range(gv[v] + 1):
                g2 = gv[v] - g1
                # v keeps `keep` with genus g1; the new vertex V takes `move` with genus g2
                if 2 * g1 - 2 + len(keep) + 1 <= 0 or 2 * g2 - 2 + len(move) + 1 <= 0:
                    continue
                gv2 = list(gv) + [g2]
                gv2[v] = g1
                ed2 = []
                for e, (u, w) in enumerate(ed):
                    a, b = u, w
                    if (e, 0) in move:
                        a = V
                    if (e, 1) in move:
                        b = V
                    ed2.append((a, b))
                ed2.append((v, V))
                out.append((tuple(gv2), ed2))
    return out


def main():
    g = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    Emax = 3 * g - 3
    raw = json.load(open(os.path.join(PKG, "data", "graphs_%d_0_%d.json" % (g, Emax))))
    flat = [(tuple(k[0]), [tuple(e) for e in k[1]]) for k, a in raw]
    levels = {}
    for gv, ed in flat:
        levels.setdefault(len(ed), []).append((gv, ed))
    keys = {E: set(key_of(gv, ed) for gv, ed in levels[E]) for E in levels}
    print("shipped list genus %d: %s" % (g, {E: len(keys[E]) for E in sorted(keys)}))
    assert len(levels[0]) == 1 and levels[0][0][0] == (g,) and not levels[0][0][1], "base case"
    print("base case E = 0: the single genus-%d vertex -- OK" % g)
    ok = True
    for E in sorted(levels):
        if E == Emax:
            break
        t0 = time.time()
        gen = set()
        for gv, ed in levels[E]:
            for gv2, ed2 in expansions(gv, ed):
                if len(ed2) != E + 1:
                    continue
                if L.genus(len(gv2), gv2, L.mat_of(len(gv2), ed2)) != g:
                    continue
                if not L.stable(len(gv2), gv2, L.mat_of(len(gv2), ed2)):
                    continue
                gen.add(key_of(gv2, ed2))
        extra = gen - keys[E + 1]
        missing = keys[E + 1] - gen
        print("  level %2d -> %2d: %5d classes generated, shipped %5d; missing from the list: %d; in the list but not"
              " generated: %d  (%.0f s)" % (E, E + 1, len(gen), len(keys[E + 1]), len(extra), len(missing),
                                            time.time() - t0), flush=True)
        if extra or missing:
            ok = False
            for k in list(extra)[:3]:
                print("      NOT IN THE SHIPPED LIST:", k)
            for k in list(missing)[:3]:
                print("      IN THE LIST BUT NOT AN EXPANSION:", k)
    print("COMPLETENESS:", "the shipped list is exactly the closure -- no class is missing" if ok else "MISMATCH")


main()
