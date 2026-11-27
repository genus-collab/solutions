# -*- coding: utf-8 -*-
"""fq3_graphs.py -- stable graphs WITH LEGS by vertex splitting (genus-3 programme), and the 42 leg-free genus-3 outer graphs.

legged_graphs(g, n, Emax): isomorphism classes of stable graphs of type (g, n) (legs labelled 0..n-1) with at most Emax edges,
generated level by level from the single vertex (g, n) by one-edge degenerations (non-separating: g_v -> g_v - 1 plus a loop;
separating: split a vertex into two stable vertices joined by a new edge, distributing its half-edges and legs), with canonical
forms under S_V (fq_mvertex._canon conventions: key = (gv, sorted edge multiset, leg -> vertex)).  Every stable graph of type (g, n)
with E edges is a one-edge degeneration of one with E - 1 edges (smooth a non-separating or separating edge), so the levels are complete.
Output: list of dicts with the fields of fq_mvertex.legged_graphs (V, gv, ed [(u, v), u <= v, loops allowed], legs, aut, H, hsAt,
legsAt, E), sorted by (E, key); |Aut| = number of (vertex perm, edge perm, loop flips) preserving the graph and fixing the legs
(fq_graphs.autos_of filtered).  Cross-checks: the leg-free lists equal the Mathematica exports StableGraph-2.bin (7) and
StableGraph-3.bin (42) as sets of canonical forms; at (g, n) with Emax = 3, the list equals fq_mvertex.legged_graphs (brute force).
outer_graphs3(): the 42 leg-free genus-3 graphs in fq_msum.outer_graphs format (idx = Mathematica order, gv, ed, H, hsAt, V, aut)
with the extra field E (number of edges).
"""
import itertools, math, os, sys
from functools import lru_cache
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_graphs import load_graphs, edges_of, autos_of, aut_count


def canon(gv, ed, legs):
    V = len(gv); best = None
    for sig in itertools.permutations(range(V)):
        gv2 = [0] * V
        for v in range(V): gv2[sig[v]] = gv[v]
        ed2 = tuple(sorted(tuple(sorted((sig[u], sig[v]))) for u, v in ed))
        legs2 = tuple(sig[l] for l in legs)
        key = (tuple(gv2), ed2, legs2)
        if best is None or key < best: best = key
    return best


def degenerations(gv, ed, legs):
    """all one-edge degenerations of the stable graph (gv, ed, legs); legs = tuple vertex of each leg"""
    V = len(gv); out = []
    for v in range(V):
        hv = [(e, s) for e, (a, b) in enumerate(ed) for s, x in enumerate((a, b)) if x == v]
        lv = [l for l, x in enumerate(legs) if x == v]
        if gv[v] >= 1:                                            # non-separating node
            gv2 = list(gv); gv2[v] -= 1
            out.append((tuple(gv2), list(ed) + [(v, v)], tuple(legs)))
        items = hv + [('L', l) for l in lv]                       # separating node: new vertex V
        for g1 in range(gv[v] + 1):
            g2 = gv[v] - g1
            for r in range(len(items) + 1):
                for S in itertools.combinations(range(len(items)), r):
                    n1 = len(items) - r + 1; n2 = r + 1
                    if 2 * g1 - 2 + n1 <= 0 or 2 * g2 - 2 + n2 <= 0: continue
                    gv2 = list(gv) + [g2]; gv2[v] = g1
                    ed2 = [list(e) for e in ed]; legs2 = list(legs)
                    for k in S:
                        it = items[k]
                        if it[0] == 'L': legs2[it[1]] = V
                        else: ed2[it[0]][it[1]] = V
                    ed2 = [tuple(sorted(e)) for e in ed2] + [(v, V)]
                    out.append((tuple(gv2), ed2, tuple(legs2)))
    return out


@lru_cache(maxsize=None)
def legged_keys(g, n, Emax):
    """list over E = 0..Emax of sorted lists of canonical keys"""
    assert 2 * g - 2 + n > 0, "unstable type (%d, %d)" % (g, n)
    start = ((g,), [], tuple([0] * n))
    levels = [[canon(*start)]]
    for E in range(1, Emax + 1):
        new = set()
        for key in levels[-1]:
            for H in degenerations(key[0], list(key[1]), key[2]):
                new.add(canon(*H))
        levels.append(sorted(new))
    return levels


def _graph_dict(key, n):
    gv, ed, legs = key; V = len(gv); ed = list(ed)
    A = [[0] * V for _ in range(V)]
    for u, v in ed:
        if u == v: A[u][u] += 1
        else: A[u][v] += 1; A[v][u] += 1
    auts = [a for a in autos_of(gv, A) if all(a[0][l] == l for l in legs)]
    H = [(e, s) for e in range(len(ed)) for s in (0, 1)]
    hvtx = [ed[e][s] for e, s in H]
    hsAt = [[h for h, x in enumerate(hvtx) if x == v] for v in range(V)]
    legsAt = [[l for l in range(n) if legs[l] == v] for v in range(V)]
    return dict(V=V, gv=gv, ed=ed, legs=legs, aut=len(auts), H=H, hsAt=hsAt, legsAt=legsAt, E=len(ed), key=key)


@lru_cache(maxsize=None)
def legged_graphs(g, n, Emax):
    """list of graph dicts (see module docstring), all E <= Emax, sorted by (E, key)"""
    out = []
    for lev in legged_keys(g, n, Emax):
        for key in lev: out.append(_graph_dict(key, n))
    return out


def counts(g, n, Emax):
    return [len(l) for l in legged_keys(g, n, Emax)]


# ------------------------------------------------------------------ the genus-3 outer graphs (Mathematica order)
def outer_graphs(g=3):
    out = []
    for gi, (gv, m, A) in enumerate(load_graphs(g)):
        V = len(gv); ed = edges_of(A)
        H = [(e, s) for e in range(len(ed)) for s in (0, 1)]
        hvtx = [ed[e][s] for e, s in H]
        hsAt = [[h for h, x in enumerate(hvtx) if x == v] for v in range(V)]
        out.append(dict(idx=gi + 1, gv=gv, ed=ed, H=H, hsAt=hsAt, V=V, aut=aut_count(gv, A), E=len(ed), key=canon(gv, ed, ())))
    return out


def outer_graphs3(): return outer_graphs(3)


def mathematica_check(g):
    """the leg-free level lists equal the Mathematica export as sets of canonical forms"""
    mine = set(k for lev in legged_keys(g, 0, 3 * g - 3) for k in lev)
    ref = set(canon(gv, edges_of(A), ()) for gv, m, A in load_graphs(g))
    return mine == ref, len(mine), len(ref)


if __name__ == "__main__":
    import time
    for g in (2, 3):
        ok, a, b = mathematica_check(g); print("genus %d leg-free: generator %d, Mathematica %d, equal %s" % (g, a, b, ok)); assert ok
    # brute-force comparison with fq_mvertex.legged_graphs at small sizes (same canonical keys)
    from fq_mvertex import legged_graphs as brute
    for (g, n, E) in [(0, 3, 0), (0, 4, 1), (0, 5, 2), (1, 1, 1), (1, 2, 2), (1, 3, 3), (2, 0, 3), (2, 1, 3), (2, 2, 3), (0, 6, 3)]:
        t = time.time(); B = brute(g, n, E); tb = time.time() - t
        t = time.time(); L = legged_graphs(g, n, E); tl = time.time() - t
        kb = sorted(canon(G['gv'], G['ed'], G['legs']) for G in B); kl = sorted(G['key'] for G in L)
        ab = sorted((canon(G['gv'], G['ed'], G['legs']), G['aut']) for G in B); al = sorted((G['key'], G['aut']) for G in L)
        assert kb == kl and ab == al, (g, n, E)
        print("(%d,%d,E<=%d): %d graphs, |Aut| agree; brute %.1f s, splitting %.2f s" % (g, n, E, len(L), tb, tl))
    for (g, n, E) in [(3, 0, 6), (3, 1, 5), (3, 2, 4), (2, 2, 5), (2, 3, 4), (1, 4, 4), (1, 5, 3), (0, 6, 3), (0, 7, 2)]:
        t = time.time(); L = legged_graphs(g, n, E)
        print("(%d,%d,E<=%d): per level %s, total %d, %.2f s" % (g, n, E, counts(g, n, E), len(L), time.time() - t))
    G3 = outer_graphs3(); assert len(G3) == 42
    print("outer genus-3 graphs: 42; edge counts:", sorted(set(G['E'] for G in G3)), "; by E:", {E: sum(1 for G in G3 if G['E'] == E) for E in range(7)})
    print("fq3_graphs ok")
