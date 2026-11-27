# -*- coding: utf-8 -*-
"""g5_graphs.py -- fast canonical form and generator of leg-free stable graphs (genus-5 package).

From the genus-5 feasibility study (hodge_vertex_quintic/natural_form/genus4_referee/genus5/g5_lib.py), unchanged in substance:
canon_count(gv, ed, legs) -> (key, number of vertex automorphisms): canonical key of the same SHAPE as fq3_graphs.canon (genera, sorted
    edge list, legs), minimised over the leaves of an individualisation-refinement search instead of all V! vertex orders.  The set of
    leaves is isomorphism invariant, so the minimum is a canonical form (in general a different representative than fq3_graphs.canon).
    The number of leaves attaining the minimum is the number of vertex automorphisms; |Aut| = that x prod_loops 2^l l! x prod_multi m!.
graphs(g, Emax): drop-in for fq3_graphs.legged_graphs(g, 0, Emax) (same dict fields, sorted by (E, key)); the lists are read from
    data/graphs_<g>_0_<Emax>.json when present (shipped for g = 2..5), otherwise generated (genus 5: about 10 s) and cached there.
Validated in the study: genus 4 = the 379 graphs of fq3_graphs with equal |Aut|; genus 5 = 4555 graphs (the known number), |Aut| equal
to the brute-force fq_graphs.aut_count for all of them.  The INDEX of a graph in these lists is used by the cost table and the task
assignment, so the shipped lists must not be regenerated in a different order.
"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
DATA = os.path.join(HERE, "data")


def _refine(V, adj, col):
    ncl = len(set(col))
    while True:
        sig = [(col[v], tuple(sorted((col[w], m) for w, m in adj[v]))) for v in range(V)]
        ranks = {s: i for i, s in enumerate(sorted(set(sig)))}
        col = [ranks[s] for s in sig]
        if len(ranks) == ncl: return col
        ncl = len(ranks)


def canon_count(gv, ed, legs):
    V = len(gv); loops = [0] * V; mult = {}
    for u, w in ed:
        if u == w: loops[u] += 1
        else:
            a, b = (u, w) if u < w else (w, u); mult[(a, b)] = mult.get((a, b), 0) + 1
    adj = [[] for _ in range(V)]
    for (a, b), m in mult.items(): adj[a].append((b, m)); adj[b].append((a, m))
    lg = [tuple(l for l, x in enumerate(legs) if x == v) for v in range(V)]
    sig0 = [(gv[v], loops[v], lg[v]) for v in range(V)]
    r0 = {s: i for i, s in enumerate(sorted(set(sig0)))}
    best = [None, 0]

    def leaf(col):
        gv2 = [0] * V
        for v in range(V): gv2[col[v]] = gv[v]
        ed2 = []
        for u, w in ed:
            a = col[u]; b = col[w]; ed2.append((a, b) if a <= b else (b, a))
        ed2.sort()
        key = (tuple(gv2), tuple(ed2), tuple(col[l] for l in legs))
        if best[0] is None or key < best[0]: best[0] = key; best[1] = 1
        elif key == best[0]: best[1] += 1

    def rec(col):
        col = _refine(V, adj, col)
        cnt = {}
        for c in col: cnt[c] = cnt.get(c, 0) + 1
        if len(cnt) == V: leaf(col); return
        c = min(k for k, n in cnt.items() if n > 1)
        for v in range(V):
            if col[v] == c:
                col2 = [2 * x for x in col]; col2[v] -= 1
                rec(col2)

    rec([r0[s] for s in sig0])
    return best[0], best[1]


def canon_fast(gv, ed, legs):
    return canon_count(gv, ed, legs)[0]


def aut_of(key, nvaut):
    gv, ed, legs = key; r = nvaut; mult = {}
    for e in ed: mult[e] = mult.get(e, 0) + 1
    for (u, w), m in mult.items():
        r *= math.factorial(m)
        if u == w: r *= 2 ** m
    return r


def levels(g, Emax):
    from fq3_graphs import degenerations
    k0, c0 = canon_count((g,), [], ()); lev = [[(k0, aut_of(k0, c0))]]
    for E in range(1, Emax + 1):
        new = {}; seen = set()
        for key, _ in lev[-1]:
            for gv2, ed2, legs2 in degenerations(key[0], list(key[1]), key[2]):
                raw = (gv2, tuple(sorted(ed2)), legs2)
                if raw in seen: continue
                seen.add(raw); k, c = canon_count(gv2, ed2, legs2)
                if k not in new: new[k] = aut_of(k, c)
        lev.append(sorted(new.items()))
    return lev


def graph_dict(key, aut):
    gv, ed, legs = key; V = len(gv); ed = [tuple(e) for e in ed]
    H = [(e, s) for e in range(len(ed)) for s in (0, 1)]
    hvtx = [ed[e][s] for e, s in H]
    hsAt = [[h for h, x in enumerate(hvtx) if x == v] for v in range(V)]
    return dict(V=V, gv=tuple(gv), ed=ed, legs=(), aut=aut, H=H, hsAt=hsAt, legsAt=[[] for _ in range(V)], E=len(ed),
                key=(tuple(gv), tuple(ed), ()))


_MEM = {}


def graphs(g, Emax):
    ck = (g, Emax)
    if ck in _MEM: return _MEM[ck]
    path = os.path.join(DATA, "graphs_%d_0_%d.json" % (g, Emax))
    if os.path.exists(path):
        raw = json.load(open(path))
        flat = [((tuple(k[0]), tuple(tuple(e) for e in k[1]), tuple(k[2])), a) for k, a in raw]
    else:
        flat = [x for lev in levels(g, Emax) for x in lev]
        tmp = path + ".tmp%d" % os.getpid()
        json.dump([[[list(k[0]), [list(e) for e in k[1]], list(k[2])], a] for k, a in flat], open(tmp, "w")); os.replace(tmp, path)
    out = [graph_dict(k, a) for k, a in flat]; _MEM[ck] = out
    return out


def legged_graphs_fast(g, n, Emax):
    """the signature of fq3_graphs.legged_graphs, for leg-free graphs only"""
    assert n == 0, "the fast list is for leg-free graphs only"
    return graphs(g, Emax)
