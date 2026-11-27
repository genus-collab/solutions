# -*- coding: utf-8 -*-
"""g5_lib.py -- fast canonical form / generator of leg-free stable graphs, and a patched fq4_split.Split for genus 5 (feasibility study).

Nothing in the parent folders is modified: fq4_split is imported and two NAMES of that module (legged_graphs, canon) are rebound at
run time to the functions of this file.

canon_fast(gv, ed, legs): canonical key (gv', sorted edges', legs') of the same SHAPE as fq3_graphs.canon, but minimised only over the
    leaves of an individualisation-refinement search (colour refinement by (genus, loops, legs, multiset of (neighbour colour,
    multiplicity)); in the first non-singleton cell every vertex is individualised in turn).  The set of leaves is isomorphism
    invariant, so the minimum is a canonical form; it is in general NOT the same representative as fq3_graphs.canon (minimum over all
    V! orders), which costs 8! = 40320 orders for the genus-5 trivalent graphs.  The number of leaves attaining the minimum is the
    number of vertex automorphisms, so |Aut| = (that number) prod_loops 2^l l! prod_multi-edges m!  (same value as fq_graphs.aut_count).
fast_graphs(g, n, Emax): drop-in replacement for fq3_graphs.legged_graphs (same dict fields, sorted by (E, key)); cached on disk in
    this folder (g5_graphs_<g>_<n>_<Emax>.json).
make_split(p, g, ab, K, variant): fq4_split.Split built on fast_graphs / canon_fast.
SplitTools: cost sampling of one inner graph (random (decoration, outer-edge subset) pairs, timed) without running the whole graph.
"""
import os, sys, json, math, time, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
if ROOT not in sys.path: sys.path.insert(0, ROOT)


def _refine(V, adj, col):
    ncl = len(set(col))
    while True:
        sig = [(col[v], tuple(sorted((col[w], m) for w, m in adj[v]))) for v in range(V)]
        ranks = {s: i for i, s in enumerate(sorted(set(sig)))}
        col = [ranks[s] for s in sig]
        if len(ranks) == ncl: return col
        ncl = len(ranks)


def canon_count(gv, ed, legs):
    """(canonical key, number of vertex automorphisms fixing the legs)"""
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


def fast_levels(g, n, Emax, log=None, deadline=None):
    """list over E = 0..Emax of sorted lists of (canonical key, |Aut|); raw degenerations are deduplicated before canonicalisation"""
    from fq3_graphs import degenerations
    assert 2 * g - 2 + n > 0
    k0, c0 = canon_count((g,), [], tuple([0] * n)); levels = [[(k0, aut_of(k0, c0))]]; t0 = time.time()
    for E in range(1, Emax + 1):
        t1 = time.time(); new = {}; seen = set(); nraw = 0
        for key, _ in levels[-1]:
            for gv2, ed2, legs2 in degenerations(key[0], list(key[1]), key[2]):
                nraw += 1; raw = (gv2, tuple(sorted(ed2)), legs2)
                if raw in seen: continue
                seen.add(raw); k, c = canon_count(gv2, ed2, legs2)
                if k not in new: new[k] = aut_of(k, c)
            if deadline is not None and time.time() - t0 > deadline: raise RuntimeError("deadline in level %d" % E)
        levels.append(sorted(new.items()))
        if log: log("   E = %2d: %5d graphs (%d degenerations, %d distinct raw, %.1f s; cumulative %.1f s)" % (E, len(new), nraw, len(seen), time.time() - t1, time.time() - t0))
    return levels


def graph_dict(key, aut, n):
    gv, ed, legs = key; V = len(gv); ed = [tuple(e) for e in ed]
    H = [(e, s) for e in range(len(ed)) for s in (0, 1)]
    hvtx = [ed[e][s] for e, s in H]
    hsAt = [[h for h, x in enumerate(hvtx) if x == v] for v in range(V)]
    legsAt = [[l for l in range(n) if legs[l] == v] for v in range(V)]
    return dict(V=V, gv=tuple(gv), ed=ed, legs=tuple(legs), aut=aut, H=H, hsAt=hsAt, legsAt=legsAt, E=len(ed), key=(tuple(gv), tuple(ed), tuple(legs)))


_MEM = {}


def fast_graphs(g, n, Emax, log=None):
    """drop-in for fq3_graphs.legged_graphs(g, n, Emax); disk cache in this folder"""
    ck = (g, n, Emax)
    if ck in _MEM: return _MEM[ck]
    path = os.path.join(HERE, "g5_graphs_%d_%d_%d.json" % ck)
    if os.path.exists(path):
        raw = json.load(open(path))
        flat = [((tuple(k[0]), tuple(tuple(e) for e in k[1]), tuple(k[2])), a) for k, a in raw]
    else:
        flat = [x for lev in fast_levels(g, n, Emax, log) for x in lev]
        json.dump([[[list(k[0]), [list(e) for e in k[1]], list(k[2])], a] for k, a in flat], open(path + ".tmp%d" % os.getpid(), "w"))
        os.replace(path + ".tmp%d" % os.getpid(), path)
    out = [graph_dict(k, a, n) for k, a in flat]; _MEM[ck] = out
    return out


def make_split(p, g, ab=(0, 0), K=None, variant='nat'):
    import fq4_split
    fq4_split.legged_graphs = fast_graphs; fq4_split.canon = canon_fast
    return fq4_split.Split(p, g, ab, K, variant)


class SplitTools(object):
    """cost sampling for one inner graph of a Split object sp (same loop as Split.graph, on random (decoration, subset) pairs)"""
    def __init__(self, sp): self.sp = sp

    def prepare(self, t):
        sp = self.sp; G = sp.inner_graphs[t]; V = G['V']; E = len(G['ed']); hsAt = G['hsAt']; t0 = time.time(); opts = []
        for v in range(V):
            tab = sp.vertex_table(G['gv'][v], len(hsAt[v])); lst = []
            for hp, vec in tab.items():
                for perm in set(itertools.permutations(hp)): lst.append((perm, vec))
            opts.append(lst)
        t1 = time.time()
        if sp.variant == 'id': subsets = [((1 << E) - 1, t)]
        else: subsets = [(S, sp.outer_pos[sp.contract(G, S)]) for S in range(1 << E) if bin(S).count("1") < sp.K]
        t2 = time.time(); steps = sp.plan(G)
        return dict(G=G, opts=opts, subsets=subsets, steps=steps, t_opts=t1 - t0, t_subsets=t2 - t1)

    def sample(self, t, budget=1.0, min_samples=4, max_samples=400, seed=0):
        """dict with ncombo, nsub, samples, skipped (valuation pruning), mean seconds per pair, estimated seconds for the whole graph"""
        sp = self.sp; P = self.prepare(t); G = P['G']; V = G['V']; E = len(G['ed']); hsAt = G['hsAt']; rng = random.Random(seed * 100003 + t)
        ncombo = 1
        for lst in P['opts']: ncombo *= len(lst)
        nsub = len(P['subsets']); out = dict(t=t, E=E, V=V, gv=list(G['gv']), val=[len(h) for h in hsAt], ncombo=ncombo, nsub=nsub,
                                             t_opts=P['t_opts'], t_subsets=P['t_subsets'])
        if ncombo == 0: out.update(samples=0, skipped=0, mean=0.0, est=P['t_opts'] + P['t_subsets'], tmax=0.0); return out
        n = 0; skipped = 0; tsum = 0.0; tmax = 0.0; t0 = time.time()
        while n < max_samples and (n < min_samples or time.time() - t0 < budget):
            combo = [rng.choice(lst) for lst in P['opts']]; S, pos = rng.choice(P['subsets']); hp = [0] * (2 * E)
            for v in range(V):
                for r, h in enumerate(hsAt[v]): hp[h] = combo[v][0][r]
            keys = [(hp[2 * e], hp[2 * e + 1]) for e in range(E)]; n += 1
            if sum(sp.valO[keys[e]] for e in range(E) if (S >> e) & 1) >= sp.K: skipped += 1; continue
            vecs = [combo[v][1] for v in range(V)]
            mats = [sp.EO[keys[e]] if (S >> e) & 1 else sp.EM[keys[e]] for e in range(E)]
            ta = time.time(); sp.colour_sum(P['steps'], vecs, mats); dt = time.time() - ta; tsum += dt; tmax = max(tmax, dt)
        mean = tsum / n
        out.update(samples=n, skipped=skipped, mean=mean, tmax=tmax, est=P['t_opts'] + P['t_subsets'] + ncombo * nsub * mean)
        return out
