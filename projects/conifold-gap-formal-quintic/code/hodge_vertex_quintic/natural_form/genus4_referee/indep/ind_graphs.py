# -*- coding: utf-8 -*-
"""ind_graphs.py -- INDEPENDENT RE-IMPLEMENTATION (genus-4 referee): leg-free stable graphs, automorphisms, contractions.

Source of the graphs: the Mathematica export  <project>/../cluster_genus4_hodge/StableGraph-g.bin  parsed here from the documented
layout (first and last byte dropped, bytes 13 dropped; per graph: V, g_1..g_V, V further bytes (not used), l_1..l_V = loops, then
the strict upper triangle of the adjacency matrix row by row).  Every graph is checked: connected, stable, genus g, and the file is
consumed exactly.  The list is then compared, as a set of canonical forms, with the 'outer' list stored in the project's result
file (which comes from a different generator: one-edge degenerations, fq3_graphs.legged_graphs); the position in that list is the
"outer index" used to compare tables.
Own code: canon (minimum over vertex permutations), n_aut (depth-first count of the half-edge bijections that commute with the
edge involution and respect vertices and genera: NOT the product formula of fq_graphs.aut_count), contract (union-find).
Graph dict: V, E, gv, ed [(u, w)], hs[v] = half-edges at v (half-edge 2e is the u-end of edge e, 2e+1 its w-end), aut, key.
"""
import itertools, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
BIN = os.path.normpath(os.path.join(PROJ, "graph_data"))


def canon(gv, ed):
    V = len(gv); best = None
    for sig in itertools.permutations(range(V)):
        g2 = [0] * V
        for v in range(V): g2[sig[v]] = gv[v]
        e2 = sorted((min(sig[u], sig[w]), max(sig[u], sig[w])) for u, w in ed)
        key = (tuple(g2), tuple(e2))
        if best is None or key < best: best = key
    return best


def n_aut(gv, ed):
    """number of automorphisms: bijections phi of the half-edges with phi(partner(h)) = partner(phi(h)) inducing a genus-preserving
    bijection of the vertices (depth-first search, vertex by vertex along the edges; the graph is connected)"""
    V = len(gv); E = len(ed)
    if E == 0: return 1
    vert = [ed[h // 2][h % 2] for h in range(2 * E)]
    at = [[h for h in range(2 * E) if vert[h] == v] for v in range(V)]
    phi = {}; sig = {}; used_h = set(); used_v = set()

    def next_free():
        for v in sorted(sig):
            for h in at[v]:
                if h not in phi: return h
        return None

    def rec():
        if len(phi) == 2 * E: return 1 if len(sig) == V else 0
        h = next_free()
        if h is None: return 0
        v = vert[h]; w = vert[h ^ 1]; cnt = 0
        for h2 in at[sig[v]]:
            if h2 in used_h or (h2 ^ 1) in used_h: continue
            w2 = vert[h2 ^ 1]; new_w = False
            if w in sig:
                if sig[w] != w2: continue
            else:
                if w2 in used_v or gv[w2] != gv[w] or len(at[w2]) != len(at[w]): continue
                new_w = True
            phi[h] = h2; phi[h ^ 1] = h2 ^ 1; used_h.add(h2); used_h.add(h2 ^ 1)
            if new_w: sig[w] = w2; used_v.add(w2)
            cnt += rec()
            del phi[h]; del phi[h ^ 1]; used_h.discard(h2); used_h.discard(h2 ^ 1)
            if new_w: del sig[w]; used_v.discard(w2)
        return cnt

    total = 0
    for v0 in range(V):
        if gv[v0] != gv[0] or len(at[v0]) != len(at[0]): continue
        sig[0] = v0; used_v.add(v0)
        total += rec()
        del sig[0]; used_v.discard(v0)
    return total


def make(gv, ed):
    gv = tuple(gv); ed = [tuple(e) for e in ed]; V = len(gv); E = len(ed)
    hs = [[2 * e + s for e in range(E) for s in (0, 1) if ed[e][s] == v] for v in range(V)]
    return dict(V=V, E=E, gv=gv, ed=ed, hs=hs, aut=n_aut(gv, ed), key=canon(gv, ed))


def parse_bin(g):
    raw = [x for x in open(os.path.join(BIN, "StableGraph-%d.bin" % g), "rb").read()[1:-1] if x != 13]
    out = []; i = 0
    while i < len(raw):
        V = raw[i]; assert 1 <= V <= 2 * g - 2, "bad vertex count"
        gv = raw[i + 1:i + 1 + V]; loops = raw[i + 1 + 2 * V:i + 1 + 3 * V]; i += 1 + 3 * V
        ed = []
        for v in range(V): ed += [(v, v)] * loops[v]
        for v in range(V):
            for w in range(v + 1, V):
                ed += [(v, w)] * raw[i]; i += 1
        E = len(ed); val = [0] * V
        for u, w in ed: val[u] += 1; val[w] += 1
        assert sum(gv) + E - V + 1 == g, "genus"
        assert all(2 * gv[v] - 2 + val[v] > 0 for v in range(V)), "stability"
        comp = list(range(V))
        for _ in range(V):
            for u, w in ed:
                m = min(comp[u], comp[w]); comp[u] = comp[w] = m
        assert len(set(comp)) == 1, "connected"
        out.append(make(gv, ed))
    assert i == len(raw), "file not consumed exactly"
    keys = [G['key'] for G in out]; assert len(set(keys)) == len(keys), "duplicate graphs"
    return out


def contract(G, S):
    """canonical key of G with the edges NOT in the bitmask S contracted (the edges in S survive)"""
    V = G['V']; ed = G['ed']; rep = list(range(V))
    def root(x):
        while rep[x] != x: x = rep[x]
        return x
    for e, (u, w) in enumerate(ed):
        if not (S >> e) & 1:
            a, b = root(u), root(w)
            if a != b: rep[max(a, b)] = min(a, b)
    roots = sorted(set(root(v) for v in range(V))); num = {r: i for i, r in enumerate(roots)}
    gsum = [0] * len(roots); nver = [0] * len(roots); nin = [0] * len(roots); keep = []
    for v in range(V): c = num[root(v)]; gsum[c] += G['gv'][v]; nver[c] += 1
    for e, (u, w) in enumerate(ed):
        if (S >> e) & 1: keep.append((num[root(u)], num[root(w)]))
        else: nin[num[root(u)]] += 1
    return canon([gsum[c] + nin[c] - nver[c] + 1 for c in range(len(roots))], keep)


def reference_outer(g, prime="89"):
    """the 'outer' list of the project's result file: list of dicts(E, gv, ed, aut) -> list of (key, aut)"""
    if g == 4: d = json.load(open(os.path.join(PROJ, "g4_split_g4_%s_id.json" % prime)))['outer']
    else: return None
    return [(canon(o['gv'], [tuple(e) for e in o['ed']]), o['aut']) for o in d]


def graphs(g):
    """the stable graphs of genus g in the order of the reference outer list when it exists (g = 4), else in the Mathematica order
    (which is the 1-based 'idx' order of the stored genus-2/3 tables: position i here <-> idx i+1)"""
    L = parse_bin(g); ref = reference_outer(g)
    if ref is None: return L
    by = {G['key']: G for G in L}
    assert set(by) == set(k for k, a in ref) and len(ref) == len(L), "Mathematica list != project list"
    out = [by[k] for k, a in ref]
    assert all(G['aut'] == a for G, (k, a) in zip(out, ref)), "automorphism counts differ from the project's"
    return out


if __name__ == "__main__":
    import math
    for g in (2, 3, 4):
        L = parse_bin(g); print("genus %d: %d graphs; by number of edges %s" % (g, len(L), [sum(1 for G in L if G['E'] == E) for E in range(3 * g - 2)]))
        assert len(L) == {2: 7, 3: 42, 4: 379}[g]
        # the product formula as a cross-check of the depth-first count (own implementation of the formula)
        for G in L:
            V = G['V']; A = [[0] * V for _ in range(V)]
            for u, w in G['ed']:
                if u == w: A[u][u] += 1
                else: A[u][w] += 1; A[w][u] += 1
            n = sum(1 for s in itertools.permutations(range(V)) if all(G['gv'][s[v]] == G['gv'][v] for v in range(V))
                    and all(A[s[i]][s[j]] == A[i][j] for i in range(V) for j in range(V)))
            for i in range(V):
                n *= 2 ** A[i][i] * math.factorial(A[i][i])
                for j in range(i + 1, V): n *= math.factorial(A[i][j])
            assert n == G['aut'], (G['gv'], G['ed'], n, G['aut'])
        print("   depth-first automorphism counts = product formula for all graphs; max |Aut| = %d" % max(G['aut'] for G in L))
    L4 = graphs(4)
    print("genus 4: Mathematica list = project outer list as sets of canonical forms (379), all |Aut| equal the project's")
    G = L4[0]; assert G['E'] == 0
    # contraction sanity: contracting everything gives the vertex graph; contracting nothing gives the graph itself
    for G in L4:
        assert contract(G, 0) == L4[0]['key'] and contract(G, (1 << G['E']) - 1) == G['key']
    print("ind_graphs ok")
