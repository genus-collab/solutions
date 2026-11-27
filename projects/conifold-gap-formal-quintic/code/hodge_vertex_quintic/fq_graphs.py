# -*- coding: utf-8 -*-
"""fq_graphs.py -- stable graphs (StableGraph-g.bin loader), automorphisms, decorations and their orbits.

Adapted from hodge_vertex_conifold/cy_torus_explore/exact_test/cylib.py (load_graphs, edges_of, autos_of,
aut_count, compositions, tuples_sum_le, hodge_poly_nonzero, vertex_options, Decorated); imports only fq_hodge.

Graph data: formalquintic/hodge_vertex_quintic/graph_data/StableGraph-g.bin (Mathematica export, one byte per
entry): for each graph  K, g_1..g_K, m_1..m_K (valences), l_1..l_K (loops), then the strict upper triangle of the
adjacency matrix row by row.  load_graphs(g) -> list of (gv, m, A) with A[i][i] = number of loops at i.

Decorated(g, hodge=True): for every graph, the labelled decorations (psi-exponent per half-edge, dilaton-leaf
exponents d >= 2 per vertex) with non-zero vertex integrals, grouped into Aut-orbits (orbs: canonical key ->
list of labelled decorations).  Vertex integrals with the Hodge class are non-zero only if the Hodge polynomial is
(hodge_poly_nonzero); with hodge=False the point class (WK) is used.  Conventions as in cylib.Decorated.
"""
import itertools, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fq_hodge import wk, hi, dmax, partitions_odd

DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "graph_data")


def load_graphs(g):
    b = open(os.path.join(DATADIR, "StableGraph-%d.bin" % g), "rb").read()
    raw = [x for x in b[1:-1] if x != 13]
    out = []; i = 0
    while i < len(raw):
        K = raw[i]; gv = raw[i + 1:i + 1 + K]; m = raw[i + 1 + K:i + 1 + 2 * K]; l = raw[i + 1 + 2 * K:i + 1 + 3 * K]
        i2 = i + 1 + 3 * K
        A = [[0] * K for _ in range(K)]
        for j in range(K): A[j][j] = l[j]
        for j in range(K):
            for ii in range(j + 1, K):
                A[ii][j] = A[j][ii] = raw[i2]; i2 += 1
        out.append((tuple(gv), tuple(m), A))
        i = i2
    return out


def edges_of(A):
    K = len(A); ed = []
    for i in range(K):
        ed += [(i, i)] * A[i][i]
    for i in range(K):
        for j in range(i + 1, K):
            ed += [(i, j)] * A[i][j]
    return ed


def autos_of(gv, A):
    """automorphisms as (vertex perm, edge map, loop flips)"""
    V = len(gv); ed = edges_of(A); out = []
    for sig in itertools.permutations(range(V)):
        if any(gv[sig[v]] != gv[v] for v in range(V)): continue
        if any(A[sig[i]][sig[j]] != A[i][j] for i in range(V) for j in range(V)): continue
        img = [tuple(sorted((sig[u], sig[v]))) for u, v in ed]
        classes = {}
        for e, im in enumerate(img): classes.setdefault(im, []).append(e)
        choices = []
        for im, cl in classes.items():
            targets = [e for e, x in enumerate(ed) if x == im]
            choices.append([dict(zip(cl, p)) for p in itertools.permutations(targets)])
        for ch in itertools.product(*choices):
            emap = {}
            for d in ch: emap.update(d)
            flipsets = [(0, 1) if ed[e][0] == ed[e][1] else (0,) for e in range(len(ed))]
            for fl in itertools.product(*flipsets):
                out.append((sig, emap, fl))
    return out


def aut_count(gv, A):
    V = len(gv); n = 0
    for sig in itertools.permutations(range(V)):
        if any(gv[sig[v]] != gv[v] for v in range(V)): continue
        if any(A[sig[i]][sig[j]] != A[i][j] for i in range(V) for j in range(V)): continue
        n += 1
    r = n
    for i in range(V):
        r *= 2 ** A[i][i] * math.factorial(A[i][i])
        for j in range(i + 1, V): r *= math.factorial(A[i][j])
    return r


def compositions(s, nd, minp):
    """ordered tuples of nd parts >= minp with sum s"""
    if nd == 0:
        if s == 0: yield ()
        return
    for f in range(minp, s - minp * (nd - 1) + 1):
        for r in compositions(s - f, nd - 1, minp):
            yield (f,) + r


def tuples_sum_le(n, budget):
    if n == 0:
        yield ()
        return
    for f in range(budget + 1):
        for r in tuples_sum_le(n - 1, budget - f):
            yield (f,) + r


_nonzero_memo = {}
def hodge_poly_nonzero(g, d):
    """is Int Hodge_w prod psi^d a non-zero polynomial of the s_k?"""
    d = tuple(sorted(d)); key = (g, d)
    r = _nonzero_memo.get(key)
    if r is None:
        n = len(d); D = 3 * g - 3 + n - sum(d)
        if D < 0 or 2 * g - 2 + n <= 0: r = False
        elif D > dmax(g): r = False
        else: r = any(hi(g, d, K) != 0 for K in partitions_odd(D, max(2 * g - 1, 1)))
        _nonzero_memo[key] = r
    return r


def vertex_options(g, hs, hodge):
    """list of (hp tuple, dils tuple) for a vertex of genus g with hs half-edges; hodge=False: point vertices"""
    res = []
    for nd in range(0, max(0, 3 * g - 3 + hs) + 1):
        dim = 3 * g - 3 + hs + nd
        if dim < 2 * nd: continue
        dilss = [()] if nd == 0 else [c for s in range(2 * nd, dim + 1) for c in compositions(s, nd, 2)]
        for dils in dilss:
            budget = dim - sum(dils)
            for hp in tuples_sum_le(hs, budget):
                if hodge:
                    if hodge_poly_nonzero(g, hp + dils): res.append((hp, dils))
                else:
                    if sum(hp) == budget and wk(g, hp + dils) != 0: res.append((hp, dils))
    return res


class Decorated(object):
    """orbits of labelled decorations for all stable graphs of genus g (see module docstring)"""
    def __init__(self, g, hodge=True):
        self.g = g; self.hodge = hodge
        self.graphs = []
        for gi, (gv, m, A) in enumerate(load_graphs(g)):
            V = len(gv); ed = edges_of(A)
            H = [(e, s) for e in range(len(ed)) for s in (0, 1)]
            hvtx = [ed[e][s] for e, s in H]
            hsAt = [[i for i, x in enumerate(hvtx) if x == v] for v in range(V)]
            opts = [vertex_options(gv[v], len(hsAt[v]), hodge) for v in range(V)]
            auts = autos_of(gv, A)
            assert len(auts) == aut_count(gv, A)
            orbs = {}; nlab = 0
            for dec in itertools.product(*opts):
                nlab += 1
                psi = [[0, 0] for _ in ed]
                for v in range(V):
                    for r, h in enumerate(hsAt[v]):
                        e, s = H[h]; psi[e][s] = dec[v][0][r]
                dil = [tuple(sorted(dec[v][1])) for v in range(V)]
                best = None
                for sig, emap, fl in auts:
                    psi2 = [[0, 0] for _ in ed]
                    for e in range(len(ed)):
                        e2 = emap[e]
                        tgt = (sig[ed[e][0]], sig[ed[e][1]])
                        ordr = (0, 1) if tgt == ed[e2] else (1, 0)
                        if ed[e2][0] == ed[e2][1] and fl[e] == 1: ordr = (ordr[1], ordr[0])
                        psi2[e2][ordr[0]] = psi[e][0]; psi2[e2][ordr[1]] = psi[e][1]
                    dil2 = [None] * V
                    for v in range(V): dil2[sig[v]] = dil[v]
                    key = (tuple(tuple(x) for x in psi2), tuple(dil2))
                    if best is None or key < best: best = key
                orbs.setdefault(best, []).append(dec)
            self.graphs.append(dict(idx=gi + 1, gv=gv, A=A, ed=ed, H=H, hsAt=hsAt, aut=len(auts), orbs=orbs,
                                    nlab=nlab, V=V, bound=sum(max(2 * h - 2, 0) for h in gv)))

    def counts(self):
        return dict(graphs=len(self.graphs), labelled=sum(G['nlab'] for G in self.graphs),
                    orbits=sum(len(G['orbs']) for G in self.graphs))


if __name__ == "__main__":
    gs = load_graphs(2)
    assert len(gs) == 7, len(gs)
    for gv, m, A in gs:
        ed = edges_of(A)
        assert sum(gv) + len(ed) - len(gv) + 1 == 2
    D = Decorated(2, hodge=True)
    print("genus 2 (Hodge vertices):", D.counts())
    D0 = Decorated(2, hodge=False)
    print("genus 2 (point vertices):", D0.counts())
    for G in D.graphs:
        print(G['idx'], 'gv', G['gv'], 'edges', G['ed'], '|Aut|', G['aut'], 'orbits', len(G['orbs']), 'bound', G['bound'])
    print("fq_graphs ok")
