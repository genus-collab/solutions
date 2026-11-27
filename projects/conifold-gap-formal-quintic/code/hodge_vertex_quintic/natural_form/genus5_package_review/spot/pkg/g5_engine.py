# -*- coding: utf-8 -*-
"""g5_engine.py -- the outer graph sum of the formal quintic by edge splitting, with ONE vertex elimination for all outer-edge subsets.

Built on fq4_split.Split (hodge_vertex_quintic/fq4_split.py, copied unchanged into this package; see its docstring for the method):
make_split rebinds the two names fq4_split.legged_graphs / fq4_split.canon to the fast graph list and canonical form of g5_graphs, so
Split.inner_graphs, Split.outer_pos and Split.contract all use the fast keys consistently.
graph_dp(sp, t, maxout): same output as Split.graph(t) ({outer index: [K coefficients of u]}, weight 1/|Aut|), but the inner/outer marks
  of the edges are carried through the vertex elimination (state = (frontier colours, set S of outer edges so far) -> series), an edge
  branching into inner (E_M) and outer (E_out) when its second end is eliminated, only while |S| < min(K - 1, maxout) before the branch
  (maxout = None: |S| <= K - 1 after the branch, i.e. exactly the subsets of Split.graph).  Work: Sum_steps 2^(closed edges) 5^(frontier+1)
  instead of 2^E Sum_steps 5^(frontier+1).  Validated in the study on all 379 genus-4 inner graphs and 18 genus-5 inner graphs (equal to
  Split.graph), and again by genus5_run.py --validate.
  With maxout = m the result is EXACT for the vertex graph and for every outer graph with at most m edges; other outer graphs are
  incomplete and must not be used.
graph_id(sp, t): the point-frame term of the inner graph t (variant 'id'), a list of K coefficients.
"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from fractions import Fraction as Fr

POPC = [bin(i).count("1") for i in range(1 << 16)]


def make_split(p, g, ab=(0, 0), K=None, variant='nat'):
    import fq4_split, g5_graphs
    fq4_split.legged_graphs = g5_graphs.legged_graphs_fast; fq4_split.canon = g5_graphs.canon_fast
    return fq4_split.Split(p, g, (Fr(ab[0]), Fr(ab[1])), K, variant)


def combo_dp(sp, steps, vecs, inner, outer, cap):
    mul = sp.mul; add = sp.add
    cur = {(): {0: [1] + [0] * (sp.K - 1)}}
    for v, conn, keep, loops in steps:
        new = {}; vec = vecs[v]
        for key, sub in cur.items():
            for c in range(5):
                vc = vec[c]
                if vc is None: continue
                facs = [(1 << e, inner[e][c][c], outer[e][c][c]) for e in loops]
                for e, pos, flip in conn:
                    k = key[pos]
                    facs.append((1 << e, inner[e][c][k], outer[e][c][k]) if flip else (1 << e, inner[e][k][c], outer[e][k][c]))
                items = []
                for S, val in sub.items():
                    x = mul(val, vc)
                    if x is not None: items.append((S, x))
                for bit, fi, fo in facs:
                    nxt = []
                    for S, x in items:
                        a = mul(x, fi)
                        if a is not None: nxt.append((S, a))
                        if POPC[S] < cap:
                            b = mul(x, fo)
                            if b is not None: nxt.append((S | bit, b))
                    items = nxt
                if not items: continue
                full = key + (c,); nk = tuple(full[k] for k in keep); dst = new.get(nk)
                if dst is None: dst = {}; new[nk] = dst
                for S, x in items:
                    y = dst.get(S); dst[S] = x if y is None else add(y, x)
        cur = new
        if not cur: return {}
    return cur.get((), {})


def decorations(sp, G):
    opts = []
    for v in range(G['V']):
        tab = sp.vertex_table(G['gv'][v], len(G['hsAt'][v])); lst = []
        for hp, vec in tab.items():
            for perm in sorted(set(itertools.permutations(hp))): lst.append((perm, vec))
        opts.append(lst)
    return opts


def graph_dp(sp, t, maxout=None):
    assert sp.variant != 'id'
    G = sp.inner_graphs[t]; V = G['V']; E = len(G['ed']); hsAt = G['hsAt']; K = sp.K; p = sp.p
    cap = K - 1 if maxout is None else min(K - 1, maxout)
    opts = decorations(sp, G)
    if any(not lst for lst in opts): return {}
    pos = {S: sp.outer_pos[sp.contract(G, S)] for S in range(1 << E) if POPC[S] <= cap}
    steps = sp.plan(G); acc = {}
    for combo in itertools.product(*opts):
        hp = [0] * (2 * E)
        for v in range(V):
            for r, h in enumerate(hsAt[v]): hp[h] = combo[v][0][r]
        keys = [(hp[2 * e], hp[2 * e + 1]) for e in range(E)]
        res = combo_dp(sp, steps, [combo[v][1] for v in range(V)], [sp.EM[k] for k in keys], [sp.EO[k] for k in keys], cap)
        for S, z in res.items():
            if z is None: continue
            q = pos[S]; a = acc.get(q)
            if a is None: acc[q] = list(z)
            else:
                for k in range(K): a[k] += z[k]
    ia = sp.R.inv(G['aut'])
    return {q: [x * ia % p for x in a] for q, a in acc.items()}


def graph_id(sp, t):
    assert sp.variant == 'id'
    r = sp.graph(t)
    return r.get(t, [0] * sp.K)


def tau_identity_ok(sp):
    """tau(u) * Sum_n 4 u^n/(5n+4) = 1 mod u^K (the rescaling used for E_M)"""
    from fq4_split import tau_pred
    R = sp.R; K = sp.K; tau = [R.of(x) for x in tau_pred(K)]; inv = [R.of(Fr(4, 5 * n + 4)) for n in range(K)]
    prod = [sum(tau[i] * inv[n - i] for i in range(n + 1)) % sp.p for n in range(K)]
    return prod == [1] + [0] * (K - 1)
