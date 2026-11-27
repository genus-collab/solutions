# -*- coding: utf-8 -*-
"""ms_graphs.py -- stable graphs with labelled legs, as fq3_graphs.legged_graphs (same one-edge degeneration step, imported read-only),
but with a canonical form by colour refinement (vertex invariants: genus, loops, degree, the labelled legs carried; refined by the
multisets of neighbour colours with edge multiplicities) and permutations only inside the final colour cells; |Aut| (vertex
permutations fixing every vertex that carries a leg, times 2^loops loops! and multi-edge factorials) counted the same way.
Output dicts have the fields used by ms_engine.graph_amp: V, gv, ed, legs, aut, H, hsAt, legsAt, E, key.
Validation (ms_graphs.py run as a script): the same number of graphs per edge level and the same Sum 1/|Aut| per level as
fq3_graphs.legged_graphs for several small (g, n).  The canonical KEYS differ from fq3_graphs (different minimisation domain).
"""
import itertools, math, os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ)
from fq3_graphs import degenerations


def _adj(V, ed):
    A = [[0] * V for _ in range(V)]
    for a, b in ed:
        if a == b: A[a][a] += 1
        else: A[a][b] += 1; A[b][a] += 1
    return A


def _cells(gv, A, legs):
    V = len(gv); n = len(legs)
    col = [(gv[v], A[v][v], sum(A[v][w] for w in range(V) if w != v), tuple(l for l in range(n) if legs[l] == v)) for v in range(V)]
    ranks = None
    while True:
        keys = sorted(set(col)); r = [keys.index(c) for c in col]
        if ranks is not None and len(set(r)) == len(set(ranks)): ranks = r; break
        ranks = r
        col = [(ranks[v], tuple(sorted((ranks[w], A[v][w]) for w in range(V) if w != v and A[v][w]))) for v in range(V)]
    cells = {}
    for v in range(V): cells.setdefault(ranks[v], []).append(v)
    return [cells[k] for k in sorted(cells)]


def _perms(cells):
    """all vertex maps sig (old -> new position) sending cell i onto the block of positions of cell i"""
    blocks = []; pos = 0
    for c in cells: blocks.append(list(range(pos, pos + len(c)))); pos += len(c)
    for choice in itertools.product(*[itertools.permutations(b) for b in blocks]):
        sig = {}
        for c, img in zip(cells, choice):
            for v, x in zip(c, img): sig[v] = x
        yield sig


def canon2(gv, ed, legs):
    V = len(gv); A = _adj(V, ed); cells = _cells(gv, A, legs); best = None
    for sig in _perms(cells):
        gv2 = [0] * V
        for v in range(V): gv2[sig[v]] = gv[v]
        ed2 = tuple(sorted(tuple(sorted((sig[a], sig[b]))) for a, b in ed))
        key = (tuple(gv2), ed2, tuple(sig[l] for l in legs))
        if best is None or key < best: best = key
    return best


def aut2(gv, ed, legs):
    V = len(gv); A = _adj(V, ed); cells = _cells(gv, A, legs); legv = set(legs); cnt = 0
    order = [v for c in cells for v in c]
    for sig in _perms(cells):
        s = [order[sig[v]] for v in range(V)]         # an automorphism candidate v -> s[v] inside the cells
        if any(s[v] != v for v in legv): continue
        if any(gv[s[v]] != gv[v] for v in range(V)): continue
        if any(A[s[i]][s[j]] != A[i][j] for i in range(V) for j in range(V)): continue
        cnt += 1
    r = cnt
    for i in range(V):
        r *= 2 ** A[i][i] * math.factorial(A[i][i])
        for j in range(i + 1, V): r *= math.factorial(A[i][j])
    return r


def legged_keys2(g, n, Emax):
    levels = [[canon2((g,), [], tuple([0] * n))]]
    for E in range(1, Emax + 1):
        new = set()
        for key in levels[-1]:
            for H in degenerations(key[0], list(key[1]), key[2]): new.add(canon2(*H))
        levels.append(sorted(new))
    return levels


def graph_dict2(key, n):
    gv, ed, legs = key; V = len(gv); ed = list(ed)
    H = [(e, s) for e in range(len(ed)) for s in (0, 1)]
    hvtx = [ed[e][s] for e, s in H]
    hsAt = [[h for h, x in enumerate(hvtx) if x == v] for v in range(V)]
    legsAt = [[l for l in range(n) if legs[l] == v] for v in range(V)]
    return dict(V=V, gv=gv, ed=ed, legs=legs, aut=aut2(gv, ed, legs), H=H, hsAt=hsAt, legsAt=legsAt, E=len(ed), key=key)


def legged_graphs2(g, n, Emax):
    return [graph_dict2(key, n) for lev in legged_keys2(g, n, Emax) for key in lev]


if __name__ == "__main__":
    from fq3_graphs import legged_graphs
    for g, n in [(0, 4), (0, 5), (1, 2), (1, 3), (2, 1), (2, 2), (1, 4), (2, 3), (3, 1), (3, 2)]:
        t0 = time.time(); A = legged_graphs2(g, n, 3 * g - 3 + n); t1 = time.time(); B = legged_graphs(g, n, 3 * g - 3 + n); t2 = time.time()
        lev = lambda L: [(E, len([x for x in L if x['E'] == E]), sum(Fr(1, x['aut']) for x in L if x['E'] == E)) for E in range(3 * g - 2 + n)]
        print("(%d,%d): %d graphs (%.1f s) vs fq3_graphs %d (%.1f s): per-level counts and Sum 1/|Aut| %s"
              % (g, n, len(A), t1 - t0, len(B), t2 - t1, "AGREE" if lev(A) == lev(B) else "DIFFER"), flush=True)
