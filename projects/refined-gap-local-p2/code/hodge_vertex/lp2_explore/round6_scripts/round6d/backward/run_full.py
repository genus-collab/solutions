# Full stable-graph sum by degeneration enumeration. Per-colouring polar test. Evidence only.
# usage: run_full.py e1 e2 g N mode nproc [chunk]
# mode 'total' : all graphs, all colourings; t_c^0 of the difference
# mode 'hot'   : only graphs with a vertex of genus >= 2, only colourings with such a vertex at colour 0,
#                one colouring out of each Galois pair (1 <-> 2)
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import time
from itertools import product, permutations, combinations
from fractions import Fraction as Fr
from math import factorial
from multiprocessing import Pool
from gs import *


def canon(genera, edges):
    nv = len(genera)
    # invariants
    val = [0] * nv
    loops = [0] * nv
    mult = {}
    for u, v in edges:
        a, b = min(u, v), max(u, v)
        mult[(a, b)] = mult.get((a, b), 0) + 1
        val[u] += 1
        val[v] += 1
        if u == v:
            loops[u] += 1
    nbr = [[] for _ in range(nv)]
    for (a, b), m in mult.items():
        if a != b:
            nbr[a].append((genera[b], val[b], m))
            nbr[b].append((genera[a], val[a], m))
    inv_ = [(genera[i], val[i], loops[i], tuple(sorted(nbr[i]))) for i in range(nv)]
    classes = {}
    for i in range(nv):
        classes.setdefault(inv_[i], []).append(i)
    keys = sorted(classes, reverse=True)
    groups = [classes[k] for k in keys]
    best = None
    count_best = 0
    import itertools
    naut = 0
    identity_key = None
    for pp in itertools.product(*[permutations(g) for g in groups]):
        order = [x for p in pp for x in p]
        pos = {v: i for i, v in enumerate(order)}
        e2 = sorted((min(pos[u], pos[v]), max(pos[u], pos[v])) for u, v in edges)
        key = tuple(e2)
        if best is None or key < best:
            best = key
            count_best = 1
        elif key == best:
            count_best += 1
    g2 = tuple(genera[v] for grp in groups for v in grp)
    # vertex automorphisms = number of orderings attaining the canonical key
    autv = count_best
    m2 = {}
    for e in best:
        m2[e] = m2.get(e, 0) + 1
    for (a, b), m in m2.items():
        if a == b:
            autv *= (2 ** m) * factorial(m)
        else:
            autv *= factorial(m)
    return (g2, best), autv


def degenerations(genera, edges):
    nv = len(genera)
    he = [[] for _ in range(nv)]
    for k, (u, v) in enumerate(edges):
        he[u].append((k, 0))
        he[v].append((k, 1))
    for v in range(nv):
        gv = genera[v]
        if gv >= 1:
            g2 = list(genera)
            g2[v] = gv - 1
            yield tuple(g2), list(edges) + [(v, v)]
        H = he[v]
        n = len(H)
        for mask in range(2 ** n):
            for g1 in range(gv + 1):
                g2_ = gv - g1
                n1 = bin(mask).count('1') + 1
                n2 = n - n1 + 1 + 1
                if 2 * g1 - 2 + n1 <= 0 or 2 * g2_ - 2 + n2 <= 0:
                    continue
                new_g = list(genera)
                new_g[v] = g1
                new_g.append(g2_)
                w = nv
                ee = [list(e) for e in edges]
                for i, (k, side) in enumerate(H):
                    if not (mask >> i) & 1:
                        ee[k][side] = w
                ee = [tuple(e) for e in ee] + [(v, w)]
                yield tuple(new_g), ee


def enumerate_graphs(g):
    start = ((g,), ())
    seen = {}
    k0, a0 = canon((g,), [])
    seen[k0] = a0
    front = [k0]
    while front:
        nxt = []
        for gen, ed in front:
            for g2, e2 in degenerations(gen, list(ed)):
                k, a = canon(g2, e2)
                if k not in seen:
                    seen[k] = a
                    nxt.append(k)
        front = nxt
    return [(k[0], list(k[1]), a) for k, a in seen.items()]


# globals for workers
_R = _Ev = None


def init(e1, e2, N, zmax):
    global _R, _Ev
    _R = Rows(e1, e2, N, zmax)
    _Ev = Evaluator(_R, Vert(e1, e2))


def work(args):
    genera, edges, aut, mode = args
    nv = len(genera)
    hot = [i for i in range(nv) if genera[i] >= 2]
    tot = None
    worst = None
    ncol = 0
    if mode == 'total':
        cols = list(product(range(3), repeat=nv))
    else:
        cols = []
        for c in product(range(3), repeat=nv):
            if not any(c[i] == 0 for i in hot):
                continue
            nz = [x for x in c if x]
            if nz and nz[0] == 2:
                continue
            cols.append(c)
    bad = []
    for c in cols:
        T = _Ev.graph(genera, edges, c)
        ncol += 1
        if T.c and T.val < 0:
            if len(edges):
                bad.append((c, T.val))
        if T.aprec < 1 and T.aprec < 10 ** 6:
            bad.append((c, 'prec', T.aprec))
        if mode == 'total':
            tot = T if tot is None else tot + T
    res = (genera, edges, aut, ncol, bad)
    if mode == 'total':
        s = tot * inv(aut)
        return res + (s.val, s.c, s.aprec)
    return res


def main():
    e1, e2, g, N, mode, nproc = Fr(sys.argv[1]), Fr(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], int(sys.argv[6])
    zmax = 3 * g + 3
    t0 = time.time()
    gr = enumerate_graphs(g)
    print('genus', g, 'graphs', len(gr), 'enumerated in', round(time.time() - t0, 1), 'sum 1/Aut (edge-free check)', sum(Fr(1, a) for _, _, a in gr), flush=True)
    if mode == 'hot':
        gr2 = [x for x in gr if any(h >= 2 for h in x[0]) and len(x[1]) >= 1]
        print('hot edge graphs', len(gr2), flush=True)
    else:
        gr2 = gr
    gr2.sort(key=lambda x: len(x[0]))
    tasks = [(a, b, c, mode) for a, b, c in gr2]
    nbad = 0
    ncol = 0
    total = None
    with Pool(nproc, initializer=init, initargs=(e1, e2, N, zmax)) as pool:
        for i, r in enumerate(pool.imap_unordered(work, tasks, chunksize=1)):
            ncol += r[3]
            if r[4]:
                nbad += 1
                print('BAD', r[0], r[1], r[4][:6], flush=True)
            if mode == 'total':
                s = LS(r[5], r[6], r[7])
                total = s if total is None else total + s
            if i % 200 == 0:
                print(' ', i, 'done, coloured terms', ncol, 'time', round(time.time() - t0), flush=True)
    print('graphs tested', len(tasks), 'coloured terms', ncol, 'graphs with a bad colouring', nbad, 'time', round(time.time() - t0))
    if mode == 'total':
        b = beta(2 * g - 2, e1, e2)
        gap = fp((-3) ** (g - 1) * b)
        B = Base(N)
        d = total - B.tc.pow(2 - 2 * g) * gap
        print('total-gap: val', d.val, 'polar', d.polar(), 'aprec', d.aprec)
        c0 = d.coef(0)
        print('t_c^0 coefficient', c0, ratrec(c0, 2 ** 63))


if __name__ == '__main__':
    main()
