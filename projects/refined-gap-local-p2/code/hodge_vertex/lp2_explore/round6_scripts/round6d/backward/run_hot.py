# Per-colouring polar test at genus g, restricted BFS to graphs with a vertex of genus >= 2. Evidence only.
# Contraction of an edge never lowers a vertex genus, so the restricted BFS reaches every such graph.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import time
from itertools import product
from fractions import Fraction as Fr
from multiprocessing import Pool
from gs import *
from run_full import canon, degenerations, init, work


def hot_graphs(g):
    k0, a0 = canon((g,), [])
    seen = {k0: a0}
    front = [k0]
    while front:
        nxt = []
        for gen, ed in front:
            for g2, e2 in degenerations(gen, list(ed)):
                if max(g2) < 2:
                    continue
                k, a = canon(g2, e2)
                if k not in seen:
                    seen[k] = a
                    nxt.append(k)
        front = nxt
    return [(k[0], list(k[1]), a) for k, a in seen.items()]


if __name__ == '__main__':
    e1, e2, g, N, nproc = Fr(sys.argv[1]), Fr(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    t0 = time.time()
    gr = hot_graphs(g)
    gr2 = [x for x in gr if x[1]]
    print('genus', g, 'eps', e1, e2, 'hot edge graphs', len(gr2), 'enumerated in', round(time.time() - t0), flush=True)
    gr2.sort(key=lambda x: len(x[0]))
    tasks = [(a, b, c, 'hot') for a, b, c in gr2]
    ncol = 0
    nbad = 0
    with Pool(nproc, initializer=init, initargs=(e1, e2, N, 3 * g + 3)) as pool:
        for i, r in enumerate(pool.imap_unordered(work, tasks, chunksize=1)):
            ncol += r[3]
            if r[4]:
                nbad += 1
                print('BAD', r[0], r[1], r[4][:6], flush=True)
            if i % 500 == 0:
                print(' ', i, 'done, coloured terms', ncol, 'time', round(time.time() - t0), flush=True)
    print('graphs tested', len(tasks), 'coloured terms', ncol, 'graphs with a bad colouring', nbad, 'time', round(time.time() - t0), flush=True)
