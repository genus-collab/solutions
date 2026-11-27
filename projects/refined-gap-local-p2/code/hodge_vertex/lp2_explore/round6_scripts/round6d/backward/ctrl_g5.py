# Controls for the genus-5 hot test: perturb F_cf at xi^k, or switch to mode A. Evidence only.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import time
from itertools import product
from fractions import Fraction as Fr
from multiprocessing import Pool
from gs import *
from run_full import enumerate_graphs

_Ev = None


def init(e1, e2, N, zmax, k, modeA):
    global _Ev
    b = Base(N)
    FF = b.XX
    if k is not None:
        FF = b.XX + LS(k, [fp(Fr(1, 9))], N + 3)
    R = Rows(e1, e2, N, zmax, base=b, FF=FF)
    E = Fr(e1) + Fr(e2)
    V = Vert(e1, e2, w=(-E, Fr(e1), Fr(e2))) if modeA else Vert(e1, e2)
    _Ev = Evaluator(R, V)


def work(args):
    genera, edges = args
    hot = [i for i in range(len(genera)) if genera[i] >= 2]
    nbad = 0
    n = 0
    worst = 0
    for c in product(range(3), repeat=len(genera)):
        if not any(c[i] == 0 for i in hot):
            continue
        nz = [x for x in c if x]
        if nz and nz[0] == 2:
            continue
        T = _Ev.graph(genera, edges, c)
        n += 1
        if T.c and T.val < 0:
            nbad += 1
            worst = min(worst, T.val)
    return n, nbad, worst


if __name__ == '__main__':
    gr = enumerate_graphs(5)
    hot = [(a, b) for a, b, c in gr if any(h >= 2 for h in a) and b]
    for label, k, modeA in [('F+xi^6/9', 6, False), ('F+xi^7/9', 7, False), ('mode A', None, True)]:
        t0 = time.time()
        with Pool(4, initializer=init, initargs=(3, 7, 14, 18, k, modeA)) as pool:
            res = pool.map(work, hot, chunksize=4)
        print(label, 'coloured terms', sum(r[0] for r in res), 'with a pole', sum(r[1] for r in res), 'worst val', min(r[2] for r in res), 'time', round(time.time() - t0), flush=True)
