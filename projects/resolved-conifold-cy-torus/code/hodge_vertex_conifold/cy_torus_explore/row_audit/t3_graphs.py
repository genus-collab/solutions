# T3: per-graph, per-colouring polar parts at generic CY weights, frame B, against controls;
# total against P_g Li_{3-2g}(Q)/(e1 e2) + constant-map term. Usage: py -I t3_graphs.py g sampleindex frame
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cylib import *
from gsum import *
from t1_rows import SAMPLES, mk

EXTRA = [
    (Fr(2, 3), Fr(5, 7), Fr(-9, 4), Fr(3, 11), None),
    (Fr(1), Fr(2), Fr(4), Fr(3), Fr(4) + 1),   # non-CY control: Ea = 7, Ee = 8
]


def get(i):
    S = SAMPLES + EXTRA
    return mk(S[i])


def main():
    g = int(sys.argv[1])
    idx = int(sys.argv[2])
    frame = sys.argv[3] if len(sys.argv) > 3 else 'B'
    Lx = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    w = get(idx)
    t0 = time.time()
    S = Sum(w, g, frame, Lx)
    print('genus', g, 'weights', [str(x) for x in w.raw], 'CY' if w.cy else 'nonCY', 'frame', frame, flush=True)
    tot, rows = S.run()
    # per colouring regularity
    nreg = nsing = 0
    edge_sing = []
    for (i, gv, A, percol, gtot) in rows:
        hasedge = sum(sum(r) for r in A) > 0
        for col, t in percol:
            pol = [k for k in range(-2 * g, 0) if lcoef(t, k)]
            if hasedge:
                if pol:
                    nsing += 1
                    edge_sing.append((i, gv, col, min(pol)))
                else:
                    nreg += 1
    print('  graphs with an edge: colourings regular %d, singular %d' % (nreg, nsing))
    if edge_sing[:6]:
        print('   first singular', edge_sing[:6])
    # single vertex graph: identify gv=[g], no edges
    pol = {}
    for (i, gv, A, percol, gtot) in rows:
        if len(gv) == 1 and sum(sum(r) for r in A) == 0:
            for col, t in percol:
                pol[col] = [(k, lcoef(t, k)) for k in range(-2 * g, 0) if lcoef(t, k)]
    pred, Pg = predicted(w, g, S.L)
    gap = lcoef(pred, 2 - 2 * g)
    print('  single vertex polar orders, colour 1:', [k for k, v in pol[(1,)]], ' colour 2:', [k for k, v in pol[(2,)]])
    if pol[(1,)]:
        print('  colour 1 leading coefficient = gap:', pol[(1,)][0][1] == gap, ' only t^(2-2g):', len(pol[(1,)]) == 1)
    print('  total polar orders:', [k for k in range(-2 * g, 0) if lcoef(tot, k)], ' leading = gap:', lcoef(tot, 2 - 2 * g) == gap)
    # total against the closed form, Taylor coefficients
    Kmax = Lx - 3
    res = []
    for k in range(-2 * g + 1, Kmax + 1):
        if k == 0:
            continue
        res.append((k, lcoef(tot, k) == lcoef(pred, k)))
    print('  total - P_g Li/(e1 e2): all non-constant coefficients vanish, t^(1-2g)..t^%d:' % Kmax, all(b for _, b in res), [k for k, b in res if not b][:8])
    c0 = (lcoef(tot, 0) - lcoef(pred, 0)) % P
    print('  [t^0] total - [t^0] closed form = constant-map term F_g(0):', c0 == constant_term(w, g))
    print('  time %.1f s' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
