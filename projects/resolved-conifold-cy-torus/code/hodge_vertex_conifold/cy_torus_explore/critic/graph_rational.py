# Critic: each frame-B stable graph (colour-summed) as a rational function of Q, by evaluation at random Q0 in F_p.
# For each graph find minimal m, n with value * Disc^m (1-Q)^n a polynomial N(Q), and deg N. Regular at Q=oo iff deg N <= 2m+n.
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'row_audit'))
from fractions import Fraction as Fr
from cylib import *
import gsum
from gsum import *
import gen_q

def newton_degree(xs, ys, itab):
    K = len(xs)
    c = list(ys)
    deg = -1
    for j in range(K):
        if c[j] % P:
            deg = j
        for i in range(K - 1, j, -1):
            c[i] = (c[i] - c[i - 1]) * itab[i][i - j - 1] % P
        # standard divided differences
    return deg

def main():
    g = int(sys.argv[1]); K = int(sys.argv[2])
    a0, a1, a2, e1 = [Fr(x) for x in sys.argv[3:7]]
    frame = sys.argv[7] if len(sys.argv) > 7 else 'B'
    e2 = a0 + a1 + a2 - e1
    w = weights(a0, a1, a2, e1, e2)
    rnd = random.Random(1234)
    xs = []; vals = None; tots = []
    t0 = time.time()
    while len(xs) < K:
        Q0 = rnd.randrange(2, P)
        disc = (w.a0 * w.a0 + (2 * w.a0 * (w.a1 + w.a2) + 4 * w.a1 * w.a2) * Q0 + (w.a1 - w.a2) ** 2 * Q0 * Q0) % P
        if pow(disc, (P - 1) // 2, P) != 1: continue
        gsum.build_rows = gen_q.rows_at_Q0(Q0)
        S = Sum(w, g, frame, 1)
        tot, rows = S.run()
        if vals is None:
            vals = [[] for _ in rows]; info = [(r[1], sum(sum(x) for x in r[2])) for r in rows]
        for i, r in enumerate(rows):
            vals[i].append(lcoef(r[4], 0))
        tots.append(lcoef(tot, 0))
        xs.append(Q0)
    print('evaluated', K, 'points %.1fs' % (time.time() - t0), flush=True)
    itab = [[0] * K for _ in range(K)]
    for i in range(K):
        for j in range(i):
            itab[i][j] = inv((xs[i] - xs[j]) % P)
    Dv = [(w.a0 * w.a0 + (2 * w.a0 * (w.a1 + w.a2) + 4 * w.a1 * w.a2) * q + (w.a1 - w.a2) ** 2 * q * q) % P for q in xs]
    Ov = [(1 - q) % P for q in xs]
    def analyse(v):
        best = None
        for s in range(0, 40):
            for n in range(0, min(s, 2 * g + 2) + 1):
                m = s - n
                ys = [x * pow(d, m, P) % P * pow(o, n, P) % P for x, d, o in zip(v, Dv, Ov)]
                deg = dd(xs, ys)
                if deg < K - 3:
                    return m, n, deg
        return None
    def dd(xs, ys):
        c = list(ys); deg = -1
        for j in range(K):
            if c[j]: deg = j
            for i in range(K - 1, j, -1):
                c[i] = (c[i] - c[i - 1]) % P * itab[i][i - j - 1] % P
        return deg
    print('g', g, 'frame', frame, 'weights', a0, a1, a2, e1, e2)
    print(' total (m,n,degN):', analyse(tots))
    for i, v in enumerate(vals):
        r = analyse(v)
        if r is None:
            print(' graph', i, info[i], 'no fit'); continue
        m, n, deg = r
        print(' graph', i, 'gv', info[i][0], 'edges', info[i][1], ' Disc^-%d (1-Q)^-%d deg N %d' % (m, n, deg), ' at Q=oo:', 'regular' if deg <= 2 * m + n else 'POLE', ' at Q=1:', 'regular' if n == 0 else 'pole', flush=True)

main()
