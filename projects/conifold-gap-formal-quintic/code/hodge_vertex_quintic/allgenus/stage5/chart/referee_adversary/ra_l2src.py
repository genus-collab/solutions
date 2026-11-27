# -*- coding: utf-8 -*-
"""ra_l2src.py -- stage5/chart/referee_adversary: exact test of the delta lemma WITH A SOURCE, (L2') of REFEREE.md s.3.2.

(L2'): variables (s, m), phase P(s) + alpha0 s m with P = Sum_(k>=2) p_k s^k (no pure m term), amplitude e^(c m) A(s, z).
Claim: N = exp(P(delta)/z) A(delta, z) EXACTLY, delta = -c z/alpha0 (for c = 0 this is the worker's (L2), CHART_PROOF_2 s.2.7).
Test: the full 2-variable Wick series (2x2 covariance, Isserlis moments, no use of the claim) on seeded random rational data.
Controls: a pure m^3 term in the phase (the multiplier is no longer linear) must break it; delta with the wrong sign must fail.
Usage: py -I -B ra_l2src.py N SEEDS   (log logs/ra_l2src_N_SEEDS.log)
"""
import os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import factorial

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
SEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 4
LOGREL = 'logs/ra_l2src_%d_%d.log' % (N, SEEDS)
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok:
        FAILS.append(name)


def bmul(a, b, D):
    """bivariate dicts {(r-power, eta-degree): coeff}, truncated at r^D"""
    o = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            if i + k <= D:
                o[(i + k, j + l)] = o.get((i + k, j + l), 0) + x * y
    return o


def bexp(V, D):
    E = {(0, 0): Fr(1)}; term = {(0, 0): Fr(1)}
    for n in range(1, D + 1):
        term = {k: v / n for k, v in bmul(term, V, D).items()}
        if not term:
            break
        for k, v in term.items():
            E[k] = E.get(k, 0) + v
    return E


def wick2(p, alpha0, c, A, q3, N):
    """N_n, n <= N: phase Sum p_k s^k + alpha0 s m + q3 m^3, amplitude e^(c m) exp(Sum A[(i, j)] z^i s^j)"""
    D = 2 * N
    V1 = {(k - 2, k): Fr(p[k]) for k in p if k >= 3 and k - 2 <= D}
    LA = {(2 * i + j, j): Fr(v) for (i, j), v in A.items() if 2 * i + j <= D}
    T = bmul(bexp(V1, D), bexp(LA, D), D)
    V2 = {(1, 3): Fr(q3)} if q3 else {}
    V2[(1, 1)] = V2.get((1, 1), 0) + Fr(c)
    Mp = bexp(V2, D)
    h11, h12, h22 = 2 * Fr(p[2]), Fr(alpha0), Fr(0)
    det = h11 * h22 - h12 * h12
    C11, C12, C22 = -h22 / det, h12 / det, -h11 / det

    def mom(a, b):
        s = Fr(0)
        for cc in range(0, min(a, b) + 1):
            if (a - cc) % 2 or (b - cc) % 2:
                continue
            i, j = (a - cc) // 2, (b - cc) // 2
            s += Fr(factorial(a) * factorial(b), factorial(cc) * factorial(i) * factorial(j) * 2 ** (i + j)) \
                * C12 ** cc * C11 ** i * C22 ** j
        return s
    res = [Fr(0)] * (N + 1)
    for (k1, a), x in T.items():
        for (k2, b), y in Mp.items():
            if (k1 + k2) % 2 == 0 and (k1 + k2) // 2 <= N and (a + b) % 2 == 0:
                res[(k1 + k2) // 2] += x * y * mom(a, b)
    return res


def predicted(p, alpha0, c, A, N, sign=1):
    """exp(P(delta)/z) A(delta, z), delta = -c z/alpha0, as a series in z"""
    dc = -sign * Fr(c) / alpha0
    L = [Fr(0)] * (N + 1)
    for k, v in p.items():
        if 1 <= k - 1 <= N:
            L[k - 1] += Fr(v) * dc ** k
    for (i, j), v in A.items():
        if 1 <= i + j <= N:
            L[i + j] += Fr(v) * dc ** j
    E = [Fr(1)] + [Fr(0)] * N
    for n in range(1, N + 1):
        E[n] = sum(i * L[i] * E[n - i] for i in range(1, n + 1)) / n
    a00 = sum(Fr(v) for (i, j), v in A.items() if i == 0 and j == 0)
    return E, a00


out("ra_l2src: (L2') delta lemma with a source, exact over Q through z^%d, %d seeds; log %s" % (N, SEEDS, LOGREL))
for seed in range(SEEDS):
    rnd = random.Random(1000 + seed)
    rq = lambda: Fr(rnd.randint(-9, 9), rnd.randint(1, 7))
    p = {k: rq() for k in range(2, 2 * N + 3)}
    p[2] = p[2] or Fr(1, 3)
    A = {(i, j): rq() for i in range(0, N + 1) for j in range(0, 2 * N + 1) if (i, j) != (0, 0)}
    alpha0 = rq() or Fr(2); c = rnd.choice([0, 1, 2, -1, -3, Fr(1, 5)])
    w = wick2(p, alpha0, c, A, 0, N); w = [x / w[0] for x in w]
    pr, _ = predicted(p, alpha0, c, A, N)
    check('seed %d: (alpha0, c) = (%s, %s): 2-variable Wick = exp(P(delta)/z) A(delta, z) through z^%d' % (seed, alpha0, c, N),
          w == pr)
    if c != 0:
        prw, _ = predicted(p, alpha0, c, A, N, sign=-1)
        check('control seed %d: delta with the wrong sign is rejected' % seed, w != prw)
    w3 = wick2(p, alpha0, c, A, Fr(1, 2), N); w3 = [x / w3[0] for x in w3]
    check('control seed %d: a pure m^3 term in the phase breaks the identity' % seed, w3 != pr)
out('RA_L2SRC %d %d: %s   (%.0f s)' % (N, SEEDS, 'ALL CHECKS PASSED' if not FAILS else 'FAILURES: %s' % FAILS, time.time() - T0))
