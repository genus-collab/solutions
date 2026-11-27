# -*- coding: utf-8 -*-
"""ra_poly.py -- stage5/chart/referee_adversary: the Laurent data of c_n at L = 1 (ra_sp.py) are those of a polynomial
c_n(L) = Sum_k a_k L^(4n - 5k), 0 <= k <= floor(4n/5) (the Stage-4 shape), and that polynomial reproduces c_n(2) computed
independently by evaluation at L = 2.  Exact over Q.  Usage: py -I -B ra_poly.py N EXTRA   (log logs/ra_poly_N_EXTRA.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from ra_core import FQ, RSym, RNum
from ra_sp import laurent
from ra_models import quintic

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
EXTRA = int(sys.argv[2]) if len(sys.argv) > 2 else 6
LOGREL = 'logs/ra_poly_%d_%d.log' % (N, EXTRA)
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


F = FQ(); K = 2 * N + EXTRA + 1                         # eps-truncation: enough Taylor data for every n <= N
out('ra_poly: polynomial shape of c_n from its Taylor data at L = 1, n <= %d, %d check coefficients; log %s' % (N, EXTRA, LOGREL))
Nn = laurent(quintic(), RSym(F, K), N)
c2 = laurent(quintic(), RNum(F, 1), N)                     # eps = 1: c_n(2) = N_n
for n in range(1, N + 1):
    tay = [Nn[n][n + e] for e in range(K - n)]           # [eps^e] c_n(1 + eps)
    ex = [4 * n - 5 * k for k in range(4 * n // 5 + 1)]; m = len(ex)
    A = [[Fr(comb(x, e)) for x in ex] + [tay[e]] for e in range(m)]
    for i in range(m):
        piv = next(r for r in range(i, m) if A[r][i] != 0); A[i], A[piv] = A[piv], A[i]
        for r in range(m):
            if r != i and A[r][i] != 0:
                f = A[r][i] / A[i][i]; A[r] = [x - f * y for x, y in zip(A[r], A[i])]
    a = [A[i][m] / A[i][i] for i in range(m)]
    rest = len(tay) - m
    ok = rest >= EXTRA and all(sum(ak * comb(x, e) for ak, x in zip(a, ex)) == tay[e] for e in range(len(tay)))
    ok2 = sum(ak * Fr(2) ** x for ak, x in zip(a, ex)) == c2[n]
    s = 'PASS' if ok and ok2 else 'FAIL'
    if s == 'FAIL':
        FAILS.append(n)
    out('%s n = %d: fit on %d Taylor coefficients, confirmed on %d more; fitted polynomial at L = 2 equals c_n(2) = %s: %s'
        % (s, n, m, rest, c2[n], ok2))
out('RA_POLY %d %d: %s   (%.0f s)' % (N, EXTRA, 'ALL CHECKS PASSED' if not FAILS else 'FAILURES: %s' % FAILS, time.time() - T0))
