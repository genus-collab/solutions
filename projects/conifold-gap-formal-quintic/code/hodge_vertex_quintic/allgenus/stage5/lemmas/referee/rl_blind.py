# -*- coding: utf-8 -*-
"""rl_blind.py -- a record of one BLIND control met while setting up rl_toy.py (new code, stdlib only; uses rl_core).
In the toy chart of rl_toy (k = 2, colours 1, -2) the control "chart amplitude without the Jacobian X_0", i.e. 1/(X_0 x_1) ~ y^(-k)
instead of 1/x_1, gives the SAME normalised chart series at the special point L = 2 (through z^K), but not at generic L.
So a control evaluated at one special parameter value can fail to fail; rl_toy therefore uses L = 3 and L = -3/2.
Usage: py -I -B rl_blind.py [K]   (default 9)
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from rl_core import sp_norm, comb_gen

K = int(sys.argv[1]) if len(sys.argv) > 1 else 9
KT, C0, C1 = 2, Fr(1), Fr(-2)
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
REL = 'logs/rl_blind_%d.log' % K; LOG = open(os.path.join(HERE, REL), 'w', encoding='utf-8')


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def chart(L, M, jac):
    b, a0, a1 = KT * L, C0 - L, C1 - L; tay = {}
    for n in range(M + 3):
        for p in range(M + 3 - n):
            if n + p < 2:
                continue
            v = Fr(0)
            if p == 0:
                v += C0 * KT * Fr((-1) ** (n + 1), n) / b ** n
            if n == 0:
                v += (C1 - C0) * Fr((-1) ** (p + 1), p) / a1 ** p
            v += -a0 * comb(KT, n) / b ** n * Fr((-1) ** p) / a1 ** p
            if v:
                tay[(n, p)] = v
    amp = {}
    for n in range(M + 1):
        for p in range(M + 1 - n):
            v = (Fr((-1) ** p) / a1 ** p if n == 0 else 0) if jac else (comb_gen(-KT, n) / b ** n if p == 0 else 0)
            if v:
                amp[((n, p), 0)] = v
    return tay, amp


out('rl_blind: toy chart, amplitude 1/x_1 versus 1/(X_0 x_1) (no Jacobian), through z^%d; log %s' % (K, REL))
for L in (Fr(2), Fr(3), Fr(7, 4), Fr(-3, 2)):
    A = sp_norm(*chart(L, 2 * K, True), 2, K); B = sp_norm(*chart(L, 2 * K, False), 2, K)
    d = next((i for i in range(K + 1) if A[i] != B[i]), None)
    out('L = %s: the two normalised chart series %s' % (L, 'AGREE through z^%d (blind control)' % K if d is None
                                                         else 'differ first at z^%d' % d))
out('RL_BLIND %d: done (an observation, not a test of the lemmas)' % K)
LOG.close()
