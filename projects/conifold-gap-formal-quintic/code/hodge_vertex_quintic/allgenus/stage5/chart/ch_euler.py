# -*- coding: utf-8 -*-
"""ch_euler.py -- Stage 5 / chart, check of Lemma (L3) (CHART_PROOF s.2): Stirling's series is the formal stationary phase of
Euler's integral, with x^m inserted, for both measures.

For a != 0 the phase a log x - x has the critical point x = a. With M = m (measure dx/x, integrand x^(a/z+m) e^(-x/z) dx/x) or
M = m + 1 (measure dx), (L3) says that the normalised formal SP series is EXACTLY
    Poch_M(z/a) * exp( Sum_k b_k (z/a)^(2k-1) ),   Poch_M(e) = prod_(i=0..M-1) (1 + i e) (M >= 0),  1/prod_(i=1..-M) (1 - i e) (M < 0),
i.e. the normalised series of z^(a/z+M) Gamma_St(a/z + M). Computed here by the one-variable Wick formula (ch_wick.sp1), exact over Q.
Controls (must FAIL): b_k -> B_2k/(2k); the Pochhammer factor dropped; exp(-Sum ...) (1/Gamma instead of Gamma).
Usage: py -I -B ch_euler.py [N] [N2]   (default N = 40 at a = 1, N2 = 20 at a = -2/7, 3/2); log logs/ch_euler_<N>.log
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from ch_core import bern, bk_list, sexp_q
from ch_wick import sp1, normalise

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
N2 = int(sys.argv[2]) if len(sys.argv) > 2 else 20
LOGREL = 'logs/ch_euler_%d.log' % N
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok: FAILS.append(name)


class QF:                                   # the field Q with the ch_core interface
    zero = Fr(0); one = Fr(1)
    def of(self, x): return Fr(x)
    def add(self, a, b): return a + b
    def sub(self, a, b): return a - b
    def neg(self, a): return -a
    def mul(self, a, b): return a * b
    def smul(self, s, a): return Fr(s) * a
    def iszero(self, a): return a == 0
    def isunit(self, a): return a != 0
    def inv(self, a): return 1 / a


Q = QF()


def lhs(a, M, N):
    """normalised Wick series of exp((a log x - x)/z) x^(M-1) dx at x = a, through z^N"""
    K = 2 * N + 2
    phi = [Fr(0), Fr(0)] + [a * Fr((-1) ** (k + 1), k) / a ** k for k in range(2, K + 1)]
    e = M - 1; amp = []; b = Fr(1)                 # (1 + t/a)^(M-1)
    for k in range(2 * N + 1):
        amp.append(b / a ** k); b = b * (e - k) / (k + 1)
    return normalise(Q, sp1(Q, phi, [amp], N))


def rhs(a, M, N, bk, mode=None):
    eps = [Fr(0), 1 / a] + [Fr(0)] * (N - 1)       # e = z/a as a series in z
    H = [Fr(0)] * (N + 1)
    for k in range(1, (N + 1) // 2 + 1):
        if 2 * k - 1 <= N: H[2 * k - 1] = bk[k] / a ** (2 * k - 1)
    if mode == 'neg': H = [-x for x in H]
    E = sexp_q(H, N); P = [Fr(1)] + [Fr(0)] * N
    if mode != 'nopoch':
        rng = range(0, M) if M >= 0 else range(1, -M + 1)
        for i in rng:
            f = [Fr(1), Fr(i) / a] + [Fr(0)] * (N - 1) if M >= 0 else [Fr(1), -Fr(i) / a] + [Fr(0)] * (N - 1)
            if M >= 0:
                P = [sum(P[j] * f[n - j] for j in range(n + 1)) for n in range(N + 1)]
            else:                                   # divide by (1 - i z/a)
                Q2 = [Fr(0)] * (N + 1)
                for n in range(N + 1): Q2[n] = P[n] + (Fr(i) / a * Q2[n - 1] if n else 0)
                P = Q2
    return [sum(P[j] * E[n - j] for j in range(n + 1)) for n in range(N + 1)]


out('ch_euler: (L3) exact over Q; a = 1 through z^%d, a = -2/7 and 3/2 through z^%d; log %s' % (N, N2, LOGREL))
BK = bk_list(N + 2); B = bern(2 * N + 4)
BKbad = [None] + [B[2 * k] / (2 * k) for k in range(1, N + 3)]
for a, NN in ((Fr(1), N), (Fr(-2, 7), N2), (Fr(3, 2), N2)):
    t = time.time()
    for m in (0, 1, 2, 3):
        for meas, M in (('dx/x', m), ('dx', m + 1)):
            L = lhs(a, M, NN)
            check('a = %s, x^%d %s: normalised SP = z^(a/z+%d) Gamma_St(a/z+%d) normalised, z^0..z^%d' % (a, m, meas, M, M, NN),
                  L == rhs(a, M, NN, BK))
    for m in (-1, -2, -3):
        for meas, M in (('dx/x', m), ('dx', m + 1)):
            check('a = %s, x^(%d) %s (negative m, used in step (e)), z^0..z^%d' % (a, m, meas, NN),
                  lhs(a, M, NN) == rhs(a, M, NN, BK))
    out('     a = %s done in %.1f s' % (a, time.time() - t))
L0 = lhs(Fr(1), 0, 12)
fb = [n for n in range(13) if L0[n] != rhs(Fr(1), 0, 12, BKbad)[n]]
check('control: b_k -> B_2k/(2k) fails (first at z^3: b_1 = B_2/2 agrees)', fb[:1] == [3], fb[:3])
L2 = lhs(Fr(1), 2, 12)
fp = [n for n in range(13) if L2[n] != rhs(Fr(1), 2, 12, BK, 'nopoch')[n]]
check('control: Pochhammer factor dropped (M = 2) fails at z^1', fp[:1] == [1], fp[:3])
fn = [n for n in range(13) if L0[n] != rhs(Fr(1), 0, 12, BK, 'neg')[n]]
check('control: exp(-Sum b_k e^(2k-1)) (1/Gamma) fails at z^1', fn[:1] == [1], fn[:3])
out('first coefficients (a = 1, M = 0): %s' % [str(x) for x in L0[:6]])
out('CH_EULER %d: %s   (%.0f s)' % (N, 'ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - T0))
LOG.close()
