# -*- coding: utf-8 -*-
"""ch_chart.py -- Stage 5 / chart: steps (a)-(d) of CHART_PROOF s.3 tested directly, with own code (ch_core, ch_wick, ch_phases).

For rational L != 0, 1 three normalised series must coincide (Theorem Gamma, step (c)):
  c(z; L)   the sigma-form: one-variable Wick formula for the full-Stirling summand at sigma = L (definition s.1.6);
  N7(z; L)  the 7-variable integrand Omega_7 at p(L) = (L, 5L, zeta^j - L): 7-variable Wick with the full (non-diagonal) 7x7
            covariance -- this is the object both orders of integration compute (Fubini (L1b));
  Nch(z; L) the chart integrand Omega_ch at p_ch(L) = (5L, zeta^j - L): 5-variable Wick (order (b), after (L2)).
At L = 1 (q = 0) the chart series must equal h_(w_0)(z) (step (d)); h by Newton sums (ch_core.h_newton).
Fields: exact Q(zeta_5) (small orders) and F_p for the two primes of ch_core.PRIMES (higher orders).
Controls: chart without the Jacobian X_0 of (L2) (must fail, first at z^1); chart with measure dy/y (must fail, first at z^2);
7-variable integrand with dy/y (Gamma(5s)) against c (must fail) and against the consistently modified sigma-form (must PASS).
Usage: py -I -B ch_chart.py NCH_EX N7_EX NCH_P N7_P [NH1]   (log logs/ch_chart_<args>.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from ch_core import QZ, FP, PRIMES, h_newton
from ch_wick import sp_series, sp1, normalise
from ch_phases import sigma_data, phase7, amp7, phasech, ampch

A = [int(x) for x in sys.argv[1:]] + [None] * 5
NCH_EX = A[0] if A[0] is not None else 2
N7_EX = A[1] if A[1] is not None else 1
NCH_P = A[2] if A[2] is not None else 3
N7_P = A[3] if A[3] is not None else 2
NH1 = A[4] if A[4] is not None else max(NCH_EX, 4)
LOGREL = 'logs/ch_chart_%d_%d_%d_%d_%d.log' % (NCH_EX, N7_EX, NCH_P, N7_P, NH1)
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok: FAILS.append(name)


def csig(F, L, N, variant=None):
    phi, amps = sigma_data(F, L, N, variant)
    return normalise(F, sp1(F, phi, amps, N))


def chart(F, L, N, variant=None):
    P, K = phasech(F, L, 2 * N + 2)
    s, H = sp_series(F, P, [ampch(F, L, 2 * N, variant)], N, 5)
    return normalise(F, s)


def seven(F, L, N, variant=None):
    s, H = sp_series(F, phase7(F, L, 2 * N + 2), [amp7(F, L, 2 * N, variant)], N, 7)
    return normalise(F, s)


def same(F, a, b): return len(a) == len(b) and all(F.eq(x, y) for x, y in zip(a, b))


def firstdiff(F, a, b): return [n for n in range(min(len(a), len(b))) if not F.eq(a[n], b[n])]


def short(x):
    s = str(x); return s if len(s) <= 34 else s[:16] + '..(%d chars)' % len(s)


out('ch_chart: exact Q(zeta5): chart through z^%d, Omega_7 through z^%d, q = 0 chart through z^%d; mod p: chart through z^%d, '
    'Omega_7 through z^%d; log %s' % (NCH_EX, N7_EX, NH1, NCH_P, N7_P, LOGREL))
Q = QZ()
t = time.time(); h = h_newton(NH1); s1 = chart(Q, Fr(1), NH1)
check('[Q(zeta)] L = 1 (q = 0): chart SP = h_(w_0)(z) exactly through z^%d (step (d))' % NH1, [Q.rat(x) for x in s1] == h,
      '%.1f s' % (time.time() - t))
for L in (Fr(2), Fr(-2, 3)):
    t = time.time(); c = csig(Q, L, max(NCH_EX, N7_EX)); cq = [Q.rat(x) for x in c]
    ch = chart(Q, L, NCH_EX)
    check('[Q(zeta)] L = %s: chart SP = c(z; L) exactly through z^%d (steps (b)-(c))' % (L, NCH_EX), same(Q, ch, c[:NCH_EX + 1]),
          '%.1f s; c_1..c_3 = %s' % (time.time() - t, [short(x) for x in cq[1:4]]))
    if L == 2 and N7_EX >= 0:
        t = time.time(); s7 = seven(Q, L, N7_EX)
        check('[Q(zeta)] L = 2: 7-variable SP of Omega_7 = c(z; L) exactly through z^%d (order (a) + Fubini)' % N7_EX,
              same(Q, s7, c[:N7_EX + 1]), '%.1f s' % (time.time() - t))
for p in PRIMES:
    F = FP(p); hp = [F.of(x) for x in h_newton(max(NCH_P, 2))]
    t = time.time()
    check('[%s] L = 1: chart SP = h through z^%d' % (F.name, max(NCH_P, 2)), same(F, chart(F, Fr(1), max(NCH_P, 2)), hp))
    for L in (Fr(3), Fr(5, 4), Fr(-1, 2)):
        c = csig(F, L, max(NCH_P, N7_P))
        check('[%s] L = %s: chart SP = c(z; L) through z^%d' % (F.name, L, NCH_P), same(F, chart(F, L, NCH_P), c[:NCH_P + 1]))
        if L == 3:
            check('[%s] L = 3: 7-variable SP of Omega_7 = c(z; L) through z^%d' % (F.name, N7_P), same(F, seven(F, L, N7_P), c[:N7_P + 1]))
    out('     [%s] done in %.1f s' % (F.name, time.time() - t))
F = FP(PRIMES[0]); L = Fr(2); NC = max(NCH_P, 2); c = csig(F, L, NC)
fj = firstdiff(F, chart(F, L, NC, 'nojac'), c)
check('control: chart without the Jacobian X_0 of (L2) fails, first at z^1', fj[:1] == [1], fj[:3])
fd = firstdiff(F, chart(F, L, NC, 'dyy'), c)
check('control: chart with measure dy/y fails, first at z^2 (z^1 is blind to it)', fd[:1] == [2], fd[:3])
N7c = max(N7_P, 2); s7d = seven(F, L, N7c, 'dyy')
f7 = firstdiff(F, s7d, csig(F, L, N7c))
check('control: Omega_7 with dy/y (Gamma(5s)) differs from c(z; L)', len(f7) > 0, f7[:3])
check('control pair: Omega_7 with dy/y equals the sigma-form built with Gamma_St(5s) (consistent change), through z^%d' % N7c,
      same(F, s7d, csig(F, L, N7c, 'g5s')))
out('CH_CHART: %s   (%.0f s)' % ('ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - T0))
LOG.close()
