# -*- coding: utf-8 -*-
"""ch_sigma.py -- Stage 5 / chart: the sigma-form c(z; L) of CHART_PROOF s.1.6 (definition checks, Z/5 covariance), exact.

c(z; L) is computed by the one-variable Wick formula (ch_wick.sp1) from the Taylor data of s.1.6, here for L an arbitrary element of
Q(zeta_5) (needed for the covariance test); the data are assembled separately from ch_phases.sigma_data (rational L only),
the Wick routine sp1 and the exponential _bexp are shared. Checked:
  (1) the two implementations agree at rational L;  (2) c_n(L) is rational at rational L;
  (3) c_1(L) = -(3/20) L^4, c_2(L) = (9/800) L^8 at several L;  (4) the stored value c_5(2) = 11495493197/31500 (GAMMA logs, SYN4 s.3 E3);
  (5) Theorem Gamma (iv), Z/5 covariance: c_n(zeta^i L) = zeta^(-in) c_n(L) exactly in Q(zeta), i = 1..4, n <= N, two values of L.
Control (must FAIL): the covariance with zeta^(+in).
Usage: py -I -B ch_sigma.py [N]   (default 8; log logs/ch_sigma_<N>.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from ch_core import QZ, bk_list
from ch_wick import sp1, normalise
from ch_phases import sigma_data, _bexp

N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
LOGREL = 'logs/ch_sigma_%d.log' % N
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()
F = QZ()


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok: FAILS.append(name)


def cgen(Lel, N):
    """normalised c(z; L) for L a Q(zeta) element (L != 0, L^5 != 1)"""
    K = 2 * N + 3; iL = F.inv(Lel); P = [F.zero] * (K + 1)
    for j in range(5):
        ic = F.inv(F.sub(F.zeta(j), Lel)); p = F.one
        for k in range(K + 1): P[k] = F.add(P[k], p); p = F.mul(p, ic)
    iLp = [F.one]
    for _ in range(K + 1): iLp.append(F.mul(iLp[-1], iL))
    phi = [F.zero] * (2 * N + 3)
    for k in range(1, 2 * N + 2):
        phi[k + 1] = F.smul(Fr(1, k * (k + 1)), F.add(F.smul(5 * (-1) ** (k + 1), iLp[k]), P[k]))
    X = [[F.zero] * (2 * N + 1) for _ in range(N + 1)]
    for k in range(1, 2 * N + 1): X[0][k] = F.smul(Fr(1, 2 * k), F.add(F.smul((-1) ** (k + 1), iLp[k]), P[k]))
    bk = bk_list(N // 2 + 2)
    for k in range(1, N // 2 + 2):
        m = 2 * k - 1
        if m > N: break
        for r in range(0, 2 * (N - m) + 1):
            v = F.add(F.smul(Fr(5) ** (1 - 2 * k), iLp[2 * k - 1 + r]), F.smul((-1) ** r, P[2 * k - 1 + r]))
            X[m][r] = F.smul(bk[k] * (-1) ** r * comb(2 * k - 2 + r, r), v)
    return normalise(F, sp1(F, phi, _bexp(F, X, N), N))


out('ch_sigma: c(z; L) exact over Q(zeta5) through z^%d; log %s' % (N, LOGREL))
LS = [Fr(2), Fr(-2, 3), Fr(7, 5), Fr(3)]
ok1 = ok2 = ok3 = True; C = {}
for L in LS:
    a = cgen(F.of(L), N); b = normalise(F, sp1(F, *sigma_data(F, L, N), N)); C[L] = a
    ok1 &= all(F.eq(x, y) for x, y in zip(a, b))
    try: [F.rat(x) for x in a]
    except AssertionError: ok2 = False
    ok3 &= F.eq(a[1], F.of(Fr(-3, 20) * L ** 4)) and F.eq(a[2], F.of(Fr(9, 800) * L ** 8))
check('two assemblies of the Taylor data (rational L in ch_phases, Q(zeta)-valued L here; shared Wick routine) agree through z^%d at L = %s' % (N, ', '.join(map(str, LS))), ok1)
check('c_n(L) is rational at rational L (n <= %d)' % N, ok2)
check('c_1(L) = -(3/20) L^4 and c_2(L) = (9/800) L^8 at these L (so c_1(1) = h_1, c_2(1) = h_2)', ok3)
if N >= 5:
    check('stored value c_5(2) = 11495493197/31500 (GAMMA logs; SYN4 s.3 E3) reproduced', F.rat(C[Fr(2)][5]) == Fr(11495493197, 31500),
          F.rat(C[Fr(2)][5]))
okc = okw = True
for L in (Fr(2), Fr(-2, 3)):
    for i in range(1, 5):
        cz = cgen(F.mul(F.zeta(i), F.of(L)), N)
        okc &= all(F.eq(cz[n], F.mul(F.zeta(-i * n), C[L][n])) for n in range(N + 1))
        okw &= all(F.eq(cz[n], F.mul(F.zeta(i * n), C[L][n])) for n in range(N + 1))
check('Z/5 covariance: c_n(zeta^i L) = zeta^(-in) c_n(L) exactly, i = 1..4, n <= %d, L = 2, -2/3' % N, okc)
check('control: the covariance with zeta^(+in) fails', not okw)
out('c_3(2), c_4(2) = %s, %s' % (F.rat(C[Fr(2)][3]), F.rat(C[Fr(2)][4])))
out('CH_SIGMA %d: %s   (%.0f s)' % (N, 'ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - T0))
LOG.close()
