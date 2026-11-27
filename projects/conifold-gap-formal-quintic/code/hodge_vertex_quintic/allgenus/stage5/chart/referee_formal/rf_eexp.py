# -*- coding: utf-8 -*-
"""rf_eexp.py -- referee check of step (e) (CHART_PROOF s.3.7; not used for Theorems Gamma, D):
    L e^(Dpsi(q)/z) c(z; L(q)) = h(z) Sum_d q^d I_d(z)   in Q((z))[[q]],
L(q) = (1 - 5^5 q)^(-1/5), Dpsi = Sum_(d>=1) [q^d](L - 1) q^d/d (so D_q Dpsi = L - 1),
I_d(z) = (1/d!) z^-d prod_(k<=5d) (5 + kz) / prod_(i<=d) [((1+iz)^5 - 1)/(iz)]   (= prod_(j>=1) (1 - zeta^j + iz)).
c_n(L) = L^(4n) ph_n(L^-5) with ph_n from the HATTED computation rf_hat (POLY's integrand), not from the sigma-form fit.
Usage: py -I -B rf_eexp.py DQ NZ      (log logs/rf_eexp_DQ_NZ.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import factorial
from rf_core import gbinom, h_series, ser_mul, ser_exp
from rf_hat import hat_series

DQ = int(sys.argv[1]) if len(sys.argv) > 1 else 4; NZ = int(sys.argv[2]) if len(sys.argv) > 2 else 8
LOGREL = 'logs/rf_eexp_%d_%d.log' % (DQ, NZ); LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8')
FAILS = []; T0 = time.time()


def out(s): print(s); LOG.write(s + '\n'); LOG.flush()


def check(ok, msg, ctrl=False):
    good = (not ok) if ctrl else ok
    out(('PASS ' if good else 'FAIL ') + ('control: ' if ctrl else '') + msg)
    if not good: FAILS.append(msg)


NN = NZ + DQ; U = (4 * NN) // 5 + 4
ph = hat_series(NN, U); ph = [[x for x in p.c] for p in ph]
qser = lambda e: [gbinom(e, d) * Fr(-3125) ** d for d in range(DQ + 1)]      # (1 - 3125 q)^e
Lq = qser(Fr(-1, 5)); Dpsi = [Fr(0)] + [Lq[d] / d for d in range(1, DQ + 1)]
cq = []                                                                       # cq[n][d] = [q^d] c_n(L(q))
for n in range(NN + 1):
    s = [Fr(0)] * (DQ + 1)
    for k, a in enumerate(ph[n]):
        if a: s = [x + a * y for x, y in zip(s, qser(Fr(5 * k - 4 * n, 5)))]
    cq.append(s)
h = h_series(NN); h5 = h_series(NN, twist=Fr(5))


def lhs(withL=True, withE=True):
    """dict (d, n) -> coefficient of q^d z^n"""
    E = {(0, 0): Fr(1)}                                                       # e^(Dpsi/z): q^d z^-m
    if withE:
        P = [Fr(1)] + [Fr(0)] * DQ
        for m in range(1, DQ + 1):
            P = ser_mul(P, Dpsi, DQ, Fr(0))
            for d in range(DQ + 1):
                if P[d]: E[(d, -m)] = E.get((d, -m), Fr(0)) + P[d] / factorial(m)
    R = {}
    for (d1, m), e in E.items():
        for n in range(NN + 1):
            for d2 in range(DQ + 1 - d1):
                if cq[n][d2]: R[(d1 + d2, n + m)] = R.get((d1 + d2, n + m), Fr(0)) + e * cq[n][d2]
    if withL:
        R2 = {}
        for (d, n), v in R.items():
            for d3 in range(DQ + 1 - d): R2[(d + d3, n)] = R2.get((d + d3, n), Fr(0)) + v * Lq[d3]
        R = R2
    return R


def rhs(hh):
    R = {}
    for d in range(DQ + 1):
        num = [Fr(1)]
        for k in range(1, 5 * d + 1): num = ser_mul(num, [Fr(5), Fr(k)], NN + d, Fr(0))
        den = [Fr(1)]
        for i in range(1, d + 1): den = ser_mul(den, [Fr(5), Fr(10 * i), Fr(10 * i * i), Fr(5 * i ** 3), Fr(i ** 4)], NN + d, Fr(0))
        inv = [1 / den[0]]
        for k in range(1, NN + d + 1): inv.append(-sum(den[i] * inv[k - i] for i in range(1, min(k, len(den) - 1) + 1)) * inv[0])
        Id = [x / factorial(d) for x in ser_mul(num, inv, NN + d, Fr(0))]        # times z^-d
        prod = ser_mul(hh, Id, NN + d, Fr(0))
        for k, v in enumerate(prod):
            if v: R[(d, k - d)] = v
    return R


def cmp(A, B):
    keys = [(d, n) for d in range(DQ + 1) for n in range(-d, NZ + 1)]
    bad = [k for k in keys if A.get(k, Fr(0)) != B.get(k, Fr(0))]
    return bad, sum(1 for k in keys if B.get(k, Fr(0)) != 0)


out('rf_eexp: q^0..q^%d, z^-d..z^%d; c_n from the hatted ph_n (rf_hat, n <= %d); log %s' % (DQ, NZ, NN, LOGREL))
bad, nz = cmp(lhs(), rhs(h))
check(not bad, 'step (e): L e^(Dpsi/z) c(z; L(q)) = h(z) Sum_d q^d I_d(z), all coefficients q^d z^n, d <= %d, -d <= n <= %d (%d non-zero)' % (DQ, NZ, nz))
for args, txt in [((False, True), 'factor L dropped'), ((True, False), 'e^(Dpsi/z) dropped')]:
    bad, _ = cmp(lhs(*args), rhs(h)); check(not bad, '%s fails, first at (q^d, z^n) = %s' % (txt, bad[:1]), ctrl=True)
bad, _ = cmp(lhs(), rhs(h5)); check(not bad, 'twist +5 in h fails, first at %s' % bad[:1], ctrl=True)
out('RF_EEXP %d %d: %s   (%.0f s)' % (DQ, NZ, 'ALL CHECKS PASSED' if not FAILS else 'FAILURES: %s' % FAILS, time.time() - T0))
LOG.close()
