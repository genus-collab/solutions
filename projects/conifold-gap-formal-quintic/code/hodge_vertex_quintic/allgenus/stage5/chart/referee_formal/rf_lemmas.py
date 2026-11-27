# -*- coding: utf-8 -*-
"""rf_lemmas.py -- referee tests of the formal lemmas of CHART_PROOF s.2 on seeded random exact data (Q), with failing controls.
Two independent SP codes are used: rf_core.sp1 (Morse + Lagrange-Buermann, one variable) and rf_wick.wickN with the
differential-operator expectation (any dimension).
  (T0) sp1 = wickN on random one-variable data;  (T1) (L3) for M = -3..5, a in Q and a in Q(zeta);
  (T2) (L2) delta lemma: N = A(X_0) exactly;     (T3) (L1b) Fubini (inner dim 1) incl. the factor beta, over Q[[w]];
  (T4) (S1);  (T5) (L1a) with a nonlinear chi.
Usage: py -I -B rf_lemmas.py N SEEDS     (log logs/rf_lemmas_N_SEEDS.log)
"""
import os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from rf_core import QZ, sp1, bk_list, bernoulli, ser_exp, ser_mul, gbinom
from rf_rings import TS
from rf_wick import wickN, make_expect_do, mat_inv, hessian, split_deg

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4; SEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
LOGREL = 'logs/rf_lemmas_%d_%d.log' % (N, SEEDS); LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8')
FAILS = []; T0 = time.time(); one = Fr(1)


def out(s): print(s); LOG.write(s + '\n'); LOG.flush()


def check(ok, msg, ctrl=False):
    good = (not ok) if ctrl else ok
    out(('PASS ' if good else 'FAIL ') + ('control: ' if ctrl else '') + msg)
    if not good: FAILS.append(msg)


def rq(R): return Fr(R.randint(-9, 9), R.randint(1, 7))


def _tpw(x, k):
    r = TS.of(1)
    for _ in range(k): r = r * x
    return r


def wick(Phi, amps, n, nv):
    C = [[-x for x in r] for r in mat_inv(hessian(Phi, nv, one), one)]
    return wickN(Phi, amps, n, nv, make_expect_do(C, one), one)


def amps1(A0, A1, n):
    """one-variable lists -> amps for sp1: A = A0 + z A1"""
    return [A0[:2 * n + 1], A1[:2 * n - 1]] + [[] for _ in range(n - 1)]


out('rf_lemmas: N = %d, %d seeds; log %s' % (N, SEEDS, LOGREL))
ok = {k: True for k in ['T0', 'T2', 'T2c1', 'T2c2', 'T3', 'T3c', 'T3c2', 'T4', 'T4c', 'T5', 'T5c']}
for seed in range(SEEDS):
    R = random.Random(1000 + seed)
    # (T0) one variable: sp1 (Lagrange-Buermann) versus wickN (differential operator)
    ph = [Fr(0), Fr(0), -abs(rq(R)) - 1] + [rq(R) for _ in range(2 * N)]
    A0 = [Fr(1)] + [rq(R) for _ in range(2 * N)]; A1 = [rq(R) for _ in range(2 * N)]
    s1 = sp1(ph, amps1(A0, A1, N), N, one)
    P1 = {(k,): c for k, c in enumerate(ph) if c}
    w1 = wick(P1, [{(k,): c for k, c in enumerate(A0)}, {(k,): c for k, c in enumerate(A1)}], N, 1)
    ok['T0'] &= s1 == w1
    # (T2) delta lemma: Phi = s G(xi) + Fh(xi), G(0) = 0, g_1 != 0, A = A(xi, z) (no s)
    g = [Fr(0), rq(R) or Fr(1)] + [rq(R) for _ in range(2 * N + 1)]; fh = [Fr(0), Fr(0)] + [rq(R) for _ in range(2 * N + 1)]
    Phi = {}
    for k in range(1, 2 * N + 2): Phi[(1, k)] = g[k]
    for k in range(2, 2 * N + 3): Phi[(0, k)] = fh[k]
    Phi = {e: c for e, c in Phi.items() if c and sum(e) <= 2 * N + 2}
    Ax = [{(0, k): c for k, c in enumerate(A0)}, {(0, k): c for k, c in enumerate(A1)}]
    r = wick(Phi, Ax, N, 2); ok['T2'] &= r == [A0[0], A1[0]] + [Fr(0)] * (N - 1)
    Ab = [{**Ax[0], (1, 1): Fr(1)}, Ax[1]]                          # control: an s-dependent amplitude term s xi
    ok['T2c1'] &= wick(Phi, Ab, N, 2) == r
    Pb = dict(Phi); Pb[(2, 1)] = Fr(1)                                     # control: s^2 xi in the phase
    ok['T2c2'] &= wick(Pb, Ax, N, 2) == r
    # (T3) Fubini, inner variable v, outer w; Phi(v, w) random with critical point 0
    D = 2 * N + 2; Phi = {}
    for a in range(D + 1):
        for b in range(D + 1 - a):
            if a + b >= 2: Phi[(a, b)] = rq(R)
    Phi[(2, 0)] = -abs(Phi[(2, 0)]) - 2; Phi[(0, 2)] = -abs(Phi[(0, 2)]) - 2; Phi[(1, 1)] = Fr(1, 3)
    Am = [{(a, b): (rq(R) if (a, b) != (0, 0) else Fr(1)) for a in range(2 * N + 1) for b in range(2 * N + 1 - a)},
          {(a, b): rq(R) for a in range(2 * N - 1) for b in range(2 * N - 1 - a)}]
    tot = wick(Phi, Am, N, 2)
    TS.n = 2 * N + 4; W = TS.var(); H = Phi[(2, 0)] * 2
    vc = TS.of(0)
    for _ in range(TS.n + 1):                                                # d_v Phi(vc, w) = 0 by fixed point
        Rv = TS.of(0)
        for (a, b), c in Phi.items():
            if a >= 1 and (a, b) != (2, 0): Rv = Rv + _tpw(vc, a - 1) * _tpw(W, b) * (c * a)
        vc = Rv * (-1 / H)

    def shift(P, vc0, kmax):
        co = [TS.of(0) for _ in range(kmax + 1)]
        for (a, b), c in P.items():
            for k in range(min(a, kmax) + 1): co[k] = co[k] + _tpw(vc0, a - k) * _tpw(W, b) * (c * comb(a, k))
        return co

    def fub(vc0, use_beta=True):
        inn = shift(Phi, vc0, D); amA = [shift(Am[0], vc0, 2 * N), shift(Am[1], vc0, 2 * N - 2)] + [[] for _ in range(N - 1)]
        B = sp1(inn, amA, N, TS.of(1))
        x = inn[2] * (1 / (inn[2].c[0])); y = x - 1; beta = TS.of(0); yp = TS.of(1)
        for k in range(TS.n): beta = beta + yp * gbinom(Fr(-1, 2), k); yp = yp * y
        if not use_beta: beta = TS.of(1)
        oph = [Fr(0), Fr(0)] + [inn[0].c[k] for k in range(2, D + 1)]
        oam = [[(beta * B[m]).c[j] for j in range(2 * (N - m) + 1)] for m in range(N + 1)]
        return sp1(oph, oam, N, one), inn[1]
    it, lin = fub(vc)
    ok['T3'] &= it == tot and lin == TS.of(0)                             # inner point critical identically in w
    ok['T3c'] &= fub(vc, False)[0] == tot                                  # control: beta dropped
    ok['T3c2'] &= fub(TS.of(0))[0] == tot                                 # control: inner critical point not moved
    # (T4) (S1): N[Phi, (d_1 Phi) h + z d_1 h] = 0
    P2 = {e: c for e, c in Phi.items()}; h0 = {(a, b): rq(R) for a in range(2 * N) for b in range(2 * N - a)}
    dP = {(a - 1, b): c * a for (a, b), c in P2.items() if a >= 1}
    def mul(A, B):
        Rr = {}
        for e1, c1 in A.items():
            for e2, c2 in B.items():
                e = (e1[0] + e2[0], e1[1] + e2[1]); Rr[e] = Rr.get(e, Fr(0)) + c1 * c2
        return {e: c for e, c in Rr.items() if sum(e) <= 2 * N}
    dh = {(a - 1, b): c * a for (a, b), c in h0.items() if a >= 1}
    ok['T4'] &= wick(P2, [mul(dP, h0), dh], N, 2) == [Fr(0)] * (N + 1)
    ok['T4c'] &= wick(P2, [mul(dP, h0), {}], N, 2) == [Fr(0)] * (N + 1)
    # (T5) (L1a): chi(t) = mu t + c2 t^2 + c3 t^3; N[ph o chi, (A o chi) chi'] = chi'(0) N[ph, A]
    mu = Fr(R.randint(1, 5), R.randint(1, 3)); c2, c3 = rq(R), rq(R); K = 2 * N + 2
    chi = [Fr(0), mu, c2, c3] + [Fr(0)] * (K - 3)
    def comp(f, K):
        res = [Fr(0)] * (K + 1); p = [Fr(1)] + [Fr(0)] * K
        for k in range(len(f)):
            if k > K: break
            res = [x + f[k] * y for x, y in zip(res, p)]; p = ser_mul(p, chi, K, Fr(0))
        return res
    phc = comp(ph, K); dchi = [Fr(1) * mu, 2 * c2, 3 * c3] + [Fr(0)] * K
    A0c = ser_mul(comp(A0, 2 * N), dchi, 2 * N, Fr(0)); A1c = ser_mul(comp(A1, 2 * N), dchi, 2 * N, Fr(0))
    lhs = sp1(phc, amps1(A0c, A1c, N), N, one)
    ok['T5'] &= lhs == [x * mu for x in s1]
    ok['T5c'] &= sp1(phc, amps1(comp(A0, 2 * N), comp(A1, 2 * N), N), N, one) == [x * mu for x in s1]

# (T1) (L3): N[a log x - x, x^(M-1); a] = a^(M-1) Poch_M(z/a) exp(Sum_k b_k (z/a)^(2k-1)) (normalised: divide by a^(M-1))
def l3(a, M, n, one, variant=None):
    zero = one - one; apw = [one]
    for k in range(2 * n + 3): apw.append(apw[-1] * a)
    ph = [zero, zero] + [one * Fr((-1) ** (k + 1), k) / apw[k - 1] for k in range(2, 2 * n + 3)]
    amp = [one * gbinom(M - 1, j) / apw[j] for j in range(2 * n + 1)]
    S = sp1(ph, [amp] + [[] for _ in range(n)], n, one); S = [x / S[0] for x in S]
    B = bernoulli(n + 2); bk = [Fr(0)] + [B[2 * k] / (2 * k * (2 * k - 1)) for k in range(1, n // 2 + 2)]
    if variant == 'B2k': bk = [Fr(0)] + [B[2 * k] / (2 * k) for k in range(1, n // 2 + 2)]
    sg = -1 if variant == 'neg' else 1
    lg = [zero] * (n + 1)
    for k in range(1, n // 2 + 2):
        if 2 * k - 1 <= n: lg[2 * k - 1] = one * (sg * bk[k]) / apw[2 * k - 1]
    rhs = ser_exp(lg, n, one)
    if variant != 'nopoch':
        if M >= 0:
            for i in range(M): rhs = ser_mul(rhs, [one, one * i / a], n, zero)
        else:
            for i in range(1, -M + 1): rhs = ser_mul(rhs, [one * Fr(i) ** k / apw[k] for k in range(n + 1)], n, zero)
    return S == rhs, [k for k in range(n + 1) if S[k] != rhs[k]][:1]
NL = max(2 * N + 4, 12); okl = True; okq = True
for a in [Fr(1), Fr(7, 3), Fr(-2, 5)]:
    for M in range(-3, 6): okl &= l3(a, M, NL if a == 1 else 12, one)[0]
for M in range(-3, 6): okq &= l3(QZ.zeta(1) - 2, M, 8, QZ.of(1))[0]
check(ok['T0'], '(T0) one-variable SP: Lagrange-Buermann (sp1) = Wick definition (wickN, differential operator), %d seeds, z^0..z^%d' % (SEEDS, N))
check(okl, '(T1) (L3) exact over Q: M = -3..5, a = 1 (through z^%d), 7/3, -2/5 (through z^12)' % NL)
check(okq, '(T1) (L3) exact over Q(zeta): M = -3..5, a = zeta - 2, through z^8')
for var, txt, M in [('B2k', 'b_k -> B_2k/(2k)', 0), ('nopoch', 'Pochhammer factor dropped', 3), ('neg', 'exp(-Sum b_k ...)', 1)]:
    r, first = l3(Fr(1), M, 10, one, var)
    check(r, '(T1) %s (M = %d) fails, first at z^%s' % (txt, M, first[0] if first else '-'), ctrl=True)
check(ok['T2'], '(T2) (L2): N[s G + Fh, A(xi, z)] = A(X_0, z) EXACTLY (no z-corrections), %d seeds, through z^%d' % (SEEDS, N))
check(ok['T2c1'], '(T2) amplitude with an s xi term gives the same series', ctrl=True)
check(ok['T2c2'], '(T2) phase with an s^2 xi term (not linear in sigma) gives the same series', ctrl=True)
check(ok['T3'], '(T3) (L1b) Fubini: N[Phi, A] = N[phi, beta B] (inner series over Q[[w]], inner point critical identically in w), %d seeds, through z^%d' % (SEEDS, N))
check(ok['T3c'], '(T3) Fubini without the factor beta', ctrl=True)
check(ok['T3c2'], '(T3) Fubini with the inner critical point frozen at v = 0', ctrl=True)
check(ok['T4'], '(T4) (S1): N[Phi, (d_v Phi) h + z d_v h] = 0, %d seeds, through z^%d' % (SEEDS, N))
check(ok['T4c'], '(T4) N[Phi, (d_v Phi) h] alone vanishes', ctrl=True)
check(ok['T5'], "(T5) (L1a): N[ph o chi, (A o chi) chi'] = chi'(0) N[ph, A] for chi = mu t + c2 t^2 + c3 t^3, %d seeds" % SEEDS)
check(ok['T5c'], "(T5) without the Jacobian chi'", ctrl=True)
out('RF_LEMMAS %d %d: %s   (%.0f s)' % (N, SEEDS, 'ALL CHECKS PASSED' if not FAILS else 'FAILURES: %s' % FAILS, time.time() - T0))
LOG.close()
