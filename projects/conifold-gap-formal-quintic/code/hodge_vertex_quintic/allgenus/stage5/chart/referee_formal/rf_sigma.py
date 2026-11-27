# -*- coding: utf-8 -*-
"""rf_sigma.py -- referee check of c(z; L), the flat one-variable full-Stirling saddle series (GAMMA s.2, CHART_PROOF s.1.6).

Integrand (measure d sigma): exp(S_0/z) exp(S_1) exp(Sum_k z^(2k-1) S_2k) at sigma = L + t, built here from first principles:
  S_0'' = 5/sigma + Sum_j 1/(zeta^j - sigma);   S_1' = (1/2)[1/sigma + Sum_j 1/(zeta^j - sigma)];
  S_2k = b_k [(5 sigma)^(1-2k) + Sum_j (zeta^j - sigma)^(1-2k)].
SP by the Lagrange-Buermann routine rf_core.sp1 (not Wick). Checks:
  (A) rational L: c_1 = -(3/20) L^4, c_2 = (9/800) L^8, stored c_5(2) = 11495493197/31500;
  (B) L = 1 + eps in truncated Laurent series (rf_rings.LS): every pole eps^-n..eps^-1 of c_n cancels and c_n(1) = h_n,
      for n <= N, with NO assumption on the shape of c_n (this is Theorem Gamma (iii), tested directly at L = 1);
  (C) controls: Stirling term of Gamma(1+5s) dropped; b_1 + 1/7 in the j = 0 term only (must leave a pole);
      Gamma(5s) for Gamma(1+5s); comparison with the twist +5 series h;
  (D) Z/5 covariance in Q(zeta): c_n(zeta^i L) = zeta^(-in) c_n(L), n <= NC; control zeta^(+in).
Usage: py -I -B rf_sigma.py N NC     (log logs/rf_sigma_N_NC.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from rf_core import QZ, bk_list, gbinom, ser_exp, ser_mul, sp1, newton_sums, h_series
from rf_rings import LS


def bexp(X, N, one):
    """A[m][j] = [z^m t^j] exp(Sum_m z^m X_m(t)), X[0][0] = 0, j <= 2(N-m)"""
    zero = one - one; T = lambda m: 2 * (N - m)
    E0 = ser_exp(X[0], T(0), one); Z = [[one] + [zero] * T(0)]
    for m in range(1, N + 1):
        acc = [zero] * (T(m) + 1)
        for i in range(1, m + 1):
            if not X[i]: continue
            pr = ser_mul(X[i], Z[m - i], T(m), zero)
            for a in range(T(m) + 1): acc[a] = acc[a] + pr[a] * i
        Z.append([x * Fr(1, m) for x in acc])
    return [ser_mul(E0, Z[m], T(m), zero) for m in range(N + 1)]


def sigma_data(L, P, Pj0, N, one, variant=None):
    """Taylor data at sigma = L. P[k] = Sum_(j=0..4) (zeta^j - L)^-k, Pj0[k] = (1 - L)^-k (only for the 'b1j0' control)."""
    zero = one - one; K = 2 * N + 2; iL = 1 / L; ipw = [one]
    for k in range(K + 2): ipw.append(ipw[-1] * iL)
    ph = [zero, zero]
    for i in range(K - 1):                                      # S_0''(L+t) = Sum_i [5(-1)^i L^(-i-1) + P_(i+1)] t^i
        ph.append((ipw[i + 1] * (5 * (-1) ** i) + P[i + 1]) * Fr(1, (i + 1) * (i + 2)))
    X = [[zero] * (2 * N + 1)] + [[] for _ in range(N)]
    for i in range(2 * N):                                      # S_1(L+t) - S_1(L)
        sg = -1 if variant == 'g5s' else 1                     # Gamma(5s): S_1 - log(5 sigma)
        X[0][i + 1] = (ipw[i + 1] * (sg * (-1) ** i) + P[i + 1]) * Fr(1, 2 * (i + 1))
    b = bk_list(N // 2 + 2)
    for k in range(1, N // 2 + 2):
        m = 2 * k - 1
        if m > N: break
        row = []
        for r in range(2 * (N - m) + 1):
            g = gbinom(1 - 2 * k, r)
            y5 = zero if variant == 'noy' else ipw[2 * k - 1 + r] * (Fr(5) ** (1 - 2 * k))
            row.append((y5 + P[2 * k - 1 + r] * (-1) ** r) * (b[k] * g))
        if variant == 'b1j0' and k == 1:                       # extra (1/7) z (1 - sigma)^-1 in the j = 0 term
            row = [row[r] + Pj0[1 + r] * Fr(1, 7) for r in range(len(row))]
        X[m] = row
    return ph, bexp(X, N, one)


def c_series(L, P, Pj0, N, one, variant=None):
    ph, amps = sigma_data(L, P, Pj0, N, one, variant)
    Nn = sp1(ph, amps, N, one); n0 = Nn[0]
    return [x / n0 for x in Nn]


def P_rational(L, K):
    """P_k(L) for rational L via Q(zeta) (asserted rational); Pj0_k = (1-L)^-k"""
    ys = [1 / (QZ.zeta(j) - QZ.of(L)) for j in range(5)]; pw = [QZ.of(1)] * 5; P = [Fr(5)]
    for k in range(1, K + 1):
        pw = [a * y for a, y in zip(pw, ys)]; s = pw[0] + pw[1] + pw[2] + pw[3] + pw[4]
        assert s.rat() is not None; P.append(s.rat())
    return P, [Fr(1) / (1 - Fr(L)) ** k for k in range(K + 1)]


def P_qz(L, K):
    ys = [1 / (QZ.zeta(j) - L) for j in range(5)]; pw = [QZ.of(1)] * 5; P = [QZ.of(5)]
    for k in range(1, K + 1):
        pw = [a * y for a, y in zip(pw, ys)]; P.append(pw[0] + pw[1] + pw[2] + pw[3] + pw[4])
    return P


def P_eps(K, R):
    """P_k(1 + eps) = (-eps)^-k + Sum_r binom(k+r-1, r) eps^r Q_(k+r), Q_m = Sum_(j>=1) (zeta^j - 1)^-m = (-1)^m p_m;
    regular parts to absolute precision eps^R; Pj0_k = (-eps)^-k exactly"""
    p = newton_sums([Fr(-2), Fr(2), Fr(-1), Fr(1, 5)], K + R + 2)
    Q = [(-1) ** m * p[m] for m in range(len(p))]
    P = [LS.of(5)]; Pj0 = [LS.of(1)]
    for k in range(1, K + 1):
        reg = LS(0, [gbinom(k + r - 1, r) * Q[k + r] for r in range(R)], R)
        pole = LS.mono(-k, (-1) ** k); P.append(reg + pole); Pj0.append(pole)
    return P, Pj0


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10; NC = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    LOGREL = 'logs/rf_sigma_%d_%d.log' % (N, NC); LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8')
    FAILS = []; T0 = time.time()

    def out(s): print(s); LOG.write(s + '\n'); LOG.flush()

    def check(ok, msg, ctrl=False):
        good = (not ok) if ctrl else ok
        out(('PASS ' if good else 'FAIL ') + ('control: ' if ctrl else '') + msg)
        if not good: FAILS.append(msg)

    out('rf_sigma: N = %d, NC = %d; log %s' % (N, NC, LOGREL))
    h = h_series(N); h5 = h_series(N, twist=Fr(5))
    # h cross-check in Q(zeta) directly (no Newton sums): T_k = 5^(1-2k) + Sum_j (zeta^j - 1)^(1-2k) = -Sum_(a in w_0) a^(1-2k)
    b = bk_list(N // 2 + 1); lg = [QZ.of(0)] * (N + 1); okT = True
    ys = [1 / (QZ.zeta(j) - 1) for j in range(1, 5)]; ws = [1 / (1 - QZ.zeta(j)) for j in range(1, 5)]

    def qpow(x, e):
        r = QZ.of(1)
        for _ in range(e): r = r * x
        return r
    for k in range(1, (N + 1) // 2 + 1):
        e = 2 * k - 1
        T = QZ.of(Fr(5) ** (-e)) + qpow(ys[0], e) + qpow(ys[1], e) + qpow(ys[2], e) + qpow(ys[3], e)
        W = QZ.of(Fr(-5) ** (-e)) + qpow(ws[0], e) + qpow(ws[1], e) + qpow(ws[2], e) + qpow(ws[3], e)
        okT = okT and T.rat() is not None and T == -W
        if e <= N: lg[e] = T * b[k]
    hq = ser_exp(lg, N, QZ.of(1))
    check(okT and all(hq[n] == QZ.of(h[n]) for n in range(N + 1)),
          'h_(w_0) from the five Stirling factors in Q(zeta) (T_k rational, T_k = -Sum_(w_0) a^(1-2k)) = Newton-sum h, z^0..z^%d' % N)
    out('h_1 = %s, h_5 = %s' % (h[1], h[5] if N >= 5 else '-'))
    check(h[1] == Fr(-3, 20) and (N < 5 or h[5] == Fr(-257071, 8064000000)), 'h_1 = -3/20 and h_5 = -257071/8064000000')
    # (A) rational L
    for L in [Fr(2), Fr(-2, 3), Fr(3), Fr(7, 5)]:
        P, Pj0 = P_rational(L, 2 * N + 4); c = c_series(L, P, Pj0, min(N, 6), Fr(1))
        check(c[1] == Fr(-3, 20) * L ** 4 and c[2] == Fr(9, 800) * L ** 8, 'L = %s: c_1 = -(3/20) L^4, c_2 = (9/800) L^8' % L)
        if L == 2: check(c[5] == Fr(11495493197, 31500), 'stored c_5(2) = 11495493197/31500 reproduced (%s)' % c[5])
    # (B) L = 1 + eps
    R = N + 4; LS.CAP = R + 2 * N + 12; t1 = time.time()
    P, Pj0 = P_eps(2 * N + 4, R); Leps = LS(0, [1, 1])
    ce = c_series(Leps, P, Pj0, N, LS.of(1))
    reg = all(x.v >= 0 for x in ce); val = all(ce[n].coef(0) == h[n] for n in range(N + 1))
    out('     eps-expansion: valuations %s; precisions %s' % ([x.v for x in ce], [x.p for x in ce]))
    check(reg, 'L = 1 + eps: every c_n, n <= %d, has no pole at eps = 0 (all eps^-n..eps^-1 cancel)   %.1f s' % (N, time.time() - t1))
    check(val, 'L = 1 + eps: c_n(1) = [eps^0] c_n = h_n for all n <= %d: (Gamma) at L = 1 directly' % N)
    check(all(ce[n].coef(0) == h5[n] for n in range(N + 1)), 'c(z; 1) equals h with twist weight +5', ctrl=True)
    # (C) controls in the eps mode (smaller order)
    Nc = min(N, 6); Pc, Pj0c = P_eps(2 * Nc + 4, Nc + 4)
    cc = c_series(Leps, Pc, Pj0c, Nc, LS.of(1), 'noy')
    bad = [n for n in range(Nc + 1) if cc[n].v < 0 or cc[n].coef(0) != h[n]]
    check(not bad, 'Stirling term of Gamma(1+5s) dropped: c(z; 1) = h fails, first at z^%s' % (bad[0] if bad else '-'), ctrl=True)
    # Gamma(5s) for Gamma(1+5s) (S_1 - log 5 sigma; the measure dy/y): BLIND at L = 1 (Gamma(1+w) = w Gamma(w), w = 5/z,
    # changes only the elementary factor; (L3): Poch_0 = Poch_1 = 1), but visible at L != 1 from z^2 on (critic's remark)
    cc = c_series(Leps, Pc, Pj0c, Nc, LS.of(1), 'g5s')
    check(all(cc[n].v >= 0 and cc[n].coef(0) == h[n] for n in range(Nc + 1)),
          'Gamma(5s) for Gamma(1+5s): still regular with value h at L = 1 through z^%d (blind there, as predicted)' % Nc)
    P2, P2j = P_rational(Fr(2), 2 * Nc + 4)
    ca = c_series(Fr(2), P2, P2j, Nc, Fr(1)); cb = c_series(Fr(2), P2, P2j, Nc, Fr(1), 'g5s')
    bad = [n for n in range(Nc + 1) if ca[n] != cb[n]]
    check(not bad, 'Gamma(5s) for Gamma(1+5s) at L = 2: c(z; 2) changes, first at z^%s' % (bad[0] if bad else '-'), ctrl=True)
    cc = c_series(Leps, Pc, Pj0c, Nc, LS.of(1), 'b1j0')
    bad = [n for n in range(Nc + 1) if cc[n].v < 0]
    check(not bad, 'b_1 + 1/7 in the j = 0 Stirling term only: c_n keeps a pole at L = 1, first n = %s (valuation %s)'
          % (bad[0] if bad else '-', cc[bad[0]].v if bad else '-'), ctrl=True)
    # (D) Z/5 covariance in Q(zeta)
    t1 = time.time(); okc = True; okbad = True
    for L in [Fr(2), Fr(-1, 3)]:
        base = c_series(QZ.of(L), P_qz(QZ.of(L), 2 * NC + 4), None, NC, QZ.of(1))
        for i in range(1, 5):
            Li = QZ.zeta(i) * L; ci = c_series(Li, P_qz(Li, 2 * NC + 4), None, NC, QZ.of(1))
            okc = okc and all(ci[n] == base[n] * QZ.zeta(-i * n) for n in range(NC + 1))
            okbad = okbad and all(ci[n] == base[n] * QZ.zeta(i * n) for n in range(NC + 1))
    check(okc, 'Z/5: c_n(zeta^i L) = zeta^(-in) c_n(L) in Q(zeta), i = 1..4, n <= %d, L = 2, -1/3   %.1f s' % (NC, time.time() - t1))
    check(okbad, 'Z/5 with zeta^(+in)', ctrl=True)
    out('RF_SIGMA %d %d: %s   (%.0f s)' % (N, NC, 'ALL CHECKS PASSED' if not FAILS else 'FAILURES: %s' % FAILS, time.time() - T0))
    LOG.close()
