# -*- coding: utf-8 -*-
"""rf_hat.py -- referee check of the hatted/flat identity (CHART_PROOF s.3.10) and of (Gamma) in hatted form (poly's "(A2)").

POLY s.5 integrand, re-derived here from the summand q^s Gamma(1+5s)/prod_j Gamma(1 + s - zeta^j/z) (1/Gamma form, q^s, NOT
GAMMA's form), hatted with sigma = L/X, z = L^-4 w, u = L^-5, at X = 1 + t, coefficients in Q[[u]]/(u^(U+1)):
  phase Ph with Ph'(X) = (1/(u X^2)) log((1 - u X^5)/(1 - u))   [= d/dX (L^4 S(L/X)), S = POLY's flat phase];
  amplitude (1 - u X^5)^(-1/2) exp(Sum_k b_k (u X w)^(2k-1) [5^(1-2k) - 5 Sum_P binom(2k-2+5P, 5P) u^P X^(5P)]) dX.
ph_n(u) := [w^n] normalised SP (rf_core.sp1 over the ring TS). Checks:
  (H1) [u^k] ph_n = 0 for floor(4n/5) < k <= U   (polynomiality, Theorem P(a), to the truncation);
  (H2) ph_n(1) = h_n                            ((Gamma) in hatted form = poly's "(A2)");
  (H3) L^(4n) ph_n(L^-5) = c_n(L) at L = 2, -2/3, 3, c_n from the FLAT sigma-form of rf_sigma (GAMMA's form) (s.3.10);
  (H4) ph_n(0) = the model's normalised SP of exp(-(X^4/4 + 1/X)/w) dX at X = 1;
  controls: twist weight +5 (5^(1-2k) -> -5^(1-2k)); amplitude factor (1 - u X^5)^(-1/2) dropped.
Usage: py -I -B rf_hat.py N      (log logs/rf_hat_N.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from rf_core import bk_list, sp1, h_series
from rf_rings import TS
from rf_sigma import bexp, c_series, P_rational


def tpoly_mul(a, b, K):
    r = [0] * (K + 1)
    for i, x in enumerate(a):
        if x and i <= K:
            for j in range(min(len(b), K + 1 - i)): r[i + j] += x * b[j]
    return r


def hat_data(N, U, variant=None):
    TS.n = U + 1; K = 2 * N + 1; one = TS.of(1); zero = TS.of(0)
    delta = [0, 5, 10, 10, 5, 1]; dm = [[1]]
    for m in range(1, K + 2): dm.append(tpoly_mul(dm[-1], delta, K))
    D = [None] + [TS([0] * (m - 1) + [Fr(comb(m + i - 1, i), m) for i in range(U + 1)]) for m in range(1, K + 2)]
    inv2 = [(-1) ** i * (i + 1) for i in range(K + 1)]            # (1+t)^-2
    S = [zero] * (K + 1)                                            # Sum_m D_m(u) delta^m(t)
    for m in range(1, min(K, U + 1) + 1):
        for i, x in enumerate(dm[m]):
            if x: S[i] = S[i] + D[m] * x
    Phd = [zero] * (K + 1)                                          # Ph'(1+t) = -(1+t)^-2 S(t)
    for i in range(K + 1):
        for j in range(K + 1 - i): Phd[i + j] = Phd[i + j] - S[j] * inv2[i]
    ph = [zero, Phd[0]] + [Phd[k] * Fr(1, k + 1) for k in range(1, K + 1)]
    assert Phd[0] == 0
    uu = TS.var(); X = [[zero] * (2 * N + 1)] + [[] for _ in range(N)]
    if variant != 'nor':                                           # X_0 = -(1/2) log(1 - u delta/(1-u)) = (u/2) S
        for i in range(1, 2 * N + 1): X[0][i] = uu * S[i] * Fr(1, 2)
    b = bk_list(N // 2 + 2)
    for k in range(1, N // 2 + 2):
        m = 2 * k - 1
        if m > N: break
        s5 = Fr(5) ** (1 - 2 * k) * (-1 if variant == 'tw' else 1)
        row = []
        for r in range(2 * (N - m) + 1):
            c = [Fr(0)] * (U + 1)
            if m <= U: c[m] += b[k] * comb(m, r) * (s5 - 5)
            for P in range(1, U + 1 - m):
                c[m + P] -= 5 * b[k] * comb(2 * k - 2 + 5 * P, 5 * P) * comb(m + 5 * P, r)
            row.append(TS(c))
        X[m] = row
    return ph, bexp(X, N, one)


def hat_series(N, U, variant=None):
    ph, amps = hat_data(N, U, variant); Nn = sp1(ph, amps, N, TS.of(1))
    return [x / Nn[0] for x in Nn]


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10; U = (4 * N) // 5 + 4
    LOGREL = 'logs/rf_hat_%d.log' % N; LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8')
    FAILS = []; T0 = time.time()

    def out(s): print(s); LOG.write(s + '\n'); LOG.flush()

    def check(ok, msg, ctrl=False):
        good = (not ok) if ctrl else ok
        out(('PASS ' if good else 'FAIL ') + ('control: ' if ctrl else '') + msg)
        if not good: FAILS.append(msg)

    out('rf_hat: hatted SP in Q[[u]]/(u^%d) through w^%d; log %s' % (U + 1, N, LOGREL))
    ph = hat_series(N, U); h = h_series(N)
    deg = [max([k for k in range(U + 1) if ph[n].c[k] != 0] or [-1]) for n in range(N + 1)]
    check(all(deg[n] <= (4 * n) // 5 for n in range(N + 1)),
          '(H1) deg ph_n <= floor(4n/5) (coefficients u^(floor(4n/5)+1)..u^%d vanish), n <= %d   %.1f s' % (U, N, time.time() - T0))
    out('     degrees: %s' % deg)
    check(all(sum(ph[n].c) == h[n] for n in range(N + 1)), '(H2) ph_n(1) = h_n for all n <= %d (poly\'s "(A2)")' % N)
    out('     ph_5(u) = %s' % ' + '.join('(%s) u^%d' % (ph[5].c[k], k) for k in range(deg[5] + 1)) if N >= 5 else '')
    okf = True; NF = min(N, 12)
    for L in [Fr(2), Fr(-2, 3), Fr(3)]:
        P, Pj0 = P_rational(L, 2 * NF + 4); c = c_series(L, P, Pj0, NF, Fr(1)); u = 1 / L ** 5
        for n in range(NF + 1):
            val = sum(ph[n].c[k] * u ** k for k in range(deg[n] + 1)) * L ** (4 * n)
            okf = okf and val == c[n]
    check(okf, '(H3) L^(4n) ph_n(L^-5) = c_n(L) (flat GAMMA form, rf_sigma) at L = 2, -2/3, 3, n <= %d (s.3.10)' % NF)
    # (H4) model at u = 0
    K = 2 * NF + 2; mph = [Fr(0), Fr(0)] + [-(Fr(comb(4, k), 4) + (-1) ** k) for k in range(2, K + 1)]
    mod = sp1(mph, [[Fr(1)] + [Fr(0)] * (2 * NF)] + [[Fr(0)] * (2 * (NF - m) + 1) for m in range(1, NF + 1)], NF, Fr(1))
    check(all(ph[n].c[0] == mod[n] / mod[0] for n in range(NF + 1)),
          '(H4) ph_n(0) = model SP of exp(-(X^4/4 + 1/X)/w) at X = 1, n <= %d; F_(0,1..3) = %s' % (NF, [str(mod[n] / mod[0]) for n in (1, 2, 3)]))
    NCt = min(N, 8); Uc = (4 * NCt) // 5 + 4
    for var, txt in [('tw', 'twist weight +5 in the hatted amplitude'), ('nor', 'factor (1 - u X^5)^(-1/2) dropped')]:
        pc = hat_series(NCt, Uc, var); hc = h_series(NCt)
        degc = [max([k for k in range(Uc + 1) if pc[n].c[k] != 0] or [-1]) for n in range(NCt + 1)]
        bad = [n for n in range(NCt + 1) if degc[n] > (4 * n) // 5 or sum(pc[n].c) != hc[n]]
        check(not bad, '%s: (H1)+(H2) fail, first at w^%s (polynomial: %s)'
              % (txt, bad[0] if bad else '-', all(degc[n] <= (4 * n) // 5 for n in range(NCt + 1))), ctrl=True)
    TS.n = U + 1
    out('RF_HAT %d: %s   (%.0f s)' % (N, 'ALL CHECKS PASSED' if not FAILS else 'FAILURES: %s' % FAILS, time.time() - T0))
    LOG.close()
