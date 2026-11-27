"""lf_lmi N1 N2: Theorem LMI (REFEREE s.3.1) on seeded random regular data, exact over Q; (L2') and the moments.

L1  the inner step of step 1: N[a nu - e^nu, e^(M nu); log a] = a^M Poch_M(z/a) exp(Sum b_k (z/a)^(2k-1)), z^N2.
L2  g = 0, F = 1: c^F(z; L) = 1 EXACTLY at every L (the s*-free content of LMI), M0 = -2..3, z^N2.
L3  random g (degree 4), random log F, (alpha0, M0) with M0 in {-1, 0, 1, 2}:
    (a) two-variable chart = one-variable c^F at L = s* + eps0, eps0 = 1/3, -2/5 (steps 1-2), z^N1;
    (b) chart at Q = 0 = exp(gt(delta)/z) F(s* + delta, z)/F(s*, 0), delta = -M0 z/alpha0 (step 3), z^N1;
    (c) Method A, sigma = L + eps tau, z = eps w over Q[eps]/(eps^(N2+1)): c^F_n(s* + eps) = c'_n(eps)/eps^n has
        no pole ([eps^k] c'_n = 0, k < n) and its value is the LMI formula, n <= N2.
L4  (L2') with the generic (non-separated) two-variable engine, c in {0, 1, -1, 2, 1/5, -3/2}; moment tables
    against the Isserlis recursion; the generating-function identity of REFEREE s.3.2.
L5  floats: the chart critical point at Q(L) = a0(L) e^(g'(L)/alpha0) is (L, -g'(L)/alpha0); branch independence of Q.
Failing controls: wrong-sign delta; singular Stirling series dropped; Poch dropped; singular b_1 + 1/7; exponent M0 for
M0 - 1/2; b_k -> B_2k/(2k) in L1; a pure m^3 term and an s^2 m term in (L2'); Q e^(2 mu).
"""
import cmath
import math
import os
import random
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import factorial
from lf_core import bk, bern, ser_exp, ser_mul, TS
from lf_wick import wick_1, wick_2sep, wick_2gen, moments_2, cov_2
from lf_models import lmi_onevar, lmi_chart, lmi_scaled, lmi_predict
from lf_log import Log, first_diff

Z, O = Fr(0), Fr(1)


def rnd(rng, den=7):
    return Fr(rng.randint(-9, 9), rng.randint(1, den))


def l3_inner(a, M, n, bern_wrong=False, poch=True):
    a = Fr(a)
    ph3 = {k: -a / factorial(k) for k in range(3, 2 * n + 3)}
    N = wick_1(ph3, {(1, 0): Fr(M)} if M else {}, 1 / a, n, Z, O)
    s = [Fr(0)] * (n + 1)
    for k in range(1, (n + 1) // 2 + 1):
        s[2 * k - 1] += (bern(2 * k) / (2 * k) if bern_wrong else bk(k)) * a ** (1 - 2 * k)
    if poch:
        for i in range(1, M):
            for p in range(1, n + 1):
                s[p] += Fr((-1) ** (p + 1) * i ** p, p) * a ** (-p)
        for i in range(1, -M + 1):
            for p in range(1, n + 1):
                s[p] += Fr(i ** p, p) * a ** (-p)
    return N, ser_exp(s, n + 1)


def chart_series(data, n):
    phs, las, phm, lam, H = data
    Css, Csm, Cmm, det = cov_2(*H)
    return wick_2sep(phs, las, phm, lam, (Css, Csm, Cmm), n, Z, O)


def method_a(G, LF, alpha0, M0, n, **ctrl):
    ph3, logA, C = lmi_scaled(G, LF, alpha0, M0, n, **ctrl)
    cp = wick_1(ph3, logA, C, n, TS([]), TS.const(1))
    poles = [(m, k) for m in range(n + 1) for k in range(m) if cp[m].c[k] != 0]
    return poles, [cp[m].c[m] for m in range(n + 1)]


def isserlis(Css, Csm, Cmm, amax, bmax):
    M = {}
    for a in range(amax + 1):
        for b in range(bmax + 1):
            if a == 0 and b == 0:
                M[a, b] = Fr(1)
            elif (a + b) % 2:
                M[a, b] = Fr(0)
            elif a > 0:      # pull one eta_s: pair it with an eta_s (a - 1 ways) or an eta_m (b ways)
                M[a, b] = (a - 1) * Css * M.get((a - 2, b), Fr(0)) + b * Csm * M.get((a - 1, b - 1), Fr(0))
            else:
                M[a, b] = (b - 1) * Cmm * M[a, b - 2]
    return M


def main():
    N1, N2 = int(sys.argv[1]), int(sys.argv[2])
    log = Log("lf_lmi_%d_%d" % (N1, N2))
    # L1
    ok = True
    for a in [1, Fr(-2, 7), Fr(3, 2), 5]:
        for M in range(-3, 4):
            N, pred = l3_inner(a, M, N2)
            ok &= (N == pred)
    log.check("L1 inner step (L3 in nu-form), a in {1, -2/7, 3/2, 5}, M = -3..3, through z^%d" % N2, ok)
    N, pred = l3_inner(Fr(3, 2), 2, N2, bern_wrong=True)
    log.control("L1 b_k -> B_2k/(2k)", N == pred, "first diff z^%s" % first_diff(N, pred))
    N, pred = l3_inner(Fr(3, 2), 2, N2, poch=False)
    log.control("L1 Poch_2 dropped (M = 2)", N == pred, "first diff z^%s" % first_diff(N, pred))
    # L2
    ok = True
    one = [O] + [Z] * N2
    for alpha0 in [1, -1, 5, Fr(-1, 3), Fr(7, 2)]:
        for M0 in range(-2, 4):
            for eps0 in [Fr(1, 2), Fr(-3)]:
                ph3, logA, C = lmi_onevar({}, {}, alpha0, M0, eps0, N2)
                ok &= (wick_1(ph3, logA, C, N2, Z, O) == one)
    log.check("L2 g = 0, F = 1: c^F(z; L) = 1 exactly, 5 alpha0 x 6 M0 x 2 L, through z^%d" % N2, ok)
    for lab, ctrl, M0 in [("singular b_1 + 1/7", dict(b1_shift=Fr(1, 7)), 0), ("exponent M0 for M0 - 1/2", dict(expo=1), 1),
                          ("Poch dropped, M0 = 2", dict(poch=False), 2), ("Poch dropped, M0 = -1", dict(poch=False), -1),
                          ("singular Stirling dropped", dict(stirling=False), 0)]:
        ph3, logA, C = lmi_onevar({}, {}, Fr(-1, 3), M0, Fr(1, 2), N2, **ctrl)
        c = wick_1(ph3, logA, C, N2, Z, O)
        log.control("L2 " + lab, c == one, "first diff z^%s" % first_diff(c, one))
    # L3
    rng = random.Random(20261114)
    cases = [(-1, 0), (-1, 1), (-1, -1), (-1, 2), (2, 1), (2, -1), (Fr(-1, 3), 2), (Fr(5, 2), 0), (5, 1), (Fr(-7, 4), -1)]
    for idx, (alpha0, M0) in enumerate(cases):
        G = {k: rnd(rng) for k in (2, 3, 4)}
        LF = {key: rnd(rng) for key in [(1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (0, 2), (2, 1)]}
        tag = "L3 case %d (alpha0 = %s, M0 = %s)" % (idx, alpha0, M0)
        for eps0 in [Fr(1, 3), Fr(-2, 5)]:
            ph3, logA, C = lmi_onevar(G, LF, alpha0, M0, eps0, N1)
            c1 = wick_1(ph3, logA, C, N1, Z, O)
            ch = chart_series(lmi_chart(G, LF, alpha0, M0, eps0, N1), N1)
            log.check("%s (a) chart = c^F at eps0 = %s, z^%d" % (tag, eps0, N1), ch == c1)
            if idx == 1 and eps0 == Fr(1, 3):
                chx = chart_series(lmi_chart(G, LF, alpha0, M0, eps0, N1, mu_scale=2), N1)
                log.control("%s chart with Q e^(2 mu)" % tag, chx == c1, "first diff z^%s" % first_diff(chx, c1))
        ch0 = chart_series(lmi_chart(G, LF, alpha0, M0, Fr(0), N1), N1)
        pred = lmi_predict(G, LF, alpha0, M0, N1)
        log.check("%s (b) chart at Q = 0 = LMI formula, z^%d" % (tag, N1), ch0 == pred)
        if M0 != 0:
            predx = lmi_predict(G, LF, alpha0, M0, N1, sign=-1)
            log.control("%s delta = +M0 z/alpha0" % tag, ch0 == predx, "first diff z^%s" % first_diff(ch0, predx))
        poles, val = method_a(G, LF, alpha0, M0, N2)
        pred2 = lmi_predict(G, LF, alpha0, M0, N2)
        log.check("%s (c) Method A: no pole, n <= %d" % (tag, N2), not poles)
        log.check("%s (c) Method A: value c^F(z; s*) = LMI formula, n <= %d" % (tag, N2), val == pred2)
        if idx in (0, 3, 5):
            for lab, ctrl in [("singular Stirling dropped", dict(stirling=False)), ("singular b_1 + 1/7", dict(b1_shift=Fr(1, 7))),
                              ("exponent M0 for M0 - 1/2", dict(expo=M0))] + ([("Poch dropped", dict(poch=False))] if M0 not in (0, 1) else []):
                px, vx = method_a(G, LF, alpha0, M0, N2, **ctrl)
                log.control("%s Method A, %s: regular" % (tag, lab), not px,
                            "first pole (n, eps^k) = %s" % (px[0],) if px else "")
    # L4
    ok = True
    for seed in range(4):
        r4 = random.Random(1000 + seed)
        Css, Csm, Cmm = rnd(r4), rnd(r4), rnd(r4)
        Mt = moments_2(Css, Csm, Cmm, 24, 24, Z, O)
        Mi = isserlis(Css, Csm, Cmm, 24, 24)
        ok &= all(Mt[a][b] == Mi[a, b] for a in range(25) for b in range(25))
    log.check("L4 moment table = Isserlis recursion, 4 random covariances, degrees <= 24", ok)
    ok = True
    Csm, Cmm = Fr(-3, 7), Fr(5, 11)
    Mt = moments_2(Fr(0), Csm, Cmm, 30, 30, Z, O)
    for a in range(31):
        for r in range(31):
            j2 = r - a
            want = Fr(0) if (j2 < 0 or j2 % 2) else Fr(factorial(r), 2 ** (j2 // 2) * factorial(j2 // 2)) * Csm ** a * Cmm ** (j2 // 2)
            ok &= (Mt[a][r] == want)
    log.check("L4 C_ss = 0: <eta_m^r eta_s^a> = r!/(2^j j!) C_sm^a C_mm^j, j = (r - a)/2 (generating-function step), r, a <= 30", ok)
    for seed, cc in enumerate([Fr(0), Fr(1), Fr(-1), Fr(2), Fr(1, 5), Fr(-3, 2)]):
        r4 = random.Random(2000 + seed)
        alpha0 = rnd(r4, 3) or Fr(2)
        P = {k: rnd(r4) for k in range(2, 6)}
        LA = {(d, m): rnd(r4) for (d, m) in [(1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (0, 2)]}
        Phi3 = {(k, 0): v for k, v in P.items() if k >= 3}
        logA = {(d, 0, m): v for (d, m), v in LA.items()}
        if cc:
            logA[(0, 1, 0)] = cc
        cov = cov_2(2 * P[2], alpha0, Fr(0))[:3]
        Ng = wick_2gen(Phi3, logA, cov, N1, Z, O)
        pred = lmi_predict(P, LA, alpha0, cc, N1)          # exp(P(delta)/z) A(delta, z)/A(0, 0), delta = -c z/alpha0
        log.check("L4 (L2') generic engine, c = %s, alpha0 = %s, z^%d" % (cc, alpha0, N1), Ng == pred)
        if seed in (1, 4):
            log.control("L4 (L2') with delta -> -delta, c = %s" % cc, Ng == lmi_predict(P, LA, alpha0, cc, N1, sign=-1))
            Ph = dict(Phi3)
            Ph[(0, 3)] = Fr(2, 3)
            log.control("L4 (L2') plus a pure m^3 term, c = %s" % cc, wick_2gen(Ph, logA, cov, N1, Z, O) == pred)
            Ph = dict(Phi3)
            Ph[(2, 1)] = Fr(-1, 4)
            log.control("L4 (L2') plus an s^2 m term, c = %s" % cc, wick_2gen(Ph, logA, cov, N1, Z, O) == pred)
    # L5 (floats)
    r5 = random.Random(555)
    worst, worst_b = 0.0, 0.0
    for trial in range(6):
        sst, alpha0 = r5.uniform(-1, 1), r5.choice([-1.0, 2.0, -0.4, 5.0])
        Gf = {k: r5.uniform(-2, 2) for k in (1, 2, 3, 4)}
        gp = lambda s: sum(k * v * (s - sst) ** (k - 1) for k, v in Gf.items())
        gpp = lambda s: sum(k * (k - 1) * v * (s - sst) ** (k - 2) for k, v in Gf.items() if k >= 2)
        for e in [1e-3, 1e-2, 5e-2]:
            L = sst + e / alpha0                             # a0(L) = e > 0
            Q = e * math.exp(gp(L) / alpha0)
            s, mu = sst, -gp(sst) / alpha0                   # Newton from the Q = 0 critical point
            for _ in range(60):
                F1, F2 = gp(s) + alpha0 * mu, alpha0 * (s - sst) - Q * math.exp(mu)
                J = [[gpp(s), alpha0], [alpha0, -Q * math.exp(mu)]]
                dt = J[0][0] * J[1][1] - J[0][1] * J[1][0]
                s, mu = s - (J[1][1] * F1 - J[0][1] * F2) / dt, mu - (-J[1][0] * F1 + J[0][0] * F2) / dt
            worst = max(worst, abs(s - L), abs(mu + gp(L) / alpha0))
            lam0 = -gp(L) - alpha0 * math.log(e)
            for kk in (-2, 1, 3):
                lamk = lam0 - alpha0 * 2j * math.pi * kk
                worst_b = max(worst_b, abs(cmath.exp(-lamk / alpha0) - Q))
    log.check("L5 critical point of Psi_Q(L) is (L, -g'(L)/alpha0) (Newton from the Q = 0 point), max error %.1e" % worst,
              worst < 1e-10)
    log.check("L5 Q = exp(-lambda/alpha0) is the same for lambda + 2 pi i k alpha0, max error %.1e" % worst_b, worst_b < 1e-12)
    log.close()


if __name__ == "__main__":
    main()
