# -*- coding: utf-8 -*-
"""lm_chart.py -- referee lmi_adversary: the two-variable chart of LMI step 2, as an independent Wick computation.

Chart (RA REFEREE.md s.3.3): Psi_Q = g(sigma) + alpha0 (sigma - s*) mu - Q e^mu, amplitude e^(M0 mu) F(sigma, z).
At the critical point belonging to L = s* + e one has mu_c = -g'(L)/alpha0 and Q e^(mu_c) = alpha0 e EXACTLY, so every
Taylor datum is a polynomial in e (no log, no Q^(1/alpha0)):
  [s^k] = G_k(L) (k >= 2),  [s m] = alpha0,  [m^k] = -alpha0 e/k! (k >= 2),  no other mixed terms;
  covariance C = -H^-1: C_ss = -e/D, C_sm = -1/D, C_mm = 2 G_2(L)/(alpha0 D), D = alpha0 + 2 e G_2(L).
Since s and m interact only through C (separable phase and amplitude), N_n = Sum X_a Y_b <eta_s^a eta_m^b>.
At e = 0 this is exactly (L2') (C_ss = 0, no pure m-vertices).  Options: mcube (a constant added to [m^3], control),
pows (amplitude factors (sigma - s*)^q are NOT allowed here: the chart needs a holomorphic amplitude).
"""
from fractions import Fraction as Fr
from math import comb, factorial
from lm_core import LS, lsc, lsm, F, finv, fr, fpow, INF
from lm_wick import gexp, vertex_table
from lm_data import reg_coef


def chart_series(N, alpha0, M0, G, f, exact=True, T=None, mcube=None, normalise=True):
    W = 2 * N
    T = T if T is not None else N + 6
    ps = {k: reg_coef(G, k, exact) for k in range(3, W + 3) if len(G) > k}
    pm = {k: lsm(1, 1).scal(fr(-alpha0 * F(Fr(1, factorial(k))))) for k in range(3, W + 3)}
    if mcube is not None:
        pm[3] = pm[3] + lsc(mcube)
    G2 = reg_coef(G, 2, exact)
    D = lsc(1).scal(alpha0) + lsm(2, 1) * G2
    iD = D.inv(rel=T)
    Css = -(lsm(1, 1) * iD)
    Csm = -iD
    Cmm = (G2 * iD).scal(fr(F(2) * finv(alpha0)))
    # amplitude in s: F(L + s, z)/F(L, 0)
    H = {}
    for m, row in f.items():
        for j in range(0, W - 2 * m + 1):
            if (m, j) != (0, 0) and len(row) > j:
                H[(m, j)] = reg_coef(row, j, exact)
    a = gexp(H, W)
    evs = vertex_table(ps, W)
    evm = vertex_table(pm, W)
    X = [dict() for _ in range(W + 1)]       # X[w][deg]: coefficient of s^w eta_s^deg
    for (m, j), am in a.items():
        for w1 in range(0, W - 2 * m - j + 1):
            for V, ev in evs[w1].items():
                w, d = w1 + 2 * m + j, w1 + 2 * V + j
                t = am * ev
                X[w][d] = X[w][d] + t if d in X[w] else t
    Y = [dict() for _ in range(W + 1)]       # e^(M0 s eta_m) times the m-vertices
    for w2 in range(W + 1):
        for V, ev in evm[w2].items():
            for l in range(0, W - w2 + 1):
                w, d = w2 + l, w2 + 2 * V + l
                t = ev.scal(fr(fpow(M0, l) * F(Fr(1, factorial(l)))))   # M0: field element
                Y[w][d] = Y[w][d] + t if d in Y[w] else t
    pw = {}

    def cpow(name, C, k):
        key = (name, k)
        if key not in pw:
            pw[key] = lsc(1) if k == 0 else cpow(name, C, k - 1) * C
        return pw[key]
    mom = {}

    def moment(da, db):
        if (da, db) in mom:
            return mom[(da, db)]
        tot = LS(0, [])
        for k in range(da % 2, min(da, db) + 1, 2):
            if (db - k) % 2:
                continue
            i, l = (da - k) // 2, (db - k) // 2
            cnt = Fr(factorial(da) * factorial(db), factorial(k) * 2 ** i * factorial(i) * 2 ** l * factorial(l))
            tot = tot + (cpow('ss', Css, i) * cpow('mm', Cmm, l) * cpow('sm', Csm, k)).scal(F(cnt))
        mom[(da, db)] = tot
        return tot
    res = []
    for n in range(N + 1):
        tot = LS(0, [])
        for w1 in range(0, 2 * n + 1):
            w2 = 2 * n - w1
            for da, xa in X[w1].items():
                for db, yb in Y[w2].items():
                    if (da + db) % 2:
                        continue
                    tot = tot + (xa * moment(da, db)) * yb
        res.append(tot)
    if normalise:
        i0 = res[0].inv(rel=T)
        res = [x * i0 for x in res]
    return res
