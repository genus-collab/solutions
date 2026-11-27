# -*- coding: utf-8 -*-
"""lm_fals.py -- referee lmi_adversary: integrands OUTSIDE the hypotheses of Theorem LMI (falsifier builders).

* nonlinear_data: singular factor Gamma(a(sigma)/z) with a = alpha0 x + alpha1 x^2 (x = sigma - s*), exact Stirling form.
  The Taylor data at L = s* + e contain ell = log a(L) (only through [t^2]: + alpha1 ell); ell is a free parameter here.
* lam_extra: amplitude exp(kappa z lambda x), lambda = -g'(L) - alpha0 ell (the critical log(-q)): data depend on log Q.
* qdep_extra: amplitude exp(z Q' f1(sigma)), Q' = alpha0 e exp((g'(L) - g'(s*))/alpha0) = Q exp(-g'(s*)/alpha0):
  data depend HOLOMORPHICALLY on Q (positive control: LMI's chart is still holomorphic at Q = 0).
"""
from fractions import Fraction as Fr
from math import comb
from lm_core import LS, lsc, lsm, F, finv, fr, fpow, lsexp, INF
from lm_wick import gexp
from lm_data import reg_coef


def tz(K):
    return [LS(0, []) for _ in range(K + 1)]


def tmul(A, B, K):
    C = tz(K)
    for i, x in enumerate(A):
        if i > K or x.ezero():
            continue
        for j, y in enumerate(B):
            if i + j > K:
                break
            if not y.ezero():
                C[i + j] = C[i + j] + x * y
    return C


def tlog1p(Y, K):
    """log(1 + Y), Y[0] = 0, t-series with LS coefficients"""
    out, P = tz(K), [lsc(1)] + tz(K)[1:]
    for i in range(1, K + 1):
        P = tmul(P, Y, K)
        c = F(Fr((-1) ** (i + 1), i))
        out = [o + x.scal(c) for o, x in zip(out, P)]
    return out


def texp(Y, K):
    """exp(Y), Y[0] = 0"""
    E = [lsc(1)] + tz(K)[1:]
    for k in range(1, K + 1):
        s = LS(0, [])
        for i in range(1, k + 1):
            if not Y[i].ezero():
                s = s + (Y[i] * E[k - i]).scal(F(i))
        E[k] = s.scal(F(Fr(1, k)))
    return E


def nonlinear_data(N, alpha0, alpha1, M0, G, f, ell, T=None):
    """(p, a, C) for wick1; M0 rational, alpha0, alpha1, ell field elements"""
    from lm_wick import beta
    W, K = 2 * N, 2 * N + 2
    T = T if T is not None else N + 8
    ae = lsm(1, 1).scal(alpha0) + lsm(1, 2).scal(alpha1)          # a(e)
    ia = ae.inv(rel=T)                                            # valuation -1
    ap = lsc(1).scal(alpha0) + lsm(2, 1).scal(alpha1)             # a'(e)
    Y = tz(K)
    Y[1], Y[2] = ap * ia, ia.scal(alpha1)                         # a(e+t)/a(e) - 1
    lg = tlog1p(Y, K)
    p = {}
    for k in range(2, K + 1):
        x = reg_coef(G, k, True) if len(G) > k else LS(0, [])
        term = ap * lg[k - 1]
        if k - 2 >= 1:
            term = term + lg[k - 2].scal(fr(2 * alpha1))
        x = x + term.scal(F(Fr(1, k)))
        if k == 2:
            x = x + lsc(1).scal(fr(alpha1 * ell))
        p[k] = x
    C = -(p[2].scal(F(2))).inv(rel=T)
    H = {}

    def add(key, val):
        H[key] = H[key] + val if key in H else val
    expo = F(Fr(M0) - Fr(1, 2))
    for j in range(1, W + 1):
        add((0, j), lg[j].scal(expo))
    ian = lsc(1)
    for n in range(1, N + 1):
        ian = ian * ia
        cn = F(beta(n, Fr(M0)))
        pw = texp([x.scal(F(-n)) for x in lg], W)                  # (a(e+t)/a(e))^(-n)
        for j in range(0, W - 2 * n + 1):
            add((n, j), (ian * pw[j]).scal(cn))
    for m, row in f.items():
        for j in range(0, W - 2 * m + 1):
            if (m, j) != (0, 0) and len(row) > j:
                add((m, j), reg_coef(row, j, True))
    a = gexp(H, W)
    return {k: v for k, v in p.items() if k >= 3}, a, C


def gprime(G):
    """g'(L) as an exact LS in e"""
    return LS(0, [fr(G[i] * i) for i in range(1, len(G))])


def lam_extra(G, alpha0, ell, kappa, f1, N):
    """log-amplitude terms of exp(kappa z lambda f1(sigma)), lambda = -g'(L) - alpha0 ell, f1(s* + x) = Sum f1[i] x^i"""
    lam = -gprime(G) - lsc(1).scal(fr(alpha0 * ell))
    lk = lam.scal(kappa)
    return {(1, j): lk * reg_coef(f1, j, True) for j in range(0, 2 * N - 1) if len(f1) > j}


def qdep_extra(G, alpha0, f1, N, T):
    """log-amplitude terms of exp(z Q' f1(sigma)); f1(s* + x) = Sum f1[i] x^i (polynomial)"""
    gp = gprime(G)
    d = gp - lsc(1).scal(gp.coef(0))                                # g'(L) - g'(s*), valuation >= 1
    Qp = lsm(1, 1).scal(alpha0) * lsexp(d.scal(finv(alpha0)), T)
    out = {}
    for j in range(0, 2 * N - 1):
        if len(f1) > j:
            out[(1, j)] = Qp * reg_coef(f1, j, True)
    return out
