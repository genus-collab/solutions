# -*- coding: utf-8 -*-
"""ch_phases.py -- Stage 5 / chart: the three formal integrands of CHART_PROOF s.3, as Taylor data at their critical points,
over a field F of ch_core (QZ exact or FP), at a parameter L (a Fraction, L != 0 and L^5 != 1 except L = 1 for the chart).

  sigma-form  (c(z; L)):  phase S_0, amplitude e^(S_1) exp(Sum_k z^(2k-1) S_2k), one variable sigma = L + t;
  7-variable  Omega_7:    phase Phi, amplitude prod_(j=0..4) 1/x_j, variables (s, eta_y, xi_0..xi_4) at p(L);
  chart       Omega_ch:   phase psi, amplitude prod_(j=1..4) 1/x_j (measure dy), variables (eta_y, xi_1..xi_4) at p_ch(L).
Variants (controls): sigma-form 'g5s' (Gamma(5s) for Gamma(1+5s): S_1 -> S_1 - log 5 sigma); 7-variable 'dyy' (measure dy/y);
chart 'nojac' (the Jacobian X_0 of (L2) dropped: amplitude 1/(X_0 prod x_j), i.e. y^-5) and 'dyy' (measure dy/y).
"""
from fractions import Fraction as Fr
from math import comb
from ch_core import log1p_lin, geom, bk_list
from ch_wick import guni, gadd, gmul, gconst, gtrunc


def psums(F, L, K):
    """P_k = Sum_(j=0..4) (zeta^j - L)^(-k), k = 0..K"""
    P = [F.zero] * (K + 1)
    for j in range(5):
        ic = F.inv(F.sub(F.zeta(j), F.of(L))); p = F.one
        for k in range(K + 1):
            P[k] = F.add(P[k], p); p = F.mul(p, ic)
    return P


def sigma_data(F, L, N, variant=None):
    """(phi, amps) for the sigma-form at sigma = L: phi[k] = [t^k](S_0(L+t) - S_0(L)), amps[m] = [z^m] amplitude (t-lists)"""
    L = Fr(L); K = 2 * N + 2; P = psums(F, L, K + 1); iL = F.of(1 / L)
    phi = [F.zero] * (K + 1)
    for k in range(1, K):                                    # S_0'(L+t) = Sum_k t^k/k [5 (-1)^(k+1) L^-k + P_k]
        phi[k + 1] = F.smul(Fr(1, k * (k + 1)), F.add(F.of(5 * Fr((-1) ** (k + 1)) / L ** k), P[k]))
    sg = -1 if variant == 'g5s' else 1
    S1 = [F.zero] * (2 * N + 1)
    for k in range(1, 2 * N + 1):
        S1[k] = F.smul(Fr(1, 2 * k), F.add(F.of(sg * Fr((-1) ** (k + 1)) / L ** k), P[k]))
    bk = bk_list(N // 2 + 2)
    # bivariate exponent X[m][r] = [z^m t^r] (S_1(L+t) - S_1(L) + Sum_k z^(2k-1) S_2k(L+t))
    X = [[F.zero] * (2 * N + 1) for _ in range(N + 1)]
    X[0] = S1
    for k in range(1, N // 2 + 2):
        m = 2 * k - 1
        if m > N: break
        for r in range(0, 2 * (N - m) + 1):
            bin_ = Fr((-1) ** r * comb(2 * k - 2 + r, r))
            v = F.add(F.of(Fr(5) ** (1 - 2 * k) * L ** (1 - 2 * k - r)), F.smul((-1) ** r, P[2 * k - 1 + r]))
            X[m][r] = F.smul(bk[k] * bin_, v)
    amps = _bexp(F, X, N)
    return phi, amps


def _bexp(F, X, N):
    """exp of a bivariate series X[m][r] (X[0][0] = 0), truncated: z^m with t-degree <= 2(N-m)"""
    T = lambda m: 2 * (N - m)
    Y = [[F.zero] * (T(m) + 1) for m in range(N + 1)]; Y[0][0] = F.one
    # exp(X_0(t)) by the t-recursion, exp(Sum_(m>=1) z^m X_m(t)) by the z-recursion m Z_m = Sum_i i X_i Z_(m-i); then multiply
    E0 = [F.one] + [F.zero] * T(0)                          # exp(X[0](t)) in t
    for r in range(1, T(0) + 1):
        s = F.zero
        for i in range(1, r + 1): s = F.add(s, F.smul(i, F.mul(X[0][i], E0[r - i])))
        E0[r] = F.smul(Fr(1, r), s)
    # Z-part: exp(Sum_(m>=1) z^m X[m](t)) via z-derivative recursion: m Z_m = Sum_i i X_i Z_(m-i)
    Zs = [[F.one] + [F.zero] * T(0)]
    for m in range(1, N + 1):
        acc = [F.zero] * (T(m) + 1)
        for i in range(1, m + 1):
            Xi = X[i]; Zm = Zs[m - i]
            for a in range(T(m) + 1):
                if a < len(Xi) and not F.iszero(Xi[a]):
                    for b in range(T(m) + 1 - a):
                        if b < len(Zm) and not F.iszero(Zm[b]):
                            acc[a + b] = F.add(acc[a + b], F.smul(i, F.mul(Xi[a], Zm[b])))
        Zs.append([F.smul(Fr(1, m), x) for x in acc])
    for m in range(N + 1):
        for a in range(T(m) + 1):
            s = F.zero
            for b in range(a + 1): s = F.add(s, F.mul(E0[b], Zs[m][a - b]))
            Y[m][a] = s
    return Y


def phase7(F, L, D):
    """Taylor polynomial (degrees 1..D) of Phi at p(L), variables (s, eta_y, xi_0..xi_4); the s-linear term with the
    transcendental coefficient lambda + 5 log 5L - Sum_j log(zeta^j - L) is omitted (it is 0 by the choice of lambda)"""
    n = 7; L = Fr(L); Lq = F.of(L); one = F.one
    lg = log1p_lin(F, F.of(5 * L), D)
    P = guni(F, n, 1, [F.zero] + [F.smul(5 * L, x) for x in lg[1:]])                # 5L log(1 + eta/5L)
    P = gadd(F, P, guni(F, n, 1, [F.zero, F.neg(one)]))                              # - eta
    P = gadd(F, P, gmul(F, guni(F, n, 0, [F.zero, F.of(5)]), guni(F, n, 1, lg), D))  # 5 s log(1 + eta/5L)
    for j in range(5):
        cj = F.sub(F.zeta(j), Lq); lj = log1p_lin(F, cj, D)
        P = gadd(F, P, guni(F, n, 2 + j, [F.zero, F.neg(one)]))
        P = gadd(F, P, guni(F, n, 2 + j, [F.mul(cj, x) for x in lj]))
        P = gadd(F, P, gmul(F, guni(F, n, 0, [F.zero, F.neg(one)]), guni(F, n, 2 + j, lj), D))
    return P


def amp7(F, L, D, variant=None):
    n = 7; A = gconst(F, n, F.one)
    for j in range(5):
        A = gmul(F, A, guni(F, n, 2 + j, geom(F, F.sub(F.zeta(j), F.of(L)), D, -1)), D)
    if variant == 'dyy': A = gmul(F, A, guni(F, n, 1, geom(F, F.of(5 * Fr(L)), D, -1)), D)
    return A


def phasech(F, L, D, qsign=1):
    """Taylor polynomial (degrees 1..D) of psi at p_ch(L) = (5L, zeta^j - L), K = q y^5/prod x_j evaluated from q = (1 - L^-5)/5^5"""
    n = 5; L = Fr(L); Lq = F.of(L); one = F.one
    c = [F.sub(F.zeta(j), Lq) for j in range(5)]
    q = Fr(qsign) * (1 - 1 / L ** 5) / 3125
    K = F.mul(F.of(q * (5 * L) ** 5), F.inv(F.mul(F.mul(c[1], c[2]), F.mul(c[3], c[4]))))
    P = guni(F, n, 0, [F.zero] + [F.smul(5, x) for x in log1p_lin(F, F.of(5 * L), D)[1:]])
    P = gadd(F, P, guni(F, n, 0, [F.zero, F.neg(one)]))
    prod = guni(F, n, 0, geom(F, F.of(5 * L), D, 5))
    for j in range(1, 5):
        P = gadd(F, P, guni(F, n, j, [F.zero, F.neg(one)]))
        P = gadd(F, P, guni(F, n, j, [F.mul(F.sub(F.zeta(j), one), x) for x in log1p_lin(F, c[j], D)]))
        prod = gmul(F, prod, guni(F, n, j, geom(F, c[j], D, -1)), D)
    if not F.iszero(K):
        P = gadd(F, P, {d: m for d, m in prod.items() if d >= 1}, K)
    return P, K


def ampch(F, L, D, variant=None):
    n = 5; A = gconst(F, n, F.one)
    if variant == 'nojac':
        return guni(F, n, 0, geom(F, F.of(5 * Fr(L)), D, -5))
    for j in range(1, 5):
        A = gmul(F, A, guni(F, n, j, geom(F, F.sub(F.zeta(j), F.of(L)), D, -1)), D)
    if variant == 'dyy': A = gmul(F, A, guni(F, n, 0, geom(F, F.of(5 * Fr(L)), D, -1)), D)
    return A
