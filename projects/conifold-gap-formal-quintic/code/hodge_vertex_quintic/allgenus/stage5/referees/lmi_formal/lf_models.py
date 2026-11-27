"""lf_models: Taylor data of the integrands tested by the lmi_formal referee (own derivation).

Quintic (colour 0): S_0 = lambda s + [5s log 5s - 5s] + Sum_(j=0..4) [(zeta^j - s) log(zeta^j - s) - (zeta^j - s)],
amplitude exp(S_1 + Sum_k z^(2k-1) S_2k) (CHART_PROOF s.1.6). For a = zeta^j - s: d^k/ds^k (a log a - a) = (k-2)!/a^(k-1)
(k >= 2); for 5s log 5s - 5s: 5 (-1)^k (k-2)!/s^(k-1).
* q_onevar(L)       the one-variable full-Stirling series c(z; L) at rational L (L != 0, L^5 != 1);
* q_chart1(L)       the two-variable chart of REFEREE s.3.3 at s* = 1 (alpha0 = -1, M0 = 0), any rational L != 0;
* q_chart0(L)       the same chart at s* = 0 (singular factor Gamma(1 + 5s/z): alpha0 = 5, M0 = 1), L^5 != 1.
General LMI data: g = Sum_(k>=2) G[k] (sigma - s*)^k, log F = Sum LF[(d, m)] (sigma - s*)^d z^m ((d, m) != (0, 0)).
* lmi_onevar(eps0)  the one-variable series c^F at L = s* + eps0 (with failing-control switches);
* lmi_chart(eps0)   the two-variable chart (Psi_Q, e^(M0 mu) F) at its critical point; eps0 = 0 is Q = 0;
* lmi_scaled()      sigma = L + eps tau, z = eps w (L1a, F3): data over TS = Q[eps]/(eps^K);
* lmi_predict()     exp(gt(delta)/z) F(s* + delta, z)/F(s*, 0), delta = -M0 z/alpha0.
"""
from fractions import Fraction as Fr
from math import comb, factorial
from lf_core import bk, gbinom, power_sums_all, power_sums_reg, ser_exp, TS


def _addto(dct, key, val):
    dct[key] = dct.get(key, 0) + val


# ---------------------------------------------------------------- quintic
def q_onevar(L, nmax, sing_stirling=True, twist_y=1):
    """Data (ph3, logA, C) of c(z; L). sing_stirling=False drops the Bernoulli part of the j = 0 factor."""
    L = Fr(L)
    R = 2 * nmax + 4
    P = power_sums_all(L, R + 2)
    ph3 = {k: (5 * (-1) ** k * L ** (1 - k) + P[k - 1]) / (k * (k - 1)) for k in range(3, 2 * nmax + 3)}
    S2 = 5 / L + P[1]
    logA = {}
    for d in range(1, 2 * nmax + 1):
        _addto(logA, (d, 0), Fr((-1) ** (d + 1), 2 * d) / L ** d + P[d] / (2 * d))
    k = 1
    while 2 * k - 1 <= nmax:
        e = 2 * k - 1
        for d in range(0, 2 * nmax - 2 * e + 1):
            Pj = P[e + d] - (0 if sing_stirling else (1 - L) ** (-(e + d)))
            val = bk(k) * (twist_y * Fr(5) ** (1 - 2 * k) * gbinom(1 - 2 * k, d) * L ** (1 - 2 * k - d)
                           + comb(e + d - 1, d) * Pj)
            _addto(logA, (d, e), val)
        k += 1
    return ph3, logA, -1 / S2


def q_chart1(L, nmax, mu_scale=1):
    """Chart at s* = 1: g = 5s log 5s - 5s + Sum_(j>=1) [..], F = (5s)^(1/2) prod_(j>=1) (..)^(-1/2) exp(Stirling).
    mu-part -Q e^(mu_scale mu); mu_scale = 2 is a failing control. Returns (phs3, logAs, phm3, logAm, H)."""
    L = Fr(L)
    R = 2 * nmax + 4
    p = power_sums_reg(L, R + 2)
    phs = {k: (5 * (-1) ** k * L ** (1 - k) + p[k - 1]) / (k * (k - 1)) for k in range(2, 2 * nmax + 3)}
    a0 = 1 - L                       # a0(L) = alpha0 (L - s*), alpha0 = -1
    c = mu_scale
    phm = {k: -a0 / c * Fr(c) ** k / factorial(k) for k in range(2, 2 * nmax + 3)}
    logAs = {}
    for d in range(1, 2 * nmax + 1):
        _addto(logAs, (d, 0), Fr((-1) ** (d + 1), 2 * d) / L ** d + p[d] / (2 * d))
    k = 1
    while 2 * k - 1 <= nmax:
        e = 2 * k - 1
        for d in range(0, 2 * nmax - 2 * e + 1):
            val = bk(k) * (Fr(5) ** (1 - 2 * k) * gbinom(1 - 2 * k, d) * L ** (1 - 2 * k - d)
                           + comb(e + d - 1, d) * p[e + d])
            _addto(logAs, (d, e), val)
        k += 1
    H = (2 * phs.pop(2), Fr(-1), 2 * phm.pop(2))
    return phs, logAs, phm, {}, H


def q_chart0(L, nmax):
    """Chart at s* = 0: g = Sum_(j=0..4) [(zeta^j - s) log(..) - (..)], F = prod_j (..)^(-1/2) exp(Stirling_j),
    a0 = 5 sigma, amplitude e^(mu) (M0 = 1)."""
    L = Fr(L)
    R = 2 * nmax + 4
    P = power_sums_all(L, R + 2)
    phs = {k: P[k - 1] / (k * (k - 1)) for k in range(2, 2 * nmax + 3)}
    a0 = 5 * L
    phm = {k: -a0 / factorial(k) for k in range(2, 2 * nmax + 3)}
    logAs = {}
    for d in range(1, 2 * nmax + 1):
        _addto(logAs, (d, 0), P[d] / (2 * d))
    k = 1
    while 2 * k - 1 <= nmax:
        e = 2 * k - 1
        for d in range(0, 2 * nmax - 2 * e + 1):
            _addto(logAs, (d, e), bk(k) * comb(e + d - 1, d) * P[e + d])
        k += 1
    H = (2 * phs.pop(2), Fr(5), 2 * phm.pop(2))
    return phs, logAs, phm, {(1, 0): Fr(1)}, H


# ---------------------------------------------------------------- general LMI data
def _g_taylor(G, eps0, j):
    return sum(G[k] * comb(k, j) * eps0 ** (k - j) for k in G if k >= j)


def _logF_taylor(LF, eps0, nmax):
    out = {}
    for (d1, m), f in LF.items():
        for d in range(0, d1 + 1):
            if (d, m) != (0, 0) and d + 2 * m <= 2 * nmax:
                _addto(out, (d, m), f * comb(d1, d) * eps0 ** (d1 - d))
    return out


def _sing_amp(alpha0, M0, eps0, nmax, stirling=True, poch=True, b1_shift=0, expo=None):
    """log of the normalised singular amplitude a0^(expo) Poch_M0(z/a0) exp(Sum b_k (z/a0)^(2k-1)) at L = s* + eps0,
    in s = sigma - L; a0 = alpha0 (eps0 + s); expo defaults to M0 - 1/2."""
    expo = Fr(M0) - Fr(1, 2) if expo is None else Fr(expo)
    out = {}
    for d in range(1, 2 * nmax + 1):
        _addto(out, (d, 0), expo * Fr((-1) ** (d + 1), d) / eps0 ** d)

    def inv_a0_pow(p, d):            # [s^d] a0^-p
        return Fr(alpha0) ** (-p) * gbinom(-p, d) * eps0 ** (-p - d)
    if poch:
        for p in range(1, nmax + 1):
            for d in range(0, 2 * nmax - 2 * p + 1):
                if M0 >= 2:
                    c = sum(Fr((-1) ** (p + 1) * i ** p, p) for i in range(1, M0))
                elif M0 <= -1:
                    c = sum(Fr(i ** p, p) for i in range(1, -M0 + 1))
                else:
                    c = 0
                if c:
                    _addto(out, (d, p), c * inv_a0_pow(p, d))
    if stirling:
        k = 1
        while 2 * k - 1 <= nmax:
            e = 2 * k - 1
            for d in range(0, 2 * nmax - 2 * e + 1):
                _addto(out, (d, e), bk(k, b1_shift if k == 1 else 0) * inv_a0_pow(e, d))
            k += 1
    return out


def lmi_onevar(G, LF, alpha0, M0, eps0, nmax, **ctrl):
    eps0, alpha0 = Fr(eps0), Fr(alpha0)
    ph3 = {}
    for j in range(3, 2 * nmax + 3):
        ph3[j] = _g_taylor(G, eps0, j) + Fr((-1) ** j) * alpha0 / (j * (j - 1) * eps0 ** (j - 1))
    S2 = 2 * _g_taylor(G, eps0, 2) + alpha0 / eps0
    logA = _sing_amp(alpha0, M0, eps0, nmax, **ctrl)
    for key, v in _logF_taylor(LF, eps0, nmax).items():
        _addto(logA, key, v)
    return ph3, logA, -1 / S2


def lmi_chart(G, LF, alpha0, M0, eps0, nmax, mu_scale=1):
    """Chart (Psi_Q, e^(M0 mu) F) at its critical point (L, mu_c(L)), L = s* + eps0; Q e^(mu_c) = a0(L) = alpha0 eps0."""
    eps0, alpha0 = Fr(eps0), Fr(alpha0)
    phs = {j: _g_taylor(G, eps0, j) for j in range(2, 2 * nmax + 3)}
    a0 = alpha0 * eps0
    c = mu_scale
    phm = {k: -a0 / c * Fr(c) ** k / factorial(k) for k in range(2, 2 * nmax + 3)}
    logAs = _logF_taylor(LF, eps0, nmax)
    logAm = {(1, 0): Fr(M0)} if M0 else {}
    H = (2 * phs.pop(2), alpha0, 2 * phm.pop(2))
    return phs, logAs, phm, logAm, H


def lmi_scaled(G, LF, alpha0, M0, nmax, **ctrl):
    """Method A over TS = Q[eps]/(eps^(nmax+1)): returns (ph3, logA, C) of c'(w; eps), c^F_n(s* + eps) = c'_n/eps^n."""
    TS.K = nmax + 1
    alpha0 = Fr(alpha0)

    def g_part(j):                   # Sum_k G[k] C(k, j) eps^(k-1)
        c = [Fr(0)] * TS.K
        for k, v in G.items():
            if k >= j and k - 1 < TS.K:
                c[k - 1] += v * comb(k, j)
        return TS(c)
    ph3 = {j: g_part(j) + Fr((-1) ** j) * alpha0 / (j * (j - 1)) for j in range(3, 2 * nmax + 3)}
    S2 = 2 * g_part(2) + alpha0      # = alpha0 + eps g''(L)
    sing = _sing_amp(alpha0, M0, Fr(1), nmax, **ctrl)     # eps0 = 1 in tau: a0 = alpha0 (1 + tau)
    logA = {key: TS.const(v) for key, v in sing.items()}
    for (d1, m), f in LF.items():    # log F(L + eps tau, eps w): eps^(d1+m) (1 + tau)^d1 w^m
        for d in range(0, d1 + 1):
            if (d, m) != (0, 0) and d + 2 * m <= 2 * nmax and d1 + m < TS.K:
                c = [Fr(0)] * TS.K
                c[d1 + m] = f * comb(d1, d)
                logA[(d, m)] = logA.get((d, m), TS([])) + TS(c)
    return ph3, logA, (-1) * S2.inv()


def lmi_predict(G, LF, alpha0, M0, nmax, sign=1):
    """Series of exp(gt(delta)/z + log F(s* + delta, z) - log F(s*, 0)); sign = -1 flips delta (control)."""
    c = sign * Fr(-M0) / Fr(alpha0)
    a = [Fr(0)] * (nmax + 1)
    for k, v in G.items():
        if 1 <= k - 1 <= nmax:
            a[k - 1] += v * c ** k
    for (d, m), f in LF.items():
        if 1 <= d + m <= nmax:
            a[d + m] += f * c ** d
    return ser_exp(a, nmax + 1)
