# -*- coding: utf-8 -*-
"""gm_fs.py -- formal saddle-point series of the FULL-STIRLING Mellin-Barnes integrand of the formal quintic.

MB integrand (s = sigma/z):  Gamma(5s + 1) prod_{j=0..4} Gamma(zeta^j/z - s) (-q)^s.  Replace EVERY Gamma by its Stirling form
Gamma_St(x) = exp((x - 1/2) log x - x + (1/2) log 2 pi + Sum_k b_k x^(1-2k)),  b_k = B_2k/(2k(2k-1))   (Gamma(y+1) = y Gamma(y)).
The log z terms cancel and the integrand is exp(S_0/z + S_1 + Sum_k z^(2k-1) S_2k) dsigma (times an elementary constant), with
   S_0' = 5 log(5 sigma) - log(1 - sigma^5) + log(-q),   S_0'' = 5/(sigma (1 - sigma^5)),
   S_1  = (1/2) log(5 sigma/(1 - sigma^5)),   S_2k = b_k [ (5 sigma)^(1-2k) + Sum_j (zeta^j - sigma)^(1-2k) ].
Saddle of colour 0: sigma = mu = L, u mu^5 = 1 (u = 1 - 5^5 q).  e^(S_1(mu)) (-S_0''(mu))^(-1/2) = const * mu, and S_0(mu(q)) = u_0 +
const (D S_0 = mu), so the formal saddle-point expansion is const * e^(u_0/z) L * c(z; L) with
   c(z; L) = exp(Sum_k z^(2k-1) S_2k(mu)) Sum_m (2m-1)!! (z/a)^m [Y^(2m)] ( e^(S_1(mu+t)-S_1(mu)) t'(Y) e^(Sum_k z^(2k-1)(S_2k(mu+t)-S_2k(mu))) ),
a = -S_0''(mu), S_0(mu+t) - S_0(mu) = -(a/2) Y^2 (t = t(Y) by Lagrange inversion).  All coefficients are in Q(mu): exact at
rational mu (over Q or mod p).  Cyclotomic sums: Sum_j (zeta^j - sigma)^(-m) = f^(m-1)(sigma)/(m-1)!, f = 5 sigma^4/(1 - sigma^5).
API: fs_series(R, mu, N) -> [c_0..c_N] (c_0 = 1).
"""
from fractions import Fraction as Fr
from math import comb
from gm_series import (smul, sinv, slog1, sexp0, spow, sint, reversion_powers, compose, bernoulli_list)


def _tay_log1(R, c, n):
    """log(1 + c t) mod t^n"""
    return slog1(R, [R.one(), R.of(c)] + [R.zero()] * (n - 2))


def _tay_g(R, mu, n):
    """g(t) = (1 - (mu+t)^5)/(1 - mu^5), constant term 1"""
    d = 1 - Fr(mu) ** 5; out = [R.zero()] * n; out[0] = R.one()
    for k in range(1, 6):
        if k < n: out[k] = R.of(-comb(5, k) * Fr(mu) ** (5 - k) / d)
    return out


def fs_series(R, mu, N, return_parts=False):
    mu = Fr(mu); nT = 2 * N + 4                                   # t-degree truncation
    l1 = _tay_log1(R, 1 / mu, nT); lg = slog1(R, _tay_g(R, mu, nT))
    dphi = [R.norm(5 * x - y) for x, y in zip(l1, lg)]            # S_0'(mu + t)
    phi = sint(R, dphi)                                            # S_0(mu+t) - S_0(mu), starts at t^2
    a = R.norm(-2 * phi[2])
    assert a == R.of(-5 / (mu * (1 - mu ** 5))), "Hessian check"
    ratio = [R.norm(x * R.inv(phi[2])) for x in phi[2:]] + [R.zero(), R.zero()]
    w = spow(R, ratio[:nT], Fr(1, 2))                             # Y = t w(t)
    nY = 2 * N + 1
    T = reversion_powers(R, w, nY + 1, nY)
    tY = T[1]; dtY = [R.norm((k + 1) * tY[k + 1]) for k in range(nY)]
    B0 = [R.norm((x - y) * R.inv(2)) for x, y in zip(l1, lg)]     # S_1(mu+t) - S_1(mu)
    E = smul(R, sexp0(R, compose(R, B0[:nY + 1], T, nY)), dtY, nY)
    # f(mu + t) = 5 (mu+t)^4/(1 - (mu+t)^5):  F[n] = [t^n]
    nF = N + nT + 2
    num = [R.of(5 * comb(4, k) * mu ** (4 - k)) if k <= 4 else R.zero() for k in range(nF)]
    den = [R.of(1 - mu ** 5)] + [R.of(-comb(5, k) * mu ** (5 - k)) if k <= 5 else R.zero() for k in range(1, nF)]
    F = smul(R, num, sinv(R, den, nF), nF)
    Bn = bernoulli_list(N + 2)
    S2k_const = {}; S2kY = {}
    for k in range(1, N // 2 + 2):
        j = 2 * k - 1
        if j > N: break
        bk = Bn[2 * k] / (2 * k * (2 * k - 1))
        ser = []
        for r in range(nT):
            # (1 + t/mu)^(1-2k) coefficient: binom(1-2k, r) mu^-r ; binom(-m, r) = (-1)^r binom(m+r-1, r), m = 2k-1
            c1 = Fr((-1) ** r * comb(2 * k - 2 + r, r)) * (5 * mu) ** (1 - 2 * k) / mu ** r
            ser.append(R.norm(R.of(bk * c1) + R.of(bk) * comb(2 * k - 2 + r, r) * F[2 * k - 2 + r]))
        S2k_const[j] = ser[0]
        ser0 = [R.zero()] + ser[1:]
        S2kY[j] = compose(R, ser0[:nY + 1], T, nY)
    # E2 = exp(Sum_j z^j S2kY[j]) as a bivariate series: E2[n] = Y-series
    E2 = [[R.one()] + [R.zero()] * (nY - 1)]
    for n in range(1, N + 1):
        acc = [R.zero()] * nY
        for j in range(1, n + 1):
            if j in S2kY:
                prod = smul(R, S2kY[j], E2[n - j], nY)
                acc = [R.norm(x + j * y) for x, y in zip(acc, prod)]
        E2.append([R.norm(x * R.inv(n)) for x in acc])
    ainv = R.inv(a); core = [R.zero()] * (N + 1)
    EE = [smul(R, E, E2[n], nY) for n in range(N + 1)]
    for n in range(N + 1):
        s = R.zero(); df = 1; ap = R.one()
        for m in range(0, n + 1):
            if m > 0: df *= (2 * m - 1); ap = R.norm(ap * ainv)
            s += R.of(df) * ap * EE[n - m][2 * m]
        core[n] = R.norm(s)
    # multiply by exp(Sum_j z^j S2k(mu))
    lc = [R.zero()] * (N + 1)
    for j, v in S2k_const.items(): lc[j] = v
    pref = sexp0(R, lc)
    c = smul(R, pref, core, N + 1)
    return (c, dict(a=a, pref=pref, core=core)) if return_parts else c
