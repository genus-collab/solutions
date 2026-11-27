# -*- coding: utf-8 -*-
"""lm_data.py -- referee lmi_adversary: LMI integrands at L = s* + e (one variable), the LMI prediction, analysis.

Integrand (RA REFEREE.md s.3.1): phase S = lambda sigma + g(sigma) + Sum_sing [a log a - a], a = alpha (sigma - s*);
amplitude A = Prod_sing a^expo exp(Sum_n stir_n (z/a)^n) * Prod (sigma - s*)^p * exp(f(sigma, z)).
The exact Stirling form of Gamma(M + a/z) is expo = M - 1/2, stir_n = beta_n(M) (DLMF 5.11.8 form; for M in Z it equals
RA's Poch_M(z/a) exp(Sum b_k (z/a)^(2k-1)), checked in lm_run E3).  lambda never enters (Taylor data of order >= 2).
Data: g(s* + x) = Sum_i G[i] x^i; f(s* + x, z) = Sum_m Sum_i f[m][i] z^m x^i (field elements); `exact` = these are
polynomials (else truncated Taylor series: then every derived coefficient carries its true precision).
"""
from fractions import Fraction as Fr
from math import comb
from lm_core import LS, lsc, lsm, F, finv, fr, fz, fpow, INF
from lm_wick import gexp, wick1, beta


class Sing:
    """Gamma(M + alpha (sigma - s*)/z), Stirling-expanded; alpha a field element, M rational.
    expo: exponent of a (default M - 1/2); stir: n -> coefficient of (z/a)^n in the log (default beta_n(M));
    phase: whether the phase a log a - a is present."""

    def __init__(self, alpha, M, expo=None, stir=None, phase=True, field_M=False):
        """M rational, or (field_M=True) a field element; expo and stir(n) are returned as field elements."""
        self.alpha = alpha
        self.M = M if field_M else F(M)
        self.expo = fr(self.M - F(Fr(1, 2))) if expo is None else F(expo)
        if stir is not None:
            self.stir = lambda n, s=stir: F(s(n))
        elif field_M:
            self.stir = lambda n, M=M: beta_field(n, M)
        else:
            self.stir = lambda n, M=Fr(M): F(beta(n, M))
        self.phase = phase


def beta_field(n, h):
    """beta_n(h) for a field element h: (-1)^(n+1) B_(n+1)(h)/(n(n+1)), B_k(h) = Sum binom(k, i) B_i h^(k-i)"""
    from lm_wick import bern
    k = n + 1
    s = 0
    for i in range(k + 1):
        s += F(comb(k, i) * bern(i)) * fpow(h, k - i)
    return fr(F(Fr((-1) ** (n + 1), n * (n + 1))) * s)


def reg_coef(G, k, exact, shift=0):
    """[t^k] of Sum_i G[i] (e + t)^i  (= Taylor coefficient at L = s* + e), as LS; truncated data keep precision."""
    c = [fr(G[i] * comb(i, k)) for i in range(k, len(G))]
    return LS(0, c, INF if exact else max(len(G) - k, 0))


def one_var(N, sings, G, f, exact=True, T=None, pows=(), extraH=None, extraP=None):
    """Taylor data of the one-variable integrand at L = s* + e: (p, a, C) for lm_core.wick1.
    extraH: dict (m, j) -> LS added to log A;  extraP: dict k -> LS added to the phase coefficients."""
    W, K = 2 * N, 2 * N + 2
    T = T if T is not None else N + 6
    p = {}
    for k in range(2, K + 1):
        x = reg_coef(G, k, exact) if len(G) > k else LS(0, [], INF if exact else 0)
        for s in sings:
            if s.phase:
                x = x + lsm(1, 1 - k).scal(fr(s.alpha * F(Fr((-1) ** k, k * (k - 1)))))
        if extraP and k in extraP:
            x = x + extraP[k]
        p[k] = x
    C = -(p[2].scal(F(2))).inv(rel=T)
    H = {}

    def add(key, val):
        H[key] = H[key] + val if key in H else val
    for s in sings:
        ia = finv(s.alpha)
        for j in range(1, W + 1):
            add((0, j), lsm(1, -j).scal(fr(F(Fr((-1) ** (j + 1), j)) * s.expo)))
        for n in range(1, N + 1):
            cn = s.stir(n)
            if fz(cn):
                continue
            base = fr(cn * fpow(ia, n))
            for j in range(0, W - 2 * n + 1):
                add((n, j), lsm(1, -n - j).scal(fr(base * F(binom_neg(n, j)))))
    for q in pows:
        for j in range(1, W + 1):
            add((0, j), lsm(1, -j).scal(F(Fr((-1) ** (j + 1), j) * Fr(q))))
    for m, row in f.items():
        for j in range(0, W - 2 * m + 1):
            if m == 0 and j == 0:
                continue
            if len(row) > j:
                add((m, j), reg_coef(row, j, exact))
            elif not exact:
                add((m, j), LS(0, [], 0))
    if extraH:
        for key, val in extraH.items():
            if 2 * key[0] + key[1] <= W:
                add(key, val)
    a = gexp(H, W)
    return {k: v for k, v in p.items() if k >= 3}, a, C


def binom_neg(n, j):
    """binom(-n, j) = (-1)^j binom(n + j - 1, j)"""
    return (-1) ** j * comb(n + j - 1, j)


def cseries(N, sings, G, f, **kw):
    p, a, C = one_var(N, sings, G, f, **kw)
    return wick1(p, a, N, C)


# ---- power series in z with field coefficients ----------------------------------------------------------------------
def zexp(h, N):
    """exp of Sum_(n>=1) h[n] z^n (h[0] ignored) up to z^N"""
    E = [F(0)] * (N + 1)
    E[0] = F(1)
    for k in range(1, N + 1):
        s = 0
        for i in range(1, k + 1):
            if i < len(h) and h[i]:
                s += i * h[i] * E[k - i]
        E[k] = fr(s * finv(F(k)))
    return E


def lmi_value(N, alpha0, M0, G, f):
    """RA's LMI formula: exp(gt(delta)/z) F(s* + delta, z)/F(s*, 0), delta = -M0 z/alpha0 (as a series in z)."""
    d1 = fr(-M0 * finv(alpha0))   # M0: field element
    h = [F(0)] * (N + 1)
    for k in range(2, len(G)):
        if k - 1 <= N:
            h[k - 1] = fr(h[k - 1] + G[k] * fpow(d1, k))
    for m, row in f.items():
        for i, c in enumerate(row):
            if (m, i) != (0, 0) and m + i <= N:
                h[m + i] = fr(h[m + i] + c * fpow(d1, i))
    return zexp(h, N)


def analyse(cs):
    """for each n: (pole order, [e^0] value or None, precision); pole order counts only known coefficients."""
    out = []
    for x in cs:
        pole = -x.v if (x.c and x.v < 0) else 0
        val = x.coef(0) if x.prec >= 1 else None
        out.append((pole, val, x.prec))
    return out


def leading(x):
    """(exponent, coefficient) of the lowest known non-zero term"""
    return (x.v, x.c[0]) if x.c else (None, None)
