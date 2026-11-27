# -*- coding: utf-8 -*-
"""lm_models.py -- referee lmi_adversary: the formal quintic summand in LMI form, and random LMI data.

Summand: Gamma(1 + 5 sigma/z) Prod_(j=0..4) Gamma((zeta^j - sigma)/z) (-q)^(sigma/z), all Gamma Stirling-expanded (CP s.1.6).
* quintic_s1(K): colour 0, s* = 1: singular factor j = 0 (alpha0 = -1, M0 = 0); g, f = the other five factors, as Taylor
  series at sigma = 1 truncated at x^K (exact rationals).  Only one rational series is needed:
  Qs(x) := Sum_(j>=1) 1/(zeta^j - 1 - x) = -(4s^3 + 3s^2 + 2s + 1)/(1 + s + s^2 + s^3 + s^4), s = 1 + x
  (partial fractions of 5 s^4/(1 - s^5) minus the j = 0 term); [x^i] Qs = Sum_(j>=1) (zeta^j - 1)^(-i-1).
* quintic_col(K, i, zeta): colour i over F_p, s* = zeta^i, explicit zeta (no symmetric-function shortcut).
* quintic_s0(K, M): s* = 0 (q = oo): singular factor Gamma(M + 5 sigma/z) (M = 1 is Gamma(1 + 5s), M = 0 is Gamma(5s));
  g, f = the five factors Gamma((zeta^j - sigma)/z): Sum_j (zeta^j - x)^(-r) = 5 [5 | r] x-free part, from 5x^4/(1 - x^5).
* quintic_fixed(L0, N): all six factors at a fixed rational point sigma = L0 (no singular factor), for engine checks.
"""
import random
from fractions import Fraction as Fr
from math import comb
from lm_core import F, finv, fr, fpow, lsc, LS
from lm_wick import bk, binom, gexp


def pdiv(num, den, K):
    """power series num/den up to x^(K-1), field coefficients"""
    inv0 = finv(den[0])
    out = []
    for k in range(K):
        s = num[k] if k < len(num) else 0
        for i in range(1, min(k, len(den) - 1) + 1):
            s -= den[i] * out[k - i]
        out.append(fr(s * inv0))
    return out


def shift_poly(cs, a):
    """coefficients of P(a + x) from those of P(s) (field elements)"""
    n = len(cs)
    return [fr(sum(cs[d] * comb(d, i) * fpow(a, d - i) for d in range(i, n))) for i in range(n)]


def _assemble(K, h1, sum_r, s_lin):
    """G and f from: h1(i) = [x^i] of g''(s*+x), sum_r(r, i) = [x^i] Sum_reg (c_j - x)^(-r) (regular Gamma((c_j-sigma)/z)),
    s_lin = (base, sign) for the factor Gamma(1 + 5 sigma/z) if regular (or None): 5 sigma = 5 (s* + x)."""
    G = [F(0), F(0)] + [fr(h1(i - 2) * F(Fr(1, i * (i - 1)))) for i in range(2, K)]
    f = {}
    row = [F(0)]
    for i in range(1, K):
        c = fr(F(Fr(1, 2)) * sum_r(1, i - 1))           # -(1/2) log(c_j - sigma): derivative (1/2)/(c_j - sigma)
        if s_lin is not None:                           # +(1/2) log(5 sigma): derivative (1/2)/sigma
            c = fr(c + F(Fr(1, 2)) * fpow(s_lin, -i) * F((-1) ** (i - 1)))
        row.append(fr(c * F(Fr(1, i))))
    f[0] = row
    k = 1
    while 2 * k - 1 < K:
        n = 2 * k - 1
        row = []
        for i in range(K):
            c = fr(F(comb(n + i - 1, i)) * sum_r(n + i, 0))
            if s_lin is not None:                       # (5 sigma)^(1-2k) = 5^(1-2k) s*^(1-2k) (1 + x/s*)^(1-2k)
                c = fr(c + fpow(F(5), -n) * fpow(s_lin, -n - i) * F(binom(-n, i)))
            row.append(fr(c * F(bk(k))))
        f[n] = row
        k += 1
    return G, f


def quintic_s1(K):
    num = shift_poly([F(1), F(2), F(3), F(4)], F(1))
    den = shift_poly([F(1)] * 5, F(1))
    qs = [fr(-x) for x in pdiv(num, den, 3 * K + 8)]

    def sum_r(r, i):                                    # [x^i] Sum_(j>=1) (e_j - x)^(-r) = binom(r+i-1, i) q[r+i-1]
        return fr(F(comb(r + i - 1, i)) * qs[r + i - 1])

    def h1(i):                                          # g'' = 5/(1+x) + Qs
        return fr(F(5 * (-1) ** i) + qs[i])
    return _assemble(K, h1, sum_r, F(1))


def quintic_col(K, i0, zeta):
    st = fpow(zeta, i0)
    cs = [fr(fpow(zeta, j) - st) for j in range(5) if j != i0]

    def sum_r(r, i):
        return fr(F(comb(r + i - 1, i)) * sum(fpow(c, -(r + i)) for c in cs))

    def h1(i):                                          # 5/(s* + x) + Sum_j 1/(c_j - x)
        return fr(5 * fpow(st, -i - 1) * F((-1) ** i) + sum_r(1, i))
    return _assemble(K, h1, sum_r, st)


def quintic_s0(K):
    def psum(r):                                        # Sum_(j=0..4) zeta^(-j r)
        return F(5) if r % 5 == 0 else F(0)

    def sum_r(r, i):
        return fr(F(comb(r + i - 1, i)) * psum(r + i))
    return _assemble(K, lambda i: sum_r(1, i), sum_r, None)        # g'' = Sum_j 1/(zeta^j - x)


def h_quintic(N):
    """h_(w_0) over the current field: exp(Sum_k b_k z^(2k-1) [5^(1-2k) + Sum_(j>=1) (zeta^j - 1)^(1-2k)])"""
    from lm_data import zexp
    num = shift_poly([F(1), F(2), F(3), F(4)], F(1))
    den = shift_poly([F(1)] * 5, F(1))
    qs = [fr(-x) for x in pdiv(num, den, N + 2)]
    h = [F(0)] * (N + 1)
    for k in range(1, (N + 1) // 2 + 1):
        n = 2 * k - 1
        h[n] = fr(F(bk(k)) * (fpow(F(5), -n) + qs[n - 1]))
    return zexp(h, N)


def ginf(N):
    """Ginf(z) = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))"""
    from lm_wick import bernpoly
    from lm_data import zexp
    h = [F(0)] * (N + 1)
    for l in range(1, N // 5 + 1):
        h[5 * l] = F((-1) ** (l + 1) * bernpoly(5 * l + 1, Fr(1, 5)) / (l * (5 * l + 1)))
    return zexp(h, N)


def quintic_fixed(L0, N):
    """(p, a, C) at sigma = L0 for all six factors (constants as LS)"""
    L0 = F(L0)
    K = 2 * N + 4
    num = [fr(5 * fpow(L0, 4 - d) * comb(4, d)) for d in range(5)]          # 5 (L0 + t)^4
    den = [fr((1 if d == 0 else 0) - comb(5, d) * fpow(L0, 5 - d)) for d in range(6)]   # 1 - (L0 + t)^5
    gen = pdiv(num, den, K + 2 * N + 4)                                          # Sum_r P_r t^(r-1)
    P = lambda r: gen[r - 1]
    p = {k: LS(0, [fr((5 * fpow(L0, 1 - k) * (-1) ** k + P(k - 1)) * F(Fr(1, k * (k - 1))))])
         for k in range(2, K)}
    C = -(p[2].scal(F(2))).inv(rel=1)
    H = {}
    for kk in range(1, 2 * N + 1):
        H[(0, kk)] = LS(0, [fr(((-1) ** (kk + 1) * fpow(L0, -kk) + P(kk)) * F(Fr(1, 2 * kk)))])
    k = 1
    while 2 * k - 1 <= N:
        n = 2 * k - 1
        for r in range(0, 2 * N - 2 * n + 1):
            v = fpow(F(5), -n) * F(binom(-n, r)) * fpow(L0, -n - r) + F(comb(n + r - 1, r)) * P(n + r)
            H[(n, r)] = LS(0, [fr(v * F(bk(k)))])
        k += 1
    a = gexp(H, 2 * N)
    return {k: v for k, v in p.items() if k >= 3}, a, C


def rand_data(rng, dg, df, mz, h=7):
    """random polynomial g (degree dg) and f (z-degree mz, x-degree df), rationals of height <= h"""
    r = lambda: Fr(rng.randint(-h, h), rng.randint(1, h))
    G = [F(0), F(r())] + [F(r()) for _ in range(dg - 1)]
    f = {m: [F(r()) for _ in range(df + 1)] for m in range(mz + 1)}
    f[0][0] = F(0)
    return G, f
