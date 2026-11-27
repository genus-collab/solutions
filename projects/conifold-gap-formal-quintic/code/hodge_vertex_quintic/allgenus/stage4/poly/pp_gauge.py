# -*- coding: utf-8 -*-
"""pp_gauge.py -- Stage 4 / poly: the first row of the regular gauge, extracted from the scalar hatted WKB function alone,
and the normalisation c(z) = P_11(0, z), frame (alpha, beta) = (0, 0).

Hatted first row (UNIFORM_REDUCTION s.4, colour 0; C0 = varpi_0 = f0 in frame (0,0)):
    (1/C0) Phi(w, u) = Sum_k w^k g_k(u) F_(-k)(tau w),   g_k = u^mu(-k,k) tau^((1-c_(-k))/4) p_(k,-k)(u),   mu(-k,k) = 4q + (0,1,2,3,3)_r (k = 5q+r).
Triangular in k (F_b(0) = (1,-1,-1,1,1)); divisibility of g_k by u^mu is a CHECK (regularity of P).
Row 0 of z P' = A P - P At gives P_1. = (z P_0.' + P_0. At)/b_0, so with b_0(0) = bt_0(0) = 1:
    c(z) := P_11(0, z) = P_00(0, z) + z d/du P_01(0, z),   [z^5q] c = p_(5q,0)(0) + p_(5q-1,1)'(0).      (PROVED, two lines)
Canonical normalisation <=> c = 1;  (N) <=> c = 1 for the QRR-normalised Phi.
"""
from fractions import Fraction as Fr
from pp_series import s_mul, s_add, s_scale, s_inv, s_pow, periods, tau_series, model_rows

CEXP = (0, -1, 3, 2, 1)
F0SIGN = (1, -1, -1, 1, 1)


def mu_k(k): return 4 * (k // 5) + (0, 1, 2, 3, 3)[k % 5]     # (4k + c_(-k))/5


def first_row(ph, NK, K, C0=None, tau=None):
    """ph: list of polynomials ph_n(u) (n = 0..NK); returns dict k -> p_(k,-k)(u) (series mod u^(K - mu)), and checks"""
    if C0 is None: C0 = periods(K)[0]
    if tau is None: tau = tau_series(K)
    F = model_rows(NK)
    iC0 = s_inv(C0, K)
    tpow = [[Fr(1)] + [Fr(0)] * (K - 1)]
    for n in range(1, NK + 1): tpow.append(s_mul(tpow[-1], tau, K))
    tfac = {b: s_pow(tau, Fr(1 - CEXP[b], 4), K) for b in range(5)}
    g = {}; p = {}
    for N in range(NK + 1):
        lhs = s_mul(iC0, list(ph[N]) + [Fr(0)] * K, K)
        for k in range(N):
            b = (-k) % 5
            lhs = s_add(lhs, s_mul(g[k], s_scale(tpow[N - k], F[b][N - k], K), K), K, -1)
        b = (-N) % 5
        gN = s_scale(lhs, Fr(1, F0SIGN[b]), K); g[N] = gN
        m = mu_k(N)
        assert all(x == 0 for x in gN[:min(m, K)]), "g_%d not divisible by u^%d" % (N, m)
        pN = s_mul(gN[m:] + [Fr(0)] * m, s_pow(tau, -Fr(1 - CEXP[b], 4), K), K)
        p[N] = pN[:K - m]
    return p, g


def cnorm(p, Q):
    """[z^5q] c(z) = p_(5q,0)(0) + p_(5q-1,1)'(0), q = 0..Q"""
    out = []
    for q in range(Q + 1):
        a = p[5 * q][0]
        b = p[5 * q - 1][1] if q >= 1 else Fr(0)
        out.append(a + b)
    return out
