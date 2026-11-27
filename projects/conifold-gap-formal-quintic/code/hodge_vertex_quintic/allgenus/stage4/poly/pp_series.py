# -*- coding: utf-8 -*-
"""pp_series.py -- Stage 4 / poly: exact power series in u (lists of Fractions, fixed length K) and the conifold-side data
over Q, re-derived (no project module imported):
  * periods(K): f0 = 1 + O(u^3), f2 = u^2 + O(u^3), analytic solutions of
    u D^4 - (1-u)(2 D^3 + 7/5 D^2 + 2/5 D + 24/625),  D = q d/dq = -(1-u) d/du  (the vanishing period g = f1 + (7/10) f2 is built in pp_res);
  * tau(K) = 1/2F1(1, 4/5; 9/5; u) = 1/Sum 4u^n/(5n+4);
  * model(N): the colour-0 model rows F_b(t), b = 0..4 (F_0(0) = 1; F_1 = T_-1 F_0, F_2 = -T_-2 F_1, F_3 = T_2 F_2,
    F_4 = -T_1 F_3, F_0 = -T_0 F_4 checked; T_k = -1 + (t/5)(k - 4 t d/dt)).
"""
from fractions import Fraction as Fr


def s_mul(a, b, K):
    out = [Fr(0)] * K
    for i, x in enumerate(a[:K]):
        if x:
            for j in range(min(len(b), K - i)): out[i + j] += x * b[j]
    return out
def s_add(a, b, K, c=1):
    out = [Fr(0)] * K
    for i in range(K):
        out[i] = (a[i] if i < len(a) else 0) + c * (b[i] if i < len(b) else 0)
    return out
def s_scale(a, c, K): return [(a[i] * c if i < len(a) else Fr(0)) for i in range(K)]
def s_inv(a, K):
    out = [Fr(0)] * K; out[0] = 1 / Fr(a[0])
    for n in range(1, K):
        out[n] = -sum((a[k] * out[n - k] for k in range(1, min(n, len(a) - 1) + 1)), Fr(0)) * out[0]
    return out
def s_pow(a, r, K):
    """a^r, a[0] = 1, r rational"""
    assert a[0] == 1
    out = [Fr(0)] * K; out[0] = Fr(1); r = Fr(r)
    for n in range(1, K):
        out[n] = sum((((r + 1) * k - n) * a[k] * out[n - k] for k in range(1, min(n, len(a) - 1) + 1)), Fr(0)) / n
    return out
def s_deriv(a, K): return [((i + 1) * a[i + 1] if i + 1 < len(a) else Fr(0)) for i in range(K)]
def s_shift(a, m, K): return [(a[i - m] if 0 <= i - m < len(a) else Fr(0)) for i in range(K)]


def _D(f, K):
    """D f = -(1-u) f'"""
    d = s_deriv(f, K)
    return [-(d[i] - (d[i - 1] if i >= 1 else 0)) for i in range(K)]


def _pf(f, K):
    d1 = _D(f, K); d2 = _D(d1, K); d3 = _D(d2, K); d4 = _D(d3, K)
    inner = [2 * d3[i] + Fr(7, 5) * d2[i] + Fr(2, 5) * d1[i] + Fr(24, 625) * f[i] for i in range(K)]
    out = [Fr(0)] * K
    for i in range(K):
        out[i] = (d4[i - 1] if i >= 1 else 0) - inner[i] + (inner[i - 1] if i >= 1 else 0)
    return out


def periods(K):
    """f0, f2 (analytic, f_k = u^k + O(u^3)), mod u^K"""
    def solve(init):
        a = [Fr(x) for x in init] + [Fr(0)] * (K + 3 - len(init))
        for n in range(3, K + 3):
            a[n] = Fr(0); c = _pf(a, K + 3)[n - 3]
            a[n] = -c / (n * (n - 1) ** 2 * (n - 2))
        assert all(x == 0 for x in _pf(a, K + 3)[:K]), "PF fails"
        return a[:K]
    f0 = solve([1, 0, 0]); f2 = solve([0, 0, 1])
    return f0, f2


def tau_series(K):
    return s_inv([Fr(4, 5 * n + 4) for n in range(K)], K)


def model_rows(N):
    """colour-0 model rows F_b(t) as lists of length N+1 (F_0 solved from the cycle at order t^(n+1), pivot computed)"""
    M = N + 2
    def T(k, F):
        out = [Fr(0)] * M
        for n in range(M):
            out[n] -= F[n]
            if n >= 1: out[n] += Fr(k - 4 * (n - 1), 5) * F[n - 1]
        return out
    def cyc(F0):
        F1 = T(-1, F0); F2 = [-x for x in T(-2, F1)]; F3 = T(2, F2); F4 = [-x for x in T(1, F3)]
        return [F0, F1, F2, F3, F4], [-x for x in T(0, F4)]
    F0 = [Fr(0)] * M; F0[0] = Fr(1)
    for n in range(1, N + 1):
        F0[n] = Fr(0); r = cyc(F0)[1][n + 1] - F0[n + 1]
        F0[n] = Fr(1); piv = cyc(F0)[1][n + 1] - F0[n + 1] - r
        F0[n] = -r / piv
    rows, back = cyc(F0)
    assert back[:N + 1] == F0[:N + 1], "model cycle fails"
    return [R[:N + 1] for R in rows]
