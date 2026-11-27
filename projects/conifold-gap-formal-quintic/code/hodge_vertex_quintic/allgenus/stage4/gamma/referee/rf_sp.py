# -*- coding: utf-8 -*-
"""rf_sp.py -- REFEREE: formal saddle-point (stationary-phase) series, written from scratch (no worker or project code).

Full-Stirling 1-form in x = 1/sigma (re-derived by the referee from the six Stirling series; REFEREE.md s.2):
   Omega = exp(W(x)/z) A(x; z) dx,  W'(x) = x^-2 log((1 - x^5)/(1 - u)),  u = s^5,  colour-0 saddle x = s = 1/L,
   A(x; z) = (1 - x^5)^(-1/2) exp( Sum_k b_k (z x)^(2k-1) [5^(1-2k) + p_(2k-1)(x)] ),  p_m(x) = Sum_{j=0..4} (zeta^j x - 1)^(-m).
c(z; L) := normalised formal Laplace series at x = s:  SP = e^(W(s)/z) (pi z/(-w2))^(1/2) A(s;0) [Sum_n c_n z^n],  c_0 = 1, where
W(s+t) - W(s) = w2 y^2 exactly (t = t(y) by Lagrange inversion) and <y^(2i)> = (2i-1)!! (z/(-2 w2))^i.
The power sums p_m are computed WITHOUT roots of unity: a_j = 1/(zeta^j x - 1) are the roots of
(1 - x^5) a^5 + 5a^4 + 10a^3 + 10a^2 + 5a + 1 = 0, hence p_m = -kappa (5p_(m-1) + 10p_(m-2) + 10p_(m-3) + 5p_(m-4) + p_(m-5)),
kappa = 1/(1 - x^5), p_0 = 5, p_(-r) = 5(-1)^r (r = 1..4).  Everything is a series in y with coefficients in Q(s): exact over Q, or mod p.
sp_generic: the same functional for an arbitrary phase/amplitude given in t = x - c (composition; small orders; cross-checks).
"""
from fractions import Fraction as Fr
from math import comb
from rf_ring import (mul, mulsum, pad, add, scal, inv, deriv, integ, log1, powr, reversion, compose, stirling_b)


def _gauss(F, w2, Fm, N):
    """c_n = Sum_{m+i=n} (2i-1)!! (-1/(2 w2))^i [y^(2i)] Fm[m]"""
    a = F.inv(F.red(-2 * w2)); c = []
    for n in range(N + 1):
        s = F.zero; df = 1; ap = F.one
        for i in range(0, n + 1):
            if i > 0: df *= 2 * i - 1; ap = F.red(ap * a)
            m = n - i
            if m < len(Fm) and 2 * i < len(Fm[m]): s += F.of(df) * ap * Fm[m][2 * i]
        c.append(F.red(s))
    return c


def _ty(F, Wt, n):
    """from W(c+t) - W(c) (t-series, starts at t^2): w2 and t(y) mod y^n"""
    w2 = Wt[2]; assert w2 != F.zero
    rho1 = [F.red(Wt[k + 2] * F.inv(w2)) for k in range(n)]
    return w2, reversion(F, powr(F, rho1, Fr(-1, 2), n), n)


def fs_series(F, s, N, drop_j0=False, bmod=None):
    """[c_0..c_N] of the full-Stirling saddle series at x = s.  Controls: drop_j0 removes the j = 0 Stirling term (the factor
    Gamma((1 - sigma)/z)); bmod = (k, delta) replaces b_k by b_k + delta in every Stirling term."""
    n = 2 * N + 2; nt = n + 3; s = F.of(s); u = F.red(s ** 5)
    X0 = [s, F.one] + [F.zero] * (nt - 2)
    X5 = mul(F, mul(F, mul(F, X0, X0, nt), mul(F, X0, X0, nt), nt), X0, nt)
    iu = F.inv(F.red(1 - u))
    ratio = [F.red((1 - X5[0]) * iu)] + [F.red(-x * iu) for x in X5[1:]]
    Wt = integ(F, mul(F, inv(F, mul(F, X0, X0, nt), nt), log1(F, ratio, nt), nt))[:nt]
    assert Wt[0] == F.zero and Wt[1] == F.zero and Wt[2] == F.red(-5 * s * s * iu * F.inv(F.of(2))), "Hessian"
    w2, ty = _ty(F, Wt, n)
    X = [F.red(s + ty[0])] + ty[1:]; dX = deriv(F, ty) + [F.zero]
    Xp = [[F.one] + [F.zero] * (n - 1), X]
    while len(Xp) <= max(N, 5): Xp.append(mul(F, Xp[-1], X, n))
    oneX5 = add(F, [F.one], scal(F, Xp[5], F.of(-1)))
    Rfac = powr(F, scal(F, oneX5, iu), Fr(-1, 2), n)
    kap = inv(F, oneX5, n)
    p = {0: [F.of(5)] + [F.zero] * (n - 1)}
    for r in range(1, 5): p[-r] = [F.of(5 * (-1) ** r)] + [F.zero] * (n - 1)
    for m in range(1, N + 1):
        comb5 = [F.red(5 * a + 10 * b + 10 * c + 5 * d + e) for a, b, c, d, e in
                 zip(p[m - 1], p[m - 2], p[m - 3], p[m - 4], p[m - 5])]
        p[m] = scal(F, mul(F, kap, comb5, n), F.of(-1))
    if drop_j0:                                                           # remove (x - 1)^(-m), i.e. the j = 0 root a_0 = 1/(x - 1)
        a0 = inv(F, add(F, X, [F.one], -1), n); a0p = [F.one] + [F.zero] * (n - 1)
        for m in range(1, N + 1):
            a0p = mul(F, a0p, a0, n); p[m] = add(F, p[m], a0p, -1)
    E = {}
    for j in range(1, N + 1, 2):
        k = (j + 1) // 2; bk = stirling_b(k)
        if bmod and bmod[0] == k: bk = bk + bmod[1]
        lj = 2 * (N - j) + 1
        th = mul(F, Xp[j], add(F, [F.of(Fr(1, 5 ** j))], p[j]), lj)
        E[j] = scal(F, th, F.of(j * bk))                                      # j b_k Theta_j
    G = [[F.one] + [F.zero] * (2 * N)]
    for m in range(1, N + 1):
        lm = 2 * (N - m) + 1
        acc = mulsum(F, [(E[j], G[m - j]) for j in range(1, m + 1, 2)], lm)
        G.append(scal(F, acc, F.inv(F.of(m))))
    RX = mul(F, Rfac, dX, n)
    Fm = [mul(F, RX, G[m], 2 * (N - m) + 1) for m in range(N + 1)]
    return _gauss(F, w2, Fm, N)


def sp_generic(F, Wt, amps, N):
    """normalised formal SP of exp(W/z) A dt at t = 0: Wt = W(c+t) - W(c) (t-series from t^2), amps = {m: t-series of [z^m] A},
    [z^0] A normalised to 1 at t = 0 by the caller.  Composition by Horner (small N only)."""
    n = 2 * N + 2; w2, ty = _ty(F, Wt, n); dty = deriv(F, ty) + [F.zero]
    Fm = []
    for m in range(N + 1):
        a = amps.get(m, [])
        Fm.append(mul(F, compose(F, pad(a, n, F), ty, n), dty, n) if any(x != F.zero for x in a) else [F.zero] * n)
    return _gauss(F, w2, Fm, N)


def bexp(F, terms, N, nt):
    """exp of Sum_j z^j terms[j] (t-series) as {m: t-series}, m <= N"""
    G = {0: [F.one] + [F.zero] * (nt - 1)}
    for m in range(1, N + 1):
        acc = [F.zero] * nt
        for j, tj in terms.items():
            if j <= m: acc = add(F, acc, scal(F, mul(F, tj, G[m - j], nt), F.of(j)))
        G[m] = scal(F, acc, F.inv(F.of(m)))
    return G


def fs_generic_modp(F, s, N, zeta):
    """the same c(z; L) by the generic route: explicit roots of unity zeta in F_p, Taylor series in t, composition (cross-check)"""
    nt = 2 * N + 6; s = F.of(s); u = F.red(s ** 5); iu = F.inv(F.red(1 - u))
    X0 = [s, F.one] + [F.zero] * (nt - 2)
    X5 = mul(F, mul(F, mul(F, X0, X0, nt), mul(F, X0, X0, nt), nt), X0, nt)
    ratio = [F.red((1 - X5[0]) * iu)] + [F.red(-x * iu) for x in X5[1:]]
    Wt = integ(F, mul(F, inv(F, mul(F, X0, X0, nt), nt), log1(F, ratio, nt), nt))[:nt]
    terms = {}; Xm = [F.one] + [F.zero] * (nt - 1)
    for m in range(1, N + 1):
        Xm = mul(F, Xm, X0, nt)
        if m % 2 == 0: continue
        ser = [F.zero] * nt; ser[0] = pow(F.inv(5), m, F.p)
        for j in range(5):
            zj = pow(zeta, j, F.p); aj = F.red(zj * s - 1); rj = F.red(zj * F.inv(aj)); am = pow(F.inv(aj), m, F.p)
            for i in range(nt):
                ser[i] = F.red(ser[i] + am * (-1) ** i * comb(m + i - 1, i) * pow(rj, i, F.p))
        terms[m] = scal(F, mul(F, Xm, ser, nt), F.of(stirling_b((m + 1) // 2)))
    G = bexp(F, terms, N, nt)
    pre = powr(F, ratio, Fr(-1, 2), nt)
    return sp_generic(F, Wt, {m: mul(F, pre, G[m], nt) for m in G}, N)
