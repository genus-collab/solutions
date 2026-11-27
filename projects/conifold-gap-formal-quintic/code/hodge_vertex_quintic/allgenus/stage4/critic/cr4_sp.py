# -*- coding: utf-8 -*-
"""cr4_sp.py -- Stage-4 critic: formal stationary phase (SP) of the full-Stirling Mellin-Barnes integrand at a FIXED u0.

New code (stdlib + cr4_core only).  This is the integrand of poly s.5 (= gamma s.2 after x = s X, z = s^4 w), written in the
hatted variable X at the colour-0 saddle X = 1:
    Phi^SP(w; u0) = SP_(X=1)[ exp(-What(X)/w) Bhat(X, w) dX ] / ( exp(-What(1)/w) sqrt(2 pi w/5) ),
    dWhat/dX = -(1/(u0 X^2)) log((1 - u0 X^5)/(1 - u0)),
    Bhat = (1 - u0 X^5)^(-1/2) exp( Sum_k b_k (u0 X w)^(2k-1) [ 5^(1-2k) - p_(2k-1) ] ),   b_k = B_2k/(2k(2k-1)),
    p_m = Sum_(j=0..4) (1 - zeta^j x)^(-m), x^5 = u0 X^5  (Newton sums of the roots of (1-x^5)Y^5 - 5Y^4 + 10Y^3 - 10Y^2 + 5Y - 1).
With What(1+t) - What(1) = (W2/2) t^2 S(t), W2 = 5/(1-u0), the Morse lemma and Lagrange-Buermann give
    Phi^SP = Sum_r ((1-u0) w/5)^r (2r-1)!! [t^2r] ( (1-Y)^(-1/2) E(t, w) S(t)^(-r-1/2) ),   Y = u0((1+t)^5 - 1)/(1-u0),
E = exp(Sum_k ...).  Everything is rational in u0, so Phi^SP(w; u0) is exact over Q for rational u0 != 1 (u0 = 0 allowed: the
model's own stationary phase).  By Theorem P / Theorem A(ii) Phi^SP is the CANONICAL WKB function; by (Gamma)/(A2) it is the QRR one.
Options (failing controls): b2shift adds to b_2; drop5 removes the Gamma(1+5s) Stirling term 5^(1-2k).
"""
from fractions import Fraction as Fr
from cr4_core import smul, sadd, sinv, spow, bernoulli


def sp_phi(F, u0, N, b2shift=0, drop5=False):
    """[w^n] Phi^SP(w; u0) for n = 0..N"""
    T = 2 * N + 3                                               # t^0 .. t^(2N+2)
    u0 = F.c(u0); one = F.c(1); om = F.r(one - u0); iom = F.inv(om)
    Z = [0, 5, 10, 10, 5, 1] + [0] * T; Z = [F.c(x) for x in Z[:T]]   # (1+t)^5 - 1
    Lg = [0] * T; Zj = [one] + [0] * (T - 1); cj = one              # Sum_j u0^(j-1) Z^j / (j (1-u0)^j) = -log(1-Y)/u0
    for j in range(1, T):
        Zj = smul(F, Zj, Z, T); cj = F.r(cj * iom) if j == 1 else F.r(cj * u0 * iom)
        Lg = sadd(F, Lg, Zj, F.r(cj * F.inv(F.c(j))))
    inv2 = [F.c((-1) ** k * (k + 1)) for k in range(T)]           # (1+t)^(-2)
    Wp = smul(F, inv2, Lg, T)                                       # What'(1+t)
    Wd = [0] + [F.r(Wp[k - 1] * F.inv(F.c(k))) for k in range(1, T)]   # What(1+t) - What(1)
    W2 = F.r(2 * Wd[2]); assert W2 == F.r(5 * iom), "W'' at the saddle"
    i2 = F.inv(F.r(W2 * F.inv(F.c(2))))
    L2 = 2 * N + 1                                                  # need S, G_r to t^(2N)
    S = [F.r(Wd[k + 2] * i2) for k in range(L2)]
    Y = [F.r(u0 * x * iom) for x in Z[:L2]]
    A0 = smul(F, spow(F, [F.r(-y) if k else one for k, y in enumerate(Y)], Fr(-1, 2), L2), spow(F, S, Fr(-1, 2), L2), L2)
    Si = sinv(F, S, L2)
    # amplitude E(t, w): p_m(x) with x^5 = y(t) = u0 (1+t)^5, D = 1 - y
    y5 = [F.r(u0 * x) for x in [1, 5, 10, 10, 5, 1] + [0] * L2][:L2]
    iD = sinv(F, [F.r(one - y5[0])] + [F.r(-v) for v in y5[1:]], L2)
    e = [None] + [[F.r(c * v) for v in iD] for c in (5, 10, 10, 5, 1)]
    p = [None]
    for m in range(1, N + 1):
        s = [0] * L2
        for i in range(1, min(m, 5) + 1):
            term = smul(F, e[i], p[m - i], L2) if i < m else [F.r(m * v) for v in e[i]]
            s = sadd(F, s, term, (-1) ** (i - 1))
        p.append(s)
    B = bernoulli(N + 2)
    c = {}
    for k in range(1, (N + 1) // 2 + 1):
        m = 2 * k - 1
        bk = F.c(B[2 * k] / (2 * k * (2 * k - 1)) + (b2shift if k == 2 else 0))
        pw = [F.c(1)] + [0] * (L2 - 1)
        for _ in range(m): pw = smul(F, pw, [u0, u0] + [0] * (L2 - 2), L2)       # (u0 (1+t))^m
        br = sadd(F, [0 if drop5 else F.c(Fr(5) ** (1 - 2 * k))] + [0] * (L2 - 1), p[m], -1)
        c[m] = [F.r(bk * v) for v in smul(F, pw, br, L2)]
    E = [[one] + [0] * (L2 - 1)]
    for n in range(1, N + 1):
        Ln = 2 * (N - n) + 1; acc = [0] * Ln
        for m in range(1, n + 1, 2): acc = sadd(F, acc, smul(F, c[m], E[n - m], Ln), m)
        E.append([F.r(v * F.inv(F.c(n))) for v in acc])
    out = [0] * (N + 1); G = A0; fac = one; df = one
    for r in range(N + 1):
        if r:
            G = smul(F, G, Si, L2); fac = F.r(fac * om * F.inv(F.c(5))); df = F.r(df * (2 * r - 1))
        for n in range(r, N + 1):
            En = E[n - r]
            s = sum(En[i] * G[2 * r - i] for i in range(min(2 * r, len(En) - 1) + 1))
            out[n] = F.r(out[n] + fac * df * s)
    return out
