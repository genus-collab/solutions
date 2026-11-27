# -*- coding: utf-8 -*-
"""lm_wick.py -- referee lmi_adversary: Bernoulli numbers/polynomials, the graded exponential and the one-variable Wick
engine (split off lm_core.py; stdlib only).
* gexp: exponential of a series in (z, t) graded by the weight 2m + j of z^m t^j (no constant term).
* wick1: N_n = Sum_r (2r-1)!! C^r [s^(2n) tau^(2r)] E,  E = exp(Sum_k p_k s^(k-2) tau^k) * Sum a_(m,j) s^(2m+j) tau^j,
  i.e. the Wick expansion of CP s.1.2 at the critical point (xi = s tau, s = z^(1/2), <tau^(2r)> = (2r-1)!! C^r).
"""
from fractions import Fraction as Fr
from lm_core import LS, lsc, F

# ---- Bernoulli numbers and polynomials (over Q; convert with F) ------------------------------------------------------
_BN = [Fr(1)]


def bern(n):
    """B_n with B_1 = -1/2 (Sum_(k<=m) binom(m+1, k) B_k = 0)"""
    while len(_BN) <= n:
        m = len(_BN)
        s = Fr(0)
        b = 1
        for k in range(m):
            s += b * _BN[k]
            b = b * (m + 1 - k) // (k + 1)
        _BN.append(-s / (m + 1))
    return _BN[n]


def binom(x, k):
    """generalised binomial coefficient, x rational, k >= 0 integer"""
    r = Fr(1)
    for i in range(k):
        r = r * (Fr(x) - i) / (i + 1)
    return r


def bernpoly(n, h):
    """B_n(h) = Sum_k binom(n, k) B_k h^(n-k)  (DLMF 24.2.5)"""
    h = Fr(h)
    return sum(binom(n, k) * bern(k) * h ** (n - k) for k in range(n + 1))


def beta(n, h):
    """beta_n(h) = (-1)^(n+1) B_(n+1)(h)/(n(n+1)): log Gamma(w+h) ~ (w+h-1/2) log w - w + (1/2) log 2pi + Sum beta_n w^-n"""
    return (-1) ** (n + 1) * bernpoly(n + 1, h) / (n * (n + 1))


def bk(k):
    return bern(2 * k) / (2 * k * (2 * k - 1))


# ---- graded exponential and the one-variable Wick engine -------------------------------------------------------------
def gexp(H, W):
    """H: dict (m, j) -> LS (z^m t^j), weight 2m+j >= 1; returns exp(H) up to weight W as dict (m, j) -> LS."""
    Hw = [[] for _ in range(W + 1)]
    for (m, j), x in H.items():
        w = 2 * m + j
        if w == 0:
            raise ValueError('gexp: constant term')
        if w <= W:
            Hw[w].append((m, j, x))
    E = [dict() for _ in range(W + 1)]
    E[0][(0, 0)] = lsc(1)
    for w in range(1, W + 1):
        acc = {}
        for i in range(1, w + 1):
            for (m1, j1, x) in Hw[i]:
                xi = x.scal(F(i))
                for (m2, j2), y in E[w - i].items():
                    key = (m1 + m2, j1 + j2)
                    t = xi * y
                    acc[key] = acc[key] + t if key in acc else t
        iw = F(Fr(1, w))
        E[w] = {k: v.scal(iw) for k, v in acc.items()}
    out = {}
    for w in range(W + 1):
        out.update(E[w])
    return out


def vertex_table(p, W):
    """ev[w][V]: coefficient of s^w tau^(w+2V) in exp(Sum_(k>=3) p_k s^(k-2) tau^k)."""
    ev = [dict() for _ in range(W + 1)]
    ev[0][0] = lsc(1)
    for w in range(1, W + 1):
        for V in range(1, w + 1):
            s = None
            for i in range(1, w - V + 2):
                pk = p.get(i + 2)
                if pk is None or (V - 1) not in ev[w - i]:
                    continue
                t = pk.scal(F(i)) * ev[w - i][V - 1]
                s = t if s is None else s + t
            if s is not None:
                ev[w][V] = s.scal(F(Fr(1, w)))
    return ev


def wick1(p, a, N, C):
    """p: dict k -> LS (k >= 3), a: dict (m, j) -> LS normalised amplitude, C = -1/(2 p_2). Returns [N_0, .., N_N]."""
    ev = vertex_table(p, 2 * N)
    CD = [lsc(1)] + [None] * (3 * N)
    Cp = lsc(1)
    df = 1
    for r in range(1, 3 * N + 1):
        Cp = Cp * C
        df *= 2 * r - 1
        CD[r] = Cp.scal(F(df))
    res = []
    for n in range(N + 1):
        tot = LS(0, [])
        for (m, j), am in a.items():
            w = 2 * n - 2 * m - j
            if w < 0:
                continue
            for V, ew in ev[w].items():
                tot = tot + (am * ew) * CD[n - m + V]
        res.append(tot)
    return res
