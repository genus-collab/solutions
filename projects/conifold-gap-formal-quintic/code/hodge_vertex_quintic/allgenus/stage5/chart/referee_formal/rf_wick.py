# -*- coding: utf-8 -*-
"""rf_wick.py -- referee: multivariable formal stationary phase, written for this report (independent of ch_wick.py).

N[Phi, A; p](z) = Sum_n N_n z^n, N_n = < [exp(Phi_3/z) A]_(z^n) >_C, C = -H^-1 (CHART_PROOF s.1.2), computed as
  E = exp(Sum_(d>=3) Phi_d s^(d-2))  (weight recursion w E_w = Sum_i i Phi_(i+2) E_(w-i)),  F_(2n) = Sum E_w1 Amp_(2n-w1),
  Amp_w = Sum_(2m+j=w) [A_m]_j,   N_n = < F_(2n) >.
Gaussian expectations, two independent algorithms (neither is the Isserlis recursion of the worker):
  expect1:   C = diag(C0) + c u u^T  (eta = eta0 + u gamma, eta0, gamma independent): closed-form binomial/double-factorial moments;
  expect_do: <P_2k> = (1/k!) (Delta_C)^k P_2k, Delta_C = (1/2) Sum C_ij d_i d_j  (differential-operator form), any C.
Polynomials are dicts {exponent tuple: coefficient}.
"""
from fractions import Fraction as Fr
from math import comb, factorial
from rf_core import dfact


def padd(R, A, f=None):
    for e, c in A.items():
        v = c if f is None else c * f
        R[e] = R[e] + v if e in R else v


def pmul(A, B):
    R = {}
    for e1, c1 in A.items():
        for e2, c2 in B.items():
            e = tuple(x + y for x, y in zip(e1, e2)); v = c1 * c2
            R[e] = R[e] + v if e in R else v
    return R


def split_deg(P):
    out = {}
    for e, c in P.items(): out.setdefault(sum(e), {})[e] = c
    return out


def wickN(Phi, Amp, N, nvar, expect, one):
    """Phi: Taylor polynomial of Phi - Phi(p) (degree-1 part must vanish; degree-2 part ignored here, it enters through C);
    Amp[m]: Taylor polynomial of A_m. Returns N_0..N_N."""
    zero = one - one; Pd = split_deg(Phi); z0 = (0,) * nvar
    assert all(c == 0 for c in Pd.get(1, {}).values()), 'not a critical point'
    E = [{z0: one}]
    for w in range(1, 2 * N + 1):
        acc = {}
        for i in range(1, w + 1):
            if Pd.get(i + 2): padd(acc, pmul(Pd[i + 2], E[w - i]), Fr(i, w))
        E.append(acc)
    Ad = [split_deg(a) for a in Amp]
    res = []
    for n in range(N + 1):
        F = {}
        for w1 in range(2 * n + 1):
            w2 = 2 * n - w1
            for m in range(w2 // 2 + 1):
                if m < len(Ad) and Ad[m].get(w2 - 2 * m): padd(F, pmul(E[w1], Ad[m][w2 - 2 * m]))
        res.append(expect(F))
    return res


def make_expect1(C0, c, u, one):
    zero = one - one; n = len(C0); D = 64
    m0 = [[(one * dfact(k - 1) * _pw(C0[i], k // 2, one)) if k % 2 == 0 else zero for k in range(D)] for i in range(n)]
    up = [[_pw(u[i], e, one) for e in range(D)] for i in range(n)]
    mg = [(one * dfact(k - 1) * _pw(c, k // 2, one)) if k % 2 == 0 else zero for k in range(D)]
    nz0 = [not (C0[i] == 0) for i in range(n)]

    def expect(P):
        tot = zero
        for a, coef in P.items():
            acc = {0: one}
            for i in range(n):
                ai = a[i]
                if ai == 0: continue
                new = {}; cs = range(0, ai + 1, 2) if nz0[i] else (0,)
                for g, val in acc.items():
                    for cc in cs:
                        e = ai - cc; f = m0[i][cc] * up[i][e] * comb(ai, cc)
                        new[g + e] = new[g + e] + val * f if (g + e) in new else val * f
                acc = new
            s = zero
            for g, val in acc.items():
                if g % 2 == 0: s = s + val * mg[g]
            tot = tot + coef * s
        return tot
    return expect


def _pw(x, k, one):
    r = one
    for _ in range(k): r = r * x
    return r


def make_expect_do(C, one):
    n = len(C); zero = one - one

    def lap(P):
        R = {}
        for e, c in P.items():
            for i in range(n):
                for j in range(n):
                    if C[i][j] == 0: continue
                    f = list(e)
                    if i == j:
                        if f[i] < 2: continue
                        k = f[i] * (f[i] - 1); f[i] -= 2
                    else:
                        if f[i] < 1 or f[j] < 1: continue
                        k = f[i] * f[j]; f[i] -= 1; f[j] -= 1
                    t = tuple(f); v = c * C[i][j] * Fr(k, 2)
                    R[t] = R[t] + v if t in R else v
        return R

    def expect(P):
        tot = zero
        for d, Q in split_deg(P).items():
            if d % 2: continue
            for _ in range(d // 2): Q = lap(Q)
            tot = tot + Q.get((0,) * n, zero) * Fr(1, factorial(d // 2))
        return tot
    return expect


def mat_inv(M, one):
    n = len(M); zero = one - one
    A = [list(r) + [one if i == j else zero for j in range(n)] for i, r in enumerate(M)]
    for col in range(n):
        piv = next(r for r in range(col, n) if not (A[r][col] == 0))
        A[col], A[piv] = A[piv], A[col]; iv = 1 / A[col][col]
        A[col] = [x * iv for x in A[col]]
        for r in range(n):
            if r != col and not (A[r][col] == 0):
                f = A[r][col]; A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    return [row[n:] for row in A]


def det(M, one):
    n = len(M); A = [list(r) for r in M]; d = one
    for col in range(n):
        piv = next((r for r in range(col, n) if not (A[r][col] == 0)), None)
        if piv is None: return one - one
        if piv != col: A[col], A[piv] = A[piv], A[col]; d = -d
        d = d * A[col][col]; iv = 1 / A[col][col]
        for r in range(col + 1, n):
            f = A[r][col] * iv; A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    return d


def hessian(Phi, nvar, one):
    zero = one - one; H = [[zero] * nvar for _ in range(nvar)]
    for e, c in Phi.items():
        if sum(e) != 2: continue
        idx = [i for i in range(nvar) for _ in range(e[i])]
        if idx[0] == idx[1]: H[idx[0]][idx[0]] = c * 2
        else: H[idx[0]][idx[1]] = c; H[idx[1]][idx[0]] = c
    return H


def uni(nvar, i, ser, kmin=0):
    """dict polynomial Sum_k ser[k] x_i^k (k >= kmin)"""
    return {tuple(k if a == i else 0 for a in range(nvar)): c for k, c in enumerate(ser) if k >= kmin and not (c == 0)}


def trunc(P, D): return {e: c for e, c in P.items() if sum(e) <= D}


def tmul(A, B, D):
    R = {}
    for e1, c1 in A.items():
        s1 = sum(e1)
        for e2, c2 in B.items():
            if s1 + sum(e2) > D: continue
            e = tuple(x + y for x, y in zip(e1, e2)); v = c1 * c2
            R[e] = R[e] + v if e in R else v
    return R
