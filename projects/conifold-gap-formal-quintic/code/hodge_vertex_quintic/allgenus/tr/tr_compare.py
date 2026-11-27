# -*- coding: utf-8 -*-
"""tr_compare.py -- TR output -> correlators of the frozen theory, and comparisons (TR_MODEL.md s.3).

For ModelCurve(kappa) (x = kappa W, u = 1) the dictionary of TR_MODEL.md s.2 gives, with Lambda = -5 kappa,
    <prod tau_{d_l}(v_l)>^frozen_g = Lambda^(1-g) kappa^(3g-3+n-Sum d) [coefficient of prod d xi_{v_l, d_l} in omega_{g,n}],
v_l = et_i (colour basis, xi = xit^i) or H^a (flat basis, xi = xif_a); descendants xi_{v,d} = (d/dx)^d xi_v.
Functions: correlators(tr, XB, g, n, kind) -> {((i or a, d), ..): value}; frozen_compare(F, corr, Frozen, g, n, kappa) -> (n_compared, n_bad).
"""
import os, sys, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)


def correlators(tr, XB, g, n, kind='c'):
    """frozen-normalised correlators from omega_{g,n}; returns (dict, number of non-zero conversion remainders)"""
    F = tr.F; kap = tr.C.kappa; Lam = -5 * kap
    T, nbad = XB.convert(tr.omega(g, n), n, kind)
    out = {}
    for key, v in T.items():
        D = 3 * g - 3 + n - sum(d for (_, d) in key)
        out[key] = v * F.of(Lam ** (1 - g) * kap ** D)
    return out, nbad


def frozen_compare(corr, Fz, g, n, F):
    """compare colour-basis correlators (field F = make_fp(p, zeta of FpRing)) with the project's MVertex6 tensors.
    Every d-tuple with psi-deficit 0 <= D <= 6 and every colour tuple is compared (zero entries included)."""
    p = Fz.p; ncmp = 0; bad = []; nnz = 0
    dmax = 3 * g - 3 + n
    for d in itertools.product(range(dmax + 1), repeat=n):
        D = 3 * g - 3 + n - sum(d)
        if D < 0 or D > 6: continue
        T = Fz.colour(g, d)
        for idx, val in T.items():
            key = tuple((i, dd) for i, dd in zip(idx, d))
            mine = corr.get(key)
            mv = mine.v if mine is not None else 0
            ncmp += 1
            if val % p: nnz += 1
            if (mv - val) % p: bad.append((key, mv, val))
    return ncmp, nnz, bad


def string_failures(corr, g, n, lower):
    """flat basis: <tau_0(1) X>_g = Sum_j <X with d_j - 1>_g for all X (n-1 points) with keys in corr / lower"""
    bad = 0; cnt = 0
    keys = set(k[1:] for k in corr if k[0] == (0, 0)) | set(lower)
    for X in keys:
        lhs = corr.get(((0, 0),) + X)
        rhs = None
        for j, (a, d) in enumerate(X):
            if d == 0: continue
            Y = X[:j] + ((a, d - 1),) + X[j + 1:]
            v = lower.get(Y)
            if v is not None: rhs = v if rhs is None else rhs + v
        z = (lhs is None or lhs.iszero()) and (rhs is None or rhs.iszero())
        if not z and (lhs is None or rhs is None or not (lhs - rhs).iszero()): bad += 1
        cnt += 1
    return cnt, bad


def dilaton_failures(corr, g, n, lower):
    """flat basis: <tau_1(1) X>_g = (2g-2+n-1) <X>_g"""
    bad = 0; cnt = 0
    keys = set(k[1:] for k in corr if k[0] == (0, 1)) | set(lower)
    for X in keys:
        lhs = corr.get(((0, 1),) + X); rhs = lower.get(X)
        ok = (lhs is None or lhs.iszero()) and (rhs is None or rhs.iszero())
        if not ok:
            ok = lhs is not None and rhs is not None and (lhs - rhs * (2 * g - 2 + n - 1)).iszero()
        bad += 0 if ok else 1; cnt += 1
    return cnt, bad


def bernoulli(n):
    """B_n (B_1 = -1/2) as a Fraction"""
    from math import comb
    B = [Fr(1)]
    for m in range(1, n + 1):
        B.append(-sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1))
    return B[n]


def Fg_predicted(g, kappa=1):
    """F_g of x = kappa W at u = 1 if F_g[x = W] = B_2g/(2g(2g-2)):  kappa^(2-2g) B_2g/(2g(2g-2))"""
    return Fr(kappa) ** (2 - 2 * g) * bernoulli(2 * g) / (2 * g * (2 * g - 2))
