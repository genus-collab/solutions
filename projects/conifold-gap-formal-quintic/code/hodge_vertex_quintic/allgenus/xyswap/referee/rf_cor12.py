# -*- coding: utf-8 -*-
"""rf_cor12.py -- F_g (g >= 2) from Corollary 1.2, eq. (1.3), of arXiv:2502.18115 (free energy for y = z, unramified y,
no logarithmic singularities):
   F_g = 1/(2-2g) Sum_{b} Res_{z->b} [hbar^(2g-1)] dy Sum_{m>=2} (d^(m-1) y/dx^(m-1)) [u^m] exp[(S(u hbar d_z) - 1) x u]/(u hbar),
b running over the singular points of x, y (here z = 0 and z = infinity), S(t) = (e^(t/2) - e^(-t/2))/t.
With y = z:  [hbar^(2g-1)] [u^m] exp(E)/(u hbar) = [hbar^(2g)] [u^(m+1)] exp(E),  E = Sum_k u^(2k+1) hbar^(2k) x^(2k)(z)/(4^k (2k+1)!),
and [hbar^(2g)] exp(E) = Sum over partitions lambda of g of u^(2g + len(lambda)) prod_i x^(2 lambda_i)/(4^lambda_i (2 lambda_i + 1)!)
/ prod (multiplicities)!.   x is a Laurent polynomial in z with numeric (rational) coefficients.  Referee's own code.
"""
from fractions import Fraction as Fr
from math import factorial
from collections import Counter
from rf_ser import pinv, lmul, coeff, lp_mul, lp_deriv, partitions


def _lp_scale(L, c):
    return {e: v * c for e, v in L.items() if v}


def hbar_u_coeffs(X, g):
    """{m+1: Laurent polynomial in z} = [hbar^(2g)][u^(m+1)] exp(E)"""
    out = {}
    for lam in partitions(g):
        term = {0: Fr(1)}
        for k in lam:
            term = lp_mul(term, _lp_scale(lp_deriv(X, 2 * k), Fr(1, 4 ** k * factorial(2 * k + 1))))
        mult = 1
        for v in Counter(lam).values():
            mult *= factorial(v)
        term = _lp_scale(term, Fr(1, mult))
        p = 2 * g + len(lam)
        acc = out.setdefault(p, {})
        for e, c in term.items():
            acc[e] = acc.get(e, Fr(0)) + c
    return out


def _laurent_from_lp(L, var_inv):
    """Laurent polynomial in z -> Laurent series (val, list) in z (var_inv False) or in zeta = 1/z (var_inv True)"""
    M = {(-e if var_inv else e): c for e, c in L.items() if c}
    if not M:
        return (0, [])
    v = min(M)
    return (v, [M.get(e, Fr(0)) for e in range(v, max(M) + 1)])


def _dlaurent(A, var_inv):
    """d/dz of a Laurent series in z, or in zeta = 1/z (then d/dz = -zeta^2 d/dzeta)"""
    v, c = A
    d = [(v + j) * a for j, a in enumerate(c)]
    if not var_inv:
        return (v - 1, d)
    return (v + 1, [-a for a in d])


def residue(X, g, at_infinity, T=None):
    """Res_{z->b} of the integrand of (1.3) (without the 1/(2-2g)), b = 0 or infinity"""
    inv = at_infinity
    dX = lp_deriv(X)
    D = _laurent_from_lp(dX, inv)                              # x'(z) as Laurent series in the local variable
    if T is None:
        T = 12 * g + 20
    n = T + 10
    lst = D[1] + [Fr(0)] * (n + 1)
    ixp = (-D[0], pinv(lst[: n + 1], n))                        # 1/x'
    coeffs = hbar_u_coeffs(X, g)
    tot = Fr(0)
    f = _laurent_from_lp({1: Fr(1)}, inv)                       # y = z
    for m in range(1, max(coeffs)):                            # f = D_x^(m-1) z at step m
        if m >= 2 and (m + 1) in coeffs:
            Lp = _laurent_from_lp(coeffs[m + 1], inv)
            prod = lmul(f, Lp, T)
            if not inv:
                tot += coeff(prod, -1)                          # Res_{z=0} (.) dz
            else:
                tot += -coeff(prod, 1)                          # Res_{z=inf} (.) dz, dz = -dzeta/zeta^2
        f = lmul(ixp, _dlaurent(f, inv), T + 4)
    return tot


def F_cor12(X, g):
    r0 = residue(X, g, False)
    ri = residue(X, g, True)
    return r0 / (2 - 2 * g), ri / (2 - 2 * g), (r0 + ri) / (2 - 2 * g)
