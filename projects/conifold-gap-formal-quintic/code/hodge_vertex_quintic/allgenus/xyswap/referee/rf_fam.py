# -*- coding: utf-8 -*-
"""rf_fam.py -- the Hurwitz family x = z^4 + c3 z^3 + c2 z^2 + c1 z + V/z, y = z, built from five RATIONAL branch points
r_1..r_5 (z^2 x' = 4 prod (z - r_j) forces e_4(r) = 0, i.e. Sum 1/r_j = 0), plus legs (functions on the curve given by
their Taylor coefficients at a point).  Referee's own code.
"""
from fractions import Fraction as Fr
from rf_ser import taylor_laurentpoly, taylor_log, lp_mul


def elem(r):
    e = [Fr(1)] + [Fr(0)] * len(r)
    for x in r:
        for k in range(len(r), 0, -1):
            e[k] += e[k - 1] * x
    return e


def family_x(r):
    """x as a Laurent polynomial and the parameters (c3, c2, c1, V) from the five branch points"""
    r = [Fr(a) for a in r]
    e = elem(r)
    assert e[4] == 0, "need Sum 1/r_j = 0"
    c3, c2, c1, V = -4 * e[1] / 3, 2 * e[2], -4 * e[3], 4 * e[5]
    X = {4: Fr(1), 3: c3, 2: c2, 1: c1, -1: V}
    return {k: v for k, v in X.items() if v}, (c3, c2, c1, V)


def phi_parts(X):
    """Phi = Int z x'(z) dz = (polynomial part) + (log coefficient) log z"""
    pol, lg = {}, Fr(0)
    for e, c in X.items():
        if e == -1:
            lg += -c                       # z * (-V z^-2) = -V/z -> -V log z
        else:
            pol[e + 1] = pol.get(e + 1, Fr(0)) + c * e / (e + 1)
    return pol, lg


def leg_lp(L):
    """leg = Laurent polynomial in z"""
    return lambda a, K: taylor_laurentpoly(L, a, K)


def leg_log(f):
    """leg = f(z) * log z, f a Laurent polynomial in z (the constant log a_i is dropped: legitimate only when f is a
    function of x, since Sum_a Res_a f(x) omega = 0; for general f the caller must not drop it)"""
    def leg(a, K):
        F = taylor_laurentpoly(f, a, K)
        G = taylor_log(a, K)
        return [sum((F[i] * G[k - i] for i in range(k + 1)), Fr(0)) for k in range(K + 1)]
    return leg


def leg_sum(*legs):
    def leg(a, K):
        out = [Fr(0)] * (K + 1)
        for lg in legs:
            for k, c in enumerate(lg(a, K)):
                out[k] += c
        return out
    return leg


def leg_scale(c, lg):
    return lambda a, K: [c * x for x in lg(a, K)]
