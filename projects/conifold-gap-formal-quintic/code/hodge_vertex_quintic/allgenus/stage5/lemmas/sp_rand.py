# -*- coding: utf-8 -*-
"""sp_rand.py -- random test data for the Stage-5 lemma checks (new code, stdlib only)."""
from fractions import Fraction as Fr
from sp_core import rand_poly


def det_q(M):
    n = len(M)
    if n == 1:
        return Fr(M[0][0])
    return sum((-1) ** b * Fr(M[0][b]) * det_q([[M[a][c] for c in range(n) if c != b] for a in range(1, n)])
               for b in range(n))


def rand_sym(rng, n, nx=None):
    """random symmetric integer matrix, det != 0 (and the leading nx-block invertible if nx given)"""
    while True:
        H = [[0] * n for _ in range(n)]
        for a in range(n):
            for b in range(a, n):
                H[a][b] = H[b][a] = (-rng.randint(2, 5) if a == b else rng.randint(-2, 2))
        if det_q(H) != 0 and (nx is None or det_q([r[:nx] for r in H[:nx]]) != 0):
            return H


def quad_form(H, vars_, nv):
    """(1/2) x^T H x as a series in the variables vars_"""
    F = {}
    for a, va in enumerate(vars_):
        for b in range(a, len(vars_)):
            e = [0] * nv; e[va] += 1; e[vars_[b]] += 1
            v = Fr(H[a][b], 2) if a == b else Fr(H[a][b])
            if v:
                F[(tuple(e), 0)] = v
    return F


def rand_phase(cb, rng, vars_, dmax, nx=None):
    H = rand_sym(rng, len(vars_), nx)
    F = cb.add(quad_form(H, vars_, cb.nv), rand_poly(cb, rng, vars_, 3, dmax, 3, 0.7))
    return H, F


def rand_amp(c, rng, vars_):
    A = c.add(c.const(1), rand_poly(c, rng, vars_, 1, 3, 2, 0.6))
    Z = rand_poly(c, rng, vars_, 0, 2, 2, 0.6)
    return c.add(A, c.zp(Z, 1))
