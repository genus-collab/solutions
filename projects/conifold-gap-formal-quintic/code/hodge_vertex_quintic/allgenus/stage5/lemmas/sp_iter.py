# -*- coding: utf-8 -*-
"""sp_iter.py -- Stage 5 / lemmas: total and iterated formal SP, built on sp_core (new code, stdlib only).

Two contexts are used: cb (W2 + 2) for phase data, c (W2) for the SP itself (phase data enter divided by z).
total_sp   : SP^norm over all integration variables at the critical point 0 (Definition 1.3 of the note).
inner_sp   : SP^norm over the variables xs at the critical point xi(y) of x -> F(x, y), as a series in the
             remaining variables y and z (the inner integral of Lemma L1b), together with G(y) = F(xi(y), y)
             and Delta(y) = det(K(y)), K = x-Hessian at xi(y).
iterated_sp: the outer SP^norm of G with amplitude (Delta(y)/Delta(0))^(-1/2) b(y; z) (the right side of L1b).
"""
from fractions import Fraction as Fr
from sp_core import Ctx, inv_q, mmul


def hessian_const(c, F, vars_):
    """constant Hessian of F at 0 in the variables vars_ (F's quadratic part, constant coefficients)"""
    n = len(vars_); H = [[Fr(0)] * n for _ in range(n)]
    for (e, j), v in F.items():
        if sum(e) == 2 and all(e[i] == 0 for i in range(c.nv) if i not in vars_) and j == 0:
            idx = [a for a, i in enumerate(vars_) for _ in range(e[i])]
            a, b = idx
            if a == b:
                H[a][a] += 2 * v
            else:
                H[a][b] += v; H[b][a] += v
    return H


def check_critical(c, F, vars_):
    """F(0) = 0 and grad F(0) = 0: no term of total degree <= 1"""
    return not any(sum(e) <= 1 and v != 0 for (e, j), v in F.items())


def total_sp(c, cb, F, amp, vars_):
    """normalised SP over vars_ (all variables of F) at 0; F z-free, F(0) = 0, grad F(0) = 0"""
    assert check_critical(c, F, vars_), "0 is not a critical point"
    H = hessian_const(c, F, vars_)
    C = [[c.const(-x) for x in row] for row in inv_q(H)]
    P3 = {k: v for k, v in F.items() if sum(k[0]) >= 3}
    return c.sp_norm(P3, c.trunc(amp), vars_, C)


def crit_point(cb, F, xs):
    """xi(y) with d_x F(xi(y), y) = 0, xi(0) = 0 (formal implicit function theorem, fixed-point iteration)"""
    H = hessian_const(cb, F, xs); Hi = inv_q(H)
    grads = [cb.deriv(F, i) for i in xs]
    xi = {i: {} for i in xs}
    for _ in range(cb.W2 + 3):
        g = [cb.compose(gr, xi) for gr in grads]
        new = {}
        for a, i in enumerate(xs):
            s = dict(xi[i])
            for b in range(len(xs)):
                s = cb.add(s, g[b], -Hi[a][b])
            new[i] = cb.trunc(s)
        if new == xi:
            break
        xi = new
    g = [cb.compose(gr, xi) for gr in grads]
    assert all(not x for x in g), "critical point not converged"
    return xi


def expand_at(cb, F, xs, xi):
    """F(xi(y) + x, y) split by x-degree: G (deg 0), K (Hessian matrix of series), R (deg >= 3)"""
    subs = {i: cb.add(xi[i], cb.var(i)) for i in xs}
    Fs = cb.compose(F, subs); parts = cb.split(Fs, xs)
    assert not parts.get(1), "linear term at the inner critical point"
    G = parts.get(0, {}); K = cb.quad_matrix(parts.get(2, {}), xs)
    R = {}
    for d, P in parts.items():
        if d >= 3:
            R = cb.add(R, P)
    return G, K, R


def inner_sp(c, cb, F, amp_fn, xs, freeze_cov=False, at_zero=False):
    """amp_fn(subs) -> amplitude as a series after the substitution x -> xi + x' (subs given).
    Controls: freeze_cov uses the inner Hessian at y = 0; at_zero expands at x = 0 and drops the linear term."""
    xi = {i: {} for i in xs} if at_zero else crit_point(cb, F, xs)
    subs = {i: cb.add(xi[i], cb.var(i)) for i in xs}
    if at_zero:
        Fs = cb.compose(F, subs); parts = cb.split(Fs, xs)
        G = parts.get(0, {}); K = cb.quad_matrix(parts.get(2, {}), xs)
        R = {}
        for d, P in parts.items():
            if d >= 3:
                R = cb.add(R, P)
    else:
        G, K, R = expand_at(cb, F, xs, xi)
    if freeze_cov:
        K = [[cb.const(cb.cpart(x)) for x in row] for row in K]
    Kc = [[c.trunc(x) for x in row] for row in K]
    Ki = c.inv(Kc); C = [[c.scal(x, -1) for x in row] for row in Ki]
    amp = c.trunc(amp_fn(subs))
    b = c.sp_norm(R, amp, xs, C)
    return G, b, c.det(Kc)


def iterated_sp(c, cb, F, amp_fn, xs, ys, nodet=False, **kw):
    G, b, D = inner_sp(c, cb, F, amp_fn, xs, **kw)
    G0 = cb.cpart(G); Gs = cb.add(G, cb.const(G0), -1)
    D0 = c.cpart(D); ratio = c.add(c.scal(D, 1 / D0), c.const(1), -1)
    A = b if nodet else c.mul(c.binom_pow(ratio, Fr(-1, 2)), b)
    return G0, total_sp(c, cb, Gs, A, ys)
