# -*- coding: utf-8 -*-
"""rf_curves.py -- curves for the referee's TR engine (rf_tr.Curve): x = P(p) + v/(p - b), y = p (one simple pole), and Airy.

pole_curve(F, Pc, v, b, c, name):  P = Sum_k Pc[k] p^k; branch points c (asserted: x'(c) = 0, x''(c) != 0).
special(F, kappa):   x = kappa (p^4/4 + 1/p), branch points c_j = w^-j (w = F.w, order 5): colour j <-> critical point zeta^-j.
rational_point(F, r, rng, b):  x = P(p) + v/(p - b) with deg P = r, leading coefficient 1/r, and r+1 RATIONAL branch points
   (p^2 P'(p) - v = Prod (p - c_i) needs e_r(c) = 0: c_1..c_r random, c_(r+1) = -e_r/e_(r-1)), pole shifted to b.
airy(F):  x = -p^2/2, y = p, branch point 0.
"""
import os, sys, itertools
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from rf_tr import Curve


def _poly_taylor(F, Pc, a, N):
    out = []
    for k in range(N + 1):
        s = F.zero
        for e in range(k, len(Pc)):
            if Pc[e] != 0:
                s = s + F.of(Pc[e]) * F.of(comb(e, k)) * _pw(F, a, e - k)
        out.append(s)
    return out


def _pw(F, a, e):
    r = F.one
    for _ in range(e): r = r * a
    return r


def pole_curve(F, Pc, v, b, c, name):
    vv = F.of(v); bb = F.of(b)
    def xt(a, N):
        out = _poly_taylor(F, Pc, a, N); d = F.one / (a - bb); pw = d
        for k in range(N + 1):
            out[k] = out[k] + vv * pw * F.of((-1) ** k); pw = pw * d
        return out
    def yxt(a, N):
        X = xt(a, N + 1); dX = [X[k + 1] * F.of(k + 1) for k in range(N + 1)]
        return [a * dX[k] + (dX[k - 1] if k >= 1 else F.zero) for k in range(N + 1)]
    for ci in c:
        X = xt(ci, 3); assert F.iszero(X[1]) and not F.iszero(X[2]), ("not a simple branch point", name)
    return Curve(F, list(c), xt, yxt, name)


def special(F, kappa):
    k = Fr(kappa)
    c = [F.pw(-j) for j in range(5)]
    return pole_curve(F, [0, 0, 0, 0, k / 4], k, 0, c, "x = %s (p^4/4 + 1/p)" % k)


def _esym(xs, k):
    tot = 0
    for cm in itertools.combinations(xs, k):
        pr = 1
        for x in cm: pr = pr * x
        tot = tot + pr
    return tot


def rational_point(F, r, rng, b=0):
    """returns (curve, Pc, v, c) with exact rational data; requires r >= 1"""
    while True:
        cs = [Fr(rng.randrange(-9, 10), rng.randrange(1, 6)) for _ in range(r)]
        if r == 1:
            if cs[0] == 0: continue
            cn = Fr(0)                                       # e_1(c1, c2) = 0 for r = 1 means c2 = -c1
            cn = -cs[0]
        else:
            er, er1 = _esym(cs, r), _esym(cs, r - 1)
            if er1 == 0: continue
            cn = -er / er1
        c = cs + [cn]
        if len(set(c)) == r + 1 and all(x != 0 for x in c): break
    # p^2 P'(p) - v = Prod (p - c_i):  coefficients of Prod in increasing degree
    poly = [Fr(1)]
    for ci in c: poly = [(poly[k - 1] if k >= 1 else 0) - ci * (poly[k] if k < len(poly) else 0) for k in range(len(poly) + 1)]
    assert poly[1] == 0
    v = -poly[0]
    dP = [poly[k + 2] for k in range(r)]                     # P'(p) = Sum_{k} poly[k+2] p^k
    Pc = [Fr(0)] + [dP[k] / (k + 1) for k in range(r)]
    assert Pc[r] == Fr(1, r)
    # shift the pole to b:  x(p) = P(p - b) + v/(p - b), branch points c + b
    if b != 0:
        Pb = [Fr(0)] * (r + 1)
        for e, ce in enumerate(Pc):
            for k in range(e + 1): Pb[k] += ce * comb(e, k) * Fr(-b) ** (e - k)
        Pc = Pb
    cb = [F.of(x + b) for x in c]
    return pole_curve(F, Pc, v, b, cb, "rational point r = %d, b = %s" % (r, b)), Pc, v, c


def airy(F):
    def xt(a, N): return [F.zero, F.zero, F.of(Fr(-1, 2))] + [F.zero] * (N - 2)
    def yxt(a, N): return [F.zero, F.zero, F.of(-1)] + [F.zero] * (N - 2)
    return Curve(F, [F.zero], xt, yxt, "Airy x = -p^2/2")
