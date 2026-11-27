# -*- coding: utf-8 -*-
"""chk_res.py -- Lemmas FL2 and FL3 for the formal residue at X = 0 (termwise in u), exact over Q.

Functional (note s.10, the poly referee's technical remark): for f0 = 1/X + X^4/4 and v = 1/w,
  Res[ e^(-f0/w) X^p ] = coefficient of X^-1 = sum_(l >= 0, k = 4l+p+1 >= 0) (-1)^k/k! (-1/4)^l/l! v^(k+l),
a series in v (entire in 1/w); applied termwise to e^(-f0/w) A, A a Laurent polynomial in X at each order u^a v^s.
FL2: Res[ d( e^(-f/w) Q ) ] = 0 for f = V(u)/X + X^4/4, V = 1 + O(u), Q a random Laurent polynomial (u-dependent).
FL3: for Y = X (1 + sum_k u^k Y_k(X^5)) (a CFU-type map), Res_Y[ e^(-f(Y)/w) b(Y) dY ] = Res_X[ pull-back ].
Controls (must fail): the residue of a non-exact form is non-zero; the pull-back without Y'(X).
Usage: py -I -B chk_res.py [NU] [NV]     (orders u^0..u^NU and v^0..v^NV, defaults 6 and 14)
"""
import os, sys, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import factorial
from chk_util import Log

NU = int(sys.argv[1]) if len(sys.argv) > 1 else 6
NV = int(sys.argv[2]) if len(sys.argv) > 2 else 14
log = Log(HERE, 'chk_res_%d_%d' % (NU, NV))
log.out('chk_res: FL2 and FL3 for the formal residue, u^0..u^%d, v^0..v^%d (v = 1/w); log %s' % (NU, NV, log.rel))


def add(A, B, cb=1):
    R = dict(A)
    for k, c in B.items():
        R[k] = R.get(k, 0) + cb * c
    return {k: c for k, c in R.items() if c}


def mul(A, B):
    R = {}
    for (p, a, s), c in A.items():
        for (p2, a2, s2), c2 in B.items():
            if a + a2 <= NU and s + s2 <= NV:
                k = (p + p2, a + a2, s + s2); R[k] = R.get(k, 0) + c * c2
    return {k: c for k, c in R.items() if c}


def mono(p, a=0, s=0, c=1):
    return {(p, a, s): Fr(c)}


def expo(A):                                   # exp(A), A = O(u) (or O(v)) termwise
    R, T, m = mono(0), mono(0), 0
    while T:
        m += 1; T = {k: c / m for k, c in mul(T, A).items()}; R = add(R, T)
    return R


def dX(A):
    return {(p - 1, a, s): c * p for (p, a, s), c in A.items() if p}


def res(A):
    """coefficient of X^-1 of e^(-f0/w) A, as {(a, s): c}"""
    out = {}
    for (p, a, s), c in A.items():
        l = 0
        while True:
            k = 4 * l + p + 1
            if s + k + l > NV:
                break
            if k >= 0:
                key = (a, s + k + l)
                out[key] = out.get(key, 0) + c * Fr((-1) ** k, factorial(k)) * Fr((-1) ** l, 4 ** l * factorial(l))
            l += 1
    return {k: c for k, c in out.items() if c}


def binom_pow(Hu, e):                          # (1 + Hu)^e, Hu = O(u), e an integer
    R, T, m = mono(0), mono(0), 0
    while T:
        T = {k: c * Fr(e - m, m + 1) for k, c in mul(T, Hu).items()}; m += 1; R = add(R, T)
    return R


rng = random.Random(86028121)
def rand_laurent(pmin, pmax, amax, smax=0, dens=0.5):
    return {(p, a, s): Fr(rng.randint(-3, 3)) for p in range(pmin, pmax + 1) for a in range(amax + 1)
            for s in range(smax + 1) if rng.random() < dens and rng.randint(-3, 3)}


V1 = {(0, a, 0): Fr(rng.randint(-3, 3), rng.randint(1, 3)) for a in range(1, NU + 1)}   # V - 1 = O(u)
V = add(mono(0), V1)
EV = expo(mul(mono(-1, 0, 1, -1), V1))          # exp(-(V - 1) v / X)
fp = add(mul(V, mono(-2, 0, 0, -1)), mono(3))  # f' = -V/X^2 + X^3
# FL2
for trial in range(3):
    Q = rand_laurent(-4, 4, 3, 1)
    A = mul(EV, add(dX(Q), mul(mul(fp, Q), mono(0, 0, 1)), -1))      # e^(f0/w) d(e^(-f/w) Q)/dX
    r = res(A)
    log.check('FL2 trial %d: Res[d(e^(-f/w) Q)] = 0 through u^%d, v^%d (Q has %d terms)' % (trial, NU, NV, len(Q)), not r,
              'non-zero entries: %s' % sorted(r)[:3])
r1 = res(mul(EV, mono(0)))
log.control('FL2: the residue of e^(-f/w) dX (not exact) is non-zero', bool(r1),
            'e.g. [u^0 v^1, v^5] = %s, %s' % (r1.get((0, 1), 0), r1.get((0, 5), 0)))
# FL3: Y = X (1 + u h(X)), h a polynomial in X^5 with u-dependence (a CFU-type map)
h = {(5 * i, a, 0): Fr(rng.randint(-3, 3), rng.randint(1, 4)) for i in range(3) for a in range(NU)}
uh = mul(mono(0, 1), h)                         # u h
beta = rand_laurent(-3, 3, 2, 1)
EVY = EV                                        # in the Y-variable: exp(-(V-1) v / Y)
lhs = res(mul(EVY, beta))
inv1 = binom_pow(uh, -1); pow4 = binom_pow(uh, 4)
D = add(add(mul(mul(V, mono(-1)), inv1), mul(mono(4, 0, 0, Fr(1, 4)), pow4)), add(mono(-1), mono(4, 0, 0, Fr(1, 4))), -1)
ED = expo(mul(D, mono(0, 0, 1, -1)))            # exp(-(f(Y(X)) - f0(X)) v)
bY = {}
for (p, a, s), c in beta.items():
    bY = add(bY, mul({(p, a, s): c}, binom_pow(uh, p)))
Yp = add(mono(0), dX(mul(mono(1), uh)))         # Y'(X) = 1 + d/dX (X u h)
rhs = res(mul(mul(ED, bY), Yp))
log.check('FL3 (residue): Res_Y[e^(-f/w) b dY] = Res_X[pull-back by Y = X(1 + u h(X^5))] through u^%d, v^%d' % (NU, NV),
          lhs == rhs, '%d coefficients; [u^1 v^3] = %s' % (len(lhs), lhs.get((1, 3), 0)))
rbad = res(mul(ED, bY))
log.control('FL3 (residue): pull-back without Y\'(X)', rbad != lhs)
log.finish('CHK_RES %d %d' % (NU, NV))
