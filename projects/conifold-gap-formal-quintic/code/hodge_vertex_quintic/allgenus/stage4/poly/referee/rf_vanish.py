# -*- coding: utf-8 -*-
"""rf_vanish.py -- REFEREE: the two vanishing solutions (own code, no project or worker module).

Model (GI Lemma V, re-derived by the referee): I^C_m(v, z) = Res_(x=0) x^m exp(-(x^4/4 + v/x)/z) dx
    = Sum_(i >= 0, j = m + 4i + 1 >= 0) (-1)^(i+j) v^j z^-(i+j) / (4^i i! j!),
Y^C = (I_0, -I_-1, -I_3, I_2, I_1)^C, pulled back by v(u) = u tau^(-5/4).
Quintic (GI Lemma Q): y^(0) = -e_1, y^(n+1)_a = Int_0^u b_a y^(n)_(a+1) du, Y^C_Q = Sum_n z^-n y^(n); b_1 y_2 = (u b_1)(y_2/u).
Both are returned as {n: u-series} = coefficient of z^-n.
"""
from fractions import Fraction as Fr
from math import factorial

ROWS = ((0, 1), (-1, -1), (3, -1), (2, 1), (1, 1))     # (m, sign) of Y^C_b


def model_vanishing(S, v, NU):
    """[Y^C_b at z^-n] for b = 0..4 (u-series mod u^K); terms with v^j, j <= NU"""
    p = S.p; vp = [S.one()]
    for j in range(1, NU + 1): vp.append(S.mul(vp[-1], v))
    out = []
    for m, sg in ROWS:
        d = {}; i = 0
        while True:
            j = m + 4 * i + 1
            if j > NU: break
            if j >= 0:
                c = Fr(sg * (-1) ** (i + j), 4 ** i * factorial(i) * factorial(j))
                n = i + j
                d[n] = S.add(d.get(n, S.z()), vp[j], c.numerator * pow(c.denominator, -1, p) % p)
            i += 1
        out.append(d)
    return out


def quintic_vanishing(S, b, NU):
    """list over n of the 5-vectors y^(n) (u-series); b[1] is u b_1; stops when y^(n) = 0 mod u^(NU+1)"""
    y = [S.z() for _ in range(5)]; y[1] = [(-1) % S.p] + [0] * (S.K - 1)
    ys = [y]
    while True:
        prev = ys[-1]; new = []
        for a in range(5):
            if a == 1: integ = S.mul(b[1], S.shift(prev[2], -1))
            else: integ = S.mul(b[a], prev[(a + 1) % 5])
            new.append(S.integ(integ))
        if all(not any(x[:NU + 1]) for x in new): break
        ys.append(new)
    return ys
