# -*- coding: utf-8 -*-
"""rf_gchk.py -- REFEREE (univ): vanishing (pure-1/z) solutions, the identity (G), symplecticity of the flat-frame gauge.

Model vanishing solution (closed formula of the brief):  I^C_m = Sum_i (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!,
   Y^C_model = (I_0, -I_-1, -I_3, I_2, I_1)(v(u));  as {n: (component a = (1-n) mod 5, LS in u)} for the coefficient of z^-n.
System vanishing solution: th Y_a = b_a Y_{a+1}, Y = Sum_n z^-n y_n e_{(1-n) mod 5}, y_0 = -1, y_n' = b_{(1-n) mod 5} y_{n-1},
   y_n(0) = 0 (no logarithm may occur; asserted).
(G): P Y^C_model = c0 Y^C_system with c0 = P_11(0, 0); checked coefficient by coefficient in z (all z^m, m in a range).
Symplecticity: P(-z)^T eta P(z) = eta, eta_{a, 3-a} = 1 (the pairing preserved by self-dual cyclic systems).
"""
import os, sys, math
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf_series import LS


def yc_model(v, nmax, prec):
    comp = {0: (0, 1), 1: (-1, -1), 2: (3, -1), 3: (2, 1), 4: (1, 1)}
    vp = [LS.const(1, prec)]
    for e in range(1, nmax + 2): vp.append(vp[-1].mul(v).trunc(prec))
    out = {}
    for a, (m, sg) in comp.items():
        i = 0
        while True:
            e = m + 4 * i + 1; n = e + i
            if n > nmax: break
            if e >= 0:
                c = Fr((-1) ** e) * Fr(-1, 4) ** i / (math.factorial(i) * math.factorial(e)) * sg
                assert (1 - n) % 5 == a and n not in out
                out[n] = (a, vp[e].scale(c))
            i += 1
    return out


def yc_system(b, nmax, prec):
    out = {0: (1, LS.const(-1, prec))}
    for n in range(1, nmax + 1):
        a = (1 - n) % 5; integrand = b[a].mul(out[n - 1][1])
        assert integrand.v() >= 0, "logarithm in y_%d" % n
        out[n] = (a, integrand.integ().trunc(prec))
    return out


def PY(p, Y, m, prec):
    """coefficient of z^m of P Y (vector of 5 LS): Sum_{k - n = m} p[k][a(n)] y_n placed in component a(n) + k"""
    vec = [LS.zero(prec) for _ in range(5)]
    for n, (a, y) in Y.items():
        k = m + n
        if 0 <= k < len(p):
            vec[(a + k) % 5] = vec[(a + k) % 5].add(p[k][a].mul(y)).trunc(prec)
    return vec


def check_G(p, b, v, prec, mmin=-6, mmax=None):
    """returns (ok, c0, list of failures, number of z-powers compared, positive part non-zero at m = ...)"""
    KG = len(p) - 1; mmax = KG - 7 if mmax is None else mmax
    nmax = KG - mmin
    Ym = yc_model(v, nmax, prec + 2); Ys = yc_system(b, nmax, prec + 2)
    c0 = p[0][1].coef(0)
    bad = []; pos = []; nchk = 0
    for m in range(mmin, mmax + 1):
        lhs = PY(p, Ym, m, prec)
        rhs = [LS.zero(prec) for _ in range(5)]
        if -m in Ys:
            a, y = Ys[-m]; rhs[a] = y.scale(c0).trunc(prec)
        for a in range(5):
            # truncation: terms with n > nmax are missing; they have u-valuation >= nmax - nmax//5 - 1 >= prec, asserted below
            d = lhs[a].sub(rhs[a]).trunc(prec)
            if not d.iszero():
                bad.append((m, a, d.v()))
                if m > 0: pos.append(m)
        nchk += 1
    nn = KG + 1 - mmax
    assert nn - nn // 5 >= prec, "z-truncation not harmless: raise KG or lower prec"
    return (not bad), c0, bad, nchk, sorted(set(pos))


def symplectic_defects(p, prec):
    """[(n, c)] where Sum_j (-1)^j p[j][c] p[n-j][(3-c-n) mod 5] != delta_{n0}  (mod u^prec)"""
    bad = []
    for n in range(len(p)):
        for c in range(5):
            cp = (3 - c - n) % 5; acc = LS.zero(prec)
            for j in range(n + 1):
                t = p[j][c].mul(p[n - j][cp]).trunc(prec)
                acc = acc.add(t if j % 2 == 0 else t.neg())
            if n == 0: acc = acc.sub(LS.const(1, prec))
            if not acc.trunc(prec).iszero(): bad.append((n, c, acc.v()))
    return bad
