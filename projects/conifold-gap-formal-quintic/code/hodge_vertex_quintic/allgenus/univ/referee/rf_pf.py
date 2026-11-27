# -*- coding: utf-8 -*-
"""rf_pf.py -- REFEREE (univ): own Frobenius solution of the conifold Picard-Fuchs equation of the formal quintic over Q, the
conifold-frame periods (varpi_0, t_c), and the cyclic systems used in E2/E4 (quintic, perturbed, model pulled back by v(u)).

PF operator in u (as documented in fq_conifold.py; only the operator is taken from there, the solution is computed here):
   P = u D^4 - 2(1-u) D^3 - (7/5)(1-u) D^2 - (2/5)(1-u) D - (24/625)(1-u),     D = -(1-u) d/du.
Analytic solutions f_k = u^k + O(u^3) (k = 0, 1, 2); log solution g log u + h, h = O(u^3); g fixed up to scale by the vanishing of
the obstructions at u^-3..u^-1 and normalised g = u + O(u^2).  varpi_0 = f_0 + alpha g + beta f_2,  t_c = g/varpi_0.
Cyclic system (UNIFORM_REDUCTION s.0): th Y_a = b_a Y_{a+1}, b = (f, B/(f^2 g^2), f, g, g), B = -1/(u(1-u)^5),
   quintic: f = t_c', g = -varpi_0/(1-u);  perturbed: f -> f F(u), g -> g Gm(u);  model at v(u): bt = (v', -v'/v, v', -v', -v').
"""
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf_series import LS, tau_series, model_v

COEF = [None, None, None, None, None]


def _ops(N):
    one_m = LS.poly([1, -1], N + 8); u = LS.poly([0, 1], N + 8)
    return one_m, u


def Dop(A, one_m):
    return one_m.mul(A.deriv()).neg()


def apply_P(A, B, N):
    """P(A log u + B) as a pair (A', B')"""
    one_m, u = _ops(N)
    dlog = one_m.mul(LS(-1, [1], N + 8)).neg()                  # D log u = -(1-u)/u
    pows = [(A, B)]
    for _ in range(4):
        a, b = pows[-1]
        pows.append((Dop(a, one_m), Dop(b, one_m).add(a.mul(dlog))))
    cf = [one_m.scale(Fr(-24, 625)), one_m.scale(Fr(-2, 5)), one_m.scale(Fr(-7, 5)), one_m.scale(-2), u]
    oa = LS.zero(N); ob = LS.zero(N)
    for j in range(5):
        oa = oa.add(cf[j].mul(pows[j][0])); ob = ob.add(cf[j].mul(pows[j][1]))
    return oa, ob


def frobenius(N):
    """(f0, f1, f2, g, h) as LS with precision N (exact rationals)"""
    W = N + 6
    # matrix of P on monomials: P(u^n) = Sum_j c[n][j] u^(n-3+j)
    c = {}
    for n in range(0, W):
        a, b = apply_P(LS.zero(W + 8), LS.poly([0] * n + [1], W + 8), W)
        assert a.iszero()
        c[n] = {j: b.coef(n - 3 + j) for j in range(0, 6) if n - 3 + j < b.prec}
    for n in range(3):
        assert c[n].get(0, 0) == 0
    def solve(rhs, seed):
        """coefficients a_n with P(Sum a_n u^n) = rhs (LS), a_0..a_2 = seed"""
        a = list(seed) + [Fr(0)] * (W - 3)
        for n in range(0, W):
            s = rhs.coef(n - 3) - sum(a[m] * c[m].get(n - m, 0) for m in range(max(0, n - 5), n))
            if n < 3:
                s -= a[n] * c[n].get(0, 0)
                assert s == 0, "obstruction at u^%d" % (n - 3)
            else:
                a[n] = s / c[n][0]
        return LS(0, a, W)
    zero = LS.zero(W + 8)
    f = [solve(zero, [1 if k == j else 0 for j in range(3)]) for k in range(3)]
    # log solution: g = x0 f0 + x1 f1 + x2 f2; the B-part R(g) must vanish at u^-3..u^-1
    R = [apply_P(fk, LS.zero(W + 8), W)[1] for fk in f]
    M = [[R[k].coef(m) for k in range(3)] for m in (-3, -2, -1)]
    # kernel of the 3x3 matrix M (expected one-dimensional, containing g = f1 + gamma f2)
    from itertools import product
    best = None
    for x0, x1 in ((0, 1),):
        # solve M (x0, x1, x2) = 0 for x2
        cands = [(-(M[m][0] * x0 + M[m][1] * x1), M[m][2]) for m in range(3)]
        x2 = None
        for num, den in cands:
            if den != 0: x2 = num / den; break
        assert x2 is not None
        assert all(num == den * x2 for num, den in cands), "obstruction not removable with g = f1 + gamma f2"
        best = (x0, x1, x2)
    gamma = best[2]
    g = f[1].add(f[2].scale(gamma))
    Rg = apply_P(g, LS.zero(W + 8), W)[1]
    h = solve(Rg.neg(), [0, 0, 0])
    a, b = apply_P(g, h, W)
    assert a.trunc(N - 4).iszero() and b.trunc(N - 4).iszero(), "log solution check"
    return [x.trunc(N) for x in f] + [g.trunc(N), h.trunc(N)], gamma


def periods(N, alpha=0, beta=0):
    (f0, f1, f2, g, h), gamma = frobenius(N)
    varpi0 = f0.add(g.scale(Fr(alpha))).add(f2.scale(Fr(beta)))
    tc = g.mul(varpi0.inv())
    return dict(varpi0=varpi0, tc=tc, gamma=gamma, f=(f0, f1, f2), g=g, h=h)


def system(N, F=(1,), Gm=(1,), alpha=0, beta=0):
    """(b, bt, tau, v, f) with b the (perturbed) quintic system and bt the model pulled back by v(u); all LS, precision ~N"""
    per = periods(N + 4, alpha, beta)
    one_m = LS.poly([1, -1], N + 4)
    Fs = LS.poly([Fr(x) for x in F], N + 4); Gs = LS.poly([Fr(x) for x in Gm], N + 4)
    f = per['tc'].deriv().mul(Fs).trunc(N + 2)
    g = per['varpi0'].mul(one_m.inv()).neg().mul(Gs).trunc(N + 2)
    B = one_m.mul(one_m).mul(one_m).mul(one_m).mul(one_m).shift(1).inv().neg().trunc(N + 1)
    b1 = B.mul(f.mul(f).mul(g).mul(g).inv())
    b = [f, b1, f, g, g]
    tau = tau_series(N + 2); v = model_v(N + 2); vp = v.deriv()
    bt = [vp, vp.mul(v.inv()).neg(), vp, vp.neg(), vp.neg()]
    return dict(b=b, bt=bt, tau=tau, v=v, f=f, per=per, B=B)


if __name__ == "__main__":
    for ab in ((0, 0), (Fr(1, 3), Fr(-2, 7))):
        per = periods(10, *ab)
        print("frame", ab, "gamma =", per['gamma'])
        print("   t_c    =", per['tc'].fmt(5))
        print("   varpi0 =", per['varpi0'].fmt(5))
