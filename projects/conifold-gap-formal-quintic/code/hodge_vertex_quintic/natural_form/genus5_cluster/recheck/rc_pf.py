# -*- coding: utf-8 -*-
"""rc_pf.py -- independent Frobenius solution of the quintic Picard-Fuchs equation at the conifold point
u = 1 - 5^5 q = 0, the conifold flat coordinate t_c(u), and the change of variable u^-(2g-2) G(u) -> t_c.

Written from the mathematical description only (the operator, the exponents, the normalisations), NOT by calling
fq_conifold / fq_hat / fq3_msum.  Exact over Fraction throughout.

  D = q d/dq = -(1-u) d/du
  P = u D^4 - 2(1-u) D^3 - (7/5)(1-u) D^2 - (2/5)(1-u) D - (24/625)(1-u)
  indicial polynomial  n (n-1)^2 (n-2)  at u^(n-3);  exponents {0,1,1,2}
  analytic solutions f_k = u^k + O(u^3), k = 0,1,2
  log solution g log u + h with g = f_1 + gamma f_2, gamma fixed by the vanishing of the u^-1 obstruction
  varpi0 = f_0 + alpha g + beta f_2,  varpi1 = g,  t_c = varpi1 / varpi0
"""
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from rc_lser import LS, LO, compose_u0, reverse_u1

NW = 48                       # window length: exponents LO .. LO + NW - 1


def _ind(n):
    return Fr(n) * Fr(n - 1) ** 2 * Fr(n - 2)


class PF(object):
    def __init__(self, nw=NW):
        self.nw = nw
        self.one_m = LS.from_u0([Fr(1), Fr(-1)], nw)         # 1 - u
        self.p = {4: LS.mono(1, 1, nw),
                  3: self.one_m.scale(-2),
                  2: self.one_m.scale(Fr(-7, 5)),
                  1: self.one_m.scale(Fr(-2, 5)),
                  0: self.one_m.scale(Fr(-24, 625))}

    def D(self, f):
        """-(1-u) f'"""
        return self.one_m.mul(f.deriv()).neg()

    def Lop(self, f):
        j = [f]
        for _ in range(4):
            j.append(self.D(j[-1]))
        acc = LS.zero(self.nw)
        for k in range(5):
            acc = acc.add(j[k].mul(self.p[k]))
        return acc

    # ---- the log-companion operator: P(A log u + B) = log u P(A) + Clog(A) + P(B)
    def Clog(self, A):
        """the non-log part produced when the derivatives hit log u (B = 0)"""
        AL = A            # the log coefficient of the current jet
        B = LS.zero(self.nw)
        jets = [(AL, B)]
        for _ in range(4):
            a, b = jets[-1]
            na = self.D(a)
            nb = self.D(b).sub(a.mul(self.one_m).shift(-1))
            jets.append((na, nb))
        acc = LS.zero(self.nw)
        for k in range(5):
            acc = acc.add(jets[k][1].mul(self.p[k]))
        return acc

    def analytic(self, k):
        """f_k = u^k + O(u^3) with zero u^0,u^1,u^2 coefficients except the k-th"""
        nmax = self.nw + LO
        b = [Fr(0)] * nmax
        b[k] = Fr(1)
        for n in range(3, nmax):
            f = LS.from_u0(b[:n], self.nw)
            r = self.Lop(f).coef(n - 3)
            b[n] = -r / _ind(n)
        f = LS.from_u0(b, self.nw)
        chk = self.Lop(f)
        bad = [e for e in range(-3, nmax - 6) if chk.coef(e) != 0]
        assert not bad, ("analytic solution %d fails at exponents %s" % (k, bad[:5]))
        return f

    def log_obstruction(self, A):
        """(u^-3, u^-2, u^-1) coefficients of Clog(A) that an analytic h cannot cancel"""
        K = self.Clog(A)
        return [K.coef(e) for e in (-3, -2, -1)], K

    def solve_h(self, A, K):
        nmax = self.nw + LO
        b = [Fr(0)] * nmax
        for n in range(3, nmax):
            f = LS.from_u0(b[:n], self.nw)
            r = self.Lop(f).coef(n - 3) + K.coef(n - 3)
            b[n] = -r / _ind(n)
        return LS.from_u0(b, self.nw)


_PFC = {}


def tc_coeffs(alpha=0, beta=0, kmax=16, nw=NW):
    """coefficients of t_c = u + ... from u^1 up to u^kmax (list, index 0 = [u^1]); plus a diagnostics dict"""
    key = (nw,)
    if key not in _PFC:
        pf = PF(nw)
        f = [pf.analytic(k) for k in range(3)]
        ob1, _ = pf.log_obstruction(f[1])
        ob2, _ = pf.log_obstruction(f[2])
        assert ob1[0] == ob1[1] == 0 and ob2[0] == ob2[1] == 0, (ob1, ob2)
        gamma = -ob1[2] / ob2[2]
        g = f[1].add(f[2].scale(gamma))
        obg, Kg = pf.log_obstruction(g)
        assert all(x == 0 for x in obg), obg
        h = pf.solve_h(g, Kg)
        _PFC[key] = (pf, f, gamma, g, h)
    pf, f, gamma, g, h = _PFC[key]
    varpi0 = f[0].add(g.scale(Fr(alpha))).add(f[2].scale(Fr(beta)))
    tc = g.mul(varpi0.inv())
    assert tc.val() == 1 and tc.coef(1) == 1, (tc.val(), tc.coef(1))
    diag = dict(gamma=gamma, f0=f[0].u0list(5), f1=f[1].u0list(6), f2=f[2].u0list(6),
                g=g.u0list(6), h=h.u0list(6), varpi0=varpi0.u0list(6), tc=[tc.coef(e) for e in range(1, 6)])
    return [tc.coef(e) for e in range(1, kmax + 1)], diag


def polar_in_tc(acoefs, alpha=0, beta=0, g=5):
    """u^-(2g-2) sum_k a_k u^k  ->  coefficients of t_c^-(2g-2) .. t_c^-1 (list of 2g-2 entries).

    u^-m G(u) = t^-m A(t)^-m G(u(t)) with u(t) the reversion of t_c and A(t) = u(t)/t, A(0) = 1, so the polar
    coefficients are exactly [t^0 .. t^(m-1)] of A^-m (G o u).  Only a_0..a_(m-1) are needed.
    """
    m = 2 * (g - 1)
    a = [Fr(x) for x in acoefs][:m]
    assert len(a) == m, "need %d u-coefficients" % m
    k = m + 2
    tcc, _ = tc_coeffs(alpha, beta, kmax=k + 2)
    w = reverse_u1(LS.from_pairs(list(enumerate(tcc, start=1)), NW), k)   # u(t) = t + ...
    nw = NW
    ut = LS.from_pairs(list(enumerate(w, start=1)), nw)
    A = ut.shift(-1)                      # A(t) = u(t)/t, A(0) = 1
    Ginv = compose_u0(a, ut, nw)
    T = A.powi(m).inv().mul(Ginv)
    return [T.coef(e) for e in range(0, m)]


def u_over_tc_pow(m, alpha=0, beta=0, kmax=8):
    """(u/t_c(u))^m as a u-series: coefficients u^0 .. u^(kmax-1)"""
    tcc, _ = tc_coeffs(alpha, beta, kmax=kmax + 2)
    tc = LS.from_pairs(list(enumerate(tcc, start=1)), NW)
    r = tc.shift(-1)                      # t_c/u, constant term 1
    return r.inv().powi(m).u0list(kmax)
