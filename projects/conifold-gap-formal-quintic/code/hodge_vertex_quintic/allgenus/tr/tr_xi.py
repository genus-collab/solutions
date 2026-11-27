# -*- coding: utf-8 -*-
"""tr_xi.py -- the dictionary TR <-> CohFT on ModelCurve (x = kappa (p^4/4 + 1/p), y = p), DNOPS (2-30)-(2-35) (TR_MODEL.md s.2).

Local coordinate at c_j: x - u_j = -w_j^2/2, w_j = sqrt(Lambda) wt_j, Lambda = -5 kappa, wt_j = c_j tau_j(t) = c_j (p - c_j) + ...
(2-30): Delta_j^(-1/2) = dy/dw_j(0) = 1/(sqrt(Lambda) c_j).  (2-34): xi^j = xit^j/sqrt(Lambda),  xit^j(p) = -1/(c_j (p - c_j)).
Normalisation (derived in TR_MODEL.md s.2): e_j = sqrt(Lambda) et_j (et_j: normalised idempotents of the frozen theory, Deltahat_j^(1/2) = c_j);
the TR CohFT has R(z) = M_0(z/kappa) (VERIFIED), hence
    omega_{g,n} = Lambda^(g-1) Sum kappa^-(3g-3+n-Sum d) <tau_d1(et_i1) .. tau_dn(et_in)>^frozen_g prod d xit^{i_l}_{d_l}(p_l), xit^i_d = (d/dx)^d xit^i,
and in the flat basis H^a = N_a Sum_i zeta^(i(a+1)) et_i (N = (1,-1,-1,1,1)) the same with <prod tau_d(H^a)> and d xif_{a,d},
    xif_{a,0} = (N_{3-a}/5) Sum_i zeta^(i(4-a)) xit^i = -N_{3-a} p^(a-1)/(p^5 - 1)  (a = 1..4),   xif_{0,0} = -p^4/(p^5 - 1).
Functions are stored as {(j, k): c} = Sum c/(p - c_j)^k (k >= 1, vanishing at infinity); differentials in the e-basis {(j, k): c}.
(2-31): Rinv(z)[j][i] = delta_ij - Sum_{n>=1} (2n-1)!! (z/Lambda)^n [wt_j^(2n-1)] xit^i;
(2-32): Sum_k Rinv[j][k]/c_k = Sum_n (2n+1)!! (z/Lambda)^n [wt_j^(2n+1)] (p - c_j).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from tr_core import smul, sinv, binom

NSIGN = (1, -1, -1, 1, 1)


def dfact(n):
    """n!! for odd n >= -1"""
    r = 1
    while n > 1: r *= n; n -= 2
    return r


class XiBasis(object):
    def __init__(self, curve, dmax):
        self.C = curve; F = curve.F; self.F = F; self.dmax = dmax; c = curve.c
        N = 2 * dmax + 4
        # q[j][n] = [t^(n-1)] of 1/x'(c_j + t), n = 0..N  (simple pole)
        self.q = []
        for j in range(5):
            l = curve.x_taylor(c[j], N + 3)
            dl = [l[k + 1] * (k + 1) for k in range(N + 2)]          # x' = dl[1] t + dl[2] t^2 + ..
            self.q.append(sinv(dl[1:N + 2], N, F))
        self.fun = {}
        for i in range(5):
            f = {(i, 1): -(F.one / c[i])}
            self.fun[('c', i, 0)] = f
            for d in range(1, dmax + 1):
                f = self.ddx(f); self.fun[('c', i, d)] = f
        for a in range(5):
            for d in range(dmax + 1):
                tot = {}
                for i in range(5):
                    w = F.z(i * (4 - a)) * F.of(Fr(NSIGN[(3 - a) % 5], 5))
                    for key, v in self.fun[('c', i, d)].items():
                        tot[key] = tot[key] + v * w if key in tot else v * w
                self.fun[('f', a, d)] = {k: v for k, v in tot.items() if not v.iszero()}
        # d(c/(p-c_j)^k) = -k c dp/(p-c_j)^(k+1) = -k c e_{j,k}
        self.vec = {key: {(j, k): v * (-k) for (j, k), v in f.items()} for key, f in self.fun.items()}

    def ddx(self, f):
        """(1/x') d/dp of Sum c/(p - c_j)^k: again of this form (1/x' = -p^2/(p^5-1) is O(p^-3) at infinity, regular at 0)"""
        F = self.F; cc = self.C.c; out = {}
        for j in range(5):
            pole = {}; r0 = F.zero
            for (jj, k), v in f.items():
                if jj == j: pole[k + 1] = pole[k + 1] + v * (-k) if (k + 1) in pole else v * (-k)
                else:
                    dinv = F.one / (cc[j] - cc[jj]); w = F.one
                    for _ in range(k + 1): w = w * dinv
                    r0 = r0 + v * (-k) * w
            qj = self.q[j]
            for m, v in pole.items():                                 # (v t^-m) * Sum_n q[n] t^(n-1): coefficient of t^-kk, kk >= 1
                for n in range(0, m + 1):
                    kk = m + 1 - n
                    out[(j, kk)] = out[(j, kk)] + v * qj[n] if (j, kk) in out else v * qj[n]
            if not r0.iszero():
                out[(j, 1)] = out[(j, 1)] + r0 * qj[0] if (j, 1) in out else r0 * qj[0]
        return {k: v for k, v in out.items() if not v.iszero()}

    def convert_vec(self, V, kind='c'):
        """V: {(j,k): c} in the e-basis -> ({(i,d): coefficient of d xit^i_d (kind 'c') or d xif_{a,d} (kind 'f')}, remainder)"""
        F = self.F; V = dict(V); out = {}
        if kind == 'c':
            for d in range(self.dmax, -1, -1):
                for i in range(5):
                    cf = V.get((i, 2 * d + 1))
                    if cf is None or cf.iszero(): continue
                    X = self.vec[('c', i, d)]; lead = X[(i, 2 * d + 1)]; w = cf / lead
                    out[(i, d)] = w
                    for key, x in X.items():
                        V[key] = V[key] - w * x if key in V else -(w * x)
        else:
            co, rem = self.convert_vec(V, 'c')
            for (i, d), v in co.items():                                 # xit^i = Sum_a N_a zeta^(i(a+1)) xif_a
                for a in range(5):
                    w = v * F.z(i * (a + 1)) * NSIGN[a]
                    out[(a, d)] = out[(a, d)] + w if (a, d) in out else w
            return {k: x for k, x in out.items() if not x.iszero()}, rem
        return out, {k: x for k, x in V.items() if not x.iszero()}

    def convert(self, W, n, kind='c'):
        """tensor in the e-basis -> tensor in the d xi basis (slot by slot); returns (tensor, number of non-zero remainders)"""
        cur = dict(W); nbad = 0
        for slot in range(n):
            groups = {}
            for k, v in cur.items():
                groups.setdefault(k[:slot] + k[slot + 1:], {})[k[slot]] = v
            new = {}
            for rest, V in groups.items():
                co, rem = self.convert_vec(V, kind)
                if rem: nbad += 1
                for key, v in co.items(): new[rest[:slot] + (key,) + rest[slot:]] = v
            cur = new
        return {k: v for k, v in cur.items() if not v.iszero()}, nbad


def wt_series(curve, j, M):
    """t = p - c_j as a series in wt_j (x - u_j = -(5/2) wt_j^2, wt_j = c_j tau_j):  t[m] = trev[m] c_j^-m"""
    F = curve.F; L = curve.local(j); ci = F.one / curve.c[j]; out = []; w = F.one
    for m in range(M + 1):
        out.append(L.trev[m] * w); w = w * ci
    return out


def xi_local(curve, i, j, M):
    """coefficients of xit^i = -1/(c_i (p - c_i)) in powers wt_j^m, m = -1..M (list index m + 1)"""
    F = curve.F; c = curve.c; ts = wt_series(curve, j, M + 2); cii = F.one / c[i]
    if i == j:
        u = sinv(ts[1:M + 3], M + 1, F)                               # 1/t = wt^-1 u(wt)
        return [-(cii * x) for x in u[:M + 2]]
    ser = sinv([c[j] - c[i] + ts[0]] + ts[1:M + 1], M, F)             # 1/(c_j - c_i + t)
    return [F.zero] + [-(cii * x) for x in ser]


def Lam(curve):
    """Lambda = -5 kappa: x - u_j = -w_j^2/2 with w_j^2 = Lambda wt_j^2 (x - u_j = (5/2) kappa wt_j^2)"""
    return -5 * curve.kappa


def laplace_Rinv(curve, nmax):
    """Rinv[n][j][i] = [z^n] of (2-31):  delta_ij - Sum_n (2n-1)!! (z/Lambda)^n [wt_j^(2n-1)] xit^i"""
    F = curve.F; X = [[xi_local(curve, i, j, 2 * nmax) for i in range(5)] for j in range(5)]; out = []; L = Lam(curve)
    for n in range(nmax + 1):
        if n == 0:
            out.append([[F.one if i == j else F.zero for i in range(5)] for j in range(5)]); continue
        w = F.of(Fr(-dfact(2 * n - 1)) / L ** n)
        out.append([[X[j][i][2 * n] * w for i in range(5)] for j in range(5)])
    return out


def laplace_unit(curve, nmax):
    """U[n][j] = [z^n] Sum_n (2n+1)!! (z/Lambda)^n [wt_j^(2n+1)] (p - c_j)   ((2-32) right-hand side times sqrt(Lambda))"""
    F = curve.F; out = []; L = Lam(curve)
    ts = [wt_series(curve, j, 2 * nmax + 2) for j in range(5)]
    for n in range(nmax + 1):
        w = F.of(Fr(dfact(2 * n + 1)) / L ** n)
        out.append([ts[j][2 * n + 1] * w for j in range(5)])
    return out


def pairing(curve, V, which):
    """Sum_j Res_{c_j} of y (which='y') or Phi (which='Phi', dPhi = y dx) against the differential V = {(j,k): c} (e-basis)"""
    F = curve.F; tot = F.zero
    for (j, k), v in V.items():
        if which == 'y':
            if k == 1: tot = tot + v
        else:
            tot = tot + v * curve.ydx_taylor(j, k)[k - 1] * F.of(Fr(1, k))
    return tot
