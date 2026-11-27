# -*- coding: utf-8 -*-
"""ref_hat.py -- REFEREE: hatted rows Qh[a][i][n] (u-series), Teleman form, symbol M_0, tau, and the factor G.

HatRows(p, NZ, K, alpha=0, beta=0, zeta=None): see ref_rows.py for the derivation.  Checks stored in .checks:
   h_auto      ph_n(1) = h_n for all n not divisible by 5 and all colours (WKB + QRR convention)
   psi_orth    Rt(zhat = 0; u) = 1 for all u (orthonormality of the normalised canonical frame, tests C0, C1, kappa)
   sympl       Rt(zhat)^T Rt(-zhat) = 1 mod (zhat^(NZ+1), u^K)
   covar       Qh_{a,n}^{(i)} = zeta^(i(a+1-n)) Qh_{a,n}^{(0)} for all u
   circ        Rt_n[j+1][i+1] = zeta^(-n) Rt_n[j][i]
tau_pred(K) (Fractions): 1 / Sum_n 4 u^n/(5n+4).   M0_of(Rt, NZ) -> (M0 list, M0inv list).
Mtau(M0list, tau, p, NZ, K) -> bivariate {(n, k)} of M_0(tau(u) zhat).   G_of(Rt, M0inv, tau, ...) = Rt . M_0(tau zhat)^(-1).
fourier(X, zeta, p) -> {(m, m'): (1/25) Sum_ij zeta^(m i) X[i][j] zeta^(m' j)} (nonzero entries).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ref_fp import (inv, of, find_zeta, smul, sinv, stheta, sscale, sadd, sone_minus, sderiv, mid, mzero, mmul, madd, mscale, mT,
                    miszero, bclean, bmul, bflip_T)
from ref_rows import periods, wkb, _A

CEXP = (0, -1, 3, 2, 1)


class HatRows(object):
    def __init__(self, p, NZ, K, alpha=0, beta=0, zeta=None):
        self.p = p; self.NZ = NZ; self.K = K; self.zeta = find_zeta(p) if zeta is None else zeta; z = self.zeta
        f0, f1, f2, g, h = periods(p, K + 3); self.f = (f0, f1, f2); self.g = g; self.h = h
        al = of(alpha, p); be = of(beta, p)
        v0 = [(f0[k] + al * g[k] + be * f2[k]) % p for k in range(K + 2)]
        tc = smul(g[:K + 2], sinv(v0, p), p); self.tc = tc
        C1 = [(-x) % p for x in sone_minus(sderiv(tc, p), p)][:K]              # D t_c = -(1-u) t_c'
        C0 = v0[:K]; self.C0 = C0; self.C1 = C1
        C0i = sinv(C0, p); C1i = sinv(C1, p)
        kh = [C1i, smul(smul(C0, C0, p), smul(C1, C1, p), p), C1i, C0i]
        i5 = inv(5, p); self.wkb_info = []; self.ph = []
        sadd_c = lambda a, b, p_, c=1: [(x + c * y) % p_ for x, y in zip(a, b)]
        self.Q = [[None] * 5 for _ in range(5)]
        for i in range(5):
            Phi, info = wkb(p, z, i, NZ); self.wkb_info.append(info); self.ph.append(Phi); zi = pow(z, i, p)
            row = [smul(((P + [0] * K)[:K]), sscale(C0i, zi, p), p) for P in Phi]
            self.Q[0][i] = row
            for a in range(4):
                nxt = _A(CEXP[a] - 1, row, zi, p, i5, stheta, sone_minus, sadd_c, sscale)
                row = [smul(x, kh[a], p) for x in nxt]; self.Q[a + 1][i] = row
        self.Rt = self.teleman()
        self.checks = self._checks()

    def teleman(self):
        p = self.p; i5 = inv(5, p); out = {}
        for n in range(self.NZ + 1):
            mats = [mzero() for _ in range(self.K)]
            for j in range(5):
                for i in range(5):
                    acc = [0] * self.K
                    for a in range(5):
                        b = (3 - a) % 5
                        t = smul(self.Q[b][j][0], self.Q[a][i][n], p)
                        acc = [(x + y) % p for x, y in zip(acc, t)]
                    for k in range(self.K): mats[k][j][i] = acc[k] * i5 % p
            for k in range(self.K): out[(n, k)] = mats[k]
        return bclean(out)

    def _checks(self):
        p = self.p; z = self.zeta; NZ = self.NZ; K = self.K; Rt = self.Rt; chk = {}
        chk['h_auto'] = all(all(info['h_auto']) for info in self.wkb_info)
        chk['psi_orth'] = all((Rt.get((0, k), mzero()) == (mid() if k == 0 else mzero())) for k in range(K))
        S = bmul(bflip_T(Rt, p), Rt, p, NZ, K)                       # Rt(-zhat)^T Rt(zhat)
        chk['sympl'] = (S == {(0, 0): mid()})
        chk['covar'] = all(self.Q[a][i][n] == sscale(self.Q[a][0][n], pow(z, (i * (a + 1 - n)) % 5, p), p)
                           for a in range(5) for i in range(5) for n in range(NZ + 1))
        chk['circ'] = all(A[(j + 1) % 5][(i + 1) % 5] == A[j][i] * pow(z, (-n) % 5, p) % p for (n, k), A in Rt.items()
                          for i in range(5) for j in range(5))
        return chk


def tau_pred(K):
    inv_ = [Fr(4, 5 * k + 4) for k in range(K)]; tau = [Fr(1)]
    for k in range(1, K): tau.append(-sum(inv_[j] * tau[k - j] for j in range(1, k + 1)))
    return tau


def M0_of(Rt, NZ, p):
    M0 = [Rt.get((n, 0), mzero()) for n in range(NZ + 1)]
    assert M0[0] == mid()
    Mi = [mid()]
    for n in range(1, NZ + 1):
        acc = mzero()
        for m in range(1, n + 1): acc = madd(acc, mmul(M0[m], Mi[n - m], p), p)
        Mi.append(mscale(acc, p - 1, p))
    return M0, Mi


def tau_powers(tau, p, NZ, K):
    """tp[n][k] = [u^k] tau^n"""
    tp = [[1] + [0] * (K - 1)]
    for n in range(1, NZ + 1): tp.append(smul(tp[-1], tau, p, K))
    return tp


def Mtau(Mlist, tau, p, NZ, K):
    tp = tau_powers(tau, p, NZ, K)
    return bclean({(n, k): mscale(Mlist[n], tp[n][k], p) for n in range(NZ + 1) for k in range(K)})


def G_of(Rt, M0inv, tau, p, NZ, K):
    return bmul(Rt, Mtau(M0inv, tau, p, NZ, K), p, NZ, K)


def degrees(G, K):
    """zhat-degree of [u^n] G for n = 0..K-1 (None if the coefficient vanishes)"""
    return [max([n for (n, k) in G if k == j] or [None], key=lambda x: -1 if x is None else x) for j in range(K)]


def fourier(X, zeta, p):
    i25 = inv(25, p); out = {}
    for m in range(5):
        for mp in range(5):
            v = sum(pow(zeta, (m * i + mp * j) % 5, p) * X[i][j] for i in range(5) for j in range(5)) % p
            if v: out[(m, mp)] = v * i25 % p
    return out


def fit_tau(Rt, M0, M0inv, p, NZ, K, use_orders=None):
    """order by order in u: tau_m from the TOP zhat-order (entry with the largest support), then the resulting degree of G_m.
    Returns (tau list, list of dicts per m: degree, number of vanishing matrix coefficients above the degree, consistency of all
    entries at the top order)."""
    tau = [1] + [0] * (K - 1); rep = []
    D = [mzero() for _ in range(NZ + 1)]                             # zhat M_0' M_0^-1
    for a in range(1, NZ + 1):
        for b in range(0, NZ + 1 - a): D[a + b] = madd(D[a + b], mmul(M0[a], M0inv[b], p), p, a)
    for m in range(1, K):
        G = G_of(Rt, M0inv, tau, p, NZ, K)
        top = NZ if use_orders is None else use_orders
        G0 = G.get((top, m), mzero()); cand = set(); nz = 0
        for i in range(5):
            for j in range(5):
                if D[top][i][j]: cand.add(G0[i][j] * inv(D[top][i][j], p) % p); nz += 1
                else: cand.add(('zeroD', G0[i][j] == 0))
        vals = [c for c in cand if not isinstance(c, tuple)]
        ok = len(vals) == 1 and all(c[1] for c in cand if isinstance(c, tuple))
        tau[m] = vals[0] if vals else 0
        G = G_of(Rt, M0inv, tau, p, NZ, K)
        deg = max([n for (n, k) in G if k == m] or [0])
        rep.append(dict(m=m, top_consistent=ok, entries_used=nz, degree=deg, vanishing_orders=NZ - deg,
                        Dnonzero_orders=sum(1 for n in range(deg + 1, NZ + 1) if not miszero(D[n]))))
    return tau, rep
