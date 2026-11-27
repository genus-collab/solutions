# -*- coding: utf-8 -*-
"""rf_gauge.py -- REFEREE (univ): the regular gauge P(u, z) between a cyclic system and the model, computed by DIRECT
undetermined coefficients in Laurent series (not the Theorem-U projection onto K[y] used by natural_form/theory/ur_gauge.py).

Gauge equation  z P' = A P - P At,  A = Sum_a b_a E_{a,a+1} (system), At = Sum_a bt_a E_{a,a+1} (model at v(u)); P maps model solutions
to system solutions.  Graded: P_k = Sum_c p[k][c] E_{c+k, c}.  Order z^k (k >= 1):
     L_k(p_k)_c := b_{c+k-1} p[k][c] - bt_{c-1} p[k][c-1] = p[k-1][c]'          (L_0(p_0) = 0 at order z^0).
L_k is singular (prod b = prod bt): kernel w^(k), w_0 = 1, w_c = w_{c-1} bt_{c-1}/b_{c+k-1}; cokernel lam^(k), lam_0 = 1,
lam_{c+1} = lam_c b_{c+k-1}/bt_c.  Write p[k] = q[k] + kappa_k w^(k) (q[k] the particular solution with q[k][0] = 0); kappa_k is fixed
by the solvability of L_{k+1}: kappa' A + kappa B + C = 0 (A = lam.w, B = lam.w', C = lam.q'), a regular-singular ODE solved in
Laurent series; at a resonance one free constant (k = 0: overall scale, fixed by symplecticity; k = 5, 10: the freedom P c(z^5)).
Normalisations at k = 5, 10, ..:  'stored'  p[k][0](0) = 0 (convention of ur_gauge);  'canon'  p[k][1](0) = 0 (P_11(0, z) constant).
"""
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf_series import LS


def solve_ode(rho, h):
    """Laurent solutions of kappa' + rho kappa = h: (particular with the resonant coefficient 0, homogeneous or None, resonance index)"""
    assert rho.v() >= -1, "irregular singular ODE"
    e = rho.coef(-1) if rho.prec > -1 else Fr(0)
    res = int(-e) if (e.denominator == 1) else None
    starts = []
    if not h.iszero(): starts.append(h.v() + 1)
    if res is not None: starts.append(res)
    if not starts: return LS.zero(min(h.prec + 1, rho.prec + 1)), None, None
    n0 = min(starts)
    top = min(h.prec, rho.prec + n0)                                   # last index computable (inclusive)
    def run(hh, free):
        kap = {}
        for n in range(n0, top + 1):
            m = n - 1
            rhs = (hh.coef(m) if m >= hh.o else Fr(0)) - sum(rho.coef(i) * kap.get(m - i, 0) for i in range(0, m - n0 + 1))
            piv = n + e
            if piv == 0:
                assert rhs == 0, "resonance inconsistent at u^%d" % n
                kap[n] = Fr(free)
            else:
                kap[n] = rhs / piv
        return LS(n0, [kap[n] for n in range(n0, top + 1)], top + 1)
    part = run(h, 0)
    hom = run(LS.zero(h.prec), 1) if (res is not None and n0 <= res <= top) else None
    return part, hom, res


def kernel(k, b, bt):
    w = [LS.const(1, min(x.prec for x in b + bt))]
    for c in range(1, 5): w.append(w[-1].mul(bt[(c - 1) % 5]).mul(b[(c + k - 1) % 5].inv()))
    return w


def cokernel(k, b, bt):
    lam = [LS.const(1, min(x.prec for x in b + bt))]
    for c in range(0, 4): lam.append(lam[-1].mul(b[(c + k - 1) % 5]).mul(bt[c % 5].inv()))
    return lam


def particular(k, r, b, bt):
    """q with q_0 = 0 and L_k(q)_c = r_c for c = 1..4; returns (q, residual of the c = 0 equation)"""
    q = [LS.zero(r[0].prec + 4)]
    for c in range(1, 5):
        q.append(r[c].add(bt[(c - 1) % 5].mul(q[c - 1])).mul(b[(c + k - 1) % 5].inv()))
    resid = bt[4].mul(q[4]).neg().sub(r[0])                            # L_k(q)_0 - r_0 = b_{k-1} q_0 - bt_4 q_4 - r_0
    return q, resid


def isqrt_fr(x):
    import math
    n, d = x.numerator, x.denominator
    rn, rd = math.isqrt(n), math.isqrt(d)
    assert rn * rn == n and rd * rd == d, "not a rational square: %s" % x
    return Fr(rn, rd)


def gauge(b, bt, KG, norm="canon", log=None):
    """p[k][c] (LS) for k = 0..KG; also returns info (resonances, residuals, precision)"""
    info = dict(res={}, resid_v={}, prec=[])
    w = kernel(0, b, bt); q = [LS.zero(w[0].prec) for _ in range(5)]
    p = []
    for k in range(1, KG + 2):
        lam = cokernel(k, b, bt)
        A = sum_ls([lam[c].mul(w[c]) for c in range(5)])
        Bc = sum_ls([lam[c].mul(w[c].deriv()) for c in range(5)])
        C = sum_ls([lam[c].mul(q[c].deriv()) for c in range(5)])
        Ai = A.inv(); rho = Bc.mul(Ai); h = C.mul(Ai).neg()
        kap, hom, res = solve_ode(rho, h); info['res'][k - 1] = res
        km1 = k - 1
        if hom is not None:
            if km1 == 0:                                                 # overall scale: P_0 symplectic, d_0(0) d_3(0) = 1
                c = isqrt_fr(1 / w[3].coef(0)) if kap.iszero() else None
                assert c is not None, "inhomogeneous at order z^0"
                kap = hom.scale(c)
            else:
                assert km1 % 5 == 0, "unexpected resonance at z^%d" % km1
                comp = 0 if norm == "stored" else 1
                a0 = q[comp].add(kap.mul(w[comp])).coef(0); a1 = hom.mul(w[comp]).coef(0)
                assert a1 != 0
                kap = kap.add(hom.scale(-a0 / a1))
        else:
            assert km1 % 5 != 0, "missing resonance at z^%d" % km1
        pk = [q[c].add(kap.mul(w[c])) for c in range(5)]
        p.append(pk); info['prec'].append(min(x.prec for x in pk))
        if k == KG + 1: break
        r = [x.deriv() for x in pk]
        q, resid = particular(k, r, b, bt); info['resid_v'][k] = resid.v()
        assert resid.iszero(), "solvability residual non-zero at z^%d: v = %d, prec %d" % (k, resid.v(), resid.prec)
        w = kernel(k, b, bt)
    return p, info


def sum_ls(lst):
    acc = lst[0]
    for x in lst[1:]: acc = acc.add(x)
    return acc


def check_equation(p, b, bt):
    """max over (k, c) of the valuation of the defect of L_k(p_k) - p_{k-1}' (should be >= its precision) and of L_0(p_0)"""
    bad = []
    for k in range(0, len(p)):
        for c in range(5):
            lhs = b[(c + k - 1) % 5].mul(p[k][c]).sub(bt[(c - 1) % 5].mul(p[k][(c - 1) % 5]))
            rhs = p[k - 1][c].deriv() if k >= 1 else LS.zero(lhs.prec)
            d = lhs.sub(rhs)
            if not d.iszero(): bad.append((k, c, d.v(), d.prec))
    return bad


def regular(p):
    return [(k, c, p[k][c].v()) for k in range(len(p)) for c in range(5) if (not p[k][c].iszero()) and p[k][c].v() < 0]


def P11(p):
    return {k: p[k][1].coef(0) for k in range(0, len(p), 5)}


def times_c(p, cz):
    """P -> P c(z), c = {k: Fraction} (k = 0 mod 5)"""
    out = []
    for k in range(len(p)):
        row = []
        for c in range(5):
            acc = None
            for j, cj in cz.items():
                if j <= k and cj != 0:
                    t = p[k - j][c].scale(cj); acc = t if acc is None else acc.add(t)
            row.append(acc if acc is not None else LS.zero(p[k][c].prec))
        out.append(row)
    return out
