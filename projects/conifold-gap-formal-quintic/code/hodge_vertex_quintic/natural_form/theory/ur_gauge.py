# -*- coding: utf-8 -*-
"""ur_gauge.py -- UNIFORM_REDUCTION.md, check (C): the REGULAR formal gauge between the formal-quintic flat system and the model system
in the natural variable, computed over Q from the periods alone (no R-matrix input).

Systems (rows Y_a, a in Z/5; th = z d/du):
    quintic:  th Y_a = b_a(u) Y_{a+1},   b = (t_c', -1/(u varpi_0^2 (1-u)^3 t_c'^2), t_c', -varpi_0/(1-u), -varpi_0/(1-u)),
              prod b_a = -1/(u (1-u)^5)                [conifold frame (0, 0): varpi_0 = varpi_0^cf, t_c = vanishing period / varpi_0]
    model:    z d/dv Y_b = bmod_b(v) Y_{b+1},  bmod = (1, -1/v, 1, -1, -1);  pulled back by v(u) = u tau(u)^(-5/4):
              bt_b(u) = v'(u) bmod_b(v(u)),   prod bt_b = prod b_a   <=>   tau(u) = 1/2F1(1, 4/5; 9/5; u).
Gauge:  Y^quintic = P(u, z) Y^model(v(u)),  P = Sum_k z^k P_k(u),  P_k = Sum_b p_{k,b}(u) E_{b+k, b}  (Z/5 grading), with
    (E_k, b):   p_{k-1,b}' = b_{b+k-1} p_{k,b} - bt_{b-1} p_{k,b-1} .
Solution: p_k = q_k + rho_k n_k (q_k particular with q_{k,0} = 0, n_k the kernel vector, n_{k,0} = 1); rho_k is fixed by the solvability
of (E_{k+1}): a first-order ODE  rho' + kappa rho = g  with a regular singular point at u = 0, solved in Laurent series.  THEOREM
(UNIFORM_REDUCTION.md s.3): all p_{k,b} are power series in u, unique up to P -> P c(z^5).  Here: checked by computation for k <= K.
Hatted form (s = u^(1/5), zhat = z/s^4, c = (0,-1,3,2,1)):
    Ghat_{ab}(u, zhat) = tau^((1-c_b)/4) Sum_k zhat^k u^(mu(b,k)) p_{k,b}(u) delta_{a,b+k},   mu(b,k) = (4k + c_b - c_{b+k})/5,
    Qhat(u, zhat) = Ghat(u, zhat) R^(tau zhat),  i.e.  R~(zhat; u) = G(u, zhat) M_0(tau zhat) with G = Psi(u)^-1 Ghat Psi(0).
API: gauge(K, NU) -> dict(p=[[Ser]*5]*(K+1), tau=Ser, info=[...]);  hatted(res) -> {(k, b): Ser in u};  first_orders(res).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ)
from fq_fp import FpRing, QRing, P2
from fq_series import Ser
from fq_conifold import Conifold
CEXP = (0, -1, 3, 2, 1)
Q = QRing()


def mu(b, k):
    v = 4 * k + CEXP[b] - CEXP[(b + k) % 5]
    assert v % 5 == 0
    return v // 5


def coefficients(NU):
    """(bq, bt, tau, v) as Ser over Q, precision about O(u^NU)"""
    C = Conifold(FpRing(P2), Nu=NU)
    N = NU
    var0 = C.varpi0.truncate(N + 2); tc = C.tc.truncate(N + 2)
    one_m = Ser.from_list(Q, [1, -1], 0, N + 2)
    tcp = tc.deriv()
    b1 = Ser.mono(Q, 1, 1, N + 2).mul(var0).mul(var0).mul(one_m.pow(3)).mul(tcp).mul(tcp).inv().neg()
    b3 = var0.mul(one_m.inv()).neg()
    bq = [tcp, b1, tcp, b3, b3]
    tau = Ser.from_list(Q, [Fr(4, 5 * n + 4) for n in range(N + 2)], 0, N + 2).inv()
    v = tau.log().scale(Fr(-5, 4)).exp().shift(1)
    vp = v.deriv()
    bt = [vp, vp.mul(v.inv()).neg(), vp, vp.neg(), vp.neg()]
    return bq, bt, tau, v, C


def prod5(lst):
    out = lst[0]
    for x in lst[1:]: out = out.mul(x)
    return out


def solve_ode(kap, g, first):
    """Laurent-series solution of rho' + kap rho = g (kap with at most a simple pole).  At a resonance (m + res kap = 0) the right-hand
    side must vanish; the free constant is set to 1 if first (k = 0, m = 0) and to 0 otherwise.  Returns (rho, residue, resonances)."""
    assert (not kap.c) or kap.v >= -1, "kappa has a pole of order %d" % (-kap.v)
    km1 = kap.coef(-1) if kap.c else Fr(0)
    if g.c: m0 = g.v + 1
    else: m0 = 0
    if first: m0 = min(m0, 0)
    mmax = min(g.prec, kap.prec + m0) - 1
    rho = {}; res = []
    for m in range(m0, mmax + 1):
        rhs = g.coef(m - 1) if g.c else Fr(0)
        for j in range(0, m - m0):
            rhs -= kap.coef(j) * rho[m - 1 - j]
        den = m + km1
        if den == 0:
            assert rhs == 0, "resonance with non-zero right-hand side at u^%d (a logarithm would be needed)" % m
            rho[m] = Fr(1) if (first and m == 0) else Fr(0); res.append(m)
        else:
            rho[m] = rhs / den
    return Ser(Q, m0, [rho[m] for m in range(m0, mmax + 1)]), km1, res


def gauge(K, NU, log=None):
    bq, bt, tau, v, C = coefficients(NU)
    lhs = prod5(bq); rhs = prod5(bt)
    d = lhs.sub(rhs)
    assert not d.c, "prod b != prod bt"
    target = Ser.mono(Q, 1, 1, NU).mul(Ser.from_list(Q, [1, -1], 0, NU).pow(5)).inv().neg()
    assert not lhs.sub(target).c, "prod b != -1/(u (1-u)^5)"

    def particular(k, r):
        """q with q_0 = 0 solving (E_k, b) for b = 1..4, and the defect S = r_0 + bt_4 q_4 of (E_k, 0)"""
        q = [Ser.zero(Q, NU + 2)]
        for b in range(1, 5):
            q.append(r[b].add(bt[b - 1].mul(q[b - 1])).mul(bq[(b + k - 1) % 5].inv()))
        return q, r[0].add(bt[4].mul(q[4]))

    def kernel(k):
        n = [Ser.const(Q, 1, NU + 2)]
        for b in range(1, 5): n.append(n[b - 1].mul(bt[b - 1]).mul(bq[(b + k - 1) % 5].inv()))
        return n

    p = []; info = []
    r = [Ser.zero(Q, NU + 2) for _ in range(5)]
    for k in range(K + 1):
        q, S = particular(k, r)
        assert not S.c, "solvability of order %d violated: %s" % (k, S.fmt(S.v + 2))
        n = kernel(k)
        A = particular(k + 1, n)[1]
        B = particular(k + 1, [x.deriv() for x in n])[1]
        Cc = particular(k + 1, [x.deriv() for x in q])[1]
        kap = B.mul(A.inv()); g = Cc.mul(A.inv()).neg()
        rho, km1, res = solve_ode(kap, g, k == 0)
        pk = [q[b].add(rho.mul(n[b])) for b in range(5)]
        p.append(pk)
        vals = [x.v if x.c else None for x in pk]
        info.append(dict(k=k, residue=km1, resonances=res, vals=vals, prec=min(x.prec for x in pk)))
        if log: log("   k = %2d: residue of kappa = %s, resonances %s, valuations of p_{k,b} (b = 0..4): %s, precision O(u^%d)" % (k, km1, res, vals, info[-1]['prec']))
        r = [x.deriv() for x in pk]
    return dict(p=p, tau=tau, v=v, info=info, bq=bq, bt=bt, C=C)


def hatted(res):
    """{(k, b): Ser in u} = the entry (a, b) = (b + k, b) of [zhat^k] Ghat"""
    tau = res['tau']; lt = tau.log(); out = {}
    for k, pk in enumerate(res['p']):
        for b in range(5):
            f = lt.scale(Fr(1 - CEXP[b], 4)).exp()
            out[(k, b)] = f.mul(pk[b]).shift(mu(b, k))
    return out


def first_orders(res):
    """{k: {b: first u-order n at which the entry (b + k, b) of [zhat^k] Ghat is non-zero}}"""
    out = {}
    for k, pk in enumerate(res['p']):
        out[k] = {b: (mu(b, k) + pk[b].v if pk[b].c else None) for b in range(5)}
    return out
