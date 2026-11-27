# -*- coding: utf-8 -*-
"""bv_sys.py -- Stage 2, falsifier (c): self-dual cyclic one-pole systems with the product B = prod b_a VARIED.

System (rows Y_a, a in Z/5, th = z d/du):  th Y_a = b_a Y_{a+1},  b = (f, B/(f^2 g^2), f, g, g)   (self-dual shape),
    f = t_c'(u) F(u),  g = -(varpi_0(u)/(1-u)) Gm(u)  (conifold frame ab),  B = -m^20 Pn(u)/Pd(u) / (u (1-u)^5).
Pn = Pd = 1, m = 1 is the formal quintic family of univ/ (un_gauge.system).  Here B is perturbed, so the model pull-back is
RECOMPUTED from B:  the model z d/dv Y_b = bmod_b(v) Y_{b+1}, bmod = (1, -1/v, 1, -1, -1) pulled back by v(u) has prod bt = -v'^5/v,
so prod bt = B  <=>  v'^5 = -B v  <=>  (with v(0) = 0, v'(0) = m^5)
    v^(4/5) = (4/5) Int_0^u s^(-1/5) beta(s)^(1/5) ds,   beta = -u B = m^20 (Pn/Pd) (1-u)^-5,
    tau := (u/v)^(4/5),   1/tau = Sum_n 4 beta5_n/(5n+4) u^n,   beta5 = beta^(1/5) = m^4 (1-u)^-1 (Pn/Pd)^(1/5),   v = u tau^(-5/4).
(Quintic: 1/tau = 2F1(1,4/5;9/5;u).)  tau(0) = m^-4, so tau^((1-c_b)/4)(0) = m^-(1-c_b) is rational.
The Theorem-U gauge is computed by a COPY of natural_form/theory/ur_gauge.gauge (whose only change is that the hard-wired check
prod b = -1/(u(1-u)^5) is replaced by prod b = prod bt); ur_gauge.solve_ode, prod5, mu are imported read-only.
API: system(NU, sysd) -> dict;  gauge(sysd, K, NU) -> res (ur_gauge format + f, T, m);  hatted(res);  rows(res, Rh, R, KZ, K);
     em_table(base, tau_list, R) (E_M = tau^(a+b+1) Eh_0 with the NEW tau);  SYSTEMS (the cases of BVAR.md).
"""
import os, sys, copy
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); AG = os.path.normpath(os.path.join(HERE, "..", ".."))
PROJ = os.path.normpath(os.path.join(AG, ".."))
for d in (PROJ, os.path.join(PROJ, "natural_form", "theory"), os.path.join(AG, "univ")):
    if d not in sys.path: sys.path.insert(0, d)
from fq_fp import FpRing, QRing, P2
from fq_series import Ser
from fq_conifold import Conifold
from ur_gauge import solve_ode, prod5, mu, CEXP
Q = QRing(); F = Fr

# name: (F, Gm, frame, Pn, Pd, m, description)
SYSTEMS = {
    "R0": ((1,), (1,), ("0", "0"), (1,), (1,), 1, "regression: formal quintic, B unperturbed"),
    "S1": ((1,), (1,), ("0", "0"), (1, F(1, 3), F(-2, 7)), (1,), 1, "quintic f, g;  B -> B (1 + u/3 - 2u^2/7)"),
    "S2": ((1, F(1, 3), F(-2, 7), F(5, 11)), (1, F(-3, 5), F(1, 2), F(2, 3)), ("0", "0"), (1, F(-5, 4), F(7, 3)), (1,), 1,
           "f (1 + u/3 - 2u^2/7 + 5u^3/11), g (1 - 3u/5 + u^2/2 + 2u^3/3);  B -> B (1 - 5u/4 + 7u^2/3)"),
    "S3": ((F(3, 2), F(1, 3)), (1, F(-3, 5), F(1, 2)), ("1/3", "-2/7"), (1, -2, 1), (1, 4, 4), 1,
           "frame (1/3,-2/7), f (3/2 + u/3), g (1 - 3u/5 + u^2/2);  B -> B (1-u)^2/(1+2u)^2  (pole at u=1 moved)"),
    "S4": ((1, F(-7, 4), 3), (1, F(11, 6)), ("0", "0"), (1, 3, 5, F(-4, 3)), (1, F(1, 2)), 1,
           "f (1 - 7u/4 + 3u^2), g (1 + 11u/6);  B -> B (1 + 3u + 5u^2 - 4u^3/3)/(1 + u/2)  (large perturbation)"),
    "S5": ((1,), (1,), ("0", "0"), (1, F(1, 3), F(-2, 7)), (1,), 2, "S1 with the residue of B multiplied by 2^20 (v'(0) = 32)"),
}


def ser_ratio(Pn, Pd, N):
    return Ser.from_list(Q, [Fr(x) for x in Pn], 0, N).mul(Ser.from_list(Q, [Fr(x) for x in Pd], 0, N).inv())


def system(NU, sysd):
    Fp, Gp, ab, Pn, Pd, m = sysd[:6]; N = NU
    C = Conifold(FpRing(P2), Nu=NU, alpha=Fr(ab[0]), beta=Fr(ab[1]))
    var0 = C.varpi0.truncate(N + 2); tc = C.tc.truncate(N + 2)
    one_m = Ser.from_list(Q, [1, -1], 0, N + 2)
    Fs = Ser.from_list(Q, [Fr(x) for x in Fp], 0, N + 2); Gs = Ser.from_list(Q, [Fr(x) for x in Gp], 0, N + 2)
    f = tc.deriv().mul(Fs); g = var0.mul(one_m.inv()).neg().mul(Gs)
    P = ser_ratio(Pn, Pd, N + 2); assert P.coef(0) == 1
    c0 = Fr(m) ** 20
    B = Ser.mono(Q, 1, 1, N + 2).mul(one_m.pow(5)).inv().neg().mul(P).scale(c0)
    b1 = B.mul(f.mul(f).mul(g).mul(g).inv())
    bq = [f, b1, f, g, g]
    beta5 = one_m.inv().mul(P.root(5)).scale(Fr(m) ** 4)              # beta^(1/5)
    tinv = Ser(Q, 0, [4 * beta5.coef(n) / (5 * n + 4) for n in range(N + 2)])
    tau = tinv.inv(); tau0 = tau.coef(0); assert tau0 == Fr(1, m ** 4)
    v = tau.scale(1 / tau0).log().scale(Fr(-5, 4)).exp().scale(Fr(m) ** 5).shift(1)   # u tau^(-5/4), tau0^(-5/4) = m^5
    vp = v.deriv()
    bt = [vp, vp.mul(v.inv()).neg(), vp, vp.neg(), vp.neg()]
    return dict(bq=bq, bt=bt, tau=tau, v=v, C=C, f=f, B=B, m=m)


def gauge(sysd, K, NU, log=None):
    """copy of ur_gauge.gauge with general coefficients (see module docstring)"""
    S = system(NU, sysd); bq, bt = S['bq'], S['bt']
    d = prod5(bq).sub(prod5(bt)); assert not d.c, "prod b != prod bt: v(u) not consistent with B"
    d = prod5(bq).sub(S['B']); assert not d.c, "prod b != B"

    def particular(k, r):
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
        q, Sd = particular(k, r)
        assert not Sd.c, "solvability of order %d violated: %s" % (k, Sd.fmt(Sd.v + 2))
        n = kernel(k)
        A = particular(k + 1, n)[1]
        Bb = particular(k + 1, [x.deriv() for x in n])[1]
        Cc = particular(k + 1, [x.deriv() for x in q])[1]
        kap = Bb.mul(A.inv()); gg = Cc.mul(A.inv()).neg()
        rho, km1, res = solve_ode(kap, gg, k == 0)
        pk = [q[b].add(rho.mul(n[b])) for b in range(5)]
        p.append(pk)
        vals = [x.v if x.c else None for x in pk]
        info.append(dict(k=k, residue=km1, resonances=res, vals=vals, prec=min(x.prec for x in pk)))
        if log: log("   k = %2d: residue of kappa = %s, resonances %s, valuations %s, O(u^%d)" % (k, km1, res, vals, info[-1]['prec']))
        r = [x.deriv() for x in pk]
    f = S['f']
    out = dict(p=p, tau=S['tau'], v=S['v'], info=info, bq=bq, bt=bt, C=S['C'], f=f, B=S['B'], m=S['m'])
    out['T'] = Ser(Q, 1, [f.coef(n) / (n + 1) for n in range(f.prec)])          # T = Int_0^u f
    return out


def hatted(res):
    """{(k, b): Ser}: entry (b+k, b) of [zhat^k] Ghat = tau^((1-c_b)/4) z^k-coefficient, shifted by u^mu(b,k) (tau(0) = m^-4 allowed)"""
    tau = res['tau']; m = res['m']; tau0 = tau.coef(0); lt = tau.scale(1 / tau0).log(); out = {}
    for k, pk in enumerate(res['p']):
        for b in range(5):
            e = 1 - CEXP[b]
            fct = lt.scale(Fr(e, 4)).exp().scale(Fr(1, m) ** e)                 # tau0^(e/4) = m^-e
            out[(k, b)] = fct.mul(pk[b]).shift(mu(b, k))
    return out


def rows(res, Rh, R, KZ, K):
    """Qh[a][i][n] = Sum_k Ghat_(k,b)(u) tau^(n-k) Rh[b][i][n-k]  (un_gauge.rows_from_gauge with the general hatted())"""
    p = R.p; hat = hatted(res); tau = res['tau']
    def red(s):
        assert s.prec >= K, "gauge precision O(u^%d) < %d" % (s.prec, K)
        assert (not s.c) or s.v >= 0, "pole in a hatted gauge entry"
        return [R.of(s.coef(n)) for n in range(K)]
    def mul(a, b):
        out = [0] * K
        for i in range(K):
            if a[i]:
                for j in range(K - i): out[i + j] += a[i] * b[j]
        return [x % p for x in out]
    tp = [[1] + [0] * (K - 1)]; t1 = red(tau)
    for mm in range(1, KZ + 1): tp.append(mul(tp[-1], t1))
    H = {key: red(s) for key, s in hat.items() if key[0] <= KZ}
    Qh = [[[[0] * K for n in range(KZ + 1)] for i in range(5)] for a in range(5)]
    for (k, b), h in H.items():
        if not any(h): continue
        a = (b + k) % 5
        for n in range(k, KZ + 1):
            ht = mul(h, tp[n - k])
            for i in range(5):
                c = Rh[b][i][n - k]
                if c:
                    q = Qh[a][i][n]
                    for j in range(K): q[j] = (q[j] + ht[j] * c) % p
    return Qh


def tau_list(res, R, K):
    return [R.of(res['tau'].coef(n)) for n in range(K)]


def em_table(base, tl):
    """E_M(a, b) = tau^(a+b+1) Eh_0(a, b) with the new tau (list of K residues); Eh_0 = u^0 part of the quintic base E_M (tau_q(0) = 1)"""
    p = base.p; K = base.K
    tp = [[1] + [0] * (K - 1)]
    for k in range(1, base.KZ + 1): tp.append(base.mul(tp[-1], tl))
    out = {}
    for (a, b), tab in base.EM.items():
        out[(a, b)] = [[(base.scale(tp[a + b + 1], x[0]) if (x is not None and x[0] % p) else None) for x in row] for row in tab]
        for row in tab:
            for x in row: assert x is None or x[0] % p, "E_M entry with zero frozen part"
    return out


def base_with_tau(base, tl):
    """shallow copy of an fq4_split.Split base whose E_M uses the new tau (nothing of base is modified)"""
    b2 = copy.copy(base); b2.EM = em_table(base, tl); b2._vt = {}
    return b2
