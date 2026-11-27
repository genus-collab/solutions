# -*- coding: utf-8 -*-
"""rf_bsys.py -- referee's independent construction of a B-perturbed self-dual cyclic system, its model pull-back and its Theorem-U gauge.
Does NOT import any bv_* module.  Project modules used read-only: fq_conifold.Conifold (periods t_c, varpi_0 only), fq_fp (FpRing for
the Conifold constructor).  Own series class rf_ser.L.

System  z d/du Y_a = b_a Y_{a+1},  b = (f, B/(f^2 g^2), f, g, g),  f = t_c' F,  g = -(varpi_0/(1-u)) Gm,  B = -m^20 Pn/(Pd u (1-u)^5).
Model   z d/dv Y_b = (1, -1/v, 1, -1, -1)_b Y_{b+1} pulled back by v(u):  bt = v' (1, -1/v, 1, -1, -1).
v(u) = u w(u) is found by a fixed-point recursion for  v'^5 = -B v  (NOT the worker's closed form):  (u w)' = beta5 w^(1/5),
beta5 = (-u B)^(1/5), w(0) = m^5;  coefficient n: (n + 4/5) w_n = [u^n](beta5 w^(1/5))|_{w_n = 0}.   tau := w^(-4/5) = (u/v)^(4/5).
Gauge P = Sum_k z^k P_k, P_k(c+k, c) = q^k_c, solving  (q^{k-1}_c)' = b_{c+k-1} q^k_c - bt_{c-1} q^k_{c-1}  (z P' = N P - P Nt).
Own solver: particular solution started at c = 2, kernel n (n_2 = 1), rho from the solvability of the next order through the LEFT null
vector l (l_{c+1} = l_c b_{c+k}/bt_c), Laurent ODE; free resonance constants 0 for k > 0, 1 at (k, m) = (0, 0).
Normalisations defined intrinsically: 'canon' P_11(0,z) constant, 'stored' P_00(0,z) constant; then sigma = C^(-1/2) > 0 with
P_0(0)^T eta P_0(0) = C eta (q0 q3 = q1 q2 = q4^2 = C).
"""
import os, sys, math
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
if PROJ not in sys.path: sys.path.insert(0, PROJ)
from rf_ser import L
F = Fr
CX = (0, -1, 3, 2, 1)

# name: F, Gm, frame, Pn, Pd, m   (S1-S5 transcribed from BVAR.md s.2; S6, S7 are the referee's own)
SYS = {
    "R0": ((1,), (1,), (0, 0), (1,), (1,), 1),
    "S1": ((1,), (1,), (0, 0), (1, F(1, 3), F(-2, 7)), (1,), 1),
    "S2": ((1, F(1, 3), F(-2, 7), F(5, 11)), (1, F(-3, 5), F(1, 2), F(2, 3)), (0, 0), (1, F(-5, 4), F(7, 3)), (1,), 1),
    "S3": ((F(3, 2), F(1, 3)), (1, F(-3, 5), F(1, 2)), (F(1, 3), F(-2, 7)), (1, -2, 1), (1, 4, 4), 1),
    "S4": ((1, F(-7, 4), 3), (1, F(11, 6)), (0, 0), (1, 3, 5, F(-4, 3)), (1, F(1, 2)), 1),
    "S5": ((1,), (1,), (0, 0), (1, F(1, 3), F(-2, 7)), (1,), 2),
    "S6": ((2, F(-1, 3), F(3, 8)), (1, F(1, 4), F(-1, 6)), (F(-1, 2), F(1, 5)), (1, F(-2, 3), F(5, 7), F(1, 4)), (1, F(1, 5)), 3),
    "S7": ((F(-5, 7), 1), (F(4, 3), F(-2, 9)), (0, F(3, 11)), (1, 0, 0, 0, 0, 6), (1, -1), 1),
}


def mu(b, k):
    e = 4 * k + CX[b % 5] - CX[(b + k) % 5]; assert e % 5 == 0; return e // 5


def qsqrt(x):
    x = Fr(x); a, b = math.isqrt(x.numerator), math.isqrt(x.denominator)
    assert a * a == x.numerator and b * b == x.denominator, "C not a rational square"
    return Fr(a, b)


def build(name, N=34):
    from fq_conifold import Conifold
    from fq_fp import FpRing, P2
    Fp, Gp, ab, Pn, Pd, m = SYS[name]
    C = Conifold(FpRing(P2), Nu=N + 4, alpha=Fr(ab[0]), beta=Fr(ab[1]))
    tc = L(0, [C.tc.coef(k) for k in range(N + 3)], N + 3); var0 = L(0, [C.varpi0.coef(k) for k in range(N + 2)], N + 2)
    om = L.poly([1, -1], N + 2)
    f = tc.der().mul(L.poly(Fp, N + 2)); g = var0.mul(om.inv()).neg().mul(L.poly(Gp, N + 2))
    rat = L.poly(Pn, N + 2).mul(L.poly(Pd, N + 2).inv())
    B = om.mul(om).mul(om).mul(om).mul(om).inv().mul(rat).sc(-Fr(m) ** 20).sh(-1)   # -m^20 Pn/Pd / (u (1-u)^5)
    B = L(B.v, B.c, N + 1)
    b = [f, B.mul(f.mul(f).mul(g).mul(g).inv()), f, g, g]
    # v'^5 = -B v by fixed point
    beta = B.sh(1).neg()                                  # -u B, beta(0) = m^20
    beta5 = beta.powfrac(Fr(1, 5), Fr(m) ** 4)
    w = [Fr(m) ** 5]
    for n in range(1, N):
        W = L(0, w + [0] * (N - n), N); r = beta5.mul(W.powfrac(Fr(1, 5), Fr(m))).co(n)
        w.append(r / (n + Fr(4, 5)))
    W = L(0, w, N); v = W.sh(1); vp = v.der()
    tau = W.powfrac(Fr(-4, 5), Fr(1, m ** 4))
    bt = [vp, vp.mul(v.inv()).neg(), vp, vp.neg(), vp.neg()]
    pb = b[0].mul(b[1]).mul(b[2]).mul(b[3]).mul(b[4]); pt = bt[0].mul(bt[1]).mul(bt[2]).mul(bt[3]).mul(bt[4])
    d1 = pb.sub(B); d2 = pt.sub(B)
    assert d1.zero() and d2.zero(), "prod b / prod bt != B"
    T = f.integ()
    return dict(name=name, b=b, bt=bt, B=B, v=v, w=W, tau=tau, f=f, T=T, m=m, C=C, N=N)


def solve_rho(kap, gg, k):
    """Laurent solution of rho' + kap rho = gg; resonance constant 1 at (k, m) = (0, 0), else 0"""
    assert kap.v >= -1, "kappa pole order > 1"
    km1 = kap.co(-1); P = min(kap.P, gg.P); m0 = min(gg.v + 1, 0) if not gg.zero() else 0
    rho = {}; res = []
    for mm in range(m0, P):
        rhs = gg.co(mm - 1) if mm - 1 >= gg.v and mm - 1 < gg.P else Fr(0)
        rhs -= sum(kap.co(j) * rho[mm - 1 - j] for j in range(0, mm - m0) if (mm - 1 - j) in rho)
        den = mm + km1
        if den == 0:
            assert rhs == 0, "log needed"
            rho[mm] = Fr(1) if (k == 0 and mm == 0) else Fr(0); res.append(mm)
        else: rho[mm] = rhs / den
    return L(m0, [rho[mm] for mm in range(m0, P)], P), km1, res


def gauge(S, K=12):
    b, bt = S['b'], S['bt']; p = []; prev = None; info = []
    for k in range(K + 1):
        r = [prev[c].der() for c in range(5)] if prev else [L(40, [], 40)] * 5
        qh = {2: L(60, [], 60)}; n = {2: L.const(1, S['N'])}
        for c in (3, 4, 0, 1):
            qh[c] = r[c].add(bt[(c - 1) % 5].mul(qh[(c - 1) % 5])).mul(b[(c + k - 1) % 5].inv())
            n[c] = n[(c - 1) % 5].mul(bt[(c - 1) % 5]).mul(b[(c + k - 1) % 5].inv())
        cons = b[(k + 1) % 5].mul(qh[2]).sub(bt[1].mul(qh[1])).sub(r[2])
        assert cons.zero(), "order %d not solvable: %r" % (k, cons)
        l = {2: L.const(1, S['N'])}
        for c in (2, 3, 4, 0):
            l[(c + 1) % 5] = l[c].mul(b[(c + k) % 5]).mul(bt[c].inv())
        al = sum_l([l[c].mul(n[c]) for c in range(5)]); be = sum_l([l[c].mul(n[c].der()) for c in range(5)])
        ga = sum_l([l[c].mul(qh[c].der()) for c in range(5)])
        rho, km1, res = solve_rho(be.mul(al.inv()), ga.mul(al.inv()).neg(), k)
        q = [qh[c].add(rho.mul(n[c])) for c in range(5)]
        info.append((k, km1, res, [x.v for x in q], min(x.P for x in q)))
        p.append(q); prev = q
    return p, info


def sum_l(lst):
    out = lst[0]
    for x in lst[1:]: out = out.add(x)
    return out


def zmul(p, mz):
    """P -> P m(z), m = {k: Fraction} scalar"""
    out = []
    for k in range(len(p)):
        row = []
        for c in range(5):
            acc = None
            for j, mj in mz.items():
                if j <= k and mj != 0:
                    t = p[k - j][c].sc(mj); acc = t if acc is None else acc.add(t)
            row.append(acc if acc is not None else L(60, [], 60))
        out.append(row)
    return out


def normalise(p, mode):
    """returns (p', C, sigma, c0, m5): mode 'canon' (P_11(0,z) const) or 'stored' (P_00(0,z) const), then sigma = C^(-1/2)"""
    K = len(p) - 1; idx = 1 if mode == "canon" else 0
    Pz = {k: p[k][idx].co(0) for k in range(0, K + 1, 5)}
    inv = {0: 1 / Pz[0]}
    for k in range(5, K + 1, 5): inv[k] = -sum(Pz[j] * inv[k - j] for j in range(5, k + 1, 5)) / Pz[0]
    pn = zmul(p, {k: Pz[0] * x for k, x in inv.items()})
    q0 = [pn[0][c].co(0) for c in range(5)]
    C = q0[0] * q0[3]; assert C == q0[1] * q0[2] == q0[4] ** 2, ("not symplectic up to a scalar", q0)
    sig = 1 / qsqrt(C); pn = zmul(pn, {0: sig})
    P11 = {k: pn[k][1].co(0) for k in range(0, K + 1, 5)}
    return pn, C, sig, P11[0], (P11.get(5, 0) / P11[0])


def ycs(bb, nmax, N):
    """pure-1/z solution: y_0 = -e_1, y_{n+1} = Int N y_n (zero constants); {n: (component, L)}"""
    y = {0: (1, L.const(-1, N))}
    for n in range(nmax):
        a, s = y[n]; a2 = (a - 1) % 5                     # (N y)_{a-1} = b_{a-1} y_a
        y[n + 1] = (a2, bb[a2].mul(s).integ())
    return y


def check_G(S, p, J=6):
    """(G): P Y^C_model(v(u)) = c0 Y^C_system mod u^J on z^j, -nJ <= j <= K - nJ (exact range); returns (ok, range, bad)"""
    K = len(p) - 1; N = J + 3
    nJ = max(n for n in range(5 * J + 5) if n - n // 5 < J)
    ym = ycs(S['bt'], nJ + 6, N); ys = ycs(S['b'], nJ + 6, N)
    for n in range(nJ + 1, nJ + 6): assert ym[n][1].v >= J, "model y_%d valuation" % n
    c0 = p[0][1].co(0); bad = []
    for j in range(-nJ, K - nJ + 1):
        acc = L(N, [], N)
        for n in range(max(0, -j), nJ + 1):
            a, y = ym[n]; acc = acc.add(p[j + n][a].mul(y))
        if j <= 0: acc = acc.sub(ys[-j][1].sc(c0))
        if any(acc.co(t) != 0 for t in range(min(acc.v, 0), J)): bad.append(j)
    return not bad, (-nJ, K - nJ), bad


def hatted(S, p, stale_tau=None):
    """{(k, b): L}: w^(-(1-c_b)/5) u^mu(b,k) p_{k,b}  (= tau^((1-c_b)/4) u^mu p); stale_tau: use that w instead (control)"""
    W = S['w'] if stale_tau is None else stale_tau; m0 = W.co(0); out = {}; cache = {}
    for k in range(len(p)):
        for b_ in range(5):
            e = 1 - CX[b_]
            if e not in cache:
                lead = root_rat(m0, Fr(-e, 5))
                cache[e] = W.powfrac(Fr(-e, 5), lead)
            out[(k, b_)] = cache[e].mul(p[k][b_]).sh(mu(b_, k))
    return out


def root_rat(x, e):
    """x^e for rational x > 0 that is an exact power (e = a/b)"""
    a, bq = e.numerator, e.denominator
    def iroot(n):
        r = round(n ** (1.0 / bq)) if n > 0 else 0
        for t in (r - 1, r, r + 1):
            if t >= 0 and t ** bq == n: return t
        raise ValueError("not an exact root")
    y = Fr(iroot(x.numerator), iroot(x.denominator))
    return y ** a
