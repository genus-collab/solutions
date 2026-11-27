# -*- coding: utf-8 -*-
"""pb_gauge.py -- Stage-5 'pb' audit: the regular gauge of Theorem U built from scratch, and Phi^can.

Quintic (conifold frame): th Y_a = b_a Y_(a+1), b = (tc', -1/(u varpi0^2 (1-u)^3 tc'^2), tc', -varpi0/(1-u), -varpi0/(1-u)).
Model pulled back by v = u tau^(-5/4): bt = (v', -v'/v, v', -v', -v').
Gauge z P' = A P - P At, P = sum_k z^k P_k, P_k graded: entry (b+k, b) = p_(k,b).  Order z^K:
    b_(b+K-1) p_(K,b) - bt_(b-1) p_(K,b-1) = p'_(K-1,b)        (b in Z/5).
The map has a one-dimensional kernel kap^(K) (kap_0 = 1, kap_b = kap_(b-1) bt_(b-1)/b_(b+K-1)) and the cokernel functional
l_b = 1/kap^(K-1)_b; so P_K = pi_K + c_K kap^(K), pi_(K,0) = 0, and c_K solves the solvability ODE of order K+1:
    5 c_K' + (sum_b kap_b'/kap_b) c_K = - sum_b pi_(K,b)'/kap_b.
Laurent series are used throughout; regularity of every p_(K,b) is ASSERTED, not imposed.  At K = 0 mod 5 the ODE has one
free constant: normalisation 'can' (P_11(0,z) = 1), 'stored' (p_(5q,0)(0) = 0, q >= 1) or a dict {K: shift}.
Model rows F_b (colour 0) from UR s.2 Prop. S: F_1 = T_-1 F_0, F_2 = -T_-2 F_1, F_3 = T_2 F_2, F_4 = -T_1 F_3, F_0 = -T_0 F_4,
T_k = -1 + (t/5)(k - 4 t d/dt).  Hatted first row and Phi^can = C0 sum_b Ghat_0b F_b(tau w).
"""
from fractions import Fraction as Fr
import pb_core as C
from pb_core import LS, red, finv, iszero

CW = (0, -1, 3, 2, 1)        # hatted weights c_b
EPS = (1, -1, -1, 1, 1)      # F_b(0)


def systems(K, N, tau_one=False):
    """quintic b_a and model bt_a to O(u^N) in the conifold frame of K (a pb_core.Conifold)."""
    one_m = LS.from_list([1, -1], 0, N + 2)
    w0 = K.varpi0.trunc(N + 2); tcp = K.tc.trunc(N + 3).deriv()
    inv1m = one_m.inv()
    den = w0.mul(w0).mul(tcp).mul(tcp).mul(one_m).mul(one_m).mul(one_m).shift(1)   # u varpi0^2 (1-u)^3 tc'^2
    b1 = den.inv().neg()
    b3 = w0.mul(inv1m).neg()
    b = [tcp, b1, tcp, b3, b3]
    if tau_one:
        v = LS.from_list([0, 1], 0, N + 3)                       # control: v = u (wrong natural variable)
    else:
        V, v = C.model_v(N + 3)
    vp = v.deriv()
    bt = [vp, vp.div(v).neg(), vp, vp.neg(), vp.neg()]
    return [x.trunc(N) for x in b], [x.trunc(N) for x in bt]


def solve_ode(rho, sig, free=None):
    """5 c' + rho c = sig (Laurent).  Returns (c with free constant 0, homogeneous h with h(0) = 1 or None)."""
    ur = rho.shift(1); us = sig.shift(1).norm()
    e = ur.coef(0)
    lo = min(us.v, 0) if us.c else 0
    top = min(us.prec, ur.prec + lo)
    res = {}
    hom = None
    for n in range(lo, top):
        s = us.coef(n) if n >= us.v else 0
        for j in range(1, n - lo + 1):
            s -= ur.coef(j) * res[n - j]
        piv = red(5 * n + e)
        if iszero(piv):
            assert iszero(s), ('resonance obstructed', n)
            res[n] = red(0); hom = n
        else:
            res[n] = red(s * finv(piv))
    c = LS(lo, [res[n] for n in range(lo, top)], top)
    h = None
    if hom is not None:
        hres = {hom: red(1)}
        htop = ur.prec + hom
        for n in range(hom + 1, htop):
            s = 0
            for j in range(1, n - hom + 1):
                s -= ur.coef(j) * hres[n - j]
            hres[n] = red(s * finv(red(5 * n + e)))
        h = LS(hom, [hres[n] for n in range(hom, htop)], htop)
    return c, h


def gauge(b, bt, Kmax, norm='can'):
    """list P[K] = [p_(K,0..4)] (LS in u), K = 0..Kmax."""
    P = []
    one = LS.const(1, b[0].prec)
    for K in range(Kmax + 1):
        r = [x.deriv() for x in P[K - 1]] if K else [LS(0, [red(0)] * b[0].prec, b[0].prec)] * 5
        kap = [one]
        for bb in range(1, 5):
            kap.append(kap[-1].mul(bt[bb - 1]).div(b[(bb + K - 1) % 5]))
        clos = kap[4].mul(bt[4]).div(b[(K + 4) % 5]).sub(one)
        assert clos.ok_zero(), 'kernel does not close (equal products?)'
        pi = [LS(0, [red(0)] * r[1].prec, r[1].prec)]
        for bb in range(1, 5):
            pi.append(r[bb].add(bt[bb - 1].mul(pi[-1])).div(b[(bb + K - 1) % 5]))
        if K:
            cons = bt[4].mul(pi[4]).neg().sub(r[0])
            assert cons.ok_zero(), ('solvability failed at order', K)
        rho = kap[0].deriv().div(kap[0])
        sig = pi[0].deriv().div(kap[0])
        for bb in range(1, 5):
            rho = rho.add(kap[bb].deriv().div(kap[bb]))
            sig = sig.add(pi[bb].deriv().div(kap[bb]))
        c, h = solve_ode(rho, sig.neg())
        if h is not None:
            p1 = pi[1].add(c.mul(kap[1]))
            h1 = h.mul(kap[1])
            if norm == 'can':
                target = red(1 if K == 0 else 0); lam = red((target - p1.coef(0)) * finv(h1.coef(0)))
            elif norm == 'stored':
                if K == 0:
                    lam = red((1 - p1.coef(0)) * finv(h1.coef(0)))
                else:
                    lam = red(-c.coef(0) * finv(h.coef(0)))
            else:
                base = red((1 if K == 0 else 0) - p1.coef(0)) * finv(h1.coef(0))
                lam = red(base + red(norm.get(K, 0)))
            c = c.add(h.scale(lam))
        else:
            assert K % 5, 'free constant expected at K = 0 mod 5'
        PK = [pi[bb].add(c.mul(kap[bb])) for bb in range(5)]
        for bb in range(5):
            s = PK[bb].norm()
            assert (not s.c) or s.v >= 0, ('irregular gauge entry', K, bb, s.v)
        P.append(PK)
    return P


def model_rows(nmax):
    """F_b (colour 0) as coefficient lists [t^0..t^nmax], from the Prop. S recursion."""
    def T(k, F):
        out = [red(0)] * len(F)
        for n in range(len(F)):
            x = -F[n]
            if n >= 1:
                x += red(Fr(k - 4 * (n - 1), 5)) * F[n - 1]
            out[n] = red(x)
        return out

    def chain(F0):
        F1 = T(-1, F0); F2 = [red(-x) for x in T(-2, F1)]; F3 = T(2, F2); F4 = [red(-x) for x in T(1, F3)]
        G0 = [red(-x) for x in T(0, F4)]
        return [F0, F1, F2, F3, F4], G0

    L = nmax + 2
    a = [red(1)] + [red(0)] * (L - 1)
    for n in range(1, nmax + 1):
        a[n] = red(0)
        _, G0 = chain(a)
        R = red(G0[n + 1] - a[n + 1])       # = 4n a_n + R_n with a_n = 0
        a[n] = red(-R * finv(4 * n))
        _, G0 = chain(a)
        assert iszero(G0[n + 1] - a[n + 1])
    rows, G0 = chain(a)
    assert all(iszero(G0[n] - a[n]) for n in range(nmax + 1)), 'model closing'
    return [r[:nmax + 1] for r in rows]


def mu0(k):
    q, r = divmod(k, 5)
    return 4 * q + (0, 1, 2, 3, 3)[r]


def phi_can(K, P, F, tau, nmax, M):
    """[w^n] Phi^can (LS in u, to O(u^M)), n = 0..nmax: C0 sum_k tau^((1-c_b)/4) u^mu p_(k,b) F_(b,n-k) tau^(n-k)."""
    C0 = K.varpi0.trunc(M)
    tw = [tau.power(Fr(1 - CW[bb], 4)).trunc(M) for bb in range(5)]
    tpow = [LS.const(1, M)]
    for j in range(1, nmax + 1):
        tpow.append(tpow[-1].mul(tau).trunc(M))
    out = []
    for n in range(nmax + 1):
        acc = LS(0, [red(0)] * M, M)
        for k in range(n + 1):
            bb = (-k) % 5; m = mu0(k)
            term = P[k][bb].trunc(M - m).shift(m)
            if term.prec < M:
                raise ValueError('gauge precision too low at k=%d: have u^%d, need u^%d' % (k, term.prec, M))
            term = term.mul(tw[bb]).mul(tpow[n - k]).scale(F[bb][n - k])
            acc = acc.add(term.trunc(M))
        out.append(C0.mul(acc).trunc(M))
    return out
