# -*- coding: utf-8 -*-
"""rf_gauge.py -- REFEREE: own solver for the regular gauge of Theorem U (UNIFORM_REDUCTION s.3), written from the statement only.
Systems (th = z d/du):  quintic  th Y_a = b_a Y_(a+1),  b_a = -c_a/(1-u),  c = (C1, 1/(u C0^2 C1^2), C1, C0, C0)  (rf_core.periods);
model pulled back by v = u tau^(-5/4):  bt_b = v'(u) bmod_b(v),  bmod = (1, -1/v, 1, -1, -1).  Gauge Y = P Yt, z P' = A P - P At,
P = Sum_k z^k P_k, P_k = Sum_b p_(k,b) E_(b+k,b):   (E_k,b)  b_(b+k-1) p_(k,b) - bt_(b-1) p_(k,b-1) = p_(k-1,b)'.
Order k: particular solution q_k (q_(k,0) = 0) + rho_k n_k (kernel, n_(k,0) = 1); rho_k from the solvability of order k+1, a first-order
ODE rho' + kappa rho = g solved in Laurent series; at k = 0 mod 5 the residue of kappa is 0 and the free constant is the normalisation:
  norm 'can'    : p_(5q,1)(0) = delta_(q,0)     (P_11(0, z) = 1, the canonical gauge)
  norm 'stored' : p_(5q,0)(0) = delta_(q,0)     (the stored normalisation of the project, as described in GAUGE_IDENTITY s.1)
Regularity of every p_(k,b) (Theorem U) is ASSERTED, not assumed.  Hatted: Ghat[(k, b)] = tau^((1-c_b)/4) u^mu(b,k) p_(k,b)(u).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_core as C
CEXP = (0, -1, 3, 2, 1)


def mu(b, k):
    v = 4 * k + CEXP[b] - CEXP[(b + k) % 5]; assert v % 5 == 0
    return v // 5


# Laurent series x = (v, c): c[j] = coefficient of u^(v+j); absolute precision v + len(c)
def Lp(x): return x[0] + len(x[1])
def Lget(x, m): return x[1][m - x[0]] if x[0] <= m < Lp(x) else 0


def Lnorm(x):
    v, c = x; i = 0
    while i < len(c) and not c[i]: i += 1
    return (v + i, c[i:])


def Ladd(F, x, y, s=1):
    P = min(Lp(x), Lp(y)); v = min(x[0], y[0])
    return Lnorm((v, [F.red(Lget(x, m) + s * Lget(y, m)) for m in range(v, P)])) if P > v else (P, [])


def Lmul(F, x, y):
    x = Lnorm(x); y = Lnorm(y); P = min(Lp(x) + y[0], Lp(y) + x[0]); v = x[0] + y[0]; n = P - v
    if n <= 0 or not x[1] or not y[1]: return (P, [])
    return Lnorm((v, C.mul(F, x[1], y[1], n)))


def Linv(F, x):
    x = Lnorm(x); assert x[1], "inverse of a series that vanishes to its precision"
    return (-x[0], C.inv(F, x[1], len(x[1])))


def Lder(F, x): return Lnorm((x[0] - 1, [F.red((x[0] + j) * c) for j, c in enumerate(x[1])]))
def Lsc(F, x, a): return Lnorm((x[0], [F.red(a * c) for c in x[1]]))


def Lode(F, kap, g, free0):
    """rho' + kap rho = g in Laurent series; returns (rho with free constant 0, homogeneous solution with rho(0) = 1 or None)"""
    kap = Lnorm(kap); g = Lnorm(g); assert kap[0] >= -1 or not kap[1], "kappa: pole of order > 1"
    km1 = Lget(kap, -1); m0 = g[0] + 1
    if free0: m0 = min(m0, 0)
    M = min(Lp(g) + 1, Lp(kap) + m0 + 1); rho = {}
    for m in range(m0, M):
        rhs = Lget(g, m - 1) - sum(Lget(kap, j) * rho[m - 1 - j] for j in range(0, m - m0))
        den = F.red(m + km1)
        if den == 0: assert F.red(rhs) == 0, "resonance with a non-zero right-hand side (log term)"; rho[m] = 0
        else: rho[m] = F.red(rhs * F.inv(den))
    ra = Lnorm((m0, [rho[m] for m in range(m0, M)])) if M > m0 else (M, [])
    if not (free0 and F.red(km1) == 0): return ra, None
    hh = {0: F.of(1)}; Mh = Lp(kap) + 1
    for m in range(1, Mh): hh[m] = F.red(-sum(Lget(kap, j) * hh[m - 1 - j] for j in range(0, m)) * F.inv(m))
    return ra, (0, [hh[m] for m in range(Mh)])


def systems(F, NU, al=0, be=0):
    C0, C1, tc, gam = C.periods(F, NU + 2, al, be); one_m = (0, [F.of(1), F.of(-1)] + [0] * NU)
    S = lambda a: (0, list(a[:NU + 1]))
    c1 = Linv(F, Lmul(F, (1, [F.of(1)] + [0] * NU), Lmul(F, Lmul(F, S(C0), S(C0)), Lmul(F, S(C1), S(C1)))))
    cs = [S(C1), c1, S(C1), S(C0), S(C0)]; om = Linv(F, one_m)
    b = [Lsc(F, Lmul(F, c, om), -1) for c in cs]
    tau = C.inv(F, C.itau(F, NU + 2), NU + 2); w = C.spow(F, tau, Fr(-5, 4), NU + 2)
    v = (1, w[:NU + 1]); vp = Lder(F, v)
    bt = [vp, Lsc(F, Lmul(F, vp, Linv(F, v)), -1), vp, Lsc(F, vp, -1), Lsc(F, vp, -1)]
    pb = b[0]; pt = bt[0]
    for j in range(1, 5): pb = Lmul(F, pb, b[j]); pt = Lmul(F, pt, bt[j])
    target = Lsc(F, Linv(F, Lmul(F, (1, [F.of(1)] + [0] * NU), Lmul(F, Lmul(F, one_m, one_m), Lmul(F, Lmul(F, one_m, one_m), one_m)))), -1)
    assert not Ladd(F, pb, target, -1)[1] and not Ladd(F, pt, target, -1)[1], "prod b = prod bt = -1/(u (1-u)^5) fails"
    return b, bt, tau, (C0, C1, tc, gam)


def gauge(F, KMAX, NU, al=0, be=0, norm="can"):
    b, bt, tau, per = systems(F, NU, al, be); binv = [Linv(F, x) for x in b]

    def particular(k, r):
        q = [(NU, [])]
        for j in range(1, 5): q.append(Lmul(F, Ladd(F, r[j], Lmul(F, bt[j - 1], q[-1])), binv[(j + k - 1) % 5]))
        return q, Ladd(F, r[0], Lmul(F, bt[4], q[4]))

    def kernel(k):
        n = [(0, [F.of(1)] + [0] * NU)]
        for j in range(1, 5): n.append(Lmul(F, Lmul(F, n[-1], bt[j - 1]), binv[(j + k - 1) % 5]))
        return n
    r = [(NU, [])] * 5; P = []; info = []
    for k in range(KMAX + 1):
        q, S = particular(k, r); assert not S[1], "solvability of order %d violated" % k
        n = kernel(k)
        A = particular(k + 1, n)[1]; B = particular(k + 1, [Lder(F, x) for x in n])[1]
        Cc = particular(k + 1, [Lder(F, x) for x in q])[1]; Ai = Linv(F, A)
        kap = Lmul(F, B, Ai); g = Lsc(F, Lmul(F, Cc, Ai), -1)
        ra, rh = Lode(F, kap, g, k % 5 == 0)
        pk = [Ladd(F, q[j], Lmul(F, ra, n[j])) for j in range(5)]
        if k % 5 == 0:
            idx = 1 if norm == "can" else 0; hn = [Lmul(F, rh, n[j]) for j in range(5)]
            X = F.red((F.of(int(k == 0)) - Lget(pk[idx], 0)) * F.inv(Lget(hn[idx], 0)))
            pk = [Ladd(F, pk[j], Lsc(F, hn[j], X)) for j in range(5)]
        assert all(x[0] >= 0 for x in pk), "Theorem U violated: p_(%d,b) has a pole: %s" % (k, [x[0] for x in pk])
        info.append(dict(k=k, res=Lget(Lnorm(kap), -1), prec=min(Lp(x) for x in pk), val=[x[0] for x in pk]))
        P.append(pk); r = [Lder(F, x) for x in pk]
    return dict(p=P, tau=tau, info=info, per=per)


def hatted(F, res, K):
    tau = res['tau']; out = {}
    tw = {b: C.spow(F, tau, Fr(1 - CEXP[b], 4), K) for b in range(5)}
    for k, pk in enumerate(res['p']):
        for b in range(5):
            m = mu(b, k)
            if m >= K: out[(k, b)] = [0] * K; continue
            x = pk[b]; assert Lp(x) + m >= K, "gauge precision too low at (k, b) = (%d, %d)" % (k, b)
            ser = [Lget(x, j) for j in range(K - m)]
            out[(k, b)] = [0] * m + C.mul(F, tw[b], ser, K - m)
    return out


def stored_m5(KM=5, NU=16):
    """[z^5] P^stored_11(0, z), exact over Q (P_stored = P_can (1 + m5 z^5 + ..))"""
    res = gauge(C.QF(), KM, NU, norm="stored")
    return Lget(res['p'][5][1], 0)


if __name__ == "__main__":
    F = C.QF(); res = gauge(F, 10, 26, norm="stored")
    print("stored: P^s_11(0,z) at z^5, z^10:", Lget(res['p'][5][1], 0), Lget(res['p'][10][1], 0))
    print("residues of kappa k = 0..10:", [str(i['res']) for i in res['info']])
    rc = gauge(F, 10, 26, norm="can")
    print("canonical: p_(5,0)(0), p_(5,1)(0), p_(5,2)(0):", [str(Lget(rc['p'][5][j], 0)) for j in range(5)])
    print("p_(0,0)(u):", [str(Lget(rc['p'][0][0], j)) for j in range(5)])
