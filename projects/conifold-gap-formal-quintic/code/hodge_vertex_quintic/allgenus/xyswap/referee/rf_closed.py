# -*- coding: utf-8 -*-
"""rf_closed.py -- referee's own implementation of the closed formula (7.9) of arXiv:2212.00320, Theorem 7.1
(y = z, x rational, initial form -y dx), for n = 1 and n = 2, SYMBOLIC in z (numeric rational parameters):
   omega^A_{g,n} = [hbar^2g] Sum_r prod_i D_{x_i}^{r_i} [w_i^{r_i}] Sum_Gamma hbar^(2 g(Gamma)) prod W_ij prod W_i,
   W_i = -(dz_i/w_i) exp(-w_i (S(w_i hbar d) - 1) x(z_i)),  W_ij = w_i w_j/((z_i - z_j)^2 - (hbar^2/4)(w_i + w_j)^2),
   D_x omega = d(omega/dx).   Expected (EO kernel): omega^EO_{g,n} = (-1)^n omega^A_{g,n}.
x' = q(z)/z^d with q a polynomial; a rational function of one variable is (N, A, K) = N(z)/(z^A q^K); of two variables
(N, A1, A2, K1, K2, L) = N(z1, z2)/(z1^A1 z2^A2 q(z1)^K1 q(z2)^K2 (z1 - z2)^L), N a dict of monomials.
"""
from fractions import Fraction as Fr
from math import factorial, comb
from collections import Counter
from rf_ser import lp_mul, lp_deriv, partitions


def setup(X):
    dX = lp_deriv({e: Fr(c) for e, c in X.items()})
    d = -min(dX)
    q = {e + d: c for e, c in dX.items()}
    return dict(X={e: Fr(c) for e, c in X.items()}, d=d, q=q, dq=lp_deriv(q))


def vertex(C, g):
    """{p: Laurent polynomial} with [hbar^2g][w^p] exp(E), E = -Sum_k w^(2k+1) hbar^(2k) x^(2k)/(4^k (2k+1)!)"""
    if g == 0:
        return {0: {0: Fr(1)}}
    out = {}
    for lam in partitions(g):
        t = {0: Fr(1)}
        for k in lam:
            t = lp_mul(t, {e: -c / (4 ** k * factorial(2 * k + 1)) for e, c in lp_deriv(C["X"], 2 * k).items()})
        m = 1
        for v in Counter(lam).values():
            m *= factorial(v)
        p = 2 * g + len(lam)
        acc = out.setdefault(p, {})
        for e, c in t.items():
            acc[e] = acc.get(e, Fr(0)) + c / m
    return out


def _lp_to_frac1(L):
    """Laurent polynomial -> (N, A, 0)"""
    A = max(0, -min(L)) if L else 0
    return ({e + A: c for e, c in L.items() if c}, A, 0)


def D1(C, F):
    N, A, K = F
    a, k = A - C["d"], K + 1
    out = {}
    for e, c in N.items():
        for f, qc in C["q"].items():                         # N' z q - a N q
            if e:
                out[e + f] = out.get(e + f, Fr(0)) + e * c * qc
            out[e + f] = out.get(e + f, Fr(0)) - a * c * qc
        for f, qc in C["dq"].items():                        # - k N z q'
            out[e + f + 1] = out.get(e + f + 1, Fr(0)) - k * c * qc
    return ({e: c for e, c in out.items() if c}, a + 1, k + 1)


def eval1(C, F, z):
    N, A, K = F
    z = Fr(z)
    qz = sum((c * z ** e for e, c in C["q"].items()), Fr(0))
    return sum((c * z ** e for e, c in N.items()), Fr(0)) / (z ** A * qz ** K)


def omega_g1_terms(C, g):
    """[(D^r f_r, T_r)] as rational functions; omega^A_{g,1} = Sum_r D^r f_r dz, Phi^A_g = Sum_r T_r"""
    out = []
    for p, L in vertex(C, g).items():
        r = p - 1
        F = _lp_to_frac1({e: -c for e, c in L.items()})     # f_r = -[w^(r+1)] e^E
        T = None
        for j in range(r):
            if j == r - 1:                                    # T_r = D^(r-1) f_r / x' = N z^d/(z^A q^(K+1))
                N, A, K = F
                T = (N, A - C["d"], K + 1)
            F = D1(C, F)
        out.append((F, T))
    return out


def omega_g1_value(C, g, z):
    """omega^EO_{g,1}/dz at z (numeric)"""
    return -sum((eval1(C, F, z) for F, _ in omega_g1_terms(C, g)), Fr(0))


def at_zero(C, F):
    """value at z = 0 of (N, A, K); raises if singular"""
    N, A, K = F
    m0 = min(e for e, c in N.items() if c)
    if m0 > A:
        return Fr(0)
    assert m0 == A, "pole at z = 0"
    return N[A] / C["q"].get(0, Fr(0)) ** K


def deg_inf(C, F):
    """degree in z at infinity of (N, A, K)"""
    N, A, K = F
    return max(e for e, c in N.items() if c) - A - max(C["q"]) * K


# ---------------- two points ----------------
def _bmul_poly1(N, P, var):
    out = {}
    for (e1, e2), c in N.items():
        for f, pc in P.items():
            k = (e1 + f, e2) if var == 0 else (e1, e2 + f)
            out[k] = out.get(k, Fr(0)) + c * pc
    return out


def _badd(acc, N, s=Fr(1)):
    for k, c in N.items():
        acc[k] = acc.get(k, Fr(0)) + s * c


def D2(C, F, var):
    """D_{x_var} on (N, A1, A2, K1, K2, L)"""
    N, A1, A2, K1, K2, L = F
    A, K = (A1, K1) if var == 0 else (A2, K2)
    a, k = A - C["d"], K + 1
    sgn = 1 if var == 0 else -1                              # d/dz_var (z1 - z2) = sgn
    dN = {}
    for (e1, e2), c in N.items():
        e = e1 if var == 0 else e2
        if e:
            dN[(e1 - 1, e2) if var == 0 else (e1, e2 - 1)] = e * c
    zq = {e + 1: c for e, c in C["q"].items()}
    zdq = {e + 1: c for e, c in C["dq"].items()}
    t = {}
    _badd(t, _bmul_poly1(dN, zq, var))
    _badd(t, _bmul_poly1(N, C["q"], var), Fr(-a))
    _badd(t, _bmul_poly1(N, zdq, var), Fr(-k))
    out = {}
    for (e1, e2), c in t.items():                            # times (z1 - z2)
        out[(e1 + 1, e2)] = out.get((e1 + 1, e2), Fr(0)) + c
        out[(e1, e2 + 1)] = out.get((e1, e2 + 1), Fr(0)) - c
    _badd(out, _bmul_poly1(N, zq, var), Fr(-L * sgn))
    out = {kk: c for kk, c in out.items() if c}
    if var == 0:
        return (out, a + 1, A2, k + 1, K2, L + 1)
    return (out, A1, a + 1, K1, k + 1, L + 1)


def omega_g2_value(C, g, b1, b2):
    """omega^EO_{g,2}/(dz1 dz2) at (b1, b2), single-edge graph, edge expanded in hbar"""
    b1, b2 = Fr(b1), Fr(b2)
    q1 = sum((c * b1 ** e for e, c in C["q"].items()), Fr(0))
    q2 = sum((c * b2 ** e for e, c in C["q"].items()), Fr(0))
    tot = Fr(0)
    for k in range(g + 1):
        for g1 in range(g - k + 1):
            g2 = g - k - g1
            V1, V2 = vertex(C, g1), vertex(C, g2)
            for a in range(2 * k + 1):
                cf = Fr(comb(2 * k, a), 4 ** k)
                for p1, L1 in V1.items():
                    for p2, L2 in V2.items():
                        r1, r2 = a + p1, 2 * k - a + p2
                        A1 = max(0, -min(L1)); A2 = max(0, -min(L2))
                        N = {}
                        for e1, c1 in L1.items():
                            for e2, c2 in L2.items():
                                N[(e1 + A1, e2 + A2)] = N.get((e1 + A1, e2 + A2), Fr(0)) + cf * c1 * c2
                        F = (N, A1, A2, 0, 0, 2 * k + 2)
                        for _ in range(r1):
                            F = D2(C, F, 0)
                        for _ in range(r2):
                            F = D2(C, F, 1)
                        N, a1, a2, k1, k2, L = F
                        num = sum((c * b1 ** e1 * b2 ** e2 for (e1, e2), c in N.items()), Fr(0))
                        tot += num / (b1 ** a1 * b2 ** a2 * q1 ** k1 * q2 ** k2 * (b1 - b2) ** L)
    return tot                                                # (-1)^2 = 1


def omega_03_value(C, b):
    """omega^EO_{0,3} = Sum_i D_{x_i}[prod_{j != i} (z_i - z_j)^-2] (the g = 0 trees of (7.9) times (-1)^3)"""
    tot = Fr(0)
    for i in range(3):
        o = [b[j] for j in range(3) if j != i]
        # f(z) = 1/((z - o0)^2 (z - o1)^2); D_x f = d/dz (f z^d / q)
        P = lp_mul(lp_mul({1: Fr(1), 0: -Fr(o[0])}, {1: Fr(1), 0: -Fr(o[0])}),
                   lp_mul({1: Fr(1), 0: -Fr(o[1])}, {1: Fr(1), 0: -Fr(o[1])}))
        num = {C["d"]: Fr(1)}
        den = lp_mul(P, C["q"])                               # g(z) = z^d/(P q); g' = (num' den - num den')/den^2
        z = Fr(b[i])
        ev = lambda L: sum((c * z ** e for e, c in L.items()), Fr(0))
        tot += (ev(lp_deriv(num)) * ev(den) - ev(num) * ev(lp_deriv(den))) / ev(den) ** 2
    return tot
