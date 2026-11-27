# -*- coding: utf-8 -*-
"""xy_closed.py -- the closed formula of arXiv:2212.00320, Theorem 7.1 (eq. (7.9)), for the curve x = z^4 + v/z, y = z.

   omega^{(g)}_{m,0} = [hbar^2g] Sum_{r_i >= 0} prod_i D_{x_i}^{r_i} [w_i^{r_i}] Sum_Gamma hbar^{2 g(Gamma)} prod_edges W_ij prod_vertices W_i,
   W_i  = -(dz_i/w_i) exp(-w_i (S(w_i hbar d_{z_i}) - 1) x(z_i)),   S(u) = sinh(u/2)/(u/2),   D_x omega = d(omega/dx),
   W_ij = w_i w_j / ((z_i - z_j)^2 - (hbar^2/4)(w_i + w_j)^2),  Gamma: connected simple graphs on m numbered vertices.
Their convention omega^{(0)}_{1,0} = -y dx gives omega^{ABDKS}_{g,n} = (-1)^n omega^{EO}_{g,n}.
For x = z^4 + v/z:  -w (S(w hbar d) - 1) x = -(1/2) w^3 hbar^2 z^2 - (1/80) w^5 hbar^4 - v Sum_{k>=1} w^{2k+1} hbar^{2k}/(4^k (2k+1) z^{2k+1}).
(A) n = 1 symbolically: a Laurent polynomial in z and v is a dict {(ez, ev): Fraction}; a rational function is (N, k) = N/(4z^5 - v)^k.
(B) n >= 2 at rational points z_i = b_i (v a number) by Taylor jets in delta_i = z_i - b_i (exact over Q).
"""
import itertools
from fractions import Fraction as Fr
from math import factorial


def exponent_terms(gmax, v=None):
    """terms of E: list of (ew, eh, ez, ev, coeff);  v = None keeps v symbolic (ev = power of v)"""
    out = [(3, 2, 2, 0, Fr(-1, 2)), (5, 4, 0, 0, Fr(-1, 80))]
    for k in range(1, gmax + 1):
        out.append((2 * k + 1, 2 * k, -(2 * k + 1), 1, Fr(-1, 4 ** k * (2 * k + 1))))
    return [t for t in out if t[1] <= 2 * gmax]


def expE_n1(g):
    """[hbar^2g] exp(E) as {ew: {(ez, ev): c}}  (symbolic v)"""
    terms = exponent_terms(g)
    cur = {(0, 0): {(0, 0): Fr(1)}}            # (ew, eh) -> zv-poly ; cur = E^n/n!
    tot = {}
    n = 0
    while cur:
        for (ew, eh), pol in cur.items():
            if eh == 2 * g:
                d = tot.setdefault(ew, {})
                for kk, c in pol.items(): d[kk] = d.get(kk, 0) + c
        n += 1; nxt = {}
        for (ew, eh), pol in cur.items():
            for (tw, th, tz, tv, tc) in terms:
                if eh + th > 2 * g: continue
                d = nxt.setdefault((ew + tw, eh + th), {})
                for (ez, ev), c in pol.items():
                    kk = (ez + tz, ev + tv); d[kk] = d.get(kk, 0) + c * tc / n
        cur = {k: p for k, p in nxt.items() if any(p.values())}
    return tot


def pmul(a, b):
    out = {}
    for (z1, v1), c1 in a.items():
        for (z2, v2), c2 in b.items():
            k = (z1 + z2, v1 + v2); out[k] = out.get(k, 0) + c1 * c2
    return {k: c for k, c in out.items() if c}


def padd(a, b, s=1):
    out = dict(a)
    for k, c in b.items(): out[k] = out.get(k, 0) + s * c
    return {k: c for k, c in out.items() if c}


def pdz(a):
    return {(z - 1, v): c * z for (z, v), c in a.items() if z}


DPOL = {(5, 0): Fr(4), (0, 1): Fr(-1)}        # D = 4 z^5 - v


def Dx(N, k):
    """D_x (N/D^k dz) = (((N z^2)' D - 20 (k+1) z^6 N) / D^(k+2)) dz"""
    Nz2 = {(z + 2, v): c for (z, v), c in N.items()}
    a = pmul(pdz(Nz2), DPOL)
    b = {(z + 6, v): 20 * (k + 1) * c for (z, v), c in N.items()}
    return padd(a, b, -1), k + 2


def omega_g1_terms(g):
    """list of (r, N_r, k_r): the r-th term D_x^r [w^r] W_1 of omega^{ABDKS}_{g,1} as N/D^k dz, and the primitive pieces
    T_r = D_x^(r-1)([w^r] W_1)/dx = (N z^2, k+1) (r >= 1)"""
    ex = expE_n1(g); out = []
    for ew, pol in sorted(ex.items()):
        r = ew - 1
        if r < 0: continue
        N, k = {kk: -c for kk, c in pol.items()}, 0      # f_r = -[w^(r+1)] e^E
        prim = None
        for j in range(r):
            if j == r - 1: prim = ({(z + 2, v): c for (z, v), c in N.items()}, k + 1)
            N, k = Dx(N, k)
        if r == 0: prim = None
        out.append((r, N, k, prim))
    return out


def common(terms):
    """sum of (N, k) rational functions over the common denominator D^K"""
    K = max(k for (_, k) in terms) if terms else 0
    tot = {}
    for N, k in terms:
        M = dict(N)
        for _ in range(K - k): M = pmul(M, DPOL)
        tot = padd(tot, M)
    return tot, K


def omega_g1(g):
    """omega^{EO}_{g,1} = -(sum_r D_x^r [w^r] W_1) as (N, K) with N/D^K dz"""
    tot, K = common([(N, k) for (_, N, k, _) in omega_g1_terms(g)])
    return {kk: -c for kk, c in tot.items()}, K


def ev_rat(N, k, z, v):
    num = sum(c * Fr(z) ** ez * Fr(v) ** ev for (ez, ev), c in N.items())
    return num / (4 * Fr(z) ** 5 - v) ** k


def bernoulli(nmax):
    B = [Fr(0)] * (nmax + 1); B[0] = Fr(1)
    for m in range(1, nmax + 1):
        B[m] = -sum(Fr(factorial(m + 1), factorial(j) * factorial(m + 1 - j)) * B[j] for j in range(m)) / (m + 1)
    return B
