# -*- coding: utf-8 -*-
"""ref_wkb.py -- referee check (key gaps/G2/referee) of lem:gau-colrows Steps 1-3 (GAPS_G2_5 s.9.4), by an independent
computation: the formal WKB amplitudes psi_i(q, z) = Sum_n z^n psi_(i,n)(q) of PF5 = (zD)^5 - 1 - q prod_(k=1..5)(5zD + kz)
(D = q d/dq) at the large-radius point, colour i (exponent u_i with D u_i = zeta^i L, L = (1 - 5^5 q)^(-1/5)), with the
constants of integration fixed at q = 0 by psi_i(0, z) = h_(w_i)(z) (colour-i form of (QRR)), computed directly from the
weight multisets w_i = {zeta^i - zeta^j}_(j != i) U {-5 zeta^i}; exact over Q(zeta) (basis 1, zeta, zeta^2, zeta^3).
Checks: transport operator T(L g) = 5 zeta^(4i) D g (computed, not assumed); no logarithm at any order; psi_(i,0) = L;
lem:set-hrot h_(w_i)(z) = h_(w_0)(zeta^-i z); the covariance psi_i(q, z) = psi_0(q, zeta^-i z) for i = 1..4.
Failing controls: colour i normalised by h_(w_0) instead of h_(w_i); PF5 perturbed by q z^2 (breaks z -> zeta z invariance).
Usage: py -I -B ref_wkb.py [nz] [nq].  Stdlib only; imports nothing from the project.
"""
import sys
from fractions import Fraction as Fr

FAIL = []


def check(name, ok, detail=""):
    print("  [%s] %s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        FAIL.append(name)


# ---- Q(zeta), zeta^4 = -1 - zeta - zeta^2 - zeta^3; elements are 4-tuples of Fractions
def zc(c):
    return (Fr(c), Fr(0), Fr(0), Fr(0))


def zpow(k):
    k %= 5
    if k == 4:
        return (Fr(-1), Fr(-1), Fr(-1), Fr(-1))
    v = [Fr(0)] * 4; v[k] = Fr(1); return tuple(v)


def zadd(x, y):
    return tuple(a + b for a, b in zip(x, y))


def zneg(x):
    return tuple(-a for a in x)


def zmul(x, y):
    p = [Fr(0)] * 7
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                p[i + j] += a * b
    for d in (6, 5, 4):                      # x^d = -x^(d-1) - x^(d-2) - x^(d-3) - x^(d-4)
        c = p[d]; p[d] = Fr(0)
        for t in range(1, 5):
            p[d - t] -= c
    return tuple(p[:4])


def zsigma(x, m):
    """Galois automorphism zeta -> zeta^m"""
    out = zc(0)
    for k, a in enumerate(x):
        if a:
            out = zadd(out, tuple(a * c for c in zpow(k * m)))
    return out


def zinv(x):
    num = zsigma(x, 2)
    num = zmul(num, zsigma(x, 3)); num = zmul(num, zsigma(x, 4))
    nrm = zmul(num, x)
    assert all(c == 0 for c in nrm[1:]) and nrm[0] != 0
    return tuple(c / nrm[0] for c in num)


def zscal(c, x):
    return tuple(Fr(c) * a for a in x)


# ---- power series in q over Q(zeta), lists of length NQ
def sadd(f, g):
    return [zadd(a, b) for a, b in zip(f, g)]


def smul(f, g):
    n = len(f); out = [zc(0)] * n
    for i, a in enumerate(f):
        if any(a):
            for j in range(n - i):
                if any(g[j]):
                    out[i + j] = zadd(out[i + j], zmul(a, g[j]))
    return out


def sD(f):
    return [zscal(k, a) for k, a in enumerate(f)]


def sq(f):
    return [zc(0)] + f[:-1]


def L_series(nq):
    out = []; c = Fr(1)
    for n in range(nq):
        out.append(zc(c)); c = c * Fr(5 * n + 1, 5) * 3125 / (n + 1)
    return out


def apply_N(F, lam):
    """(zD + lam) on F = [f_0, f_1, ..] (coefficients of z^m, each a q-series): (N F)_m = lam f_m + D f_(m-1)"""
    out = []
    for m in range(len(F)):
        t = smul(lam, F[m])
        if m >= 1:
            t = sadd(t, sD(F[m - 1]))
        out.append(t)
    return out


def pf5(F, lam, pert=False):
    """PF5 conjugated by e^(u/z), Du = lam, applied to F (truncated in z at len(F)); pert adds the term q z^2"""
    G = F
    for _ in range(5):
        G = apply_N(G, lam)
    H = F
    for k in range(1, 6):                    # (5N + kz) H: 5 lam h_m + 5 D h_(m-1) + k h_(m-1)
        NH = apply_N(H, lam)
        H = [sadd([zscal(5, a) for a in NH[m]], [zscal(k, a) for a in H[m - 1]] if m >= 1 else [zc(0)] * len(H[m]))
             for m in range(len(H))]
    out = [sadd(sadd(G[m], [zneg(a) for a in F[m]]), [zneg(a) for a in sq(H[m])]) for m in range(len(F))]
    if pert:
        out = [sadd(out[m], sq(F[m - 2])) if m >= 2 else out[m] for m in range(len(F))]
    return out
