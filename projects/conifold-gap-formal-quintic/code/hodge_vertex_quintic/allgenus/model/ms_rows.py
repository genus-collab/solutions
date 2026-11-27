# -*- coding: utf-8 -*-
"""ms_rows.py -- flat-basis data of the frozen (model) vertex theory V = M_0.omega^ over Q.

Input (read only): natural_form/lg/lg_symbol_rows_Q.json, g[a][n] = [zhat^n] J_{c_a}(-4 zhat), c = (0,-1,3,2,1), n <= 40.
Symbol rows: R^_a^(i)(zhat) = N_a zeta^(i(a+1)) g^a(zeta^-i zhat), N = (1,-1,-1,1,1).

Basis used internally: h^a := N_a H^a (a in Z/5), i.e. h = (1, -H, -H^2, H^3, H^4).  In this basis (derived by Fourier
transform of the project's canonical-frame formulas, fq3_mvertex docstring; regression in ms_regress.py):
  TFT:     omega_{g,n}(h^{b_1}, .., h^{b_n}) = 5 [ Sum b_i = 3 - 3g mod 5 ],      eta(h^a, h^b) = 5 [a + b = 3 mod 5];
  legs:    R^-1_k h^a = (-1)^k g^a_k h^(a-k);
  leaves:  T(e) = -R^-1_(e-1) h^0 = (-1)^e g^0_(e-1) h^(1-e)      (e >= 2, carries psi^e);
  edges:   E(k,l) = [psi^k psi'^l] (eta^-1 - R^-1(psi) eta^-1 R^-1(psi')^T)/(psi + psi')
                  = Sum_{t=0..l} (-1)^t Gh(k+1+t, l-t),   Gh(m,n) = -(1/5)(-1)^(m+n) Sum_a g^a_m g^(3-a)_n h^(a-m) (x) h^(3-a-n).
CHARGE of a half-edge state (index b, psi-exponent e): b + e - 1 mod 5.  Then: every vertex has total charge 0 (leaves are
neutral), the two ends of an edge carry opposite charges (q, -q), a leg tau_d(h^a) injects the charge a + d - 1 whatever the
order k of R^-1 it uses.  In charge variables
  ECOEF(k, l, q) = (1/5)(-1)^(k+l) Sum_{t=0..l} (-1)^t g^(q+2+t)_(k+1+t) g^(1-q-t)_(l-t)     (first end: charge q, exponent k).
"""
import json, os, sys
from fractions import Fraction as Fr
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
ROWS_JSON = os.path.join(PROJ, "natural_form", "lg", "lg_symbol_rows_Q.json")
SCRATCH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scratch")

C_EXP = (0, -1, 3, 2, 1)          # c_a
N_SIGN = (1, -1, -1, 1, 1)        # N_a : H^a = N_a h^a
NMAX = 40

with open(ROWS_JSON) as _f:
    _d = json.load(_f)
assert _d["c"] == list(C_EXP) and _d["N"] == list(N_SIGN)
G = [[Fr(x) for x in row] for row in _d["g"]]       # G[a][n]
assert all(len(row) >= NMAX + 1 for row in G) and all(row[0] == 1 for row in G)


def g(a, n):
    return G[a % 5][n]


@lru_cache(maxsize=None)
def ecoef(k, l, q):
    """edge coefficient: first half-edge (charge q, exponent k), second half-edge (charge -q, exponent l)"""
    s = Fr(0)
    for t in range(l + 1):
        s += (-1) ** t * g(q + 2 + t, k + 1 + t) * g(1 - q - t, l - t)
    return s * (-1) ** (k + l) / 5


@lru_cache(maxsize=None)
def eloop(k, l):
    """loop: the charge is free, summed here"""
    return sum(ecoef(k, l, q) for q in range(5))


@lru_cache(maxsize=None)
def tleaf(e):
    return (-1) ** e * g(0, e - 1)


def legw(a, k):
    """weight of the order-k term of R^-1 on a leg tau_d(h^a): the leg enters the vertex with exponent d + k"""
    return (-1) ** k * g(a, k)


def symplectic_ok(nmax=20):
    """Sum_{m+n'=n, a-m=r} (-1)^n' g^a_m g^(3-a)_n' = delta_{n,0} for every r mod 5"""
    for n in range(nmax + 1):
        for r in range(5):
            s = Fr(0)
            for m in range(n + 1):
                a = (r + m) % 5
                s += (-1) ** (n - m) * g(a, m) * g(3 - a, n - m)
            if s != (1 if n == 0 else 0): return False
    return True


def edge_symmetric_ok(kmax=8):
    return all(ecoef(k, l, q) == ecoef(l, k, (-q) % 5) for k in range(kmax) for l in range(kmax) for q in range(5))


def edge_gf_ok(kmax=8):
    """(psi + psi') E = eta^-1 - R^-1 eta^-1 R^-1^T coefficientwise: E(m-1, n) + E(m, n-1) = Gh(m, n) in charge variables"""
    for m in range(kmax):
        for n in range(kmax):
            if m == 0 and n == 0: continue
            for q in range(5):
                # Gh(m, n) with first index b, charge of the E(m-1, n) term: q = b + (m-1) - 1
                b = (q + 2 - m) % 5; a = (b + m) % 5
                gh = -Fr(1, 5) * (-1) ** (m + n) * g(a, m) * g(3 - a, n)
                lhs = (ecoef(m - 1, n, q) if m >= 1 else 0) + (ecoef(m, n - 1, (b + m - 1) % 5) if n >= 1 else 0)
                if lhs != gh: return False
    return True


if __name__ == "__main__":
    print("rows loaded: 5 x %d rationals" % (NMAX + 1))
    print("symplectic through zhat^20:", symplectic_ok(20))
    print("edge symmetry E(k,l;q) = E(l,k;-q), k,l < 8:", edge_symmetric_ok())
    print("edge generating function identity, m,n < 8:", edge_gf_ok())
    print("T(2), T(3), T(4) =", tleaf(2), tleaf(3), tleaf(4))
    print("E(0,0;q) =", [str(ecoef(0, 0, q)) for q in range(5)])
