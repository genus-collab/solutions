# -*- coding: utf-8 -*-
"""rf_flat.py -- referee's flat layer: exact correlators (two primes + CRT + rational reconstruction), the calibration S, the
vanishing covectors Pi_i and the K-generators, all re-derived here (no import from allgenus/model/ms_*).

Calibration (exact over Q, v = 1 after homogeneity): z d/dv S = (H *) S, S = 1 + S_1/z + ..., homogeneous, integration constants 0;
column b != 1 is log-free:  S_k H^b = s(k,b) v^x(k,b) H^(b+k),  x(k+1,b) = x(k,b) + 1 - [b+k = 1 mod 5],  s(k+1,b) = s(k,b) beta(b+k)/x(k+1,b).
Pi_i(tau_k(H^b)) = (1/5) (-1)^(k-i) eta(S_(k-i) H^2, H^b)  (= (1/5)[z^i] eta(S(-z)H^2, H^b z^k));  PiC = Pi_0.
K-generators: ('H',) = tau_0(H) and kappa(X, k) = [S(z) X z^k]_+ = Sum_j S_j X z^(k-j), X in {1, H^2, H^3, H^4}.
Weight of tau_k(H^a): w = (4 - c_a - 4k)/5; v-exponent of a genus-h correlator: e = 2 - 2h - Sum w.
"""
import os, sys
from fractions import Fraction as Fr
from functools import lru_cache
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_core as rc
import rf_gsum as gsm

C = rc.C_EXP; N = rc.N_SIGN
BETA = (Fr(1), Fr(-1), Fr(1), Fr(-1), Fr(-1))          # beta_b at v = 1 (the -1/v sits at b = 1)
_M = {}


def models(K=16):
    if K not in _M: _M[K] = (rc.Model(rc.P61, K), rc.Model(rc.P89, K))
    return _M[K]


def ratrec(a, Mod):
    """rational reconstruction of a mod Mod with |num|, den < sqrt(Mod/2)"""
    import math
    B = math.isqrt(Mod // 2); r0, r1 = Mod, a % Mod; s0, s1 = 0, 1
    while r1 > B:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B: raise ValueError("no reconstruction")
    return Fr(r1, s1)


def corrQ(h, ins, K=16):
    """exact <prod tau_k(H^a)>_h at v = 1 (ins = ((a,k),..)) via 2^61-1, 2^89-1, CRT"""
    m1, m2 = models(K)
    x1 = gsm.corr(m1, h, ins); x2 = gsm.corr(m2, h, ins)
    Mod = m1.p * m2.p
    x = (x1 + m1.p * ((x2 - x1) * pow(m1.p, m2.p - 2, m2.p) % m2.p)) % Mod
    return ratrec(x, Mod)


def weight(a, k): return Fr(4 - C[a % 5] - 4 * k, 5)


def vexp(h, ins):
    e = Fr(2 - 2 * h) - sum(weight(a, k) for a, k in ins)
    return e


def charge_ok(ins): return sum(a + k - 1 for a, k in ins) % 5 == 0


@lru_cache(maxsize=None)
def S(k, b):
    """S_k H^b = (s, x, target index): coefficient at v = 1, power of v, b + k mod 5 (b != 1)"""
    b %= 5; assert b != 1
    if k == 0: return Fr(1), 0, b
    s, x, t = S(k - 1, b)
    beta_pow = -1 if t == 1 else 0
    nx = x + beta_pow + 1; assert nx != 0
    return s * BETA[t] / nx, nx, (t + 1) % 5


def eta(a, b):
    return 5 * N[a % 5] * N[b % 5] if (a + b) % 5 == 3 else 0


def pi(i, a, k):
    """Pi_i(tau_k(H^a)) at v = 1 (coefficient); 0 if k < i"""
    if k < i: return Fr(0)
    m = k - i; s, x, t = S(m, 2)
    return Fr(1, 5) * (-1) ** m * s * eta(t, a)


def kappa(b, k):
    """kappa(H^b, k) = [S(z) H^b z^k]_+ at v = 1: list of (coef, (a, psi-degree))"""
    out = []
    for j in range(k + 1):
        s, x, t = S(j, b); out.append((s, (t, k - j)))
    return out


def in_K(a, k):
    return all(pi(i, a, k) == 0 for i in range(1, k + 1))


if __name__ == "__main__":
    # S symplectic: Sum_{j} (-1)^j eta(S_j H^a, S_(m-j) H^b) = 0 for m >= 1 (a, b != 1)
    bad = 0
    for a in (0, 2, 3, 4):
        for b in (0, 2, 3, 4):
            for m in range(1, 16):
                tot = Fr(0)
                for j in range(m + 1):
                    s1, x1, t1 = S(j, a); s2, x2, t2 = S(m - j, b); tot += (-1) ** j * s1 * s2 * eta(t1, t2)
                bad += tot != 0
    print("S symplectic on the log-free columns through z^-15:", bad == 0)
    bad = sum(1 for b in (0, 2, 3, 4) for k in range(14) for i in range(0, k + 1) if sum(c * pi(i, a, kk) for c, (a, kk) in kappa(b, k)) != 0)
    print("Pi_i(kappa(X,k)) = 0 for all i >= 0, k < 14:", bad == 0)
    print("PiC on the vanishing tower:", [(a, k, str(pi(0, a, k))) for k in range(8) for a in range(5) if pi(0, a, k)])
    print("basis insertions in K (k <= 8):", [(a, k) for k in range(9) for a in range(5) if in_K(a, k)])
