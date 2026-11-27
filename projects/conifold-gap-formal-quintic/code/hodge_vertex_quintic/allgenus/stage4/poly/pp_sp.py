# -*- coding: utf-8 -*-
"""pp_sp.py -- Stage 4 / poly: formal stationary phase (SP) of the Mellin-Barnes summand of the formal quintic, in hatted variables.

Colour-0 I-function summand q^s Gamma(1+5s)/prod_j Gamma(1+s-zeta^j/z), s = sigma/z; Stirling Gamma(1+y) = sqrt(2 pi) y^(y+1/2) e^-y G*(y),
log G*(y) ~ Sum_k b_k y^(1-2k), b_k = B_2k/(2k(2k-1)).  With x = 1/sigma = s X (hatted), z = s^4 w, 3125 q = 1 - u:
    phase  -W^(X;u)/w,   W^(X;u) = Sum_(m>=1) (u^(m-1)/m) [1/X + X^(5m-1)/(5m-1)]    (critical point X = 1, W^''(1) = 5/(1-u))
    amplitude B^(X) = (1 - u X^5)^(-1/2) G*(5/(u X w)) / prod_j G*((1 - zeta^j s X)/(u X w))
        log B^ = -(1/2) log(1 - u X^5) + Sum_k b_k (u X w)^(2k-1) [5^(1-2k) - 5 Sum_(p>=0) C(2k-2+5p, 5p) (u X^5)^p].
Phi^SP(w, u) := SP_(X=1)[e^(-W^/w) B^ dX] / sqrt(2 pi w / 5)  (= 1 + O(w)).  Computed mod a prime: X = 1 + e tau, e = sqrt(w),
E = exp(Psi) as a series in e, Wick moments int tau^(2r) e^(-a tau^2/2) = (2r-1)!! a^-r sqrt(2 pi/a).
API: phi_sp(p, N, K) -> list over n <= N of u-series (length K) of [w^n] Phi^SP.
"""
from fractions import Fraction as Fr
from math import comb


def _of(x, p): x = Fr(x); return x.numerator % p * pow(x.denominator % p, p - 2, p) % p


def _bern(N):
    B = [Fr(0)] * (N + 1); B[0] = Fr(1)
    for m in range(1, N + 1): B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


def phi_sp(p, N, K):
    E2 = 2 * N; T = 3 * E2 + 2                                     # e-orders 0..2N, tau-degree bound
    Bn = _bern(2 * N + 2)
    # Psi[e][tdeg] = u-series (list length K), e = 1..2N
    Psi = [[None] * (T + 1) for _ in range(E2 + 1)]
    def add(e, d, useries):
        if e > E2 or d > T: return
        cur = Psi[e][d]
        if cur is None: Psi[e][d] = [x % p for x in useries]
        else: Psi[e][d] = [(a + b) % p for a, b in zip(cur, useries)]
    # (i) phase: -Sum_(k>=3) [t^k] W^ * e^(k-2) tau^k,  [t^k] W^ = Sum_m u^(m-1) [(-1)^k + C(5m-1,k)/(5m-1)]/m
    for k in range(3, E2 + 3):
        ser = [(-_of((Fr((-1) ** k) + Fr(comb(5 * m - 1, k), 5 * m - 1)) / m, p)) % p for m in range(1, K + 1)]
        add(k - 2, k, ser)
    # (ii) -(1/2) log(1 - u (1+t)^5) minus its t = 0 value: Sum_(P>=1) u^P [(1+t)^(5P) - 1]/(2P), t = e tau
    for P in range(1, K):
        for i in range(1, min(5 * P, E2) + 1):
            ser = [0] * K; ser[P] = _of(Fr(comb(5 * P, i), 2 * P), p); add(i, i, ser)
    # (iii) Stirling: Sum_k b_k u^(2k-1) w^(2k-1) (1+t)^(2k-1) [5^(1-2k) - 5 Sum_P C(2k-2+5P,5P) u^P (1+t)^(5P)]
    for k in range(1, N + 2):
        bk = Bn[2 * k] / (2 * k * (2 * k - 1)); ew = 2 * (2 * k - 1)
        if ew > E2 or 2 * k - 1 >= K: continue
        for P in range(0, K - (2 * k - 1)):
            coef = (Fr(1, 5 ** (2 * k - 1)) if P == 0 else 0) - 5 * comb(2 * k - 2 + 5 * P, 5 * P)
            deg = 2 * k - 1 + 5 * P
            for i in range(0, min(deg, E2 - ew) + 1):
                ser = [0] * K; ser[2 * k - 1 + P] = _of(bk * coef * comb(deg, i), p); add(ew + i, i, ser)
    # E = exp(Psi):  e E_e = Sum_j j Psi_j E_(e-j)
    def smul(a, b):
        out = [0] * K
        for i, x in enumerate(a):
            if x:
                for j in range(K - i): out[i + j] = (out[i + j] + x * b[j]) % p
        return out
    Ee = [{0: [1] + [0] * (K - 1)}]
    for e in range(1, E2 + 1):
        acc = {}
        for j in range(1, e + 1):
            for d1, s1 in enumerate(Psi[j]):
                if s1 is None or not any(s1): continue
                for d2, s2 in Ee[e - j].items():
                    s = smul(s1, s2); d = d1 + d2
                    if d in acc: acc[d] = [(x + j * y) % p for x, y in zip(acc[d], s)]
                    else: acc[d] = [(j * y) % p for y in s]
        ie = pow(e, p - 2, p)
        Ee.append({d: [x * ie % p for x in s] for d, s in acc.items()})
    # Wick: [w^n] = Sum_r (2r-1)!! ((1-u)/5)^r [tau^(2r)] E_(2n)
    i5 = pow(5, p - 2, p); oneminus = [1, p - 1] + [0] * (K - 2)
    powa = [[1] + [0] * (K - 1)]
    for r in range(1, 3 * N + 2): powa.append([x * i5 % p for x in smul(powa[-1], oneminus)])
    out = []
    for n in range(N + 1):
        tot = [0] * K
        for d, s in Ee[2 * n].items():
            if d % 2: assert not any(s), "odd moment"; continue
            r = d // 2; df = 1
            for k in range(1, 2 * r, 2): df = df * k % p
            tot = [(x + df * y) % p for x, y in zip(tot, smul(powa[r], s))]
        out.append(tot)
    return out
