# -*- coding: utf-8 -*-
"""rf_sp.py -- REFEREE: formal stationary phase of the hatted Mellin-Barnes/Stirling integrand by the LAGRANGE-BUERMANN (Morse
lemma) method, independent of the worker's Wick-contraction code (pp_sp) and of (HPF).

Integrand (POLY s.5, re-derived by the referee): exp(-W^(X;u)/w) B^(X;w) dX, X = 1 + t,
  W^ = Sum_m (u^(m-1)/m) [1/X + X^(5m-1)/(5m-1)],  W^(1+t) - W^(1) = om_2 t^2 S(t),  om_2 = 5/(2(1-u)),
  B^ = (1 - u X^5)^(-1/2) exp( Sum_k b_k (u X w)^(2k-1) [5^(1-2k) - Sum_j (1 - zeta^j y)^(1-2k)]_(y^(5P) -> u^P X^(5P)) ),
  the zeta-sum computed in F_p (series inverse of (1 - zeta^j y)^(2k-1)); its non-multiples of y^5 are checked to vanish.
With y = t S(t)^(1/2):  [y^2r] g(t(y)) t'(y) = [t^2r] g(t) S(t)^(-r-1/2)  (residue substitution), hence
  Phi^SP := SP/sqrt(2 pi w/5) = Sum_r (2r-1)!! ((1-u) w/5)^r [t^2r] ( B_1(t) X(t,w) S(t)^(-r-1/2) ),
  B_1 = sqrt(1-u) (1 - u X^5)^(-1/2) = (1 - u((1+t)^5 - 1)/(1-u))^(-1/2),  X = exp(Stirling part).
Tests: (S1) [w^n] Phi^SP == ph_n^QRR (referee's large-radius polynomials, rf_qexp) in all u^0..u^(K-1), n <= N, two primes;
(S2) the model rows: SP of Y^m exp(-(1/Y + Y^4/4)/w) dY at Y = 1 (normalised by sqrt(2 pi w/5) e^(-5/(4w))) equals eps_b F_b(w),
m = c_b, eps = (+,-,-,+,+) (Theorem P step 4 with signs); control: the sign pattern of I_m alone (all +) fails.
Run: py -I -B rf_sp.py N
"""
import os, sys, time
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_core import P61, P89, CEXP, md, Ser, bernoulli, zeta5, model_rows
from rf_qexp import hatted_polys


def tmul(S, A, B, deg):
    out = [S.z() for _ in range(deg + 1)]
    for i, a in enumerate(A[:deg + 1]):
        if not any(a): continue
        for j in range(min(len(B), deg + 1 - i)):
            if any(B[j]): out[i + j] = S.add(out[i + j], S.mul(a, B[j]))
    return out


def tpow(S, A, alpha, deg):
    """A(t)^alpha for A_0 = 1 (u-series coefficients): n g_n = Sum_k (alpha k - (n-k)) A_k g_(n-k)"""
    p = S.p; al = md(alpha, p); g = [S.one()] + [S.z() for _ in range(deg)]
    for n in range(1, deg + 1):
        acc = S.z()
        for k in range(1, min(n, len(A) - 1) + 1):
            c = (al * k - (n - k)) % p
            if c and any(A[k]): acc = S.add(acc, S.mul(A[k], g[n - k]), c)
        g[n] = S.scale(acc, pow(n, -1, p))
    return g


def stirling_sigma(p, kmax, Pmax):
    """sig[k][P] = [y^(5P)] Sum_j (1 - zeta^j y)^(1-2k) in F_p; also checks that other powers vanish"""
    zt = zeta5(p); D = 5 * Pmax + 1; sig = {}
    for k in range(1, kmax + 1):
        tot = [0] * D
        for j in range(5):
            zj = pow(zt, j, p); poly = [1] + [0] * (D - 1)       # (1 - zj y)^(2k-1)
            for _ in range(2 * k - 1):
                poly = [(poly[i] - (zj * poly[i - 1] if i else 0)) % p for i in range(D)]
            inv = [pow(poly[0], -1, p)] + [0] * (D - 1)
            for n in range(1, D): inv[n] = -sum(poly[i] * inv[n - i] for i in range(1, n + 1)) % p * inv[0] % p
            tot = [(a + b) % p for a, b in zip(tot, inv)]
        assert all(tot[i] == 0 for i in range(D) if i % 5), "cyclotomic sum has non-multiples of y^5"
        sig[k] = [tot[5 * P] for P in range(Pmax + 1)]
    return sig


def phi_sp(p, N, K, ctrl=''):
    S = Ser(p, K); T = 2 * N; Bn = bernoulli(2 * N + 4)
    if ctrl == 'nostirling': Bn = [Fr(0)] * len(Bn)
    if ctrl == 'flipstirling': Bn = [-x for x in Bn]                 # R versus R^-1 convention
    one_m_inv = S.inv([1, p - 1] + [0] * (K - 2))
    def om(k):                                                       # [t^k] W^(1+t), u-series
        return [md((Fr((-1) ** k) + Fr(comb(5 * m - 1, k), 5 * m - 1)) / m, p) for m in range(1, K + 1)]
    om2 = om(2); assert om2 == S.scale(one_m_inv, md(Fr(5, 2), p))
    iom2 = S.inv(om2)
    Sser = [S.one()] + [S.mul(om(k + 2), iom2) for k in range(1, T + 1)]
    eps = [S.z() for _ in range(T + 1)]                              # u((1+t)^5 - 1)/(1-u)
    ufac = S.shift(one_m_inv, 1)
    for i in range(1, 6):
        if i <= T: eps[i] = S.scale(ufac, comb(5, i))
    B1 = tpow(S, [S.one()] + [S.scale(e, p - 1) for e in eps[1:]], Fr(-1, 2), T)
    kmax = (N + 1) // 2; sig = stirling_sigma(p, kmax, K)
    ell = {}
    for k in range(1, kmax + 1):
        bk = Bn[2 * k] / (2 * k * (2 * k - 1)); e = 2 * k - 1; lk = [S.z() for _ in range(T + 1)]
        for P in range(0, K - e):
            c = ((md(Fr(1, 5 ** e), p) if P == 0 else 0) - sig[k][P]) * md(bk, p) % p
            if not c: continue
            for i in range(0, min(e + 5 * P, T) + 1):
                lk[i][e + P] = (lk[i][e + P] + c * comb(e + 5 * P, i)) % p
        ell[k] = lk
    X = [[S.one()] + [S.z() for _ in range(T)]]                       # X_n(t) = [w^n] exp(Sum_k w^(2k-1) ell_k)
    for n in range(1, N + 1):
        deg = 2 * (N - n); acc = [S.z() for _ in range(deg + 1)]
        for k in range(1, kmax + 1):
            e = 2 * k - 1
            if e > n: break
            prod = tmul(S, ell[k], X[n - e], deg)
            acc = [S.add(a, b, e) for a, b in zip(acc, prod)]
        X.append([S.scale(a, pow(n, -1, p)) for a in acc])
    E = [tmul(S, B1, X[j], 2 * (N - j)) for j in range(N + 1)]
    Sinv = tpow(S, Sser, Fr(-1), T); Tcur = tpow(S, Sser, Fr(-1, 2), T)
    Tr = [Tcur[:1]]
    for r in range(1, N + 1):
        Tcur = tmul(S, Tcur, Sinv, T); Tr.append(Tcur[:2 * r + 1])
    af = S.scale([1, p - 1] + [0] * (K - 2), pow(5, -1, p)); afp = [S.one()]
    for r in range(1, N + 1): afp.append(S.mul(afp[-1], af))
    out = []
    for n in range(N + 1):
        tot = S.z()
        for j in range(n + 1):
            r = n - j; acc = S.z()
            for i in range(0, 2 * r + 1):
                if i < len(E[j]) and any(E[j][i]): acc = S.add(acc, S.mul(E[j][i], Tr[r][2 * r - i]))
            df = 1
            for odd in range(1, 2 * r, 2): df = df * odd % p
            tot = S.add(tot, S.mul(afp[r], acc), df)
        out.append(tot)
    return out


def model_sp(p, N):
    """SP of Y^m exp(-(1/Y + Y^4/4)/w) dY at Y = 1, normalised; dict m -> list over w^n"""
    S = Ser(p, 1); T = 2 * N
    phi = lambda k: md(Fr((-1) ** k) + Fr(comb(4, k), 4), p)
    i2 = pow(phi(2), -1, p); Sser = [[1]] + [[phi(k + 2) * i2 % p] for k in range(1, T + 1)]
    Sinv = tpow(S, Sser, Fr(-1), T); Tcur = tpow(S, Sser, Fr(-1, 2), T); Tr = [Tcur[:1]]
    for r in range(1, N + 1):
        Tcur = tmul(S, Tcur, Sinv, T); Tr.append(Tcur[:2 * r + 1])
    res = {}
    for m in (0, -1, 3, 2, 1):
        amp = tpow(S, [[1], [1]] + [[0]] * (T - 1), Fr(m), T)        # (1+t)^m
        ser = []
        for r in range(N + 1):
            acc = sum(amp[i][0] * Tr[r][2 * r - i][0] for i in range(2 * r + 1)) % p
            df = 1
            for odd in range(1, 2 * r, 2): df = df * odd % p
            ser.append(acc * df % p * pow(pow(5, r, p), -1, p) % p)
        res[m] = ser
    return res


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10; t0 = time.time(); ok = True
    K = 4 * N // 5 + 3
    print("rf_sp: Lagrange-Buermann stationary phase, w^0..w^%d, u^0..u^%d, primes 2^61-1, 2^89-1" % (N, K - 1))
    for p in (P61, P89):
        sp = phi_sp(p, N, K); polys, _ = hatted_polys(p, N)
        same = [n for n in range(N + 1) if sp[n] == (polys[n] + [0] * K)[:K]]
        poly = all(all(x == 0 for x in sp[n][4 * n // 5 + 1:]) for n in range(N + 1))
        print("p=%d: (S1) [w^n] Phi^SP == ph_n^QRR (all u-coefficients) for %d/%d orders; vanishing above u^floor(4n/5): %s  [%.1f s]"
              % (p, len(same), N + 1, poly, time.time() - t0))
        ok = ok and len(same) == N + 1 and poly
        for ct in ('nostirling', 'flipstirling'):
            sc = phi_sp(p, min(N, 10), 4 * min(N, 10) // 5 + 3, ct); Kc = 4 * min(N, 10) // 5 + 3
            bad = [n for n in range(len(sc)) if sc[n] != (polys[n] + [0] * Kc)[:Kc]]
            print("   CONTROL '%s' (w^0..w^%d): mismatching orders %s -> %s" % (ct, len(sc) - 1, bad[:6],
                  'fails as it should' if bad else 'UNEXPECTED'))
            ok = ok and bool(bad)
        F = model_rows(p, N); ms = model_sp(p, N); eps = (1, -1, -1, 1, 1)
        s2 = all(ms[CEXP[b]] == [(eps[b] * x) % p for x in F[b]] for b in range(5))
        ctrl = all(ms[CEXP[b]] == F[b] for b in range(5))
        print("   (S2) model: SP(Y^c_b e^(-f/w) dY) = eps_b F_b(w), eps = (+,-,-,+,+), b = 0..4, w^0..w^%d: %s; control eps = +1: %s"
              % (N, s2, "fails as it should" if not ctrl else "UNEXPECTED"))
        ok = ok and s2 and not ctrl
    print("RESULT rf_sp: %s  [%.1f s]" % ("PASS" if ok else "FAIL", time.time() - t0))
