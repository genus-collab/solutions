# -*- coding: utf-8 -*-
"""pp_res.py -- Stage 4 / poly: the formal RESIDUE at x = 0 of the Mellin-Barnes/Stirling integrand versus the quintic's vanishing
solution, exact over Q, frame (0, 0).

Integrand (original variables, x = 1/sigma):  exp(-W(x;u)/z) at(x,z) dx,  W = Lam/x + G(x),  Lam = -log(1-u),
G(x) = Sum_(m>=1) x^(5m-1)/(m(5m-1)),  at(x,z) = (1-x^5)^(-1/2) G*(5/(xz)) / prod_j G*((1-zeta^j x)/(xz))  (at(0,z) = 1).
V(u,z) := Res_(x=0) = Sum_(k>=1) (-Lam/z)^k/k! [x^(k-1)] (e^(-G/z) at)   (finitely many z-powers at each u-order).
Vanishing solution (gauge s.2, Lemma Q): y^(0) = -e_1, y^(n+1)_a = Int_0^u b_a y^(n)_(a+1), b = (t_c', -1/(u w0^2 (1-u)^3 t_c'^2), t_c',
-w0/(1-u), -w0/(1-u)), w0 = f0, t_c = g/f0, g = f1 + (7/10) f2;  psi_C := w0 (Y^C_Q)_0 = w0 Sum_n z^-n y^(n)_0.
Prediction (POLY.md s.5): V = psi_C exactly (all positive powers of z produced by the Stirling factors cancel).
Run: py -I -B pp_res.py N   (u^1..u^N)
"""
import os, sys, time
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from pp_series import s_mul, s_inv, s_deriv, _pf, periods
from pp_core import bernoulli_list


def lmul(A, B):                                     # Laurent polys in z: dict power -> Fr
    out = {}
    for i, x in A.items():
        for j, y in B.items(): out[i + j] = out.get(i + j, 0) + x * y
    return {k: v for k, v in out.items() if v}


def residue_V(N):
    """V as dict u-power -> Laurent dict in z, u^1..u^N"""
    Bn = bernoulli_list(2 * N + 4)
    # log(e^(-G/z) at) as a series in x (coeffs Laurent in z), degree <= N
    L = [dict() for _ in range(N + 1)]
    for m in range(1, N // 5 + 2):
        if 5 * m - 1 <= N: L[5 * m - 1][-1] = L[5 * m - 1].get(-1, 0) - Fr(1, m * (5 * m - 1))
    for P in range(1, N // 5 + 1): L[5 * P][0] = L[5 * P].get(0, 0) + Fr(1, 2 * P)       # -(1/2) log(1-x^5)
    for k in range(1, N + 2):
        bk = Bn[2 * k] / (2 * k * (2 * k - 1)); e = 2 * k - 1
        for P in range(0, N + 1):
            deg = e + 5 * P
            if deg > N: break
            c = (Fr(1, 5 ** e) if P == 0 else 0) - 5 * comb(2 * k - 2 + 5 * P, 5 * P)
            L[deg][e] = L[deg].get(e, 0) + bk * c
    E = [dict() for _ in range(N + 1)]; E[0] = {0: Fr(1)}                                 # exp
    for n in range(1, N + 1):
        acc = {}
        for k in range(1, n + 1):
            if L[k] and E[n - k]:
                t = lmul(L[k], E[n - k])
                for a, b in t.items(): acc[a] = acc.get(a, 0) + k * b
        E[n] = {a: b / n for a, b in acc.items() if b}
    Lam = [Fr(0)] + [Fr(1, m) for m in range(1, N + 1)]
    V = {}; Lk = [Fr(1)] + [Fr(0)] * N; fact = 1
    for k in range(1, N + 1):
        Lk = s_mul(Lk, Lam, N + 1); fact *= k
        for i in range(N + 1):
            if Lk[i]:
                for zp, c in E[k - 1].items():
                    V.setdefault(i, {}); V[i][zp - k] = V[i].get(zp - k, 0) + (-1) ** k * Lk[i] * c / fact
    return {i: {a: b for a, b in d.items() if b} for i, d in V.items()}


def vanishing(N):
    K = N + 2; f0, f2 = periods(K)
    f1 = [Fr(0), Fr(1)] + [Fr(0)] * (K - 2)                                              # f1 = u + O(u^3): solve PF
    for n in range(3, K):
        f1[n] = Fr(0); c = _pf(f1 + [Fr(0)] * 3, K + 3)[n - 3]; f1[n] = -c / (n * (n - 1) ** 2 * (n - 2))
    g = [f1[i] + Fr(7, 10) * f2[i] for i in range(K)]
    assert all(x == 0 for x in _pf(g + [Fr(0)] * 3, K + 3)[:K - 3])
    tc = s_mul(g, s_inv(f0, K), K); tcp = s_deriv(tc, K); w0 = f0
    om = s_inv([Fr(1), Fr(-1)] + [Fr(0)] * (K - 2), K)
    b = [tcp, None, tcp, [-x for x in s_mul(w0, om, K)], [-x for x in s_mul(w0, om, K)]]
    one_m3 = [Fr(comb(3, i)) * (-1) ** i if i <= 3 else Fr(0) for i in range(K)]
    b1u = s_mul([-x for x in s_inv(s_mul(s_mul(w0, w0, K), s_mul(tcp, tcp, K), K), K)], s_inv(one_m3, K), K)  # u b_1
    y = [[Fr(0)] * K for _ in range(5)]; y[1][0] = Fr(-1); out = {0: [Fr(0)] * K}
    ys = [y]
    for n in range(1, 5 * N // 4 + 6):
        prev = ys[-1]; new = []
        for a in range(5):
            if a == 1: integ = s_mul(b1u, [prev[2][i + 1] if i + 1 < K else Fr(0) for i in range(K)], K)
            else: integ = s_mul(b[a], prev[(a + 1) % 5], K)
            new.append([Fr(0)] + [integ[i - 1] / i for i in range(1, K)])
        ys.append(new)
    psi = {}
    for n, yn in enumerate(ys):
        t = s_mul(w0, yn[0], K)
        for i in range(K):
            if t[i]: psi.setdefault(i, {}); psi[i][-n] = psi[i].get(-n, 0) + t[i]
    return psi


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10; t0 = time.time()
    V = residue_V(N); psi = vanishing(N)
    print("pp_res: u^1..u^%d, exact over Q, frame (0,0)" % N); ok = True
    for i in range(1, N + 1):
        v = V.get(i, {}); s = psi.get(i, {}); pos = [k for k in v if k > 0]
        same = (v == s); ok = ok and same and not pos
        print("u^%d: V == w0 (Y^C_Q)_0: %s; positive z-powers in V: %s; z-powers %s..%s"
              % (i, same, pos if pos else "none", min(v) if v else "-", max(v) if v else "-"))
    print("RESULT residue = vanishing solution (no Stirling remnant) through u^%d: %s  [%.1f s]" % (N, "PASS" if ok else "FAIL", time.time() - t0))
