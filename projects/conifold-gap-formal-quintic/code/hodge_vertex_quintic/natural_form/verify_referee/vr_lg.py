# -*- coding: utf-8 -*-
"""vr_lg.py -- VERIFIER: exact (Fraction) cross-checks that use NO F_p arithmetic, no QRR constants and no code of the referee.

 (A) identity row of the symbol from the u = 0 slice of the hatted Picard-Fuchs equation (verifier's own derivation):
         Phi_j = A_{-j} Phi_{j-1},  (A_e F)[n] = F[n] + ((4(n-1) - e)/5) F[n-1]   (colour 0, u = 0),
         Phi_5 - Phi_0 - 3 zhat Phi_4 = 0,   coefficient of zhat^(n+1):  4n ph_n + (terms in ph_0..ph_(n-1)) = 0.
     (no resonance at u = 0: the u^0 coefficient of ph_n is fixed for every n >= 1, no normalisation constant enters.)
 (B) stationary phase of Int exp((x^4 + u/x)/z) dx by LAGRANGE INVERSION (not by Gaussian moments of the exponent expansion):
     x = x_c y, phase (u/x_c)(y^4/4 + 1/y) = x_c^4 (y^4 + 4/y);  W(y) = y^4 + 4/y = 5 + 10 t^2 + 5 t^4 + Sum_{k>=5} 4 (-1)^k t^k, y = 1+t;
     put 10 T^2 = W - 5, t = Sum_j b_j T^j;  Int e^{W/eps} dy / (leading) = Sum_m (2m+1) b_(2m+1) (2m-1)!! (-eps/20)^m.
     With eps = -4 zhat: f_m = (2m+1) b_(2m+1) (2m-1)!! / 5^m.
 (C) LG differential equation re-derived: with I_k = Int x^k e^{W/z} dx, z d/du I_k = I_(k-1) and (k+1) z I_k + 4 I_(k+4) - u I_(k-1) = 0;
     k = -4:  u (z d/du)^5 I + 3 z (z d/du)^4 I - 4 I = 0;  indicial polynomial at u = 0: r (r-1)^2 (r-2) (r-3).  (symbolic check below)
 (D) tau_n = [u^n] 1/Sum 4u^n/(5n+4), n <= 6.
 (E) hand derivation of G_1 from the G-equation in a general conifold frame (alpha): x_m = Ghat_1(m, 1-m),
         x_(m+1) - x_m = t_(-m)/5,  x_3 = 0,   t = (1 + tau_1) (0,1,2,-2,-1)/5 - b,  b = (0, alpha, alpha + 2 gamma - 1 - 2 alpha .. )
     (see the formula in the code), gamma = 7/10.
Run: py -I -B vr_lg.py
"""
import json, os
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))


def identity_row(N):
    ph = [Fr(1)]
    def A(e, F):
        return [F[n] + (Fr(4 * (n - 1) - e, 5) * F[n - 1] if n else 0) for n in range(len(F))]
    for n in range(1, N + 1):
        Phi = [ph + [Fr(0), Fr(0)]]
        for j in range(1, 6): Phi.append(A(-j, Phi[-1]))
        r = Phi[5][n + 1] - Phi[0][n + 1] - 3 * Phi[4][n]             # with ph_n = 0
        ph.append(-r / (4 * n))
        Phi = [ph + [Fr(0)]]
        for j in range(1, 6): Phi.append(A(-j, Phi[-1]))
        assert Phi[5][n + 1] - Phi[0][n + 1] - 3 * Phi[4][n] == 0
    return ph


def smul(a, b, N):
    out = [Fr(0)] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(min(len(b), N + 1 - i)): out[i + j] += x * b[j]
    return out


def lagrange(N):
    """f_m, m = 0..N, by inversion of T = t sqrt(1 + t^2/2 + Sum_{k>=5} (2/5)(-1)^k t^(k-2))"""
    M = 2 * N + 2
    h = [Fr(0)] * (M + 1); h[0] = Fr(1); h[2] = Fr(1, 2)
    for k in range(5, M + 3):
        if k - 2 <= M: h[k - 2] += Fr(2, 5) * (-1) ** k
    r = [Fr(0)] * (M + 1); r[0] = Fr(1)                               # r = sqrt(h): r_n = (h_n - Sum_{0<i<n} r_i r_(n-i))/2
    for n in range(1, M + 1): r[n] = (h[n] - sum(r[i] * r[n - i] for i in range(1, n))) / 2
    assert smul(r, r, M) == h
    Tser = [Fr(0)] + r[:M]                                           # T(t) = t r(t), mod t^(M+1)
    # inversion t(T) = Sum b_j T^j: fixed point  t = T - (T(t) - t)  coefficient by coefficient
    b = [Fr(0), Fr(1)] + [Fr(0)] * (M - 1)
    for it in range(M):
        pw = [Fr(1)] + [Fr(0)] * M; comp = [Fr(0)] * (M + 1)          # comp = T(t(T))
        for j in range(1, M + 1):
            pw = smul(pw, b, M)
            if Tser[j]: comp = [c + Tser[j] * x for c, x in zip(comp, pw)]
        err = [comp[j] - (1 if j == 1 else 0) for j in range(M + 1)]
        if not any(err): break
        b = [b[j] - err[j] for j in range(M + 1)]
    assert not any(err), "inversion did not converge"
    f = []
    for m in range(N + 1):
        df = 1
        for q in range(1, 2 * m, 2): df *= q
        f.append((2 * m + 1) * b[2 * m + 1] * df / Fr(5) ** m)
    return f


def lg_ode_check():
    """the operator u d^5 + 3 d^4 (z = 1) on u^r: indicial polynomial; and the relation (k+1) I_k + 4 I_(k+4) - u I_(k-1) = 0 at k = -4"""
    def indicial(r): return r * (r - 1) * (r - 2) * (r - 3) * (r - 4) + 3 * r * (r - 1) * (r - 2) * (r - 3)
    roots = [r for r in range(-5, 10) if indicial(r) == 0]
    dbl = [r for r in roots if (indicial(Fr(r) + Fr(1, 10 ** 6)) / Fr(1, 10 ** 6)) ** 2 < Fr(1, 10 ** 6)]
    return roots, dbl


def g1_hand(alpha, gamma=Fr(7, 10)):
    tau1 = Fr(-4, 9); mb = (0, 1, 2, -2, -1)
    # theta_u log of Psi_a at u^1: a = 0: -alpha; a = 1: -alpha - c1; a = 2: +alpha + c1; a = 3: +alpha; a = 4: 0,  C1 = -(1 + c1 u), c1 = 2 gamma - 1 - 2 alpha
    c1 = 2 * gamma - 1 - 2 * alpha
    bflat = (-alpha, -alpha - c1, alpha + c1, alpha, Fr(0))          # index a;  m = -(a+1) mod 5
    b = [bflat[(-m - 1) % 5] for m in range(5)]
    t = [(1 + tau1) * Fr(mb[m], 5) - b[m] for m in range(5)]
    assert sum(t) == 0
    x = {3: Fr(0)}
    x[4] = x[3] + t[(-3) % 5] / 5; x[0] = x[4] + t[(-4) % 5] / 5; x[1] = x[0] + t[0] / 5; x[2] = x[1] + t[(-1) % 5] / 5
    assert x[2] + t[(-2) % 5] / 5 == x[3]
    return {"(0,1)": x[0], "(1,0)": x[1], "(2,4)": x[2], "(3,3)": x[3], "(4,2)": x[4]}


if __name__ == "__main__":
    N = 12
    fa = identity_row(N); fb = lagrange(N)
    print("(A) identity row from the u = 0 hatted PF recursion:", ", ".join(str(x) for x in fa[:7]))
    print("(B) Lagrange-inversion stationary phase (eps = -4 zhat):", ", ".join(str(x) for x in fb[:7]))
    print("    (A) == (B) for zhat^0..zhat^%d: %s" % (N, fa == fb))
    d = json.load(open(os.path.join(os.path.dirname(HERE), "referee", "G_fourier_NZ16_K12.json")))
    ref = [Fr(x) for x in d['identity_row']]
    print("    (A) == referee's identity row (rational reconstruction at 2^521-1) for zhat^0..zhat^%d: %s" % (N, fa == ref[:N + 1]))
    print("    in eps = -4 zhat: 1, %s, %s, %s  (hand check quoted: 3/80, 9/12800)" % (fa[1] / -4, fa[2] / 16, fa[3] / -64))
    print("(C) indicial roots of u d^5 + 3 d^4 at u = 0: %s, double: %s" % lg_ode_check())
    F = [Fr(4, 5 * k + 4) for k in range(8)]; tau = [Fr(1)]
    for k in range(1, 8): tau.append(-sum(F[j] * tau[k - j] for j in range(1, k + 1)))
    print("(D) tau_1..tau_6:", ", ".join(str(x) for x in tau[1:7]))
    for al in (Fr(0), Fr(1, 3), Fr(8, 45)):
        print("(E) hand-derived G_1 Fourier components, frame alpha = %s: %s" % (al, {k: str(v) for k, v in g1_hand(al).items()}))
    # (F) counts of testable (component, k) pairs for rule W
    Wt = (0, 4, 3, 2, 1)
    for NZ, K in ((16, 12), (24, 18)):
        tot = 0; test = 0
        for aa in range(5):
            for bb in range(5):
                for k in range(1, NZ + 1):
                    if (aa + bb + 2 + k) % 5: continue
                    tot += 1; test += max(1, -((-(4 * k - 2 + Wt[aa] + Wt[bb])) // 5)) < K
        print("(F) zhat <= %d, u < %d: Z/5-allowed (component, k) pairs %d, of which testable (bound < K): %d" % (NZ, K, tot, test))
