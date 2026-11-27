# rf4_check_a.py -- referee check A (item 1, prop:mod-Rtr; lem:set-M0sympl). Exact over Q(zeta5). Usage: py -I -B rf4_check_a.py N
# A1 T-system closes; A2 symplecticity M_0(-t)^T M_0(t) = 1 and the referee's own route: Rhat(-t)^T eta^-1 Rhat(t) = 1;
# A3 DNOPS (2-31) by residue extraction (finite part at i = j) equals M_0(t)_ij, t = -1/(kappa lam), kappa = 1, 4, -1, 1/2;
# A4 (2-32) LHS = RHS (RHS = Laplace transform of dy = dp, computed by residue extraction, independent of A3);
# A5 the IBP identity (1) of GAPS_G4_2 with It_m from residues; controls that must fail.
import os, sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sys
from fractions import Fraction as F
from rf4_field import K, ZERO, ONE, smul, spow, dfact
from rf4_model import M0_coeffs, Rinv_dnops, W_local, residue, EPS

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
ok = True
Fs, closing, Rh, M = M0_coeffs(N)
print('A1 T-system closing relation through t^%d: %s; F_0 = %s' % (N + 1, closing, [str(c) for c in Fs[0][:5]]))
ok &= closing and Fs[0][1] == F(-3, 20) and Fs[0][2] == F(9, 800)

# A2: symplecticity, and eta-form for Rhat (eta^-1_ab = (1/5) delta(a+b = 3 mod 5))
bad = 0; bad_eta = 0
for n in range(N + 1):
    for i in range(5):
        for j in range(5):
            s = ZERO; s2 = ZERO
            for m in range(n + 1):
                sg = (-1) ** m
                for k in range(5):
                    s = s + M[m][k][i] * M[n - m][k][j] * sg
                for a in range(5):
                    s2 = s2 + Rh[m][a][i] * Rh[n - m][(3 - a) % 5][j] * F(sg, 5)
            target = ONE if (n == 0 and i == j) else ZERO
            bad += 0 if s == target else 1
            bad_eta += 0 if s2 == target else 1
print('A2 M0(-t)^T M0(t) = 1: %d failures; Rhat(-t)^T eta^-1 Rhat(t) = 1: %d failures (orders 0..%d)' % (bad, bad_eta, N))
ok &= bad == 0 and bad_eta == 0
# control: with eta^-1 replaced by the identity/5 the identity must fail
ctrl = 0
for i in range(5):
    for j in range(5):
        s2 = ZERO
        for a in range(5):
            s2 = s2 + Rh[0][a][i] * Rh[0][a][j] * F(1, 5)
        ctrl += 0 if s2 == (ONE if i == j else ZERO) else 1
print('   control (eta^-1 -> 1/5): %d failures at order 0 (must be > 0)' % ctrl)
ok &= ctrl > 0

for kappa in (F(1), F(4), F(-1), F(1, 2)):
    X = Rinv_dnops(kappa, N)
    mis = 0; ctrl_noT = 0; ctrl_sign = 0
    for n in range(N + 1):
        sc = (F(-1) / kappa) ** n                     # [lam^-n] of t^n, t = -1/(kappa lam)
        for j in range(5):
            for i in range(5):
                mis += 0 if X[n][j][i] == M[n][i][j] * sc else 1
                ctrl_noT += 0 if X[n][j][i] == M[n][j][i] * sc else 1
                ctrl_sign += 0 if X[n][j][i] == M[n][i][j] * sc * (-1) ** n else 1
    # A4 (2-32): Sum_k (R^-1)^i_k Delta_k^(-1/2) vs sqrt(lam/2pi) Int dy e^.. = L_i[dp]; both carry 1/sqrt(Lambda):
    # LHS*sqrt(L) = Sum_k X^i_k / c_k ; RHS*sqrt(L) = Sum_n (2n-1)!! L^-n lam^-n [w'^2n] (dp/dw'), w' = W
    Lam = -5 * kappa
    cs = [K.z(-j) for j in range(5)]
    bad32 = 0
    for i in range(5):
        Wt = W_local(kappa, cs[i], 2 * N + 6)
        h = [Wt[k + 1] / cs[i] for k in range(2 * N + 5)]
        for n in range(N + 1):
            lhs = ZERO
            for k in range(5):
                lhs = lhs + X[n][i][k] / cs[k]
            e = -2 * n - 1
            Wp = [c * cs[i] ** e for c in spow(h, e, 2 * N + 4)]
            rhs = residue(Wp, e) * F(dfact(2 * n - 1)) / Lam ** n   # Res_t W^(-2n-1) dt = [w'^2n](dp/dw')
            bad32 += 0 if lhs == rhs else 1
    print('A3 kappa=%s: X = M0(t)^T, t = -1/(kappa lam): %d mismatches of %d; controls: no transpose %d, '
          'sign-flipped t %d | A4 (2-32): %d mismatches' % (kappa, mis, 25 * (N + 1), ctrl_noT, ctrl_sign, bad32))
    ok &= mis == 0 and ctrl_noT > 0 and ctrl_sign > 0 and bad32 == 0

# A5: the IBP identity (1): R^-1^j_i = (1/5) Sum_k c_i^(k-1) It^(j)_(2-k), It_m := sqrt(L) L_j[p^m dp] (kappa = 1)
kappa = F(1); Lam = F(-5)
X = Rinv_dnops(kappa, N)
cs = [K.z(-j) for j in range(5)]
bad = 0
for j in range(5):
    Wt = W_local(kappa, cs[j], 2 * N + 6)
    h = [Wt[k + 1] / cs[j] for k in range(2 * N + 5)]
    It = {}
    for m in (-2, -1, 0, 1, 2):
        base = [cs[j], ONE] + [ZERO] * (2 * N + 3)
        pm = spow(base, m, 2 * N + 4)
        for n in range(N + 1):
            e = -2 * n - 1
            Wp = [c * cs[j] ** e for c in spow(h, e, 2 * N + 4)]
            It[(m, n)] = residue(smul(pm, Wp, 2 * N + 4), e) * F(dfact(2 * n - 1)) / Lam ** n
    for i in range(5):
        for n in range(N + 1):
            s = ZERO
            for k in range(5):
                s = s + cs[i] ** (k - 1) * It[(2 - k, n)] * F(1, 5)
            bad += 0 if s == X[n][j][i] else 1
print('A5 IBP form (1) of GAPS_G4_2 (kappa = 1): %d mismatches of %d' % (bad, 25 * (N + 1)))
ok &= bad == 0
print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
