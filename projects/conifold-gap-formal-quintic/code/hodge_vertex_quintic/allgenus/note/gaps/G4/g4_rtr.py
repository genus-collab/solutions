# g4_rtr.py -- evidence for prop:mod-Rtr (gap G4 item 1). Exact over Q(zeta5). Usage: py -I -B g4_rtr.py N
# (1) M_0 from the T-system (pivot 4n), colour-colour form M_0 = Psi_0^-1 Rhat.
# (2) X^j_i(zeta) := DNOPS (2-31) on x = kappa W_v, y = p, read with the finite part at i = j (local expansion in w_j).
# (3) checks: X(zeta) = M_0(t)^T, t = -1/(kappa s^4 zeta) [the claim]; the integration-by-parts form
#     X^j_i = (1/5) Sum_k c_i^(k-1) I~^(j)_(2-k); symplecticity of M_0; (2-32). Controls that must fail.
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g4_field import Qz, ZERO, ONE, smul, sinv, ssqrt1, srevert, scompose, dfact

EPS = (1, -1, -1, 1, 1)
NU = (0, -1, 3, 2, 1)


def tsystem(N):
    """colour-0 solution F_b(t) of F1 = T_-1 F0, F2 = -T_-2 F1, F3 = T_2 F2, F4 = -T_1 F3, F0 = -T_0 F4, F0(0) = 1"""
    def T(k, F):
        out = [-x for x in F] + [Fr(0)]
        for n in range(1, len(out)):
            out[n] += Fr(k - 4 * (n - 1), 5) * F[n - 1]
        return out[:N + 2]

    def chain(F0):
        F1 = T(-1, F0); F2 = [-x for x in T(-2, F1)]; F3 = T(2, F2); F4 = [-x for x in T(1, F3)]
        return [F0, F1, F2, F3, F4], [-x for x in T(0, F4)]
    f = [Fr(1)] + [Fr(0)] * (N + 1)
    for n in range(1, N + 1):
        f[n] = Fr(0)
        Fs, back = chain(f[:])
        r = back[n + 1] - f[n + 1]              # residual at order n+1 with f_n = 0 (only f_{n-4..n-1} enter)
        f[n] = -r / (4 * n)                     # (back - F0) from f_n t^n starts with +4n f_n t^(n+1): pivot 4n
    Fs, back = chain(f[:])
    for n in range(N + 1):
        assert back[n] == f[n], ('T-system does not close', n)
    return [F[:N + 1] for F in Fs]


def M0_colour(N):
    """M0[n][i][j] = [t^n] M_0(t)_{ij} in Q(zeta); Rh[n][b][j] = [t^n] Rhat_b^(j)"""
    F = tsystem(N)
    Rh = [[[Qz.zeta(j * (b + 1) - j * n) * F[b][n] for j in range(5)] for b in range(5)] for n in range(N + 1)]
    M = [[[ZERO] * 5 for _ in range(5)] for _ in range(N + 1)]
    for n in range(N + 1):
        for i in range(5):
            for j in range(5):
                s = ZERO
                for b in range(5):
                    s = s + Qz.zeta(-i * (b + 1)) * Rh[n][b][j] * Fr(EPS[b], 5)
                M[n][i][j] = s
    return F, Rh, M


def local_data(kappa, sig, j, N):
    """at c_j = zeta^-j sig on x = kappa (p^4/4 + v/p), v = sig^5: series t(wt) = p - c_j, with x - x(c_j) = (5 kappa/2) wt^2"""
    c = Qz.zeta(-j) * sig; v = Fr(sig) ** 5
    from math import comb
    xk = [ZERO] * (N + 4)
    for k in range(2, N + 4):
        a = Qz.of(Fr(comb(4, k), 4)) * (c ** (4 - k) if k <= 4 else ZERO) if k <= 4 else ZERO
        b = (Qz.of(v) / c) * Qz.of((-1) ** k) * (ONE / c) ** k
        xk[k] = (a + b) * kappa
    x2 = xk[2]
    u = [ONE] + [xk[m + 2] / x2 for m in range(1, N + 2)]
    S = ssqrt1(u, N + 1)
    wt = [ZERO] + [c * s for s in S[:N + 1]]          # wt as a series in t
    t_of_wt = srevert(wt, N + 1)
    return c, t_of_wt


def __pow(self, k):
    r = ONE
    for _ in range(k): r = r * self
    return r
Qz.__pow__ = __pow


def X_matrix(kappa, sig, N):
    """X[n][j][i] = [zeta^-n] of (2-31), n = 0..N"""
    Lam = Fr(-5 * kappa)
    M = 2 * N + 2
    loc = [local_data(kappa, sig, j, M) for j in range(5)]
    cs = [l[0] for l in loc]
    X = [[[ZERO] * 5 for _ in range(5)] for _ in range(N + 1)]
    for j in range(5):
        c_j, t = loc[j]
        for i in range(5):
            ci = cs[i]
            if i == j:                                  # xit = -1/(c_i t(wt)) = -(1/c_i) wt^-1 (t/wt)^-1
                q = sinv(t[1:M + 1], M - 1)
                a = {m - 1: -(q[m] / ci) for m in range(M)}
            else:
                ser = [c_j - ci + t[0]] + t[1:M + 1]
                q = sinv(ser, M)
                a = {m: -(q[m] / ci) for m in range(M + 1)}
            X[0][j][i] = ONE if i == j else ZERO
            for n in range(1, N + 1):
                X[n][j][i] = -(a.get(2 * n - 1, ZERO) * Fr(dfact(2 * n - 1)) / Lam ** n)
    return X, loc, Lam


def Itilde(loc, Lam, j, m, N):
    """[zeta^-n] sqrt(Lam) L_j[p^m dp] = Sum_n (2n-1)!! Lam^-n [wt^(2n)] (p^m dp/dwt)"""
    c_j, t = loc[j]; M = 2 * N + 2
    p = [c_j + t[0]] + t[1:M + 1]
    dp = [t[k + 1] * (k + 1) for k in range(M)]
    if m >= 0:
        pm = [ONE] + [ZERO] * M
        for _ in range(m): pm = smul(pm, p, M)
    else:
        pinv = sinv(p, M); pm = [ONE] + [ZERO] * M
        for _ in range(-m): pm = smul(pm, pinv, M)
    f = smul(pm, dp, M)
    return [f[2 * n] * Fr(dfact(2 * n - 1)) / Lam ** n for n in range(N + 1)]


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    F, Rh, M = M0_colour(N)
    print('T-system colour 0, F_0 =', [str(x) for x in F[0][:6]])
    assert F[0][1] == Fr(-3, 20) and F[0][2] == Fr(9, 800)
    ok = True
    # symplecticity M_0(-t)^T M_0(t) = 1
    bad = 0
    for n in range(N + 1):
        for i in range(5):
            for j in range(5):
                s = ZERO
                for m in range(n + 1):
                    for k in range(5):
                        s = s + M[m][k][i] * M[n - m][k][j] * ((-1) ** m)
                if not (s == (ONE if (n == 0 and i == j) else ZERO)): bad += 1
    print('symplecticity M0(-t)^T M0(t) = 1 through t^%d: failures %d' % (N, bad)); ok &= bad == 0
    for kappa, sig in ((1, 1), (4, 1), (-1, 1), (1, 2), (4, 2)):
        X, loc, Lam = X_matrix(kappa, sig, N)
        sc = Fr(-1, kappa * sig ** 4)                   # t = sc / zeta
        res = {'claim X = M0(t)^T': 0, 'ctrl X = M0(t)': 0, 'ctrl X = M0(-t)^T': 0, 'ctrl kappa ignored': 0}
        for n in range(N + 1):
            for j in range(5):
                for i in range(5):
                    if not (X[n][j][i] == M[n][i][j] * sc ** n): res['claim X = M0(t)^T'] += 1
                    if not (X[n][j][i] == M[n][j][i] * sc ** n): res['ctrl X = M0(t)'] += 1
                    if not (X[n][j][i] == M[n][i][j] * (-sc) ** n): res['ctrl X = M0(-t)^T'] += 1
                    if not (X[n][j][i] == M[n][i][j] * Fr(-1, sig ** 4) ** n): res['ctrl kappa ignored'] += 1
        ibp = 0
        I = {(j, m): Itilde(loc, Lam, j, m, N) for j in range(5) for m in (-2, -1, 0, 1, 2)}
        cs = [l[0] for l in loc]
        for n in range(N + 1):
            for j in range(5):
                for i in range(5):
                    s = ZERO
                    for k in range(5):
                        s = s + (cs[i] ** k / cs[i]) * I[(j, 2 - k)][n] * Fr(1, 5)
                    if not (s == X[n][j][i]): ibp += 1
        u32 = 0
        for n in range(N + 1):
            for i in range(5):
                lhs = ZERO
                for k in range(5): lhs = lhs + X[n][i][k] / cs[k]
                rhs = Rh[n][0][i] * sc ** n * Fr(1, sig)
                if not (lhs == rhs) or not (lhs == I[(i, 0)][n]): u32 += 1
        print('kappa=%2d sig=%d:' % (kappa, sig), ', '.join('%s: %d mismatches' % kv for kv in res.items()),
              '| IBP form: %d | (2-32): %d' % (ibp, u32))
        ok &= res['claim X = M0(t)^T'] == 0 and ibp == 0 and u32 == 0
        ok &= res['ctrl X = M0(t)'] > 0 and res['ctrl X = M0(-t)^T'] > 0
        if kappa != 1: ok &= res['ctrl kappa ignored'] > 0
    print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')


main()
