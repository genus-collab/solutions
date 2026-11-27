# rf4_model.py -- referee's own model objects for the G4 checks (independent of the g4_*.py scripts).
# T-system and M_0 (outline prop:set-modrows, def:set-M0); DNOPS (2-31) by residue extraction (no series reversion);
# local TR data at a branch point of x = kappa (p^4/4 + v/p), v = 1, in the coordinate t = p - a.
import os, sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from rf4_field import K, ZERO, ONE, smul, sinv, spow, ssqrt, sderiv, dfact

EPS = (1, -1, -1, 1, 1)
NU = (0, -1, 3, 2, 1)


def tsystem(N):
    """colour-0 F_b(t), b = 0..4, to order t^N: F1 = T_-1 F0, F2 = -T_-2 F1, F3 = T_2 F2, F4 = -T_1 F3, F0 = -T_0 F4,
    T_k = -1 + (t/5)(k - 4 t d/dt), F0(0) = 1. Solved order by order through the closing relation (pivot 4n)."""
    M = N + 2

    def T(k, G):
        out = [-g for g in G] + [F(0)]
        for n in range(len(G)):
            out[n + 1] += F(k - 4 * n, 5) * G[n]
        return out[:M + 1]

    def run(F0):
        F1 = T(-1, F0)
        F2 = [-g for g in T(-2, F1)]
        F3 = T(2, F2)
        F4 = [-g for g in T(1, F3)]
        back = [-g for g in T(0, F4)]
        return (F0, F1, F2, F3, F4), back
    f = [F(1)] + [F(0)] * M
    for n in range(1, N + 1):
        f[n] = F(0)
        _, back = run(f)
        res = back[n + 1] - f[n + 1]
        f[n] = -res / (4 * n)
    Fs, back = run(f)
    closing = all(back[n] == f[n] for n in range(N + 2))
    return [list(Fb[:N + 1]) for Fb in Fs], closing


def M0_coeffs(N):
    """M[n][i][j] = [t^n] M_0(t)_ij ; Rh[n][b][j] = [t^n] Rhat_b^(j)(t) = zeta^(j(b+1) - j n) F_b,n"""
    Fs, closing = tsystem(N)
    Rh = [[[K.z(j * (b + 1) - j * n) * Fs[b][n] for j in range(5)] for b in range(5)] for n in range(N + 1)]
    M = [[[ZERO] * 5 for _ in range(5)] for _ in range(N + 1)]
    for n in range(N + 1):
        for i in range(5):
            for j in range(5):
                acc = ZERO
                for b in range(5):
                    acc = acc + K.z(-i * (b + 1)) * Rh[n][b][j] * F(EPS[b], 5)
                M[n][i][j] = acc
    return Fs, closing, Rh, M


def x_taylor(kappa, a, N):
    """[t^k] of x(a+t) = kappa ((a+t)^4/4 + 1/(a+t)), v = 1, k = 0..N"""
    from math import comb
    out = []
    for k in range(N + 1):
        c = ZERO
        if k <= 4:
            c = c + a ** (4 - k) * F(comb(4, k), 4)
        c = c + a ** (-1 - k) * ((-1) ** k)
        out.append(c * kappa)
    return out


def W_local(kappa, a, N):
    """W(t) with w = sqrt(Lambda) W, x = x(a) - w^2/2, Lambda = -5 kappa, W = a t sqrt(g), g(0) = 1 (orientation dw/dp = sqrt(L) a)"""
    xt = x_taylor(kappa, a, N + 2)
    x2 = xt[2]
    g = [xt[k + 2] / x2 for k in range(N + 1)]
    h = ssqrt(g, N)
    # -2 (x - x(a)) = -2 x2 t^2 g = Lambda a^2 t^2 g   (x2 = (5/2) kappa a^2 at v = 1)
    assert x2 == a * a * F(5, 2) * kappa
    return [ZERO] + [a * c for c in h[:N]]          # W(t) = a t h(t)


def residue(series_coeffs, shift):
    """Res_{t=0} of t^shift * Sum c_n t^n dt  = c_(-1-shift)"""
    k = -1 - shift
    return series_coeffs[k] if 0 <= k < len(series_coeffs) else ZERO


def Rinv_dnops(kappa, Nmax):
    """X[N][j][i] = [lam^-N] R^-1(lam^-1)^j_i from DNOPS (2-31), finite part at i = j; xi^i(q) = -1/(sqrt(L) c_i (q - c_i)).
    [lam^-N] = -L^-N (2N-3)!! (1/c_i) Res_t (c_j - c_i + t)^-2 W(t)^(-2N+1), N >= 1; [lam^0] = delta_ij."""
    Lam = F(-5 * kappa)
    cs = [K.z(-j) for j in range(5)]
    Nser = 2 * Nmax + 4
    X = [[[ZERO] * 5 for _ in range(5)] for _ in range(Nmax + 1)]
    for j in range(5):
        Wt = W_local(kappa, cs[j], Nser + 2)
        hser = [Wt[k + 1] / cs[j] for k in range(Nser + 1)]          # W = c_j t h(t)
        for i in range(5):
            ci = cs[i]
            X[0][j][i] = ONE if i == j else ZERO
            for N in range(1, Nmax + 1):
                e = -2 * N + 1
                Wp = [c * cs[j] ** e for c in spow(hser, e, Nser)]     # W^e = c_j^e t^e h^e
                if i == j:
                    ser = Wp; shift = e - 2                            # t^-2 * t^e * h^e
                else:
                    base = [cs[j] - ci, ONE] + [ZERO] * (Nser - 1)
                    ser = smul(spow(base, -2, Nser), Wp, Nser); shift = e
                r = residue(ser, shift)
                X[N][j][i] = -(r / ci) * F(dfact(2 * N - 3)) / Lam ** N
    return X


def sigma_local(a, N, kappa=1):
    """local involution sigma(t) at the branch point a (x(a + sigma) = x(a + t)), by fixed point s = -u(t)/h(s)"""
    xt = x_taylor(kappa, a, N + 3)
    x2 = xt[2]
    g = [xt[k + 2] / x2 for k in range(N + 2)]
    h = ssqrt(g, N + 1)
    u = [ZERO] + h[:N]                     # u(t) = t h(t)
    s = [ZERO, -ONE] + [ZERO] * (N - 1)
    for _ in range(N + 1):
        hs = compose(h, s, N)
        s = [-c for c in smul(u, sinv(hs, N), N)]
    return s, xt


def compose(f, g, N):
    """f(g(t)) with g[0] = 0"""
    out = [ZERO] * (N + 1); pw = [ONE] + [ZERO] * N
    for k in range(min(len(f), N + 1)):
        if k > 0:
            pw = smul(pw, g, N)
        if not f[k].is0():
            out = [o + f[k] * q for o, q in zip(out, pw)]
    return out
