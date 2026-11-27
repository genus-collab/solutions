"""rps_s2.py -- referee pb_steps: the model rows of UR Prop. S and the stationary phase (S2) of PB s.2.6.

Own code (stdlib only).  Prop. S (colour 0): F_1 = T_-1 F_0, F_2 = -T_-2 F_1, F_3 = T_2 F_2, F_4 = -T_1 F_3, F_0 = -T_0 F_4,
T_k = -1 + (t/5)(k - 4 t d/dt);  on t^(n-1): (t/5)(k - 4 theta) t^(n-1) = ((k - 4(n-1))/5) t^n.
F_0 is solved order by order from the closing relation; the slope (pivot) of [t^(n+1)] in F_0[n] is recorded.
(S2): SP^norm_(eta=1)[eta^m e^(-(1/eta + eta^4/4)/t) d eta] = Sum_r (2r-1)!! (t/5)^r [e^(2r)] (1+e)^m S(e)^(-r-1/2),
S(e) = (f(1+e) - 5/4)/((5/2) e^2)  (Lagrange-Buermann in the Morse coordinate; a method different from Wick contraction).
"""
from fractions import Fraction as Fr
from math import comb

C = (0, -1, 3, 2, 1)
EPS = (1, -1, -1, 1, 1)

def T_op(k, F, M):
    out = [Fr(0)] * (M + 1)
    for n in range(M + 1):
        out[n] = -F[n]
        if n >= 1:
            out[n] += Fr(k - 4 * (n - 1), 5) * F[n - 1]
    return out

def chain(F0, M, shift=0):
    """the rows F_1..F_4 from F_0 and the closing value -T_0 F_4 (shift: control, T_k -> T_(k+shift))."""
    F1 = T_op(-1 + shift, F0, M)
    F2 = [-x for x in T_op(-2 + shift, F1, M)]
    F3 = T_op(2 + shift, F2, M)
    F4 = [-x for x in T_op(1 + shift, F3, M)]
    cl = [-x for x in T_op(0 + shift, F4, M)]
    return [F0, F1, F2, F3, F4], cl

def solve_F(M, shift=0):
    """F_b through t^M with F_0(0) = 1; returns (F, pivots) where pivots[n] = d[t^(n+1)](closing - F_0)/d F_0[n]."""
    F0 = [Fr(1)] + [Fr(0)] * (M + 1)
    piv = {}
    for n in range(1, M + 1):
        F0[n] = Fr(0)
        _, cl = chain(F0, M + 1, shift)
        g0 = cl[n + 1] - F0[n + 1]
        F0[n] = Fr(1)
        _, cl = chain(F0, M + 1, shift)
        g1 = cl[n + 1] - F0[n + 1]
        piv[n] = g1 - g0
        F0[n] = -g0 / (g1 - g0)
    F, cl = chain(F0, M + 1, shift)
    resid = [cl[n] - F0[n] for n in range(M + 1)]
    return [f[:M + 1] for f in F], piv, resid

def ser_log1(a, K):
    """log(a) for a series a with a[0] = 1, mod e^(K+1)."""
    da = [n * a[n] for n in range(1, K + 1)] + [Fr(0)]
    inv = [Fr(1)] + [Fr(0)] * K
    for n in range(1, K + 1):
        inv[n] = -sum(a[i] * inv[n - i] for i in range(1, n + 1))
    q = [sum(da[i] * inv[n - i] for i in range(n + 1)) for n in range(K)]
    return [Fr(0)] + [q[n - 1] / n for n in range(1, K + 1)]

def ser_exp(L, K):
    F = [Fr(1)] + [Fr(0)] * K
    for n in range(1, K + 1):
        F[n] = sum(i * L[i] * F[n - i] for i in range(1, n + 1)) / n
    return F

def sp_model(m, M):
    """SP^norm of eta^m e^(-(1/eta + eta^4/4)/t) at eta = 1, through t^M (exact)."""
    K = 2 * M
    cfs = [Fr((-1) ** n) + Fr(comb(4, n), 4) for n in range(K + 3)]   # f(1+e) = Sum_n cfs[n] e^n
    S = [cfs[n + 2] / Fr(5, 2) for n in range(K + 1)]
    logS = ser_log1(S, K)
    if m >= 0:
        g = [Fr(comb(m, n)) for n in range(K + 1)]
    else:
        g = [Fr((-1) ** n * comb(-m + n - 1, n)) for n in range(K + 1)]
    out = []
    dfac = Fr(1)
    for r in range(M + 1):
        if r >= 1:
            dfac *= (2 * r - 1)
        Sp = ser_exp([-(r + Fr(1, 2)) * x for x in logS], 2 * r)
        coef = sum(g[i] * Sp[2 * r - i] for i in range(2 * r + 1))
        out.append(dfac * Fr(1, 5 ** r) * coef)
    return out
