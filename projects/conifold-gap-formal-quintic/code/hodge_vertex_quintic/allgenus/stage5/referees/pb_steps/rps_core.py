"""rps_core.py -- referee pb_steps (Stage 5): exact series over Q for the residue of Theorem P(b), step 6.

Written from the formulas of PB_AUDIT s.2.1 and s.2.7, re-derived by hand; imports no project module (stdlib only).

Original variables:  V = Res_(x=0) e^(-W/z) at dx,  W = Lam/x + G(x),  Lam = -log(1-u),
    G = Sum_m x^(5m-1)/(m(5m-1)),
    log at = -(1/2) log(1-x^5) + Sum_k b_k (xz)^(2k-1) [5^(1-2k) - 5 Sum_P C(2k-2+5P,5P) x^(5P)].
    V is computed as Sum_k (-Lam/z)^k/k! [x^(k-1)] (e^(-G/z) at): dict {(E, beta): coeff of u^E z^beta}.
Hatted variables (x = sX, z = s^4 w, s^5 = u):  e^(-What/w) B,  What = Sum_m (u^(m-1)/m)[1/X + X^(5m-1)/(5m-1)],
    B(X) = at(sX, s^4 w).  The integrand at u^N is e^(-f0/w) P_N(X, w), f0 = 1/X + X^4/4, P_N a Laurent polynomial;
    H = Res_X[e^(-What/w) B dX] is computed DIRECTLY at each hatted order, keeping e^(-f0/w) whole:
    Res_X[X^p w^q e^(-f0/w)] = w^q Sum_l (-1)^(a+l) w^(-a-l)/(a! 4^l l!),  a = 4l + p + 1 >= 0.
"""
from fractions import Fraction as Fr
from math import comb, factorial

def add(d, k, v):
    if v:
        nv = d.get(k, 0) + v
        if nv:
            d[k] = nv
        else:
            d.pop(k, None)

_B = [Fr(1)]
def bernoulli(n):
    while len(_B) <= n:
        m = len(_B)
        _B.append(-sum(comb(m + 1, k) * _B[k] for k in range(m)) / (m + 1))
    return _B[n]

def bk(k):
    """b_k = B_2k/(2k(2k-1)); b_1 = 1/12, b_2 = -1/360."""
    return bernoulli(2 * k) / (2 * k * (2 * k - 1))

def mul_u(a, b, N):
    """product of two u-series (lists of Fractions) mod u^(N+1)."""
    out = [Fr(0)] * (N + 1)
    for i, x in enumerate(a[:N + 1]):
        if x:
            for j in range(0, N + 1 - i):
                if j < len(b) and b[j]:
                    out[i + j] += x * b[j]
    return out

# ---------------------------------------------------------------- original variables
def orig_log_amp(A, flip=False, nosqrt=False, gfirst=False):
    """log(e^(-G/z) at) as a list over x^a (a = 0..A) of dicts {z-power: coeff}.
    Controls: flip (b_k -> -b_k), nosqrt (drop (1-x^5)^(-1/2)), gfirst (keep only x^4/4 in G)."""
    L = [dict() for _ in range(A + 1)]
    if not nosqrt:
        for n in range(1, A // 5 + 1):
            add(L[5 * n], 0, Fr(1, 2 * n))
    for m in range(1, (A + 1) // 5 + 1):
        if 5 * m - 1 <= A and (m == 1 or not gfirst):
            add(L[5 * m - 1], -1, -Fr(1, m * (5 * m - 1)))
    sg = -1 if flip else 1
    k = 1
    while 2 * k - 1 <= A:
        n = 2 * k - 1
        b = sg * bk(k)
        for P in range(0, (A - n) // 5 + 1):
            c = Fr(-5 * comb(n - 1 + 5 * P, 5 * P))
            if P == 0:
                c += Fr(1, 5 ** n)
            add(L[n + 5 * P], n, b * c)
        k += 1
    return L

def exp_x(L, A):
    """exp of a series in x without constant term; coefficients are dicts (Laurent polynomials in z)."""
    F = [dict() for _ in range(A + 1)]
    F[0] = {0: Fr(1)}
    for a in range(1, A + 1):
        acc = {}
        for i in range(1, a + 1):
            if L[i]:
                for p, c in L[i].items():
                    ic = i * c
                    for q, d in F[a - i].items():
                        add(acc, p + q, ic * d)
        F[a] = {p: c / a for p, c in acc.items()}
    return F

def orig_V(E, **kw):
    """V through u^E: dict {(N, beta): coeff of u^N z^beta}."""
    F = exp_x(orig_log_amp(E - 1, **kw), E - 1)
    lam = [Fr(0)] + [Fr(1, n) for n in range(1, E + 1)]
    pw = [Fr(1)] + [Fr(0)] * E
    V = {}
    for k in range(1, E + 1):
        pw = mul_u(pw, lam, E)
        fk = Fr((-1) ** k, factorial(k))
        for N in range(k, E + 1):
            if pw[N]:
                c = fk * pw[N]
                for p, d in F[k - 1].items():
                    add(V, (N, p - k), c * d)
    return V

# ---------------------------------------------------------------- hatted variables
def hat_logs(N, flip=False, noregroup=False, nosqrt=False, wrongw=False):
    """log E1 = -(What - f0)/w and log B, as u-series (lists over u^n) of dicts {(X-power, w-power): coeff}.
    Controls: noregroup (Lam/u -> 1, i.e. drop e^(-(Lam/u-1)/(Xw))), flip, nosqrt,
    wrongw (hatting with z = s^5 w in the amplitude: (xz)^n -> u^(6n/5)..., realised as an extra u per Stirling term)."""
    LE = [dict() for _ in range(N + 1)]
    LB = [dict() for _ in range(N + 1)]
    for n in range(1, N + 1):
        if not noregroup:
            add(LE[n], (-1, -1), -Fr(1, n + 1))
        add(LE[n], (5 * n + 4, -1), -Fr(1, (n + 1) * (5 * n + 4)))
        if not nosqrt:
            add(LB[n], (5 * n, 0), Fr(1, 2 * n))
    sg = -1 if flip else 1
    k = 1
    while 2 * k - 1 <= N:
        n = 2 * k - 1
        b = sg * bk(k)
        for P in range(0, N - n + 1):
            c = Fr(-5 * comb(n - 1 + 5 * P, 5 * P))
            if P == 0:
                c += Fr(1, 5 ** n)
            uo = n + P + (1 if wrongw else 0)
            if uo <= N:
                add(LB[uo], (n + 5 * P, n), b * c)
        k += 1
    return LE, LB

def exp_u(L, N):
    """exp of a u-series without constant term; coefficients dicts keyed by tuples (added componentwise)."""
    F = [dict() for _ in range(N + 1)]
    F[0] = {(0, 0): Fr(1)}
    for n in range(1, N + 1):
        acc = {}
        for i in range(1, n + 1):
            for k1, c in L[i].items():
                ic = i * c
                for k2, d in F[n - i].items():
                    add(acc, (k1[0] + k2[0], k1[1] + k2[1]), ic * d)
        F[n] = {k: c / n for k, c in acc.items()}
    return F

_FACT = [1]
def fact(n):
    while len(_FACT) <= n:
        _FACT.append(_FACT[-1] * len(_FACT))
    return _FACT[n]

def res_f0(Pn, jmin, jmax):
    """Res_X[e^(-f0/w) Sum c X^p w^q] restricted to w^j, jmin <= j <= jmax: dict {j: coeff}."""
    out = {}
    for (p, q), c in Pn.items():
        lmax = (q - p - 1 - jmin) // 5
        for l in range(0, lmax + 1):
            a = 4 * l + p + 1
            if a < 0:
                continue
            j = q - a - l
            if j > jmax:
                continue
            add(out, j, c * Fr((-1) ** (a + l), fact(a) * 4 ** l * fact(l)))
    return out

def hat_integrand(N, **kw):
    LE, LB = hat_logs(N, **kw)
    L = []
    for n in range(N + 1):
        d = dict(LE[n])
        for k, v in LB[n].items():
            add(d, k, v)
        L.append(d)
    return exp_u(L, N)

def hat_B(N, **kw):
    """B(X) alone (u-series of dicts {(X-power, w-power)}), for the CFU / reduction side."""
    _, LB = hat_logs(N, **kw)
    return exp_u(LB, N)

def jmin_for(N, E):
    """smallest w-power j whose original image u^(N + (1-4j)/5) z^j lies within u^E."""
    return -((5 * (E - N) - 1) // 4)

def hat_H(Nh, E, **kw):
    """direct hatted residue H at hatted orders 0..Nh, all w-powers j in [jmin_for(N, E), 2N+5]."""
    P = hat_integrand(Nh, **kw)
    return [res_f0(P[n], jmin_for(n, E), 2 * n + 5) for n in range(Nh + 1)], P
