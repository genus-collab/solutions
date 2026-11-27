"""rps_sp.py -- referee pb_steps: the hatted stationary phase Phi^SP over A = Q[[u]] (exact), for PB s.2.6.

Phi^SP := SP_(X=1)[e^(-What/w) B dX] / (sqrt(2 pi w/5) e^(-What(1)/w))
        = (5/a)^(1/2) Sum_r (2r-1)!! (w/a)^r [e^(2r)] B(1+e) S(e)^(-r-1/2),   a = What''(1) = 5/(1-u),
S(e) = (What(1+e) - What(1))/((a/2) e^2),  [e^k u^(m-1)] What(1+e) = (1/m)[(-1)^k + C(5m-1,k)/(5m-1)]  (Lagrange-Buermann).
Own code, stdlib only.  u-series are lists of Fractions; e-series of u-series are lists of u-series.
"""
from fractions import Fraction as Fr
from math import comb
from rps_core import mul_u

def inv_u(a, N):
    out = [Fr(0)] * (N + 1)
    out[0] = 1 / a[0]
    for n in range(1, N + 1):
        out[n] = -sum(a[i] * out[n - i] for i in range(1, n + 1) if i < len(a)) / a[0]
    return out

def log_u(a, N):
    """log a for a u-series with a[0] = 1."""
    assert a[0] == 1
    L = [Fr(0)] * (N + 1)
    for n in range(1, N + 1):
        L[n] = (n * a[n] - sum(k * L[k] * a[n - k] for k in range(1, n))) / n
    return L

def exp_u0(L, N):
    """exp of a u-series with L[0] = 0."""
    F = [Fr(1)] + [Fr(0)] * N
    for n in range(1, N + 1):
        F[n] = sum(k * L[k] * F[n - k] for k in range(1, n + 1)) / n
    return F

def pow_u(a, alpha, N):
    """a^alpha for a = 1 + O(u), alpha rational."""
    return exp_u0([alpha * x for x in log_u(a, N)], N)

def zero(N):
    return [Fr(0)] * (N + 1)

def uadd(a, b, c=1):
    return [x + c * y for x, y in zip(a, b)]

def log_eu(S, K, N):
    """log S, S an e-series (length >= K+1) of u-series with S[0] = 1."""
    L = [zero(N) for _ in range(K + 1)]
    for n in range(1, K + 1):
        acc = [n * x for x in S[n]]
        for k in range(1, n):
            acc = uadd(acc, mul_u(L[k], S[n - k], N), -k)
        L[n] = [x / n for x in acc]
    return L

def exp_eu(L, K, N):
    F = [zero(N) for _ in range(K + 1)]
    F[0][0] = Fr(1)
    for n in range(1, K + 1):
        acc = zero(N)
        for k in range(1, n + 1):
            acc = uadd(acc, mul_u(L[k], F[n - k], N), k)
        F[n] = [x / n for x in acc]
    return F

def hat_SP(Bhat, nmax, N):
    """[w^n] Phi^SP as u-series mod u^(N+1), n = 0..nmax; Bhat = rps_core.hat_B(N)."""
    K = 2 * nmax
    Wk = [[Fr(0)] * (N + 1) for _ in range(K + 3)]
    for k in range(K + 3):
        for i in range(N + 1):
            m = i + 1
            Wk[k][i] = Fr(1, m) * (Fr((-1) ** k) + Fr(comb(5 * m - 1, k), 5 * m - 1))
    inv_w2 = inv_u(Wk[2], N)
    S = [mul_u(Wk[k + 2], inv_w2, N) for k in range(K + 1)]
    logS = log_eu(S, K, N)
    one_minus_u = [Fr(1), Fr(-1)] + [Fr(0)] * (N - 1)
    sq = pow_u(one_minus_u, Fr(1, 2), N)                       # (5/a)^(1/2)
    wa = [x / 5 for x in one_minus_u]                          # 1/a = (1-u)/5
    # B(1+e) per w-power q: e-series of u-series
    Bq = {}
    for i in range(N + 1):
        for (p, q), c in Bhat[i].items():
            d = Bq.setdefault(q, [zero(N) for _ in range(K + 1)])
            for k in range(0, min(p, K) + 1):
                d[k][i] += c * comb(p, k)
    out = []
    Spow = {}
    for n in range(nmax + 1):
        tot = zero(N)
        for r in range(0, n + 1):
            q = n - r
            if q not in Bq:
                continue
            if r not in Spow:
                Spow[r] = exp_eu([[-(r + Fr(1, 2)) * x for x in s] for s in logS[:2 * r + 1]], 2 * r, N)
            Sp = Spow[r]
            co = zero(N)
            for k in range(0, 2 * r + 1):
                co = uadd(co, mul_u(Bq[q][k], Sp[2 * r - k], N))
            dfac = 1
            for t in range(1, 2 * r, 2):
                dfac *= t
            war = [Fr(1)] + [Fr(0)] * N
            for _ in range(r):
                war = mul_u(war, wa, N)
            tot = uadd(tot, [dfac * x for x in mul_u(war, co, N)])
        out.append(mul_u(sq, tot, N))
    return out

def rhs_s26(pi, F, eps, cb, nmax, N, use_tau=True, use_eps=True, tau_in_F=True):
    """[w^n] Sum_b eps_b pi_(c_b)(u,w) tau^((1-c_b)/4) F_b(tau w), n = 0..nmax, mod u^(N+1)."""
    tinv = [Fr(4, 4 + 5 * n) for n in range(N + 1)]            # 1/tau = 2F1(1,4/5;9/5;u)
    out = [zero(N) for _ in range(nmax + 1)]
    cache = {}
    def tpow(al):
        if al not in cache:
            cache[al] = pow_u(tinv, -al, N)
        return cache[al]
    for b in range(5):
        m = cb[b]
        sg = eps[b] if use_eps else 1
        for (i, j), c in pi[m].items():
            if i > N or j > nmax:
                continue
            for n in range(j, nmax + 1):
                k = n - j
                al = (Fr(1 - m, 4) if use_tau else 0) + (k if tau_in_F else 0)
                tw = tpow(al)
                coef = sg * c * F[b][k]
                for i2 in range(0, N - i + 1):
                    out[n][i + i2] += coef * tw[i2]
    return out
