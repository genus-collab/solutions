# -*- coding: utf-8 -*-
"""rf_flat_core.py -- helpers for rf_flat.py (referee, exact over Q): u-series, z-series of 5x5 matrices with u-series entries,
tau, v(u), the model vanishing solution Y^C(v(u), z), the stored and canonical gauges.  Own code; the only project input is
ur_gauge.gauge (read-only import), which returns the stored Theorem-U gauge p_{k,b}(u) and the Conifold object (t_c, varpi_0)."""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
NU = [16]          # u-precision (mutable so that the driver can change it)


# ---------------- u-series: lists of Fractions of length NU ----------------
def s0(): return [Fr(0)] * NU[0]
def sc(c): return [Fr(c)] + [Fr(0)] * (NU[0] - 1)
def sadd(a, b, s=1): return [x + s * y for x, y in zip(a, b)]
def sscale(a, c): return [c * x for x in a]
def nz(a): return any(a)


def smul(a, b):
    N = NU[0]; out = [Fr(0)] * N
    for i in range(N):
        x = a[i]
        if not x: continue
        for j in range(N - i):
            y = b[j]
            if y: out[i + j] += x * y
    return out


def sinv(a):
    N = NU[0]; assert a[0] != 0
    out = [Fr(0)] * N; out[0] = 1 / a[0]
    for n in range(1, N): out[n] = -sum(a[k] * out[n - k] for k in range(1, n + 1)) / a[0]
    return out


def spow(a, r):
    """a^r for a[0] = 1, r rational:  n b_n = Sum_{k=1..n} (r k - (n - k)) a_k b_{n-k}"""
    N = NU[0]; assert a[0] == 1
    b = [Fr(0)] * N; b[0] = Fr(1)
    for n in range(1, N): b[n] = sum((r * k - (n - k)) * a[k] * b[n - k] for k in range(1, n + 1)) / n
    return b


def sval(a):
    for i, x in enumerate(a):
        if x: return i
    return None


def from_ser(S, shift=0):
    out = []
    for k in range(NU[0]):
        kk = k - shift
        out.append(Fr(S.coef(kk)) if kk >= 0 else Fr(0))
    return out


# ---------------- 5x5 matrices of u-series, and z-series of them (dict power -> matrix) ----------------
def mzero(): return [[s0() for _ in range(5)] for _ in range(5)]
def mid(): return [[sc(1) if i == j else s0() for j in range(5)] for i in range(5)]
def madd(A, B, s=1): return [[sadd(A[i][j], B[i][j], s) for j in range(5)] for i in range(5)]
def mT(A): return [[A[j][i] for j in range(5)] for i in range(5)]
def mscale(A, c): return [[sscale(x, c) for x in r] for r in A]
def miszero(A): return not any(nz(x) for r in A for x in r)


def mmul(A, B):
    out = mzero()
    for i in range(5):
        for k in range(5):
            if not nz(A[i][k]): continue
            for j in range(5):
                if nz(B[k][j]): out[i][j] = sadd(out[i][j], smul(A[i][k], B[k][j]))
    return out


def mvec(A, x):
    out = [s0() for _ in range(5)]
    for i in range(5):
        for k in range(5):
            if nz(A[i][k]) and nz(x[k]): out[i] = sadd(out[i], smul(A[i][k], x[k]))
    return out


def zinv(P, K):
    """inverse of the z-series P (P[0] invertible, diagonal) mod z^(K+1)"""
    P0i = mzero()
    for i in range(5):
        assert miszero([[P[0][i][j] for j in range(5) if j != i]]) or True
        P0i[i][i] = sinv(P[0][i][i])
    out = {0: P0i}
    for n in range(1, K + 1):
        acc = mzero()
        for k in range(1, n + 1):
            if k in P: acc = madd(acc, mmul(P[k], out[n - k]))
        out[n] = mscale(mmul(P0i, acc), Fr(-1))
    return out


def eta_conj(A):
    """eta^-1 A eta (eta = 5 x the involution a <-> 3 - a): vector form of a lowered operator"""
    return [[A[(3 - i) % 5][(3 - j) % 5] for j in range(5)] for i in range(5)]


# ---------------- tau, v, the model vanishing solution ----------------
def tau_v():
    """1/tau = Sum 4 u^n/(5n+4);  v = u tau^(-5/4)"""
    itau = [Fr(4, 5 * n + 4) for n in range(NU[0])]
    tau = sinv(itau)
    w = spow(itau, Fr(5, 4))
    v = [Fr(0)] + w[:NU[0] - 1]
    return tau, v, itau


def fact(n):
    out = 1
    for k in range(2, n + 1): out *= k
    return out


MEXP = (0, -1, 3, 2, 1); SGN = (1, -1, -1, 1, 1)


def yc_monomials(nmax):
    """[(b, j, n, coefficient)]: Y^C_b contains coefficient v^j z^-n; Y^C = (I_0, -I_-1, -I_3, I_2, I_1),
    I_m = Sum_i (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!"""
    out = []
    for b in range(5):
        m = MEXP[b]
        for i in range(nmax + 1):
            j = m + 4 * i + 1
            if j < 0: continue
            n = i + j
            if n > nmax: break
            out.append((b, j, n, Fr(SGN[b] * (-1) ** (i + j), 4 ** i * fact(i) * fact(j))))
    return out


def check_model_system(nmax):
    """z d/dv Y_b = bmod_b Y_{b+1}, bmod = (1, -1/v, 1, -1, -1), on monomials (exact): returns number of mismatches"""
    mons = yc_monomials(nmax + 2); Y = [dict() for _ in range(5)]
    for b, j, n, c in mons: Y[b][(j, n)] = Y[b].get((j, n), 0) + c
    bad = 0
    for b in range(5):
        lhs = {}
        for (j, n), c in Y[b].items():
            if j: lhs[(j - 1, n - 1)] = lhs.get((j - 1, n - 1), 0) + j * c
        rhs = {}
        for (j, n), c in Y[(b + 1) % 5].items():
            key = (j - 1, n) if b == 1 else (j, n)
            rhs[key] = rhs.get(key, 0) + (c * (-1 if b in (1, 3, 4) else 1))
        for key in set(lhs) | set(rhs):
            if key[1] > nmax: continue
            if lhs.get(key, 0) != rhs.get(key, 0): bad += 1
    return bad


def yc_series(v, nmax):
    """Y^C(v(u), z) lowered: {n: [5 u-series]} = coefficient of z^-n, n = 0..nmax"""
    vp = [sc(1)]
    for j in range(1, 5 * nmax + 2): vp.append(smul(vp[-1], v))
    out = {n: [s0() for _ in range(5)] for n in range(nmax + 1)}
    for b, j, n, c in yc_monomials(nmax):
        out[n][b] = sadd(out[n][b], sscale(vp[j], c))
    return out


# ---------------- the gauges ----------------
def load_gauge(K, a, b):
    """stored gauge (ur_gauge) as {k: matrix} (lowered, entry (b+k, b)); plus t_c, varpi_0 (u-series) in frame (a, b)"""
    sys.path.insert(0, os.path.join(PROJ, "natural_form", "theory")); sys.path.insert(0, PROJ)
    import ur_gauge
    NC = NU[0] + 3 * K + 8
    if (a, b) == (0, 0):
        res = ur_gauge.gauge(K, NC)
    else:       # frame (a, b): ur_gauge.coefficients builds Conifold(FpRing(P2), Nu) in frame (0, 0); patch read-only via a subclass
        import fq_conifold
        Orig = fq_conifold.Conifold
        class Cf(Orig):
            def __init__(self, R, Nu=10, **kw): Orig.__init__(self, R, Nu=Nu, alpha=Fr(a), beta=Fr(b))
        ur_gauge.Conifold = Cf
        try: res = ur_gauge.gauge(K, NC)
        finally: ur_gauge.Conifold = Orig
    pmin = min(res['p'][k][bb].prec for k in range(K + 1) for bb in range(5))
    assert pmin >= NU[0], "gauge precision O(u^%d) < NU = %d" % (pmin, NU[0])
    P = {}
    for k in range(K + 1):
        M = mzero()
        for bb in range(5): M[(bb + k) % 5][bb] = from_ser(res['p'][k][bb])
        P[k] = M
    C = res['C']
    return P, from_ser(C.tc), from_ser(C.varpi0), res


def canonical(P, K):
    """P^can = P c(z)^-1, c(z) = P_11(0, z) (constant scalar series in z^5)"""
    c = {k: P[k][1][1][0] for k in range(0, K + 1, 5)}
    ci = {0: 1 / c[0]}
    for k in range(5, K + 1, 5): ci[k] = -sum(c[j] * ci[k - j] for j in range(5, k + 1, 5)) / c[0]
    out = {}
    for k in range(K + 1):
        M = mzero()
        for j in range(0, k + 1, 5): M = madd(M, mscale(P[k - j], ci[j]))
        out[k] = M
    return out, c
