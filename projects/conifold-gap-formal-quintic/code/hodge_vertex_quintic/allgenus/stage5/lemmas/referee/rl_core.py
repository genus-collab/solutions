# -*- coding: utf-8 -*-
"""rl_core.py -- Stage 5 / lemmas / referee: an independent exact engine for the formal stationary phase (SP).

New code, Python stdlib only.  Imports no project module and no module of the worker (sp_core, sp_iter, ...).
Definition (FORMAL_SP_LEMMAS s.1.4 = CHART_PROOF s.1.2): phase F with nondegenerate critical point c, Hessian H,
C = -H^-1, amplitude a(t; z), t = x - c:
    SP^norm_c[F; a] = sum_k z^k N_k,  N_k = sum_alpha [s^(2k) eta^alpha]( exp(F_3(s eta)/s^2) a(s eta; s^2) ) <eta^alpha>_C.
Algorithm, deliberately different from the worker's sp_core:
  * scaled variables t = s eta, z = s^2; the exponential is computed by the power-series recurrence in s
    (m P_m = sum_i i Q_i P_(m-i)), not by repeated multiplication with weight truncation;
  * Gaussian moments from the generating function, <eta^alpha>_C = alpha! [xi^alpha] exp(xi^T C xi/2), not from Wick's recursion;
  * one-variable SPs also by the Lagrange-Buermann (Morse) formula, a third method.
Coefficients: Fraction, or Ser (truncated Laurent series in one parameter with absolute cut-off Ser.T).
"""
from fractions import Fraction as Fr
from math import comb, factorial


class Ser:
    """sum_{p <= T} c_p d^p, a truncated Laurent series in one parameter d; T = Ser.T (class attribute)."""
    T = 10

    def __init__(self, d=None):
        self.c = {p: Fr(v) for p, v in (d or {}).items() if v != 0 and p <= Ser.T}

    @staticmethod
    def const(v):
        return Ser({0: v})

    def _co(self, o):
        return o if isinstance(o, Ser) else Ser.const(o)

    def __add__(self, o):
        o = self._co(o); r = dict(self.c)
        for p, v in o.c.items():
            r[p] = r.get(p, 0) + v
        return Ser(r)
    __radd__ = __add__

    def __neg__(self):
        return Ser({p: -v for p, v in self.c.items()})

    def __sub__(self, o):
        return self + (-self._co(o))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Ser):
            o = Fr(o)
            return Ser({p: v * o for p, v in self.c.items()})
        r, T = {}, Ser.T
        for p, v in self.c.items():
            for p2, v2 in o.c.items():
                if p + p2 <= T:
                    r[p + p2] = r.get(p + p2, 0) + v * v2
        return Ser(r)
    __rmul__ = __mul__

    def __truediv__(self, o):
        return self * o.inv() if isinstance(o, Ser) else self * (1 / Fr(o))

    def __rtruediv__(self, o):
        return self.inv() * o

    def __eq__(self, o):
        return self.c == self._co(o).c

    def __ne__(self, o):
        return not self.__eq__(o)

    __hash__ = None

    def val(self):
        return min(self.c) if self.c else None

    def coef(self, p):
        return self.c.get(p, Fr(0))

    def inv(self):
        v = self.val(); a0 = self.c[v]; N = Ser.T - abs(v); w = []
        for n in range(N + 1):
            s = Fr(1) if n == 0 else Fr(0)
            for i in range(1, n + 1):
                s -= self.c.get(v + i, 0) * w[n - i]
            w.append(s / a0)
        return Ser({-v + n: w[n] for n in range(N + 1)})

    def deriv(self):
        return Ser({p - 1: p * v for p, v in self.c.items() if p})

    def __pow__(self, k):
        assert isinstance(k, int) and k >= 0
        R = Ser.const(1)
        for _ in range(k):
            R = R * self
        return R

    def __repr__(self):
        return 'Ser(%s)' % {p: str(v) for p, v in sorted(self.c.items())}


def ser_geom(c0, r, n):
    """(c0 + r d)^n as Ser, c0 != 0 a rational, r rational, n an integer (binomial series in d)"""
    c0, r = Fr(c0), Fr(r)
    return Ser({p: Fr(comb_gen(n, p)) * c0 ** (n - p) * r ** p for p in range(Ser.T + 1)}) if r else Ser.const(c0 ** n)


def comb_gen(n, p):
    """generalised binomial coefficient C(n, p), n any rational, p >= 0"""
    out = Fr(1)
    for j in range(p):
        out = out * (Fr(n) - j) / (j + 1)
    return out


# ---- polynomials in eta (dict: exponent tuple -> coefficient; coefficients Fraction or Ser) ----------------------------
def padd(A, B, cb=1):
    R = dict(A)
    for k, v in B.items():
        R[k] = R.get(k, 0) + (v * cb if cb != 1 else v)
    return {k: v for k, v in R.items() if v != 0}


def pmul(A, B):
    R = {}
    for ka, va in A.items():
        for kb, vb in B.items():
            k = tuple(x + y for x, y in zip(ka, kb)); R[k] = R.get(k, 0) + va * vb
    return {k: v for k, v in R.items() if v != 0}


def pscal(A, c):
    return {k: v * c for k, v in A.items() if v * c != 0}


def exp_scaled(Q, M, n):
    """P = exp(sum_{m>=1} s^m Q[m]) through s^M; Q[m] polynomials in n variables (Q[0] ignored). Returns list P[0..M]."""
    one = {(0,) * n: Fr(1)}; P = [one]
    for m in range(1, M + 1):
        acc = {}
        for i in range(1, m + 1):
            if i < len(Q) and Q[i]:
                acc = padd(acc, pmul(pscal(Q[i], i), P[m - i]))
        P.append(pscal(acc, Fr(1, m)))
    return P


def gf_moments(C, n, maxdeg):
    """<eta^alpha>_C for |alpha| <= maxdeg (even), from exp(xi^T C xi/2): alpha! [xi^alpha] G_r, G_r = (xi^T C xi/2)^r/r!"""
    Qf = {}
    for i in range(n):
        for j in range(i, n):
            e = [0] * n; e[i] += 1; e[j] += 1
            v = C[i][j] if i != j else C[i][i] * Fr(1, 2)
            if v != 0:
                Qf[tuple(e)] = v
    G = {(0,) * n: Fr(1)}; mom = dict(G)
    for r in range(1, maxdeg // 2 + 1):
        G = pscal(pmul(Qf, G), Fr(1, r))
        for al, v in G.items():
            f = 1
            for x in al:
                f *= factorial(x)
            mom[al] = v * f
    return mom


def sp_scaled(n, Q, A, C, K):
    """[N_0, .., N_K] of SP^norm: Q[m] (m >= 1) = s^m-part of the exponent F_3(s eta)/s^2, A[m] = s^m-part of the scaled
    amplitude a(s eta; s^2) (z-powers >= 0), C the covariance (n x n, entries Fraction or Ser)."""
    P = exp_scaled(Q, 2 * K, n); Tk = []
    for k in range(K + 1):
        acc = {}
        for m in range(2 * k + 1):
            if 2 * k - m < len(A) and A[2 * k - m]:
                acc = padd(acc, pmul(P[m], A[2 * k - m]))
        Tk.append(acc)
    maxdeg = max([sum(al) for T in Tk for al in T] + [0])
    mom = gf_moments(C, n, maxdeg + (maxdeg % 2))
    out = []
    for T in Tk:
        s = Fr(0)
        for al, v in T.items():
            if sum(al) % 2 == 0:
                m_ = mom.get(al, 0)
                if m_ != 0:
                    s = s + v * m_
        out.append(s)
    return out


# ---- data builders ------------------------------------------------------------------------------------------------
def scaled_phase(taylor, n, M):
    """taylor: dict alpha -> coefficient of t^alpha in F(c+t) - F(c) (|alpha| >= 2).  Returns (Q, H): Q[m] = degree-(m+2) part
    (m = 1..M), H the Hessian (from the degree-2 part)."""
    Q = [dict() for _ in range(M + 1)]
    H = [[Fr(0)] * n for _ in range(n)]
    for al, v in taylor.items():
        d = sum(al)
        if d == 2:
            idx = [i for i in range(n) for _ in range(al[i])]
            if idx[0] == idx[1]:
                H[idx[0]][idx[0]] = H[idx[0]][idx[0]] + 2 * v
            else:
                H[idx[0]][idx[1]] = H[idx[0]][idx[1]] + v; H[idx[1]][idx[0]] = H[idx[1]][idx[0]] + v
        elif 3 <= d <= M + 2:
            Q[d - 2][al] = v
        elif d < 2 and v != 0:
            raise ValueError('linear or constant term in the phase data')
    return Q, H


def scaled_amp(amp, n, M):
    """amp: dict (alpha, j) -> coefficient of t^alpha z^j (j >= 0).  Returns A[m] = sum over |alpha| + 2j = m."""
    A = [dict() for _ in range(M + 1)]
    for (al, j), v in amp.items():
        m = sum(al) + 2 * j
        if m <= M and v != 0:
            A[m][al] = A[m].get(al, 0) + v
    return A


def inv_mat(H):
    """inverse of a square matrix over Q (Gauss-Jordan); for 1 x 1 also over Ser"""
    n = len(H)
    if n == 1:
        return [[1 / H[0][0]]] if not isinstance(H[0][0], Ser) else [[H[0][0].inv()]]
    A = [[Fr(x) for x in r] + [Fr(int(i == j)) for j in range(n)] for i, r in enumerate(H)]
    for c in range(n):
        r = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[r] = A[r], A[c]
        iv = 1 / A[c][c]; A[c] = [x * iv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]


def sp_norm(taylor, amp, n, K):
    """SP^norm through z^K from Taylor data (phase: alpha -> coeff; amplitude: (alpha, j) -> coeff)"""
    Q, H = scaled_phase(taylor, n, 2 * K)
    Hi = inv_mat(H); C = [[-x for x in r] for r in Hi]
    return sp_scaled(n, Q, scaled_amp(amp, n, 2 * K), C, K)


# ---- Bernoulli numbers, Stirling series ------------------------------------------------------------------------------
def bernoulli(N):
    B = [Fr(1)]
    for k in range(1, N + 1):
        B.append(-sum(comb(k + 1, j) * B[j] for j in range(k)) / (k + 1))
    return B


def bk_list(K):
    B = bernoulli(2 * K + 2)
    return [None] + [B[2 * k] / (2 * k * (2 * k - 1)) for k in range(1, K + 1)]


def ser_exp(L, N):
    """exp of a power series L[0..N] with L[0] = 0 (recurrence E' = L' E)"""
    E = [Fr(1)] + [Fr(0)] * N
    for k in range(1, N + 1):
        E[k] = sum(j * L[j] * E[k - j] for j in range(1, k + 1)) / k
    return E


def stirling_m(m, N, bk):
    """St_m(w) through w^N by its DEFINITION (FORMAL_SP_LEMMAS s.9.1):
       exp((1/w + m - 1/2) log(1 + m w) - m + sum_k b_k w^(2k-1) (1 + m w)^(1-2k))"""
    m = Fr(m); lg = [Fr(0)] + [(-1) ** (j + 1) * m ** j / j for j in range(1, N + 2)]
    L = [lg[k + 1] + (m - Fr(1, 2)) * lg[k] for k in range(N + 1)]
    L[0] -= m
    for k in range(1, (N + 1) // 2 + 1):
        for i in range(N + 2 - 2 * k):
            L[2 * k - 1 + i] += bk[k] * comb_gen(1 - 2 * k, i) * m ** i
    assert L[0] == 0
    return ser_exp(L, N)


# ---- one-variable SP by Lagrange-Buermann (Morse coordinate), a third method ------------------------------------------
def sp_lb(f, g, K):
    """f[m] = [t^m](F(c+t) - F(c)) for m = 2..2K+2 (f[0] = f[1] = 0); g[j][m] = [t^m z^j] amplitude.
    N_k = sum_{j+r=k} (2r-1)!! C^r [t^(2r)] g_j S^(-r-1/2),  S = (F(c+t)-F(c))/(f_2 t^2),  C = -1/(2 f_2)."""
    f2 = f[2]; C = -1 / (2 * f2); M = 2 * K
    S = [Fr(1) if p == 0 else (f[p + 2] / f2 if p + 2 < len(f) else 0) for p in range(M + 1)]
    out = []
    pows = {}
    for r in range(K + 1):
        beta = -Fr(2 * r + 1, 2); P = [Fr(1)]
        for mm in range(1, 2 * r + 1):
            P.append(sum(((beta + 1) * i - mm) * S[i] * P[mm - i] for i in range(1, mm + 1)) / mm)
        pows[r] = P
    for k in range(K + 1):
        s = Fr(0) * C if not isinstance(C, Fr) else Fr(0)
        for r in range(k + 1):
            j = k - r
            if j >= len(g):
                continue
            P = pows[r]; gj = g[j]
            coef = sum((gj[i] if i < len(gj) else 0) * P[2 * r - i] for i in range(2 * r + 1))
            dfact = 1
            for q in range(1, 2 * r, 2):
                dfact *= q
            s = s + coef * dfact * C ** r if r else s + coef
        out.append(s)
    return out
