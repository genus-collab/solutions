# -*- coding: utf-8 -*-
"""ra_sp.py -- stage5/chart/referee_adversary: the normalised one-variable formal stationary-phase series c(z; L) of a
Mellin-Barnes-type summand, expanded in eps = L - sigma* at a colour point sigma* (a zero of the singular factor).

Integrand in sigma (fiat Stirling forms: log(a/z) read as log a; every log z is sigma-independent or absent):
  e^((lambda sigma + g(sigma))/z) * Gamma_St(M0 + a0(sigma)/z)^count * prod_groups prod_i Gamma_St(M + a_i(sigma)/z) * F(sigma, z),
  a0 = alpha0 (sigma - sigma*), lambda chosen so that sigma = L is critical (it never enters the normalised series).
Method (new, not used by the worker): sigma = L + eps tau, z = eps z', L = sigma* + eps.  In (tau, z') every Taylor datum is
a power series in eps (no negative powers), the singular factor becomes eps-free, and the Wick series in z' has coefficients
cp_n(eps) in F[[eps]].  Hence c_n(sigma* + eps) = eps^(-n) cp_n(eps):  c_n is regular at sigma* iff [eps^k] cp_n = 0 for k < n,
and then c_n(sigma*) = [eps^n] cp_n.  Ring R: RSym (truncated eps-series) or RNum (one numeric L).
Prediction (the residue picture, = the chart argument's conclusion): the normalised value at sigma* is
  exp( Sum_(k>=2) g*_k delta^k / z + log A_reg(sigma* + delta, z) - log A_reg(sigma*, 0) ),  delta = -M0 z / alpha0,
g* the Taylor coefficients at sigma* of the whole regular phase, A_reg the whole regular amplitude.
"""
from fractions import Fraction as Fr
from math import comb
from ra_core import stirling_b, gbinom, RSym

BSTD = stirling_b(40)


def _bk(spec, N):
    """resolve a Stirling specification into {k: b} for 2k - 1 <= N"""
    if spec is None or spec == 'drop':
        return {}
    out = {k: BSTD[k] for k in range(1, N // 2 + 2) if 2 * k - 1 <= N}
    if isinstance(spec, dict):
        for k, v in spec.items():
            if 2 * k - 1 <= N:
                out[k] = out.get(k, Fr(0)) + Fr(v)
    return out


class Group:
    """regular factors a_i(sigma) = alpha sigma + beta_i (a_i(sigma*) != 0) with common alpha, Poch shift M, Stirling spec bk.
    psum_eps(R, r): R-element Sum_i a_i(sigma* + eps)^(-r);  psum_star(F, r): F-element Sum_i a_i(sigma*)^(-r)."""

    def __init__(self, label, alpha, psum_eps, psum_star, M=0, bk='std', phase=True, expo=Fr(-1, 2)):
        self.label, self.alpha, self.pe, self.ps = label, Fr(alpha), psum_eps, psum_star
        self.M, self.bk, self.phase, self.expo = M, bk, phase, Fr(expo)

    def pochs(self):
        if self.M >= 0:
            return [(l, 1) for l in range(self.M)]
        return [(-l, -1) for l in range(1, -self.M + 1)]


class Sing:
    """the singular factor Gamma_St(M0 + alpha0 (sigma - sigma*)/z)^count"""

    def __init__(self, alpha0=-1, M0=0, bk='std', expo=Fr(-1, 2), count=1):
        self.alpha0, self.M0, self.bk, self.expo, self.count = Fr(alpha0), M0, bk, Fr(expo), Fr(count)


def build(spec, R, N):
    """Taylor data in (tau, z'): Fk[k] (k = 2..2N+2) and the log-amplitude LA[m][j] (j <= 2(N-m))"""
    F = R.F; KF = 2 * N + 2; Kr = getattr(R, 'K')
    Fk = [R.zero() for _ in range(KF + 1)]
    LA = [[R.zero() for _ in range(2 * (N - m) + 1)] for m in range(N + 1)]

    def addLA(m, j, x):
        if m <= N and j <= 2 * (N - m) and (m, j) != (0, 0):
            LA[m][j] = R.add(LA[m][j], x)

    S = spec['sing']; a0 = S.alpha0; c0 = S.count
    for k in range(2, KF + 1):
        Fk[k] = R.add(Fk[k], R.const(c0 * a0 * (-1) ** k / (k * (k - 1))))
    for r in range(1, 2 * N + 1):
        addLA(0, r, R.const(c0 * S.expo * Fr((-1) ** (r + 1), r)))
    for k, bv in _bk(S.bk, N).items():
        m = 2 * k - 1
        for r in range(0, 2 * (N - m) + 1):
            addLA(m, r, R.const(c0 * bv * a0 ** (1 - 2 * k) * gbinom(1 - 2 * k, r)))
    pl = [(l, 1) for l in range(S.M0)] if S.M0 >= 0 else [(-l, -1) for l in range(1, -S.M0 + 1)]
    for (l, sg) in pl:                                  # log(1 + tau + (l/alpha0) z')
        for p in range(1, 2 * N + 1):
            for i in range(0, p + 1):
                cf = c0 * sg * Fr((-1) ** (p + 1), p) * comb(p, i) * (Fr(l) / a0) ** i
                if cf != 0:
                    addLA(i, p - i, R.const(cf))
    for G in spec.get('groups', []):
        cache = {}

        def P(r, G=G, cache=cache):
            if r not in cache:
                cache[r] = G.pe(R, r)
            return cache[r]
        al = G.alpha
        if G.phase:
            for k in range(2, KF + 1):
                if k - 1 >= Kr:
                    break
                Fk[k] = R.add(Fk[k], R.mul(R.epsp(k - 1, al ** k * (-1) ** k / (k * (k - 1))), P(k - 1)))
        if G.expo != 0:
            for r in range(1, 2 * N + 1):
                if r >= Kr:
                    break
                addLA(0, r, R.mul(R.epsp(r, G.expo * Fr((-1) ** (r + 1), r) * al ** r), P(r)))
        for k, bv in _bk(G.bk, N).items():
            m = 2 * k - 1
            for r in range(0, 2 * (N - m) + 1):
                if m + r >= Kr:
                    break
                addLA(m, r, R.mul(R.epsp(m + r, bv * gbinom(1 - 2 * k, r) * al ** r), P(m + r)))
        for (l, sg) in G.pochs():                       # log(1 + (alpha tau + l z') eps / a)
            for p in range(1, 2 * N + 1):
                if p >= Kr:
                    break
                for i in range(0, p + 1):
                    cf = sg * Fr((-1) ** (p + 1), p) * comb(p, i) * al ** (p - i) * Fr(l) ** i
                    if cf != 0:
                        addLA(i, p - i, R.mul(R.epsp(p, cf), P(p)))
    for k0, g in spec.get('gstar', {}).items():          # extra phase g(sigma* + x) = Sum g_k x^k
        for k in range(2, min(k0, KF) + 1):
            if k0 - 1 >= Kr:
                break
            Fk[k] = R.add(Fk[k], R.epsp(k0 - 1, Fr(g) * comb(k0, k)))
    for (m, i), f in spec.get('fextra', {}).items():      # extra log-amplitude f z^m (sigma - sigma*)^i
        if m + i >= Kr:
            continue
        for r in range(0, i + 1):
            addLA(m, r, R.epsp(m + i, Fr(f) * comb(i, r)))
    return Fk, LA


def bexp(R, LA, N):
    """A = exp(LA) (LA[0][0] = 0), A[m] a tau-list of length 2(N-m)+1"""
    F = R.F
    T = [2 * (N - m) for m in range(N + 1)]
    E0 = [R.const(1)] + [R.zero() for _ in range(T[0])]
    for j in range(1, T[0] + 1):
        s = R.zero()
        for i in range(1, j + 1):
            s = R.add(s, R.scal(i, R.mul(LA[0][i], E0[j - i])))
        E0[j] = R.scal(Fr(1, j), s)
    A = [E0]
    for m in range(1, N + 1):
        acc = [R.zero() for _ in range(T[m] + 1)]
        for i in range(1, m + 1):
            Li, Ami = LA[i], A[m - i]
            for a in range(min(len(Li), T[m] + 1)):
                if R.iszero(Li[a]):
                    continue
                x = R.scal(i, Li[a])
                for b in range(min(len(Ami), T[m] + 1 - a)):
                    acc[a + b] = R.add(acc[a + b], R.mul(x, Ami[b]))
        A.append([R.scal(Fr(1, m), y) for y in acc])
    return A


def wick(R, Fk, A, N):
    """N_n, n = 0..N, of e^(F(tau)/z') A(tau, z') at tau = 0 (F = Sum_(k>=2) Fk tau^k)"""
    C = R.scal(Fr(-1, 2), R.inv(Fk[2]))
    F3 = [R.zero(), R.zero(), R.zero()] + Fk[3:]
    PV = [[R.const(1)]]
    for V in range(1, 2 * N + 1):
        D = 2 * (N + V); prev = PV[-1]; new = [R.zero() for _ in range(D + 1)]
        for a, pa in enumerate(prev):
            if R.iszero(pa):
                continue
            for b in range(3, min(len(F3), D + 1 - a)):
                new[a + b] = R.add(new[a + b], R.mul(pa, F3[b]))
        PV.append([R.scal(Fr(1, V), y) for y in new])
    Cp = [R.const(1)]
    df = [Fr(1)]
    for r in range(1, 3 * N + 2):
        Cp.append(R.mul(Cp[-1], C)); df.append(df[-1] * (2 * r - 1))
    out = []
    for n in range(N + 1):
        tot = R.zero()
        for m in range(n + 1):
            for V in range(0, 2 * (n - m) + 1):
                r = n - m + V; s = R.zero(); P = PV[V]
                for j in range(0, min(2 * (N - m), 2 * (n - m) - V) + 1):
                    d = 2 * r - j
                    if d < len(P) and not R.iszero(A[m][j]):
                        s = R.add(s, R.mul(A[m][j], P[d]))
                tot = R.add(tot, R.mul(R.scal(df[r], s), Cp[r]))
        out.append(tot)
    return out


def laurent(spec, R, N):
    """RSym: list over n of (negative part [eps^-n..eps^-1], value c_n(sigma*), next Taylor coefficients)"""
    Fk, LA = build(spec, R, N)
    Nn = wick(R, Fk, bexp(R, LA, N), N)
    return Nn


def predict(spec, F, N):
    """the residue-picture value at sigma*: list of F-elements [z^0..z^N]"""
    S = spec['sing']; dc = -Fr(S.M0) / S.alpha0
    PL = [F.zero] * (N + 1); gt = {}
    for G in spec.get('groups', []):
        al = G.alpha
        if G.phase:
            for k in range(2, N + 2):
                gt[k] = F.red(gt.get(k, 0) + F.c(al ** k * (-1) ** k / (k * (k - 1))) * G.ps(F, k - 1))
        for p in range(1, N + 1):
            if G.expo != 0 and dc != 0:
                PL[p] = F.red(PL[p] + F.c(G.expo * Fr((-1) ** (p + 1), p) * (al * dc) ** p) * G.ps(F, p))
            for (l, sg) in G.pochs():
                PL[p] = F.red(PL[p] + F.c(sg * Fr((-1) ** (p + 1), p) * (al * dc + l) ** p) * G.ps(F, p))
        for k, bv in _bk(G.bk, N).items():
            for r in range(0, N - (2 * k - 1) + 1):
                PL[2 * k - 1 + r] = F.red(PL[2 * k - 1 + r] + F.c(bv * gbinom(1 - 2 * k, r) * (al * dc) ** r)
                                          * G.ps(F, 2 * k - 1 + r))
    for k, g in spec.get('gstar', {}).items():
        gt[k] = F.red(gt.get(k, 0) + F.c(g))
    for (m, i), f in spec.get('fextra', {}).items():
        if 1 <= m + i <= N and (m, i) != (0, 0):
            PL[m + i] = F.red(PL[m + i] + F.c(Fr(f) * dc ** i))
    for k, g in gt.items():
        if 1 <= k - 1 <= N:
            PL[k - 1] = F.red(PL[k - 1] + g * F.c(dc ** k))
    return sexp(F, PL)


def sexp(F, L):
    """exp of a power series with L[0] = 0"""
    N = len(L) - 1; E = [F.one] + [F.zero] * N
    for n in range(1, N + 1):
        s = 0
        for i in range(1, n + 1):
            s += i * L[i] * E[n - i]
        E[n] = F.red(s * F.c(Fr(1, n)))
    return E
