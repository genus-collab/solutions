# -*- coding: utf-8 -*-
"""ch_wick.py -- Stage 5 / chart: graded multivariate polynomials and the formal stationary-phase functional (Wick).

A graded polynomial is {deg: {mono (tuple of exponents): coeff}} over a field/ring F of ch_core (QZ, FP or RS).
N[Phi, A; p](z) = Sum_n N_n z^n with  N_n = Sum_m Sum_V < [Phi_3^V A_m / V!]_(2(n-m+V)) >_C,  C = -H^-1  (CHART_PROOF s.1.2):
phase Phi given by its Taylor expansion at the critical point (degree-1 part must vanish, degree-2 part gives H),
amplitude A = Sum_m z^m A_m. Moments <eta^g>_C by the Isserlis recursion <eta_i eta^b> = Sum_j C_ij b_j <eta^(b - e_j)>.
"""
from fractions import Fraction as Fr


def gadd(F, A, B, cb=None):
    R = {d: dict(m) for d, m in A.items()}
    for d, m in B.items():
        od = R.setdefault(d, {})
        for e, c in m.items():
            v = c if cb is None else F.mul(cb, c)
            od[e] = F.add(od[e], v) if e in od else v
    return R


def gscale(F, A, c):
    return {d: {e: F.mul(c, v) for e, v in m.items()} for d, m in A.items()}


def gmul(F, A, B, dmax, dmin=0):
    R = {}
    for da, Ad in A.items():
        for db, Bd in B.items():
            d = da + db
            if d > dmax or d < dmin: continue
            od = R.setdefault(d, {})
            for e1, c1 in Ad.items():
                for e2, c2 in Bd.items():
                    e = tuple(x + y for x, y in zip(e1, e2)); v = F.mul(c1, c2)
                    od[e] = F.add(od[e], v) if e in od else v
    return R


def guni(F, n, i, ser, dmin=0):
    """univariate list ser (coefficients of t^k) in variable i, degrees >= dmin"""
    out = {}
    for k, c in enumerate(ser):
        if k >= dmin and not F.iszero(c):
            out[k] = {tuple(k if a == i else 0 for a in range(n)): c}
    return out


def gconst(F, n, c): return {0: {(0,) * n: c}}


def gtrunc(A, dmax, dmin=0): return {d: m for d, m in A.items() if dmin <= d <= dmax}


def gclean(F, A):
    return {d: {e: c for e, c in m.items() if not F.iszero(c)} for d, m in A.items()}


def linear_is_zero(F, Phi):
    return all(F.iszero(c) for c in Phi.get(1, {}).values())


def hessian(F, Phi, n):
    H = [[F.zero] * n for _ in range(n)]
    for e, c in Phi.get(2, {}).items():
        idx = [a for a in range(n) for _ in range(e[a])]
        if idx[0] == idx[1]: H[idx[0]][idx[0]] = F.smul(2, c)
        else: H[idx[0]][idx[1]] = c; H[idx[1]][idx[0]] = c
    return H


def matinv(F, M):
    n = len(M); A = [list(r) + [F.one if i == j else F.zero for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        r = next(r for r in range(c, n) if F.isunit(A[r][c])); A[c], A[r] = A[r], A[c]
        iv = F.inv(A[c][c]); A[c] = [F.mul(x, iv) for x in A[c]]
        for r in range(n):
            if r != c and not F.iszero(A[r][c]):
                f = A[r][c]; A[r] = [F.sub(x, F.mul(f, y)) for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]


def det(F, M):
    n = len(M); A = [list(r) for r in M]; d = F.one
    for c in range(n):
        piv = [r for r in range(c, n) if F.isunit(A[r][c])]
        if not piv: return F.zero
        r = piv[0]
        if r != c: A[c], A[r] = A[r], A[c]; d = F.neg(d)
        d = F.mul(d, A[c][c]); iv = F.inv(A[c][c])
        for r in range(c + 1, n):
            if not F.iszero(A[r][c]):
                f = F.mul(A[r][c], iv); A[r] = [F.sub(x, F.mul(f, y)) for x, y in zip(A[r], A[c])]
    return d


class Wick:
    def __init__(self, F, C):
        self.F, self.C, self.n = F, C, len(C)
        self.memo = {(0,) * self.n: F.one}

    def M(self, g):
        r = self.memo.get(g)
        if r is not None: return r
        F, C = self.F, self.C
        i = next(k for k in range(self.n) if g[k]); b = list(g); b[i] -= 1; s = F.zero
        for j in range(self.n):
            if b[j] and not F.iszero(C[i][j]):
                bb = list(b); bb[j] -= 1
                s = F.add(s, F.smul(b[j], F.mul(C[i][j], self.M(tuple(bb)))))
        self.memo[g] = s
        return s


def sp_series(F, Phi, amps, N, n):
    """Unnormalised [N_0..N_N] at the critical point 0 of the graded Taylor polynomial Phi (n variables).
    amps[m] = graded polynomial A_m. Returns (series, H)."""
    assert linear_is_zero(F, Phi), 'not a critical point'
    H = hessian(F, Phi, n); Cm = matinv(F, H); C = [[F.neg(x) for x in r] for r in Cm]
    W = Wick(F, C); P3 = gtrunc(Phi, 2 * N + 2, 3)
    out = [F.zero] * (N + 1)
    for m, A in enumerate(amps[:N + 1]):
        if not A: continue
        E = gtrunc(A, 2 * (N - m)); V = 0
        while E:
            for d, md in E.items():
                if d % 2: continue
                k = d // 2 - V + m
                if m <= k <= N:
                    s = F.zero
                    for e, c in md.items(): s = F.add(s, F.mul(c, W.M(e)))
                    out[k] = F.add(out[k], s)
            V += 1
            if V > 2 * (N - m): break
            E = gscale(F, gmul(F, E, P3, 2 * (N - m + V), 3 * V), F.of(Fr(1, V)))
    return out, H


def sp1(F, phi, amps, N):
    """one variable: phi = Taylor coefficients t^0..t^(2N+2) of the phase at the critical point (phi[1] must be 0),
    amps[m] = list of t-coefficients of A_m. Unnormalised [N_0..N_N]."""
    assert F.iszero(phi[1]), 'not a critical point'
    C = F.neg(F.inv(F.smul(2, phi[2])))
    mom = [F.one]
    for k in range(1, 3 * N + 2): mom.append(F.mul(mom[-1], F.smul(2 * k - 1, C)))
    P3 = [F.zero] * 3 + list(phi[3:2 * N + 3])
    EV = [[F.one]]                                   # EV[V] = P3^V/V!, needed up to degree 2(N+V)
    for V in range(1, 2 * N + 1):
        prev = EV[-1]; K = 2 * (N + V); c = [F.zero] * (K + 1)
        for i, x in enumerate(prev):
            if F.iszero(x): continue
            for j in range(3, min(len(P3), K + 1 - i)):
                if not F.iszero(P3[j]): c[i + j] = F.add(c[i + j], F.mul(x, P3[j]))
        EV.append([F.smul(Fr(1, V), x) for x in c])
    out = [F.zero] * (N + 1)
    for m, A in enumerate(amps[:N + 1]):
        for k in range(m, N + 1):
            s = F.zero
            for V in range(0, 2 * (k - m) + 1):
                d = 2 * (k - m + V); E = EV[V]
                for i in range(3 * V, min(d, len(E) - 1) + 1):
                    if d - i < len(A) and not F.iszero(E[i]) and not F.iszero(A[d - i]):
                        s = F.add(s, F.mul(F.mul(E[i], A[d - i]), mom[k - m + V]))
            out[k] = F.add(out[k], s)
    return out


def normalise(F, s):
    i0 = F.inv(s[0]); return [F.mul(x, i0) for x in s]
