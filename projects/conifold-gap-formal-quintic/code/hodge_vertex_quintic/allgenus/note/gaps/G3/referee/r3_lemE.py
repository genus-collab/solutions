# -*- coding: utf-8 -*-
"""r3_lemE.py -- referee check (independent code) of Lemma E in the ring AA = Q((x))[[u, psi']] (G3-3, O-8).

Different construction from g3_lemE.py: R is the CAYLEY transform R = (1 - X/2)^-1 (1 + X/2) of a random infinitesimally
symplectic X = Sum_(k=1..4) z^k X_k(u) (X(-z)^T eta + eta X(z) = 0), so R(0, z) has every power of z; the edge E_R is computed
by the DIFFERENCE QUOTIENT E_R = -Sum_N r_N Sum_(k<N) psi^(N-1-k) (-psi')^k eta^-1 R^-1(psi')^T (R^-1 = Sum r_N z^N), which uses
(Sp) in the form eta^-1 = R^-1(-psi') eta^-1 R^-1(psi')^T, and it is compared with the defining quotient N/(psi + psi').
Claim (vector form): Sum_(m,a) alpha_(a,m+1) [psi^m] E_R^(ab)(psi') = -(R^-1(psi') (R f)_(>=0)(psi'))_b,  alpha_(m+1) = [x^(-m-1)] eta f(-x).
f: z^-1 J^C(v(u)), z^-3 J^C(v(u)) (J^C = eta^-1 Y^C, v = u tau^(-5/4)), and a random f with [z^-n] f = O(u^floor(n/3)).
Controls: (C1) f with constant coefficients (no decay): the u^0 psi'^0 value of the contraction depends on the cut L;
(C2) R non-symplectic: the difference quotient differs from N/(psi + psi') and N(psi, -psi) != 0; (C3) f with an extra z^0 term
(f not in HH_-): the identity fails; (C4) sign psi' -> -psi' on the right: fails; (P) R^-1(0, x) is not a polynomial (all
coefficients up to the cut non-zero), so R^-1(x) x^(-j-1) is not in Q((x^-1))[[psi']].
Usage (working directory: the scratch folder):  py -I -B r3_lemE.py N M   -- log on stdout, nothing written, no path printed.
"""
import sys
import random
from fractions import Fraction as Fr
from math import factorial

NU = (0, -1, 3, 2, 1)
EPS = (1, -1, -1, 1, 1)
D5 = range(5)
ETA = [[Fr(5) if (a + b) % 5 == 3 else Fr(0) for b in D5] for a in D5]
ETAI = [[ETA[a][b] / 25 for b in D5] for a in D5]


def smul(a, b):
    N = len(a)
    out = [Fr(0)] * N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                if b[j]:
                    out[i + j] += x * b[j]
    return out


def sadd(a, b, c=1):
    return [x + c * y for x, y in zip(a, b)]


def Z(N):
    return [Fr(0)] * N


def mz(N):
    return [[Z(N) for _ in D5] for _ in D5]


def mconst(C, N):
    return [[[C[i][j]] + [Fr(0)] * (N - 1) for j in D5] for i in D5]


def mmul(A, B):
    N = len(A[0][0])
    out = mz(N)
    for i in D5:
        for k in D5:
            if any(A[i][k]):
                for j in D5:
                    if any(B[k][j]):
                        out[i][j] = sadd(out[i][j], smul(A[i][k], B[k][j]))
    return out


def madd(A, B, c=1):
    return [[sadd(A[i][j], B[i][j], c) for j in D5] for i in D5]


def mT(A):
    return [[A[j][i] for j in D5] for i in D5]


def mv(A, x):
    N = len(x[0])
    out = [Z(N) for _ in D5]
    for i in D5:
        for k in D5:
            if any(A[i][k]) and any(x[k]):
                out[i] = sadd(out[i], smul(A[i][k], x[k]))
    return out


def zmul(P, Q, KZ):
    N = len(P[0][0][0])
    out = [mz(N) for _ in range(KZ + 1)]
    for i, A in enumerate(P):
        for j, B in enumerate(Q):
            if i + j <= KZ:
                out[i + j] = madd(out[i + j], mmul(A, B))
    return out


def cayley(X, KZ, sign=1):
    """R = (1 - sX/2)^-1 (1 + sX/2) mod z^(KZ+1), s = sign, from (1 - sX/2) R = 1 + sX/2:
    R_0 = 1, R_k = s X_k/2 + Sum_(i=1..4) (s X_i/2) R_(k-i)."""
    N = len(X[1][0][0])
    half = [[[[sign * e / 2 for e in s] for s in row] for row in A] for A in X]
    R = [mconst([[Fr(int(i == j)) for j in D5] for i in D5], N)]
    for k in range(1, KZ + 1):
        acc = half[k] if k < len(half) else mz(N)
        for i in range(1, min(len(half) - 1, k) + 1):
            acc = madd(acc, mmul(half[i], R[k - i]))
        R.append(acc)
    return R


def rand_X(N, rng):
    X = [mz(N)]
    for k in range(1, 5):
        S = [[None] * 5 for _ in D5]
        for a in D5:
            for b in range(a, 5):
                s = [Fr(rng.randint(-3, 3)) for _ in range(N)]
                if a == b and k % 2 == 0:
                    s = Z(N)
                S[a][b] = s
                S[b][a] = s if k % 2 == 1 else [-e for e in s]
        X.append(mmul(mconst(ETAI, N), S))
    return X


def model_v(N):
    F = [Fr(4, 5 * n + 4) for n in range(N + 1)]
    g = [Fr(1)]
    for n in range(N - 1):
        rhs = Fr(5, 4) * sum((i + 1) * F[i + 1] * g[n - i] for i in range(n + 1))
        rhs -= sum(F[i] * (n - i + 1) * g[n - i + 1] for i in range(1, n + 1))
        g.append(rhs / (n + 1))
    return [Fr(0)] + g[:N - 1]


def yc1(b, n):
    m, tot, vp = NU[b], Fr(0), None
    for i in range(n + 1):
        l = m + 4 * i + 1
        if l >= 0 and i + l == n:
            tot += Fr((-1) ** i, 4 ** i * factorial(i)) * Fr((-1) ** l, factorial(l))
            vp = l
    return EPS[b] * tot, vp


def jc(N, shift, v):
    """f = z^-shift J^C(v(u)) as {n: vector}; J^C_a = Y^C_(3-a)/5."""
    pw = [[Fr(1)] + [Fr(0)] * (N - 1)]
    f = {}
    for n in range(0, 2 * N + 2):
        vec = []
        for a in D5:
            c, e = yc1((3 - a) % 5, n)
            while e is not None and len(pw) <= e:
                pw.append(smul(pw[-1], v))
            vec.append([c / 5 * x for x in pw[e]] if c else Z(N))
        if any(any(x) for x in vec):
            f[n + shift] = vec
    return f
