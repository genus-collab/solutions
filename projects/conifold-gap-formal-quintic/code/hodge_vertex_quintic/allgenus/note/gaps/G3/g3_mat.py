# -*- coding: utf-8 -*-
"""g3_mat.py -- 5x5 matrices over truncated u-series (exact over Q), matrix polynomials in z, exp, random infinitesimally
symplectic X(z) (X(-z)^T eta + eta X(z) = 0: X_k = eta^-1 S_k, S_k symmetric for k odd, antisymmetric for k even).
Imported by g3_lemE.py (gap cluster G3).  Nothing here writes a file.
"""
import os
import sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g3_series as S

D5 = range(5)


def mzero(N):
    return [[S.zero(N) for _ in D5] for _ in D5]


def mid(N):
    m = mzero(N)
    for a in D5:
        m[a][a] = S.one(N)
    return m


def mmul(A, B):
    N = len(A[0][0])
    out = mzero(N)
    for i in D5:
        for k in D5:
            if any(A[i][k]):
                for j in D5:
                    if any(B[k][j]):
                        out[i][j] = S.add(out[i][j], S.mul(A[i][k], B[k][j]))
    return out


def madd(A, B, c=1):
    return [[S.add(A[i][j], S.scal(c, B[i][j])) for j in D5] for i in D5]


def mT(A):
    return [[A[j][i] for j in D5] for i in D5]


def mvec(A, x):
    N = len(x[0])
    out = [S.zero(N) for _ in D5]
    for i in D5:
        for k in D5:
            if any(A[i][k]) and any(x[k]):
                out[i] = S.add(out[i], S.mul(A[i][k], x[k]))
    return out


def const_mat(rows, N):
    m = mzero(N)
    for i in D5:
        for j in D5:
            m[i][j][0] = Fr(rows[i][j])
    return m


ETA = [[5 if (a + b) % 5 == 3 else 0 for b in D5] for a in D5]
ETAI = [[Fr(ETA[a][b], 25) for b in D5] for a in D5]


def pz_mul(P, Q, KZ):
    """product of matrix polynomials in z (lists of matrices), truncated at z^KZ."""
    N = len(P[0][0][0])
    out = [mzero(N) for _ in range(KZ + 1)]
    for i, A in enumerate(P):
        for j, B in enumerate(Q):
            if i + j <= KZ:
                out[i + j] = madd(out[i + j], mmul(A, B))
    return out


def pz_exp(X, KZ):
    N = len(X[1][0][0])
    out = [mid(N)] + [mzero(N) for _ in range(KZ)]
    term = [mid(N)] + [mzero(N) for _ in range(KZ)]
    for m in range(1, KZ + 1):
        term = pz_mul(term, X, KZ)
        term = [[[S.scal(Fr(1, m), e) for e in row] for row in A] for A in term]
        out = [madd(a, b) for a, b in zip(out, term)]
    return out


def rand_X(N, KX, rng, sign=1):
    X = [mzero(N)]
    for k in range(1, KX + 1):
        Sm = [[None] * 5 for _ in D5]
        for a in D5:
            for b in D5:
                if b < a:
                    continue
                s = [Fr(rng.randint(-2, 2)) for _ in range(N)]
                Sm[a][b] = s
                Sm[b][a] = s if k % 2 == 1 else S.scal(-1, s)
                if a == b and k % 2 == 0:
                    Sm[a][a] = S.zero(N)
        Xk = mmul(const_mat(ETAI, N), Sm)
        X.append([[S.scal(sign, e) for e in row] for row in Xk])
    return X
