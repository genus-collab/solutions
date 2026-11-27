# -*- coding: utf-8 -*-
"""g3_lemE.py -- exact check over Q (mod u^N) of Lemma E in the ring Q((x))[[u, psi']] (gap cluster G3, item 3, O-8).

R = exp(X(z)), X = Sum_(k=1..3) z^k X_k(u) infinitesimally symplectic (X(-z)^T eta + eta X(z) = 0) with random integer
entries, so that R(0, z) has EVERY power of z (the flat situation, R not a polynomial at any u-order); f in
z^-1 Q[z^-1][[u]] (x) H infinite in z^-1 with u-adic decay: (i) f = z^-1 J^C(v(u)), J^C = eta^-1 Y^C, v(u) = u tau^(-5/4);
(ii) f = z^-2 J^C; (iii) a random decaying f, [z^-n] f = O(u^ceil(n/2)).  Claim (Lemma E, vector form):
    (Omega_s(f, .) (x) 1) E_R(psi') = - R(psi')^-1 (R f)_(>=0)(psi'),
E_R = (eta^-1 - R^-1(psi) eta^-1 R^-1(psi')^T)/(psi + psi'), contraction of psi^m with eta([z^(-m-1)] f(-z), .).
Both sides are computed mod u^N and mod psi'^(M+1) with every z-truncation large enough to be exact (s.4 of GAPS_G3_4).
Controls: (C1) a NON-decaying f (constant coefficients, z^-1..z^-L): the contraction at u^0 psi'^0 changes with L (the
sum over psi-powers has no limit; the hypothesis f in HH^u_- is needed); (C2) a non-symplectic R: the numerator of E_R is
not divisible by psi + psi'; (C3) the sign psi' -> -psi' on the right side: the identity fails.
Usage: py -I -B g3_lemE.py N M      (log to stdout; nothing written)
"""
import os
import sys
import random
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g3_series as S
from g3_frames import yc_closed
from g3_mat import (D5, mzero, mid, mmul, madd, mT, mvec, const_mat, ETA, ETAI, pz_mul, pz_exp, rand_X)


def lemE_sides(R, Ri, f, N, M):
    """LHS_n and RHS_n (n <= M) of Lemma E; f = {n: vector}, n >= 1."""
    nf = max(f)
    EI = const_mat(ETAI, N)
    Nn = {}
    for m in range(nf + M + 2):
        for n in range(M + 1):
            t = mmul(mmul(Ri[m], EI), mT(Ri[n]))
            t = [[S.scal(-1, e) for e in row] for row in t]
            if m == 0 and n == 0:
                t = madd(t, EI)
            Nn[(m, n)] = t
    lhs, rhs = [], []
    for n in range(M + 1):
        tot = [S.zero(N) for _ in D5]
        for m in range(nf):
            E = mzero(N)
            for j in range(n + 1):
                E = madd(E, Nn[(m + 1 + j, n - j)], (-1) ** j)
            fm = f.get(m + 1)
            if fm is None:
                continue
            alpha = mvec(const_mat(ETA, N), [S.scal((-1) ** (m + 1), x) for x in fm])
            for b in D5:
                for a in D5:
                    tot[b] = S.add(tot[b], S.mul(alpha[a], E[a][b]))
        lhs.append(tot)
    Rf = []
    for j in range(M + 1):
        acc = [S.zero(N) for _ in D5]
        for n, vec in f.items():
            if j + n < len(R):
                w = mvec(R[j + n], vec)
                acc = [S.add(x, y) for x, y in zip(acc, w)]
        Rf.append(acc)
    for n in range(M + 1):
        acc = [S.zero(N) for _ in D5]
        for l in range(n + 1):
            w = mvec(Ri[n - l], Rf[l])
            acc = [S.add(x, S.scal(-1, y)) for x, y in zip(acc, w)]
        rhs.append(acc)
    return lhs, rhs, Rf, Nn


def jc_f(N, shift, VP):
    """f = z^-shift J^C(v(u)) as {n: vector}: J^C,a = Y^C_(3-a)/5."""
    f = {}
    for n in range(0, 4 * N):
        vec = []
        for a in D5:
            b = (3 - a) % 5
            c = yc_closed(b, n)
            vec.append(S.scal(c / 5, VP.trunc(n - n // 5, N)) if c else S.zero(N))
        if any(any(x) for x in vec):
            f[n + shift] = vec
    return f


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    rng = random.Random(1115)
    VP = S.VPowers(N + 2)
    print("g3_lemE.py  N = %d (mod u^N), M = %d (psi' order); exact over Q" % (N, M))
    fs = {"z^-1 J^C": jc_f(N, 1, VP), "z^-2 J^C": jc_f(N, 2, VP)}
    fr = {}
    for n in range(1, 2 * N):
        vec = [S.zero(N) for _ in D5]
        for a in D5:
            for k in range((n + 1) // 2, N):
                vec[a][k] = Fr(rng.randint(-3, 3))
        fr[n] = vec
    fs["random decaying"] = fr
    nf = max(max(f) for f in fs.values())
    KZ = nf + M + 2
    X = rand_X(N, 3, rng)
    R = pz_exp(X, KZ)
    Ri = pz_exp([[[S.scal(-1, e) for e in row] for row in A] for A in X], KZ)
    Rm = [mT(A) if k % 2 == 0 else [[S.scal(-1, e) for e in row] for row in mT(A)] for k, A in enumerate(R)]
    sp = pz_mul(pz_mul(Rm, [const_mat(ETA, N)], KZ), R, KZ)
    sp_ok = all(sp[k] == (const_mat(ETA, N) if k == 0 else mzero(N)) for k in range(KZ + 1))
    nz0 = sum(1 for k in range(1, KZ + 1) if any(R[k][i][j][0] for i in D5 for j in D5))
    print("[Sp] R(-z)^T eta R(z) = eta through z^%d, mod u^%d: %s;  z-powers of R(0,z) that are non-zero: %d of %d"
          % (KZ, N, sp_ok, nz0, KZ))
    ok = sp_ok and nz0 == KZ
    for name, f in fs.items():
        lhs, rhs, Rf, Nn = lemE_sides(R, Ri, f, N, M)
        good = lhs == rhs
        nonzero = sum(1 for n in range(M + 1) for b in D5 if any(lhs[n][b]))
        print("[E] f = %-16s (terms z^-1..z^-%d): LHS = RHS on all coefficients psi'^0..psi'^%d mod u^%d: %s (non-zero "
              "components %d of %d)" % (name, max(f), M, N, good, nonzero, 5 * (M + 1)))
        ok &= good and nonzero > 0
        rhs_c = []
        for n in range(M + 1):
            acc = [S.zero(N) for _ in D5]
            for l in range(n + 1):
                w = mvec(Ri[n - l], [S.scal((-1) ** l, x) for x in Rf[l]])
                acc = [S.add(x, S.scal(-1, y)) for x, y in zip(acc, w)]
            rhs_c.append(acc)
        print("[C3] same f, sign psi' -> -psi' in (R f)_(>=0): identity holds: %s (expected False)" % (lhs == rhs_c))
        ok &= lhs != rhs_c
    divis = all(not any(any(e) for row in madd_sum(Nn, D) for e in row) for D in range(M + 1))
    print("[DIV] numerator of E_R divisible by psi + psi' (N(psi,-psi) = 0 through degree %d): %s" % (M, divis))
    ok &= divis
    ok &= controls(rng)
    print("RESULT:", "ALL CHECKS PASSED" if ok else "FAILURE")


def madd_sum(Nn, D):
    out = mzero(len(Nn[(0, 0)][0][0]))
    for m in range(D + 1):
        if (m, D - m) in Nn:
            out = madd(out, Nn[(m, D - m)], (-1) ** (D - m))
    return out


def controls(rng):
    N, M = 1, 1
    X = rand_X(N, 3, rng)
    ok = True
    vals = []
    for L in (8, 13, 18, 23, 28):
        KZ = L + M + 2
        R = pz_exp(X, KZ)
        Ri = pz_exp([[[S.scal(-1, e) for e in row] for row in A] for A in X], KZ)
        c = [[Fr(1)], [Fr(2)], [Fr(-1)], [Fr(3)], [Fr(1)]]
        f = {n: [x[:] for x in c] for n in range(1, L + 1)}
        lhs, rhs, _, _ = lemE_sides(R, Ri, f, N, M)
        vals.append(tuple(lhs[0][b][0] for b in D5))
        ok &= lhs == rhs
    distinct = len(set(vals))
    print("[C1] NON-decaying f = Sum_(n=1..L) z^-n c: [u^0 psi'^0] of the contraction for L = 8, 13, 18, 23, 28 takes %d "
          "distinct values (first two: %s, %s): no limit as L grows (expected); each truncation satisfies the finite identity: %s"
          % (distinct, [str(x) for x in vals[0][:2]], [str(x) for x in vals[1][:2]], ok))
    good = distinct == len(vals) and ok
    X2 = rand_X(2, 2, rng)
    R = pz_exp(X2, 8)
    E1 = mzero(2)
    E1[0][1][0] = Fr(1)
    Rbad = pz_mul(R, [mid(2), E1], 8)
    inv = pz_exp([[[S.scal(-1, e) for e in row] for row in A] for A in X2], 8)
    E1n = [mid(2), [[S.scal(-1, e) for e in row] for row in E1]]
    Ri = pz_mul(E1n, inv, 8)
    # (1 + z E)^-1 = 1 - z E exactly, since E^2 = 0 for E = E_01
    Nn = {}
    EI = const_mat(ETAI, 2)
    for m in range(9):
        for n in range(9 - m):
            t = mmul(mmul(Ri[m], EI), mT(Ri[n]))
            t = [[S.scal(-1, e) for e in row] for row in t]
            if m == 0 and n == 0:
                t = madd(t, EI)
            Nn[(m, n)] = t
    nondiv = sum(1 for D in range(1, 8) if any(any(e) for row in madd_sum(Nn, D) for e in row))
    print("[C2] non-symplectic R (1 + z E_01) exp(X): N(psi, -psi) != 0 in %d of 7 degrees 1..7 (expected > 0)" % nondiv)
    return good and nondiv > 0


if __name__ == "__main__":
    main()
