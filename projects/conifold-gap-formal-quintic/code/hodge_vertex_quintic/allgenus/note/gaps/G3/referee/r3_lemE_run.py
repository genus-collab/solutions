# -*- coding: utf-8 -*-
"""r3_lemE_run.py -- driver of the referee check of Lemma E (helpers and the description of the checks: r3_lemE.py).
Usage (working directory: the scratch folder):  py -I -B r3_lemE_run.py N M   -- log on stdout, nothing written, no path printed.
"""
import os
import sys
import random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r3_lemE import (Fr, D5, ETA, ETAI, Z, mz, mconst, mmul, madd, mT, mv, sadd, smul, cayley, rand_X, model_v, jc)


def nmat(Ri, p, q, N):
    """[psi^p psi'^q] of eta^-1 - R^-1(psi) eta^-1 R^-1(psi')^T."""
    t = mmul(mmul(Ri[p], mconst(ETAI, N)), mT(Ri[q]))
    t = [[[-e for e in s] for s in row] for row in t]
    return madd(t, mconst(ETAI, N)) if p == q == 0 else t


def edges(Ri, nm, M, N):
    """three expansions of E_R: difference quotient (uses (Sp)), N/(psi+psi') expanded in psi'/psi and in psi/psi'."""
    W = [mmul(mconst(ETAI, N), mT(Ri[l])) for l in range(M + 2)]
    E1, E2, E3 = {}, {}, {}
    for m in range(nm):
        for n in range(M + 1):
            acc = mz(N)
            for k in range(n + 1):
                acc = madd(acc, mmul(Ri[m + k + 1], W[n - k]), -((-1) ** k))
            E1[(m, n)] = acc
            E2[(m, n)] = mz(N)
            for j in range(n + 1):
                E2[(m, n)] = madd(E2[(m, n)], nmat(Ri, m + 1 + j, n - j, N), (-1) ** j)
            E3[(m, n)] = mz(N)
            for j in range(m + 1):
                E3[(m, n)] = madd(E3[(m, n)], nmat(Ri, m - j, n + 1 + j, N), (-1) ** j)
    return E1, E2, E3


def sides(R, Ri, E, f, M, N, sgn=1):
    """LHS: contraction of alpha_(m+1) = eta (-1)^(m+1) f_(m+1) with [psi^m] E; RHS: -R^-1(psi')(R f)_(>=0)(sgn psi');
    STEP4: -R^-1(psi') eta^-1 beta_(>=0)(-psi'), beta_j = Sum_m (r_(j+m))^T alpha_m (the intermediate formula of the proof)."""
    lhs, rhs, st4 = [], [], []
    alpha = {m: mv(mconst(ETA, N), [[((-1) ** m) * e for e in x] for x in vec]) for m, vec in f.items() if m >= 1}
    Rf = []
    for j in range(M + 1):
        acc = [Z(N) for _ in D5]
        for n, vec in f.items():
            if 0 <= j + n < len(R):
                acc = [sadd(x, y) for x, y in zip(acc, mv(R[j + n], vec))]
        Rf.append(acc)
    beta = []
    for j in range(M + 1):
        acc = [Z(N) for _ in D5]
        for m, al in alpha.items():
            acc = [sadd(x, y) for x, y in zip(acc, mv(mT(Ri[j + m]), al))]
        beta.append(acc)
    for n in range(M + 1):
        tot = [Z(N) for _ in D5]
        for (m, nn), Em in E.items():
            if nn == n and (m + 1) in alpha:
                tot = [sadd(t, y) for t, y in zip(tot, mv(mT(Em), alpha[m + 1]))]
        lhs.append(tot)
        r = [Z(N) for _ in D5]
        s4 = [Z(N) for _ in D5]
        for l in range(n + 1):
            r = [sadd(x, y, -1) for x, y in zip(r, mv(Ri[n - l], [[(sgn ** l) * e for e in w] for w in Rf[l]]))]
            s4 = [sadd(x, y, -1) for x, y in zip(s4, mv(mmul(Ri[n - l], mconst(ETAI, N)), [[((-1) ** l) * e for e in w]
                                                                                          for w in beta[l]]))]
        rhs.append(r)
        st4.append(s4)
    return lhs, rhs, st4


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    rng = random.Random(77)
    print("r3_lemE_run.py  N = %d (mod u^N), M = %d (psi' order); exact over Q" % (N, M))
    v = model_v(N)
    fs = {"z^-1 J^C": jc(N, 1, v), "z^-3 J^C": jc(N, 3, v)}
    fr = {n: [[Fr(rng.randint(-4, 4)) if k >= n // 3 else Fr(0) for k in range(N)] for _ in D5] for n in range(1, 3 * N)}
    fs["random decaying"] = fr
    nf = max(max(f) for f in fs.values())
    KZ = nf + M + 3
    X = rand_X(N, rng)
    R, Ri = cayley(X, KZ, 1), cayley(X, KZ, -1)
    eta = mconst(ETA, N)
    defect = 0
    for k in range(KZ + 1):
        acc = mz(N)
        for i in range(k + 1):
            acc = madd(acc, mmul(mmul(mT(R[i]), eta), R[k - i]), (-1) ** i)
        if acc != (eta if k == 0 else mz(N)):
            defect += 1
    prod_ok = True
    for k in range(KZ + 1):
        acc = mz(N)
        for i in range(k + 1):
            acc = madd(acc, mmul(R[i], Ri[k - i]))
        prod_ok &= acc == (mconst([[Fr(int(i == j)) for j in D5] for i in D5], N) if k == 0 else mz(N))
    nzx = sum(1 for k in range(1, KZ + 1) if any(Ri[k][i][j][0] for i in D5 for j in D5))
    print("[Sp] R(-z)^T eta R(z) = eta through z^%d mod u^%d: defects %d; R R^-1 = 1: %s;  [P] non-zero x^d coefficients of "
          "R^-1(0, x), 1 <= d <= %d: %d" % (KZ, N, defect, prod_ok, KZ, nzx))
    ok = defect == 0 and prod_ok and nzx == KZ
    E1, E2, E3 = edges(Ri, nf, M, N)
    print("[DIV] N/(psi+psi') expanded in psi'/psi = expanded in psi/psi' (divisibility, (Sp)): %s; = difference quotient: %s"
          % (E2 == E3, E1 == E3))
    ok &= E1 == E2 == E3
    for name, f in fs.items():
        lhs, rhs, st4 = sides(R, Ri, E1, f, M, N)
        nz = sum(1 for n in range(M + 1) for b in D5 if any(lhs[n][b]))
        lhs4, rhs4, _ = sides(R, Ri, E1, f, M, N, sgn=-1)
        f0 = dict(f)
        f0[0] = [[Fr(1)] + [Fr(0)] * (N - 1) for _ in D5]
        l0, r0, _ = sides(R, Ri, E1, f0, M, N)
        print("[E] f = %-15s LHS = RHS: %s (non-zero components %d of %d); Step-4 form = LHS: %s; [C4] sign flip holds: %s; "
              "[C3] extra z^0 term holds: %s" % (name, lhs == rhs, nz, 5 * (M + 1), st4 == lhs, lhs4 == rhs4, l0 == r0))
        ok &= lhs == rhs and nz > 0 and st4 == lhs and lhs4 != rhs4 and l0 != r0
    vals = []
    X1 = rand_X(1, rng)
    for L in (6, 9, 12, 15):
        R1, Ri1 = cayley(X1, L + 3, 1), cayley(X1, L + 3, -1)
        fc = {n: [[Fr(c)] for c in (1, 2, -1, 3, 1)] for n in range(1, L + 1)}
        Ec, _, _ = edges(Ri1, L, 0, 1)
        lc, rc, _ = sides(R1, Ri1, Ec, fc, 0, 1)
        vals.append(tuple(lc[0][b][0] for b in D5))
        ok &= lc == rc
    print("[C1] non-decaying f, cut L = 6, 9, 12, 15: %d distinct values of the u^0 psi'^0 contraction (expected 4); each cut "
          "satisfies the finite identity" % len(set(vals)))
    ok &= len(set(vals)) == 4
    Xb = rand_X(2, rng)
    Rib = cayley(Xb, 12, -1)
    Rib[1] = madd(Rib[1], mconst([[Fr(int(i == 0 and j == 1)) for j in D5] for i in D5], 2))   # break (Sp)
    B1, B2, B3 = edges(Rib, 4, 3, 2)
    print("[C2] non-symplectic R^-1: the two expansions of N/(psi+psi') agree (divisibility): %s (expected False); difference "
          "quotient = psi'/psi expansion: %s (True for every R: the psi-free defect eta^-1 - R^-1(-psi') eta^-1 R^-1(psi')^T "
          "contributes only negative powers of psi), = psi/psi' expansion: %s (expected False)" % (B2 == B3, B1 == B2, B1 == B3))
    ok &= B2 != B3 and B1 == B2 and B1 != B3
    print("RESULT:", "ALL CHECKS PASSED" if ok else "FAILURE")


if __name__ == "__main__":
    main()
