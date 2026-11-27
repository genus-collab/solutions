# -*- coding: utf-8 -*-
"""rf_lemE.py -- referee check of Lemma E (REDUCTION_DRAFT_2.md s.3.1), exact over Q, own code (no worker module, no project module).

Lemma E: R(z) = 1 + O(z) symplectic (R(-z)^T eta R(z) = eta, vector form), E_R(psi, psi1) = (eta^-1 - R^-1(psi) eta^-1 R^-1(psi1)^T)/(psi + psi1),
f in z^-1 H[[z^-1]].  Contracting Omega(f, .) (Omega(f, z^m e) = eta([z^(-m-1)] f(-z), e)) into slot 1 of E_R gives
        C(psi1) = - R(psi1)^-1 (R f)_{>=0}(psi1).
Test: random R = exp(X), X(z) = Sum_{k>=1} eta^-1 S_k z^k, S_k symmetric (k odd) / antisymmetric (k even); random f with F terms; all
coefficients psi1^n that are complete in the truncation are compared.  Also the Z/5-graded case, the slot-2 contraction, and two
controls: (C1) a NON-symplectic R (the numerator does not vanish on psi = -psi1, so E_R is not a power series: the divisibility test
fails), (C2) the formula with the sign of the psi1-power flipped (must NOT agree).
Run: py -I -B rf_lemE.py [seed]   Log: rf_lemE.log (appended).
"""
import os, sys, random, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
ETA = [[Fr(5) if (a + b) % 5 == 3 else Fr(0) for b in range(5)] for a in range(5)]
ETAI = [[Fr(1, 5) if (a + b) % 5 == 3 else Fr(0) for b in range(5)] for a in range(5)]


def mz(): return [[Fr(0)] * 5 for _ in range(5)]
def mid(): return [[Fr(int(a == b)) for b in range(5)] for a in range(5)]
def madd(A, B, s=1): return [[A[i][j] + s * B[i][j] for j in range(5)] for i in range(5)]
def mmul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
def mT(A): return [[A[j][i] for j in range(5)] for i in range(5)]
def msc(A, c): return [[c * x for x in r] for r in A]
def mv(A, x): return [sum(A[i][k] * x[k] for k in range(5)) for i in range(5)]
def iszero(A): return all(x == 0 for r in A for x in r)


def ser_mul(A, B, N):
    """matrix series (lists of matrices, index = power of z) product mod z^(N+1)"""
    out = [mz() for _ in range(N + 1)]
    for i, a in enumerate(A):
        if i > N or iszero(a): continue
        for j, b in enumerate(B):
            if i + j > N: break
            out[i + j] = madd(out[i + j], mmul(a, b))
    return out


def ser_exp(X, N):
    """exp of a matrix series with X[0] = 0, mod z^(N+1)"""
    out = [mid()] + [mz() for _ in range(N)]
    term = [mid()] + [mz() for _ in range(N)]
    for j in range(1, N + 1):
        term = ser_mul(term, X, N); term = [msc(t, Fr(1, j)) for t in term]
        out = [madd(o, t) for o, t in zip(out, term)]
    return out


def ser_inv(A, N):
    """inverse of a matrix series with A[0] = 1, mod z^(N+1)"""
    out = [mid()] + [mz() for _ in range(N)]
    for n in range(1, N + 1):
        acc = mz()
        for k in range(1, n + 1):
            if k < len(A): acc = madd(acc, mmul(A[k], out[n - k]))
        out[n] = msc(acc, -1)
    return out


def rand_X(rng, N, graded):
    X = [mz()]
    for k in range(1, N + 1):
        S = mz()
        for a in range(5):
            for b in range(a, 5):
                x = Fr(rng.randint(-9, 9), rng.randint(1, 7))
                if k % 2 == 0 and a == b: x = Fr(0)
                S[a][b] = x; S[b][a] = x if k % 2 == 1 else -x
        if graded:      # vector form X_k = eta^-1 S_k supported on (a, b), a = b + k: S supported on (3 - a, b), i.e. a' + b = 3 - k
            S = [[S[a][b] if (a + b + k) % 5 == 3 else Fr(0) for b in range(5)] for a in range(5)]
        X.append(mmul(ETAI, S))
    return X


def sympl_defect(R, N):
    """first order k <= N at which R(-z)^T eta R(z) - eta is non-zero (None if zero)"""
    Rm = [msc(mT(R[k]), (-1) ** k) for k in range(len(R))]
    S = ser_mul([mmul(r, ETA) for r in Rm], R, N)
    S[0] = madd(S[0], ETA, -1)
    bad = [k for k in range(N + 1) if not iszero(S[k])]
    return bad[0] if bad else None


def edge(R, N):
    """E_{m,n} (matrices) for m + n <= N - 1, R^-1, and the divisibility defect"""
    Ri = ser_inv(R, N)
    Num = {}
    for m in range(N + 1):
        for n in range(N + 1 - m):
            Num[(m, n)] = msc(mmul(mmul(Ri[m], ETAI), mT(Ri[n])), -1)
    Num[(0, 0)] = madd(Num[(0, 0)], ETAI)
    E = {}
    for m in range(N):
        for n in range(N - m):
            acc = mz()
            for t in range(n + 1): acc = madd(acc, Num[(m + 1 + t, n - t)], (-1) ** t)
            E[(m, n)] = acc
    # divisibility: [psi^0 psi1^n] of (psi + psi1) E = E_{0,n-1} must equal Num_{0,n}; and Num_{0,0} = 0
    div_bad = [n for n in range(1, N) if not iszero(madd(Num[(0, n)], E[(0, n - 1)], -1))]
    if not iszero(Num[(0, 0)]): div_bad = [0] + div_bad
    return E, Ri, div_bad


def contract(E, f, N, slot):
    """C(psi1)_n = Sum_{m,a} Omega(f, z^m e_a) E_{mn}^{ab} (slot 1), resp. with E_{nm}^{ba} (slot 2)"""
    F = max(f); out = {}
    for n in range(N - F):
        vec = [Fr(0)] * 5
        for m in range(F):
            if (m + 1) not in f: continue
            cov = [(-1) ** (m + 1) * x for x in mv(ETA, f[m + 1])]
            M = E[(m, n)] if slot == 1 else mT(E[(n, m)])
            for b in range(5): vec[b] += sum(cov[a] * M[a][b] for a in range(5))
        out[n] = vec
    return out


def predicted(R, Ri, f, N, flip=False):
    F = max(f); Rf = {}
    for j in range(0, N - F):
        acc = [Fr(0)] * 5
        for k in range(j + 1, min(len(R), j + F + 1)):
            if (k - j) in f: acc = [x + y for x, y in zip(acc, mv(R[k], f[k - j]))]
        Rf[j] = acc
    out = {}
    for n in range(N - F):
        vec = [Fr(0)] * 5
        for j in range(n + 1):
            sg = (-1) ** (n - j) if flip else 1
            vec = [x - sg * y for x, y in zip(vec, mv(Ri[n - j], Rf[j]))]
        out[n] = vec
    return out


def one_test(rng, N, F, graded, log, symplectic=True):
    if symplectic: R = ser_exp(rand_X(rng, N, graded), N)
    else: R = [mid()] + [[[Fr(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(5)] for _ in range(5)] for _ in range(N)]
    sd = sympl_defect(R, N)
    E, Ri, div_bad = edge(R, N)
    f = {n: [Fr(rng.randint(-9, 9), rng.randint(1, 9)) for _ in range(5)] for n in range(1, F + 1)}
    c1 = contract(E, f, N, 1); c2 = contract(E, f, N, 2); pr = predicted(R, Ri, f, N); pf = predicted(R, Ri, f, N, flip=True)
    ok1 = all(c1[n] == pr[n] for n in c1); ok2 = all(c2[n] == pr[n] for n in c2)
    nz = sum(1 for n in c1 if any(c1[n])); okf = all(c1[n] == pf[n] for n in c1)
    log("   graded=%s symplectic=%s N=%d F=%d: sympl. defect first order %s, divisibility defect %s; slot1 = formula: %s, slot2 = formula: %s "
        "(%d of %d compared coefficients non-zero); control C2 (flipped sign) agrees: %s" % (graded, symplectic, N, F, sd, div_bad[:3], ok1, ok2, nz, len(c1), okf))
    return ok1, ok2, okf, div_bad, sd


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261113
    rng = random.Random(seed); LOG = open(os.path.join(HERE, "rf_lemE.log"), "a", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    t0 = time.time(); log("=" * 100); log("rf_lemE seed %d [%s]" % (seed, time.strftime("%Y-%m-%dT%H:%M:%S")))
    allok = True
    for graded in (False, True):
        for (N, F) in ((8, 3), (10, 4), (9, 6)):
            ok1, ok2, okf, db, sd = one_test(rng, N, F, graded, log)
            allok = allok and ok1 and ok2 and (not okf) and not db and sd is None
    ok1, ok2, okf, db, sd = one_test(rng, 8, 3, False, log, symplectic=False)
    log("   control C1 (non-symplectic R): divisibility defect present: %s (expected True)" % bool(db))
    allok = allok and bool(db)
    log("RF_LEME: %s  (%.1f s)" % ("ALL PASS" if allok else "FAIL", time.time() - t0)); LOG.close()


if __name__ == "__main__":
    main()
