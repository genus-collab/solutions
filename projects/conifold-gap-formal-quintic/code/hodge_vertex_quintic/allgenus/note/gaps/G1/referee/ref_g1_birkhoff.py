# -*- coding: utf-8 -*-
"""ref_g1_birkhoff.py -- referee evidence (not proof) for prop:set-cyclicLR (gap cluster G1, item 2).

Independent re-implementation (the prover's g1_cyclic_check.py uses a w-shift trick; here the formal log ell = log q is kept explicitly).
Exact rationals, stdlib only.  Run:  py -I -B ref_g1_birkhoff.py [N] [K]   (N = q-order, default 16; K = number of periods, default 12).
Writes ref_g1_birkhoff.log next to itself (relative name printed).

(1) Periods c_k in Q[[q]][ell]: I/z = Sum_k c_k H^(k mod 5) z^-k at the zeta point, c_k = Sum_j ell^j/j! [w^(k-j)] Fhat,
    Fhat = Sum_d q^d prod_(m=1..5d)(5w+m)/prod_(m=1..d)((w+m)^5 - w^5).
(2) Birkhoff recursion of GAPS_G1_2 Lemma C, with ell explicit: phi^(0)_k = c_k/C_0, kappa_a = d_T phi^(a)_(a+1),
    phi^(a+1)_k = d_T phi^(a)_k/kappa_a, d_T = C_1^-1 D, D ell = 1.  Checks: every kappa_a is ell-free and
    kappa = (1, L^5/(C_0^2 C_1^3), 1, C_0/C_1, C_0/C_1); closing phi^(5)_(k+5) = phi^(0)_k (H^5 = 1).
(3) Operator identity of step 3 of GAPS_G1_3 s.2.8: d_T^2 o kappa_1^-1 o d_T^2 o C_0^-1 = (C_0/C_1) PF4 as differential operators with
    Q[[q]] coefficients (kappa_1 = L^5/(C_0^2 C_1^3) inserted in closed form), coefficient by coefficient in D^0..D^4.
Controls that MUST fail: (c1) non-equivariant denominators (w+m)^5: closing fails; (c2) twist O(4) (numerator prod_(m<=4d)(4w+m)):
kappa_1 = L^5/(C_0^2 C_1^3) fails; (c3) operator identity with kappa_1 replaced by 1/(C_0^2 C_1^3) fails; (c4) kappa_3 = C_1/C_0 fails.
"""
import os, sys
from fractions import Fraction as Fr
from math import factorial

HERE = os.path.dirname(os.path.abspath(__file__))
LOGNAME = "ref_g1_birkhoff.log"
N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
K = int(sys.argv[2]) if len(sys.argv) > 2 else 12

# ---- truncated power series in q (lists of length N+1) ----
def S(c=0):
    return [Fr(c)] + [Fr(0)] * N
def add(a, b):
    return [x + y for x, y in zip(a, b)]
def sub(a, b):
    return [x - y for x, y in zip(a, b)]
def mul(a, b):
    r = [Fr(0)] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                if b[j]:
                    r[i + j] += x * b[j]
    return r
def inv(a):
    r = [Fr(0)] * (N + 1)
    r[0] = 1 / a[0]
    for n in range(1, N + 1):
        r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1) if a[k]) / a[0]
    return r
def Dq(a):
    return [n * x for n, x in enumerate(a)]
def iszero(a):
    return all(x == 0 for x in a)

# ---- Q[[q]][ell]: dict power -> series ----
def Lmul_s(f, s):
    return {j: mul(a, s) for j, a in f.items()}
def Ladd(f, g):
    r = dict(f)
    for j, a in g.items():
        r[j] = add(r[j], a) if j in r else a
    return r
def LD(f):
    r = {}
    for j, a in f.items():
        r[j] = add(r[j], Dq(a)) if j in r else Dq(a)
        if j:
            t = [j * x for x in a]
            r[j - 1] = add(r[j - 1], t) if (j - 1) in r else t
    return r
def Lfree(f):
    """ell-free part if f is ell-free, else None"""
    for j, a in f.items():
        if j and not iszero(a):
            return None
    return f.get(0, S(0))
def Leq(f, g):
    for j in set(f) | set(g):
        if not iszero(sub(f.get(j, S(0)), g.get(j, S(0)))):
            return False
    return True

# ---- polynomials / series in w truncated at degree K ----
def wmul(p, r):
    o = [Fr(0)] * (K + 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(r):
                if i + j > K:
                    break
                o[i + j] += x * y
    return o
def winv(p):
    o = [Fr(0)] * (K + 1)
    o[0] = 1 / p[0]
    for n in range(1, K + 1):
        o[n] = -sum(p[k] * o[n - k] for k in range(1, min(n, len(p) - 1) + 1)) / p[0]
    return o

def periods(twist=5, equiv=True):
    """ctilde[m] = [w^m] Fhat as q-series, m = 0..K."""
    ct = [[Fr(0)] * (N + 1) for _ in range(K + 1)]
    for d in range(N + 1):
        num = [Fr(1)]
        for m in range(1, twist * d + 1):
            num = wmul(num, [Fr(m), Fr(twist)])
        den = [Fr(1)]
        for m in range(1, d + 1):
            pm = [Fr(comb5) * Fr(m) ** (5 - e) for e, comb5 in enumerate((1, 5, 10, 10, 5, 1))]   # (w+m)^5
            if equiv:
                pm[5] -= 1                                                                 # minus w^5
            den = wmul(den, pm)
        A = wmul(num, winv(den))
        for m in range(K + 1):
            ct[m][d] += A[m]
    return ct

def birkhoff(ct):
    C0 = ct[0]
    C0i = inv(C0)
    C1 = add(S(1), Dq(mul(ct[1], C0i)))          # D(log q + I~_1/C_0)
    C1i = inv(C1)
    c = [{j: [x / factorial(j) for x in ct[k - j]] for j in range(k + 1)} for k in range(K + 1)]
    phi = [Lmul_s(ck, C0i) for ck in c]
    dT = lambda f: Lmul_s(LD(f), C1i)
    kap, allphi = [], [phi]
    for a in range(5):
        ka = Lfree(dT(phi[a + 1]))
        if ka is None:
            return C0, C1, kap, allphi, False
        kap.append(ka)
        kai = inv(ka)
        new = [None] * (K + 1)
        for k in range(a + 1, K + 1):
            new[k] = Lmul_s(dT(phi[k]), kai)
        phi = new
        allphi.append(phi)
    return C0, C1, kap, allphi, True

# ---- differential operators Sum_j a_j D^j (left coefficients a_j in Q[[q]]) ----
def binom(n, k):
    return factorial(n) // (factorial(k) * factorial(n - k))
def ocomp(A, B):
    r = [S(0) for _ in range(len(A) + len(B) - 1)]
    for j, a in enumerate(A):
        for m, b in enumerate(B):
            Dl = b
            for l in range(j + 1):
                if l:
                    Dl = Dq(Dl)
                t = mul(a, [binom(j, l) * x for x in Dl])
                r[j - l + m] = add(r[j - l + m], t)
    return r
def oeq(A, B):
    n = max(len(A), len(B))
    A = A + [S(0)] * (n - len(A)); B = B + [S(0)] * (n - len(B))
    return all(iszero(sub(a, b)) for a, b in zip(A, B))


def main():
    lines = []
    def log(s):
        print(s); lines.append(s)
    bad = [0]
    def check(name, ok, must_fail=False):
        good = (not ok) if must_fail else ok
        if not good:
            bad[0] += 1
        tag = ("PASS (control fails as required)" if good else "FAIL (control did not fail)") if must_fail else ("PASS" if ok else "FAIL")
        log("%s  %s" % (tag, name))
    log("ref_g1_birkhoff: N = %d (q-order), K = %d (periods c_0..c_K); exact rationals, ell explicit" % (N, K))
    q = [Fr(0), Fr(1)] + [Fr(0)] * (N - 1)
    L5 = [Fr(3125) ** n for n in range(N + 1)]
    ct = periods()
    C0, C1, kap, allphi, ellfree = birkhoff(ct)
    check("every kappa_a (a = 0..4) is ell-free", ellfree and len(kap) == 5)
    C0i, C1i = inv(C0), inv(C1)
    k1 = mul(L5, inv(mul(mul(C0, C0), mul(C1, mul(C1, C1)))))
    k3 = mul(C0, C1i)
    check("kappa_0 = 1", kap[0] == S(1))
    check("kappa_1 = L^5/(C_0^2 C_1^3)", kap[1] == k1)
    check("kappa_2 = 1", kap[2] == S(1))
    check("kappa_3 = C_0/C_1", kap[3] == k3)
    check("kappa_4 = C_0/C_1", kap[4] == k3)
    phi0, phi5 = allphi[0], allphi[5]
    check("closing phi^(5)_(k+5) = phi^(0)_k, k = 0..%d (H^5 = 1)" % (K - 5), all(Leq(phi5[k + 5], phi0[k]) for k in range(K - 4)))
    check("product kappa_0..kappa_4 = (L/C_1)^5", mul(mul(mul(kap[0], kap[1]), mul(kap[2], kap[3])), kap[4]) == mul(L5, inv(mul(mul(C1, C1), mul(C1, mul(C1, C1))))))
    # operator identity
    dT = [S(0), C1i]
    PF4 = [mul(q, S(-120)), mul(q, S(-1250)), mul(q, S(-4375)), mul(q, S(-6250)), sub(S(1), mul(q, S(3125)))]
    def Mop(k1x):
        M = ocomp(dT, dT)
        M = ocomp(M, [inv(k1x)])
        M = ocomp(M, dT)
        M = ocomp(M, dT)
        return ocomp(M, [C0i])
    check("operator identity d_T^2 kappa_1^-1 d_T^2 C_0^-1 = (C_0/C_1) PF4 (all D-coefficients)", oeq(Mop(k1), [mul(k3, a) for a in PF4]))
    # PF4 c_k: constants (Lemma D)
    c = [{j: [x / factorial(j) for x in ct[kk - j]] for j in range(kk + 1)} for kk in range(5)]
    def applyop(A, f):
        r, g = {}, f
        for j, a in enumerate(A):
            if j:
                g = LD(g)
            r = Ladd(r, Lmul_s(g, a))
        return r
    pf = [applyop(PF4, c[kk]) for kk in range(5)]
    check("PF4 c_k = 0 (k <= 3), PF4 c_4 = 1", all(Leq(pf[kk], {0: S(0)}) for kk in range(4)) and Leq(pf[4], {0: S(1)}))
    # controls
    ctn = periods(equiv=False)
    C0n, C1n, kapn, allphin, okn = birkhoff(ctn)
    check("c1 non-equivariant: kappa_0..kappa_3 unchanged", okn and all(kapn[a] == kap[a] for a in range(4)))
    check("c1 non-equivariant: closing", okn and all(Leq(allphin[5][k + 5], allphin[0][k]) for k in range(K - 4)), must_fail=True)
    ct4 = periods(twist=4)
    C04, C14, kap4, _, ok4 = birkhoff(ct4)
    check("c2 twist O(4): kappa_1 = L^5/(C_0^2 C_1^3)", ok4 and kap4[1] == mul(L5, inv(mul(mul(C04, C04), mul(C14, mul(C14, C14))))), must_fail=True)
    check("c3 operator identity with kappa_1 -> 1/(C_0^2 C_1^3)", oeq(Mop(inv(mul(mul(C0, C0), mul(C1, mul(C1, C1))))), [mul(k3, a) for a in PF4]), must_fail=True)
    check("c4 kappa_3 = C_1/C_0", kap[3] == mul(C1, C0i), must_fail=True)
    log("first coefficients: kappa_1 = %s ; kappa_3 = %s" % ([str(x) for x in kap[1][:4]], [str(x) for x in kap[3][:4]]))
    log("TOTAL: %d unexpected result(s)" % bad[0])
    with open(os.path.join(HERE, LOGNAME), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("log written to " + LOGNAME)


if __name__ == "__main__":
    main()
