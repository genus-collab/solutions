"""Stage 2, D(iii): exact verification of every step of DICT_PROOF.md (independent of critic/cr_*.py).
Usage: py -I -B dict_verify.py WMAX [WSMALL]   (WMAX for T1, T4-T7; WSMALL for the monomial-level T2, T3)
  T1 dimension counts            (Lemma 0)   T2 L cap Q[v][p^+-1] = Q[p,q], both directions (Lemma 1)
  T3 I(M) in M, and route A = closed formulas (Lemma 2)
  T4 all legs in M                (Cor. 2)    T5 lambda mod v bijective (Lemma 3)   T6 lambda bijective (Thm, via Nakayama)
  T7 K_0 = log-free: dim, every leg in L (independent L_test), K_0 -> L/Q(v)[x] bijective (Theorem D(iii))
  C1-C3 failing controls.
"""
import os
import sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dict_core import (ZERO, mul, I_A, legs_A, PHI0, WT0, split, in_Qpq, in_M, xpow, weight_set)  # noqa: E402
from dict_nf import NF, rank_Q, kernel_Q, L_test  # noqa: E402

WMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
WSMALL = int(sys.argv[2]) if len(sys.argv) > 2 else 40
KMAX = WMAX // 4 + 1
FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL", msg)


def triples(w):
    return [(b, k, (w - WT0[b] - 4 * k) // 5) for b in range(5) for k in range(w // 4 + 1)
            if w - WT0[b] - 4 * k >= 0 and (w - WT0[b] - 4 * k) % 5 == 0]


def nvx(w):
    return sum(1 for n in range(w // 4 + 1) if (w - 4 * n) % 5 == 0)


def shift(A, m):
    return {(j, mm + m, l): c for (j, mm, l), c in A.items()}


def t1():
    for w in range(WMAX + 1):
        S = triples(w)
        n1 = sum(1 for t in S if t[0] == 1)
        nn = len(S) - n1
        f = w // 5 + 1 - (1 if w % 5 == 0 else 0)
        npq = w // 4 + 1
        check(nn == f, "T1 w=%d #b!=1 triples %d != %d" % (w, nn, f))
        check(n1 == nvx(w), "T1 w=%d #b=1 triples != dim Q[v,x]_w" % w)
        check(npq - nvx(w) == f, "T1 w=%d dim Q[p,q]_w - dim Q[v,x]_w != formula" % w)
    print("T1 counts, w <= %d: done" % WMAX)


def t2():
    n_in = n_out = 0
    for w in range(WSMALL + 1):
        for m in range(0, w // 5 + w // 4 + 4):
            j = w - 5 * m
            ok = L_test({(j, m): Fr(1)}, w)
            check(ok == (j >= -m), "T2 w=%d v^%d p^%d: L_test %s" % (w, m, j, ok))
            n_in += ok
            n_out += not ok
    print("T2 monomials v^m p^j, w <= %d: %d in L (all with j >= -m), %d not in L (all with j < -m)" % (WSMALL, n_in, n_out))


def closed_I(a, c):
    """Route B closed formula for I(p^a q^c) as Laurent-log dict."""
    out = {}
    e1, e2 = a - c + 4, a - c - 1
    out[(e1, c, 0) if e1 else (0, c, 1)] = -Fr(1, e1) if e1 else Fr(-1)
    out[(e2, c + 1, 0) if e2 else (0, c + 1, 1)] = Fr(1, e2) if e2 else Fr(1)
    return out


def t3():
    nb = 0
    for w in range(WSMALL + 1):
        for c in range(w // 4 + 1):
            a = w - 4 * c
            mono = {(a - c, c, 0): Fr(1)}
            r = I_A(mono)
            ok, why = in_M(r)
            check(ok, "T3 I(p^%d q^%d) not in M: %s" % (a, c, why))
            check(r == closed_I(a, c), "T3 I(p^%d q^%d) != closed formula" % (a, c))
            nb += 1
        for n in range(w // 4 + 1):
            if (w - 4 * n) % 5:
                continue
            m = (w - 4 * n) // 5
            el = mul({(0, m, 1): Fr(1)}, xpow(n))
            r = I_A(el)
            ok, why = in_M(r)
            check(ok, "T3 I(log p v^%d x^%d) not in M: %s" % (m, n, why))
            lg = split(r).get(1, {})
            want = {(j, mm + m): c / (n + 1) for (j, mm, l), c in xpow(n + 1).items()}
            # expected: log coefficient = -v^m x^(n+1)/(n+1) + (a constant in p)
            diff = {k: lg.get(k, ZERO) + want.get(k, ZERO) for k in set(lg) | set(want)}
            diff = {k: c for k, c in diff.items() if c != 0}
            check(all(j == 0 for (j, mm) in diff), "T3 log coefficient of I(log p x^%d) wrong" % n)
            nb += 1
    print("T3 I on the monomial basis of M_w, w <= %d: %d basis elements, all images in M" % (WSMALL, nb))


def t4(T):
    for (b, k), A in T.items():
        ok, why = in_M(A)
        check(ok, "T4 leg b=%d k=%d not in M: %s" % (b, k, why))
        check(weight_set(A) <= {WT0[b] + 4 * k}, "T4 leg b=%d k=%d inhomogeneous" % (b, k))
    print("T4 legs Phi_k^b, k <= %d, all b: in M, homogeneous" % KMAX)


def lam(T, w):
    S = triples(w)
    return S, [NF(shift(T[(b, k)], m)) for (b, k, m) in S]


def t5_t6_t7(T, label="", verbose=True):
    bad = []
    tot = 0
    for w in range(WMAX + 1):
        S, rows = lam(T, w)
        # lambda mod v: F/vF (triples with m = 0) -> N/vN (NF coordinates at v = 0); dim N/vN = [w>=1] + [4|w]
        rows0 = [{k: c for k, c in r.items() if (k[2] if k[0] == 'p' else k[1]) == 0}
                 for r, t in zip(rows, S) if t[2] == 0]
        r0 = rank_Q(rows0)
        d0 = (1 if w >= 1 else 0) + (1 if w % 4 == 0 else 0)
        r1 = rank_Q(rows)
        logrows = [{k: c for k, c in r.items() if k[0] == 'L'} for r in rows]
        kern = kernel_Q(logrows)
        dimK = len(kern)
        f = w // 5 + 1 - (1 if w % 5 == 0 else 0)
        legsK, nfK = [], []
        for y in kern:
            A = {}
            for i, c in y.items():
                b, k, m = S[i]
                for key, d in shift(T[(b, k)], m).items():
                    A[key] = A.get(key, ZERO) + c * d
            A = {key: c for key, c in A.items() if c != 0}
            sp = split(A)
            assert set(sp) <= {0}, "kernel vector not log-free"
            legsK.append(sp.get(0, {}))
            nfK.append(NF(A))
        nL = sum(1 for F in legsK if L_test(F, w))
        nQ = sum(1 for F in legsK if in_Qpq(F))
        rK = rank_Q(nfK)
        ok = (r0 == len(rows0) == d0 and r1 == len(S) and dimK == f and nL == dimK and nQ == dimK and rK == f)
        tot += dimK
        if not ok:
            bad.append(w)
        if verbose and (w <= 20 or w % 40 == 0 or not ok):
            print(" %s w=%3d |S|=%2d rank(lambda mod v)=%d/%d rank(lambda)=%2d dim K_0=%2d (formula %2d) in L %d/%d "
                  "in Q[p,q] %d rank mod Q(v)[x]=%2d" % (label, w, len(S), r0, d0, r1, dimK, f, nL, dimK, nQ, rK))
    return bad, tot


if __name__ == "__main__":
    print("dict_verify.py WMAX=%d WSMALL=%d KMAX=%d (exact over Q)" % (WMAX, WSMALL, KMAX))
    t1()
    t2()
    t3()
    T = legs_A(KMAX)
    t4(T)
    bad, tot = t5_t6_t7(T, "main")
    check(not bad, "T5-T7 fail at weights %s" % bad)
    print("T5-T7 w <= %d: %d K_0 basis vectors; failing weights: %s" % (WMAX, tot, bad))
    print("== Controls (must FAIL)")
    T2c = legs_A(min(KMAX, 12), {**PHI0, 2: {(4, 0, 0): Fr(1, 4), (-1, 1, 0): Fr(1)}})   # Phi_0^2 := x
    W0 = WMAX
    WMAX = min(WMAX, 40)
    badc, _ = t5_t6_t7(T2c, "C2", verbose=False)
    print(" C2 (Phi_0^2 := x, a function of x): first failing weights", badc[:8], "(expected: non-empty, contains 4)")
    check(4 in badc, "C2 control did not fail at w = 4")
    T3c = legs_A(min(KMAX, 12), {**PHI0, 2: {(-6, 2, 0): Fr(1)}})                     # Phi_0^2 := v^2 p^-6
    out = [(b, k) for (b, k), A in T3c.items() if b == 2 and not in_M(A)[0]]
    print(" C3 (Phi_0^2 := v^2/p^6, outside Q[p,q]): legs outside M:", len(out), "of", len([1 for (b, k) in T3c if b == 2]))
    check(len(out) > 0, "C3: no leg outside M")
    nfail = 0
    for w in range(4, 41):
        S = triples(w)
        for (b, k, m) in S:
            if b == 2:
                A = split(T3c[(b, k)])
                if 1 not in A:
                    nfail += not L_test({(j, mm + m): c for (j, mm), c in A[0].items()}, w)
    print(" C3: log-free single legs v^m Phi_k^2 failing L_test (w <= 40):", nfail, "(expected > 0)")
    check(nfail > 0, "C3 L_test did not fail")
    WMAX = W0
    print(" C1 = the 'not in L' half of T2 (bad monomials v^m p^j, j < -m, all rejected by the linear-algebra test)")
    print("RESULT:", "ALL CHECKS PASSED" if not FAILS else "%d FAILURES" % len(FAILS))
